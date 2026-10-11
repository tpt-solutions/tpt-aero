//! Horseshoe VLM solver: AIC assembly, native linear solve, and force
//! reconstruction.
//!
//! Conventions (all verified against the lumped-vortex flat-plate analytic
//! result Γ = π c V α):
//!
//! * Body frame: *x* aft (LE→TE), *y* right span, *z* up.
//! * Nose-up α ⇒ freestream in body frame `V∞ = V·(cos α, 0, sin α)` — the
//!   relative wind comes *from below* and moves *up* through the body frame.
//! * Horseshoe with circulation Γ > 0 flows `a → b` (left-to-right bound
//!   vortex, trailing legs to +∞ along body +x), which produces downwash at
//!   its own collocation point and upward Kutta–Joukowski force.
//! * Panel force `F = ρ Γ (V_loc × (b − a))`; lift is the component along
//!   `l̂ = (−sin α, 0, cos α)`, drag along `V̂∞`; pitching moment
//!   `M_y = (r × F)_y` is positive nose-up.

use alloc::{vec, vec::Vec};

use tpt_aero_core::{Real, Vec3};

use crate::aircraft::{ControlRegion, ControlSymmetry};
use crate::lattice::Lattice;
use crate::linalg::{lu_solve, lu_solve_generic, LuFactors, SolveError};
use crate::vortex::{horseshoe_velocity, semi_infinite_velocity};

/// Errors from a VLM analysis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VlmError {
    /// Lattice construction failed.
    Lattice(crate::lattice::LatticeError),
    /// The AIC matrix was numerically singular.
    SingularAic,
    /// Dimension mismatch (e.g. twist vector length ≠ element count).
    DimensionMismatch,
    /// Non-finite angle of attack, or non-positive/non-finite speed or density.
    InvalidFlightCondition,
    /// A cached (prepared) solve was requested on a model with a
    /// [`WakeModel::Freestream`] wake. The cache assumes the body-axis wake,
    /// so solve each angle with [`VlmModel::solve`] instead.
    FreestreamWakeNotPrepared,
    /// The ground-plane height is not finite and positive, or the geometry
    /// reaches down to or below the ground plane.
    InvalidGroundHeight,
    /// The reference area, chord or span is not finite and positive, so the
    /// coefficients would be meaningless.
    InvalidReference,
}

/// Direction of the horseshoe trailing legs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum WakeModel {
    /// Trailing legs along body `+x` for every angle of attack (the default).
    /// Keeps the AIC independent of the flow state, so the model is exactly
    /// linear in V∞ and the prepared/cached solve paths are available.
    #[default]
    BodyAxis,
    /// Trailing legs along the freestream direction `(cos α, 0, sin α)`. The
    /// AIC then depends on α, so only [`VlmModel::solve`] (and
    /// [`VlmModel::solve_generic`]) support it. Identical to
    /// [`WakeModel::BodyAxis`] at α = 0.
    Freestream,
}

impl core::fmt::Display for VlmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            VlmError::Lattice(e) => write!(f, "lattice: {e}"),
            VlmError::SingularAic => write!(f, "AIC matrix is numerically singular"),
            VlmError::DimensionMismatch => {
                write!(f, "twist vector length must equal element count")
            }
            VlmError::InvalidFlightCondition => write!(
                f,
                "angle of attack must be finite; speed and density must be finite and positive"
            ),
            VlmError::FreestreamWakeNotPrepared => write!(
                f,
                "prepared solves need the body-axis wake; solve freestream-wake models per angle"
            ),
            VlmError::InvalidReference => {
                write!(
                    f,
                    "reference area, chord and span must be finite and positive"
                )
            }
            VlmError::InvalidGroundHeight => {
                write!(f, "ground height must be finite and positive")
            }
        }
    }
}

impl core::error::Error for VlmError {}

impl From<SolveError> for VlmError {
    fn from(e: SolveError) -> Self {
        match e {
            SolveError::Singular | SolveError::NotPositiveDefinite => VlmError::SingularAic,
            SolveError::DimensionMismatch => VlmError::DimensionMismatch,
        }
    }
}

/// Reference quantities used to normalize forces into coefficients.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ReferenceGeometry {
    /// Reference planform area (m²).
    pub area: f64,
    /// Reference chord (m) — mean aerodynamic chord by convention.
    pub chord: f64,
    /// Reference span (m).
    pub span: f64,
    /// Moment reference point (body frame, m).
    pub moment_center: Vec3<f64>,
}

impl ReferenceGeometry {
    /// Derive the standard references from a planform (area, span, MAC;
    /// moment center at the MAC quarter chord on the MAC span station).
    #[must_use]
    pub fn from_planform(planform: &tpt_aero_core::WingPlanform) -> Self {
        let mut num_x = 0.0;
        let mut num_y = 0.0;
        let mut den = 0.0;
        for w in planform.sections().windows(2) {
            let dy = w[1].le.y - w[0].le.y;
            let (c0, c1) = (w[0].chord, w[1].chord);
            let mac = (2.0 / 3.0) * (c0 + c1 - c0 * c1 / (c0 + c1));
            let area = 0.5 * (c0 + c1) * dy;
            let y_mac = dy * (c0 + 2.0 * c1) / (3.0 * (c0 + c1));
            let x_le_mac = w[0].le.x + (y_mac / dy) * (w[1].le.x - w[0].le.x);
            num_x += (x_le_mac + 0.25 * mac) * area;
            num_y += (w[0].le.y + y_mac) * area;
            den += area;
        }
        let (x_qc, y_mac) = if den.abs() < 1.0e-12 {
            (0.0, 0.0)
        } else {
            (num_x / den, num_y / den)
        };
        Self {
            area: planform.area(),
            chord: planform.mac(),
            span: planform.span(),
            moment_center: Vec3::new(x_qc, y_mac, 0.0),
        }
    }
}

/// Non-dimensional aerodynamic result of a VLM analysis (primal `f64`).
///
/// # Near-zero-lift caveat
///
/// `cdi` is a small residual of the discrete force, not a separately modeled
/// quantity. Oswald-style ratios built from it, `C_L² / (π·AR·C_Di)`, are
/// ill-conditioned as `C_L → 0`, and should be reported only away from zero
/// lift (|C_L| ≳ 0.1 for the lattices used in the examples). Likewise, a
/// cambered section's zero-lift angle is offset from the thin-airfoil value by
/// roughly `1/AR` (about −0.05° at AR 10 for NACA m = 0.02, p = 0.4, against
/// −2.077° thin-airfoil), so small-|C_L| results near that angle carry that
/// offset.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct VlmSolution {
    /// Total lift coefficient.
    pub cl: f64,
    /// Induced drag coefficient (near-field, vortex drag only). See the
    /// near-zero-lift caveat on this type.
    pub cdi: f64,
    /// Pitching-moment coefficient about the reference point (nose-up +).
    pub cm: f64,
    /// Rolling-moment coefficient about the reference point.
    pub cl_roll: f64,
    /// Yawing-moment coefficient about the reference point.
    pub cn: f64,
    /// Horseshoe strengths Γ (m²/s), element-ordered.
    pub gamma: Vec<f64>,
    /// Section lift coefficients per spanwise strip (from the −y tip to the +y tip).
    pub strip_cl: Vec<f64>,
    /// Total dimensional force (N) in body axes.
    pub force: Vec3<f64>,
}

impl VlmSolution {
    /// Roll and yaw moment coefficients in aircraft-standard signs (roll
    /// positive right-wing-down, yaw positive nose-right) **about the
    /// stability axes** at angle of attack `alpha_rad`, from the body-axis
    /// `cl_roll` and `cn` of this solution. Use these, not the raw body-axis
    /// fields, when forming stability or control derivatives.
    #[must_use]
    pub fn stability_roll_yaw(&self, alpha_rad: f64) -> (f64, f64) {
        let (s, c) = alpha_rad.sin_cos();
        (
            -(self.cl_roll * c + self.cn * s),
            self.cl_roll * s - self.cn * c,
        )
    }
}

/// Generic (`Real`) result fields — the differentiable path returns these so
/// gradients flow through every quantity.
#[derive(Debug, Clone)]
pub struct VlmSolutionFields<R: Real> {
    /// Total lift coefficient.
    pub cl: R,
    /// Induced drag coefficient.
    pub cdi: R,
    /// Pitching-moment coefficient about the reference point.
    pub cm: R,
    /// Rolling-moment coefficient.
    pub cl_roll: R,
    /// Yawing-moment coefficient.
    pub cn: R,
    /// Horseshoe strengths Γ (m²/s).
    pub gamma: Vec<R>,
}

/// Attitude and body rates for [`PreparedVlm::solve_state`], in aircraft-standard
/// conventions: `alpha` nose-up positive, `beta` positive when the relative
/// wind comes from the right of the nose, `p` right wing down, `q` nose up,
/// `r` nose right (rad/s, about the moment reference point). Rates are about
/// the **stability axes** (body axes rotated by `alpha` about `y`), the
/// convention of stability derivatives.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FlightState {
    /// Angle of attack (rad).
    pub alpha: f64,
    /// Sideslip angle (rad).
    pub beta: f64,
    /// Roll rate (rad/s).
    pub p: f64,
    /// Pitch rate (rad/s).
    pub q: f64,
    /// Yaw rate (rad/s).
    pub r: f64,
}

/// Wind-axis force and stability-axis moment coefficients from
/// [`PreparedVlm::solve_state`], all in aircraft-standard signs: side force
/// positive to the right, roll moment positive right-wing-down, pitch
/// positive nose-up, yaw positive nose-right.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct StateLoads {
    /// Lift coefficient (normal to the relative wind, in the plane of symmetry).
    pub c_lift: f64,
    /// Drag coefficient along the relative wind (near-field, induced only).
    pub c_drag: f64,
    /// Side-force coefficient.
    pub c_side: f64,
    /// Pitching-moment coefficient (reference chord).
    pub c_pitch: f64,
    /// Rolling-moment coefficient (reference span).
    pub c_roll: f64,
    /// Yawing-moment coefficient (reference span).
    pub c_yaw: f64,
    /// Dimensional force on the surface (N), model frame (x aft, y right,
    /// z up), before any wind-axis projection.
    pub force_n: [f64; 3],
    /// Dimensional moment (N·m) about the model's moment reference point, model frame.
    pub moment_nm: [f64; 3],
}

/// A VLM model with its geometry-only work done up front: the LU factors of
/// the AIC matrix and the Gamma-to-local-velocity influence coefficients.
///
/// Both depend only on the lattice (the trailing legs are fixed along body
/// +x), not on alpha, V or rho, so an angle-of-attack sweep pays the O(N^3)
/// factorization and O(N^2) kernel evaluations once and each angle costs only
/// O(N^2) arithmetic. Results match [`VlmModel::solve`] with zero twist.
/// Obtain one with [`VlmModel::prepare`].
/// Twist-independent horseshoe geometry of one model, reusable across twists.
///
/// Built by [`VlmModel::influence_kernel`]; consumed by
/// [`VlmModel::prepare_with_kernel`]. Holds `N²` velocity vectors for the AIC
/// and `N²` for the load influence, so it costs memory proportional to `N²`.
#[derive(Debug, Clone)]
pub struct InfluenceKernel {
    /// Number of elements `N`.
    n: usize,
    /// Row-major `N x N`: velocity at collocation point `i` per unit Γ_j (1/m).
    aic: Vec<Vec3<f64>>,
    /// Row-major `N x N`: velocity at element `j`'s bound midpoint per unit Γ_k.
    load: Vec<Vec3<f64>>,
}

impl InfluenceKernel {
    /// Number of elements the kernel was built for.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.n
    }

    /// True for an empty lattice (never produced by a constructed model).
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.n == 0
    }
}

#[derive(Debug)]
pub struct PreparedVlm<'a> {
    model: &'a VlmModel,
    lu: LuFactors,
    /// Control-point normal of each element.
    normals: Vec<Vec3<f64>>,
    /// Row-major `N x N`: velocity at element `j`'s bound midpoint per unit Gamma_k.
    influence: Vec<Vec3<f64>>,
}

impl PreparedVlm<'_> {
    /// Solve one flight condition (zero twist).
    ///
    /// # Errors
    /// [`VlmError::InvalidFlightCondition`] for a non-finite angle or a
    /// non-positive/non-finite speed or density.
    pub fn solve(&self, alpha_rad: f64, v_inf: f64, density: f64) -> Result<VlmSolution, VlmError> {
        if !(alpha_rad.is_finite()
            && v_inf.is_finite()
            && v_inf > 0.0
            && density.is_finite()
            && density > 0.0)
        {
            return Err(VlmError::InvalidFlightCondition);
        }
        let n = self.normals.len();
        let (s, c) = alpha_rad.sin_cos();
        let v = Vec3::new(v_inf * c, 0.0, v_inf * s);
        let rhs: Vec<f64> = self.normals.iter().map(|ni| -v.dot(ni)).collect();
        let gamma = self.lu.solve(&rhs)?;
        let vlocs: Vec<Vec3<f64>> = (0..n)
            .map(|j| {
                let row = &self.influence[j * n..(j + 1) * n];
                row.iter()
                    .zip(&gamma)
                    .fold(v, |acc, (u, g)| acc + u.scaled(*g))
            })
            .collect();
        let m = self.model;
        let force = m.total_force::<f64>(&gamma, &vlocs, density);
        let moment = m.total_moment::<f64>(&gamma, &vlocs, density);
        let (cl, cdi, cm, cl_roll, cn) =
            m.coefficients::<f64>(&force, &moment, &alpha_rad, v_inf, density);
        let strip_cl = m.strip_coefficients(&gamma, &vlocs, v_inf, density, alpha_rad);
        Ok(VlmSolution {
            cl,
            cdi,
            cm,
            cl_roll,
            cn,
            gamma,
            strip_cl,
            force,
        })
    }

    /// Solve a general flight state: attitude (`alpha`, `beta`) and body rates
    /// about the model's moment reference point, zero twist.
    ///
    /// Rotation adds `-omega x r` to the air velocity at every control point
    /// and bound-vortex midpoint. Trailing legs stay along body +x (as in
    /// [`VlmModel::solve`]), which is accurate for small sideslip.
    ///
    /// # Errors
    /// [`VlmError::InvalidFlightCondition`] for non-finite state or a
    /// non-positive/non-finite speed or density.
    pub fn solve_state(
        &self,
        state: &FlightState,
        v_inf: f64,
        density: f64,
    ) -> Result<StateLoads, VlmError> {
        let FlightState {
            alpha,
            beta,
            p,
            q,
            r,
        } = *state;
        if !([alpha, beta, p, q, r].iter().all(|x| x.is_finite())
            && v_inf.is_finite()
            && v_inf > 0.0
            && density.is_finite()
            && density > 0.0)
        {
            return Err(VlmError::InvalidFlightCondition);
        }
        let m = self.model;
        m.check_ground()?;
        let refg = m.reference;
        let (sa, ca) = alpha.sin_cos();
        let (sb, cb) = beta.sin_cos();
        // Model frame: x aft, y right, z up. A wind from the right (beta > 0)
        // moves toward -y relative to the body.
        let v = Vec3::new(v_inf * ca * cb, -v_inf * sb, v_inf * sa * cb);
        // Stability axes: x_s forward along the relative wind projected on the
        // symmetry plane, z_s down. In the model frame (x aft, z up) that is
        // forward = (-cos a, 0, -sin a) and down = (sin a, 0, -cos a), so the
        // standard (p, q, r) become an angular velocity in the model frame.
        let forward = Vec3::new(-ca, 0.0, -sa);
        let down = Vec3::new(sa, 0.0, -ca);
        let omega = forward.scaled(p) + Vec3::new(0.0, q, 0.0) + down.scaled(r);
        let air_at = |point: Vec3<f64>| v - omega.cross(&(point - refg.moment_center));

        let elements = m.lattice.elements();
        let n = elements.len();
        let rhs: Vec<f64> = elements
            .iter()
            .zip(&self.normals)
            .map(|(e, ni)| -air_at(e.collocation).dot(ni))
            .collect();
        let gamma = self.lu.solve(&rhs)?;
        let vlocs: Vec<Vec3<f64>> = (0..n)
            .map(|j| {
                let (p0, p1) = elements[j].bounds::<f64>();
                let row = &self.influence[j * n..(j + 1) * n];
                row.iter()
                    .zip(&gamma)
                    .fold(air_at((p0 + p1).scaled(0.5)), |acc, (u, g)| {
                        acc + u.scaled(*g)
                    })
            })
            .collect();
        let force = m.total_force::<f64>(&gamma, &vlocs, density);
        let moment = m.total_moment::<f64>(&gamma, &vlocs, density);

        // Wind axes: drag along the relative wind, lift normal to it in the
        // plane of symmetry, side force completing the right-handed triad.
        let drag_dir = v.normalized();
        // Lift is drag × ŷ. At |β| = 90° the wind is along ŷ and the cross
        // product vanishes, so take its limit from the symmetry plane instead:
        // the unit vector normal to the relative wind in the x–z plane.
        let lift_cross = drag_dir.cross(&Vec3::new(0.0, 1.0, 0.0));
        let lift_dir = if lift_cross.dot(&lift_cross) > 1.0e-24 {
            lift_cross.normalized()
        } else {
            Vec3::new(-sa, 0.0, ca)
        };
        let side_dir = lift_dir.cross(&drag_dir);
        let q_s = 0.5 * density * v_inf * v_inf * refg.area * m.beta;
        Ok(StateLoads {
            c_lift: force.dot(&lift_dir) / q_s,
            c_drag: force.dot(&drag_dir) / q_s,
            c_side: force.dot(&side_dir) / q_s,
            c_pitch: moment.y / (q_s * refg.chord),
            // Moments about the stability axes (see above).
            c_roll: moment.dot(&forward) / (q_s * refg.span),
            c_yaw: moment.dot(&down) / (q_s * refg.span),
            force_n: [force.x, force.y, force.z],
            moment_nm: [moment.x, moment.y, moment.z],
        })
    }
}

/// Everything one solve produces that downstream post-processing reuses.
struct SolveCore<R: Real> {
    fields: VlmSolutionFields<R>,
    /// Local velocity at each bound-vortex midpoint.
    vlocs: Vec<Vec3<R>>,
    /// Total dimensional body-frame force (N).
    force: Vec3<R>,
}

/// A horseshoe-vortex VLM model: lattice + reference geometry.
#[derive(Debug, Clone)]
pub struct VlmModel {
    lattice: Lattice,
    reference: ReferenceGeometry,
    controls: Vec<ControlRegion>,
    wake: WakeModel,
    ground_height: Option<f64>,
    /// Göthert/Prandtl–Glauert factor `β = √(1 − M²)`; `1` when incompressible.
    beta: f64,
}

impl VlmModel {
    /// Build a VLM model from a planform with the given discretization.
    ///
    /// # Errors
    /// Propagates [`crate::LatticeError`].
    pub fn from_planform(
        planform: &tpt_aero_core::WingPlanform,
        n_span: usize,
        n_chord: usize,
        camber: crate::lattice::CamberModel,
    ) -> Result<Self, VlmError> {
        Self::from_planform_spaced(
            planform,
            n_span,
            n_chord,
            camber,
            crate::lattice::SpanSpacing::Uniform,
        )
    }

    /// Like [`VlmModel::from_planform`], with an explicit spanwise strip
    /// distribution (see [`crate::SpanSpacing`]).
    ///
    /// # Errors
    /// Propagates [`crate::LatticeError`].
    pub fn from_planform_spaced(
        planform: &tpt_aero_core::WingPlanform,
        n_span: usize,
        n_chord: usize,
        camber: crate::lattice::CamberModel,
        spacing: crate::lattice::SpanSpacing,
    ) -> Result<Self, VlmError> {
        let lattice = Lattice::build_spaced(planform, n_span, n_chord, camber, spacing)
            .map_err(VlmError::Lattice)?;
        Ok(Self {
            lattice,
            reference: ReferenceGeometry::from_planform(planform),
            controls: Vec::new(),
            wake: WakeModel::BodyAxis,
            ground_height: None,
            beta: 1.0,
        })
    }

    /// Build a subsonic compressible model at free-stream Mach number `mach`
    /// by Göthert's rule: the planform is stretched spanwise by
    /// `β = √(1 − M²)`, solved incompressibly, and every force and moment
    /// coefficient is divided by `β`.
    ///
    /// This reproduces the Helmbold–Diederich lift slope
    /// `C_Lα = 2πAR / (2 + √(4 + (βAR)²))` for a straight wing. Coefficients are
    /// referenced to the stretched planform (`S' = βS`, `b' = βb`), which is
    /// what the Göthert transformation requires. Use `V∞` as the true speed.
    ///
    /// # Errors
    /// [`VlmError::InvalidFlightCondition`] unless `0 ≤ mach < 1`;
    /// [`VlmError::Lattice`] if the stretched planform is invalid.
    pub fn from_planform_compressible(
        planform: &tpt_aero_core::WingPlanform,
        n_span: usize,
        n_chord: usize,
        camber: crate::lattice::CamberModel,
        mach: f64,
    ) -> Result<Self, VlmError> {
        if !(mach.is_finite() && (0.0..1.0).contains(&mach)) {
            return Err(VlmError::InvalidFlightCondition);
        }
        let beta = (1.0 - mach * mach).sqrt();
        let sections = planform
            .sections()
            .iter()
            .map(|sec| tpt_aero_core::WingSection {
                le: tpt_aero_core::Vec3::new(sec.le.x, sec.le.y * beta, sec.le.z),
                chord: sec.chord,
                twist: sec.twist,
            })
            .collect();
        let stretched = tpt_aero_core::WingPlanform::from_sections(sections)
            .map_err(|_| VlmError::Lattice(crate::lattice::LatticeError::InvalidPlanform))?;
        let mut model = Self::from_planform_spaced(
            &stretched,
            n_span,
            n_chord,
            camber,
            crate::lattice::SpanSpacing::Uniform,
        )?;
        model.beta = beta;
        Ok(model)
    }

    /// Göthert compressible model of a multi-surface aircraft. Horizontal
    /// surfaces (roll 0) are stretched in body *y* by `β`; fins (roll π/2) lie
    /// in the x–z plane and are unchanged by the transformation. The reference
    /// comes from the first surface, which should be the wing.
    ///
    /// # Errors
    /// [`VlmError::InvalidFlightCondition`] unless `0 ≤ mach < 1` or a
    /// surface is rolled to an angle other than 0 or π/2;
    /// [`VlmError::Lattice`] if a stretched planform is invalid.
    pub fn from_surfaces_compressible(
        surfaces: &[crate::aircraft::Surface],
        mach: f64,
    ) -> Result<Self, VlmError> {
        if !(mach.is_finite() && (0.0..1.0).contains(&mach)) {
            return Err(VlmError::InvalidFlightCondition);
        }
        let beta = (1.0 - mach * mach).sqrt();
        let half_pi = core::f64::consts::FRAC_PI_2;
        let stretched = surfaces
            .iter()
            .map(|s| {
                let factor = if s.roll_rad.abs() < 1.0e-12 {
                    beta
                } else if (s.roll_rad - half_pi).abs() < 1.0e-12 {
                    1.0
                } else {
                    return Err(VlmError::InvalidFlightCondition);
                };
                let sections: Vec<tpt_aero_core::WingSection> = s
                    .planform
                    .sections()
                    .iter()
                    .map(|sec| tpt_aero_core::WingSection {
                        le: Vec3::new(sec.le.x, sec.le.y * factor, sec.le.z),
                        chord: sec.chord,
                        twist: sec.twist,
                    })
                    .collect();
                let planform =
                    tpt_aero_core::WingPlanform::from_sections(sections).map_err(|_| {
                        VlmError::Lattice(crate::lattice::LatticeError::InvalidPlanform)
                    })?;
                let mut out = s.clone();
                out.planform = planform;
                out.offset = Vec3::new(s.offset.x, s.offset.y * factor, s.offset.z);
                Ok(out)
            })
            .collect::<Result<Vec<_>, VlmError>>()?;
        let mut model = Self::from_surfaces(&stretched)?;
        model.beta = beta;
        Ok(model)
    }

    /// The Göthert/Prandtl–Glauert factor `β` (`1` for an incompressible model).
    #[must_use]
    pub const fn compressibility_beta(&self) -> f64 {
        self.beta
    }

    /// Choose the trailing-leg direction (see [`WakeModel`]). Defaults to
    /// [`WakeModel::BodyAxis`].
    #[must_use]
    pub fn with_wake(mut self, wake: WakeModel) -> Self {
        self.wake = wake;
        self
    }

    /// Place the model over a flat ground plane `height` metres below its
    /// origin (body `z = −height`). The ground is modelled by image horseshoe
    /// vortices mirrored through the plane with opposite circulation, which
    /// reduces induced downwash and so the induced drag, and increases lift
    /// at a given angle of attack (ground effect).
    #[must_use]
    pub fn with_ground_height(mut self, height: f64) -> Self {
        self.ground_height = Some(height);
        self
    }

    /// The ground-plane height below the origin, if the model is in ground
    /// effect.
    #[must_use]
    pub const fn ground_height(&self) -> Option<f64> {
        self.ground_height
    }

    /// Image-vortex contribution of the ground plane at point `p` from a
    /// horseshoe `(a, b, d)` of circulation `g`; zero when out of ground
    /// effect. The image mirrors `z → −2h − z`, flips the direction's `z`
    /// component, and carries circulation `−g`.
    fn image_horseshoe<R: Real>(
        &self,
        a: &Vec3<R>,
        b: &Vec3<R>,
        d: &Vec3<R>,
        p: &Vec3<R>,
        g: R,
    ) -> Vec3<R> {
        match self.ground_height {
            None => Vec3::zero(),
            Some(h) => {
                let two_h = R::from_f64(-2.0 * h);
                let mirror = |v: &Vec3<R>| Vec3::new(v.x, v.y, two_h - v.z);
                let dm = Vec3::new(d.x, d.y, d.z.neg());
                horseshoe_velocity(&mirror(a), &mirror(b), &dm, p, g.neg())
            }
        }
    }

    /// Reject a ground height that is not finite and positive.
    fn check_ground(&self) -> Result<(), VlmError> {
        let r = self.reference;
        if !([r.area, r.chord, r.span]
            .iter()
            .all(|v| v.is_finite() && *v > 0.0))
        {
            return Err(VlmError::InvalidReference);
        }
        let Some(h) = self.ground_height else {
            return Ok(());
        };
        if !(h.is_finite() && h > 0.0) {
            return Err(VlmError::InvalidGroundHeight);
        }
        // The ground is the plane z = −h. Geometry at or below it would sit in
        // the image system's singularity, so require every bound point and
        // collocation point to lie strictly above it.
        let lowest = self
            .lattice
            .elements()
            .iter()
            .flat_map(|e| [e.bound_a.z, e.bound_b.z, e.collocation.z])
            .fold(f64::INFINITY, f64::min);
        if lowest.is_finite() && lowest <= -h {
            return Err(VlmError::InvalidGroundHeight);
        }
        Ok(())
    }

    /// The trailing-leg direction this model uses.
    #[must_use]
    pub const fn wake(&self) -> WakeModel {
        self.wake
    }

    /// Build from a pre-made lattice. The reference geometry is estimated from
    /// the **first surface** only (surface index 0), as in
    /// [`VlmModel::from_surfaces`]: its area, and its span taken as the sum of
    /// its strip widths (each strip has `n_chord` panels, so its width is
    /// counted once). Use [`VlmModel::with_reference`] to override.
    #[must_use]
    pub fn from_lattice(lattice: Lattice) -> Self {
        let n_chord = lattice.n_chord().max(1) as f64;
        let first = lattice.elements().iter().filter(|e| e.surface == 0);
        let area: f64 = first.clone().map(|e| e.area).sum();
        let span: f64 = first.map(|e| e.strip_dy).sum::<f64>() / n_chord;
        let chord = if area.abs() < 1.0e-12 {
            1.0
        } else {
            area / span.max(1.0e-12)
        };
        Self {
            controls: Vec::new(),
            wake: WakeModel::BodyAxis,
            ground_height: None,
            beta: 1.0,
            lattice,
            reference: ReferenceGeometry {
                area,
                chord,
                span,
                moment_center: Vec3::zero(),
            },
        }
    }

    /// The underlying lattice.
    #[must_use]
    pub const fn lattice(&self) -> &Lattice {
        &self.lattice
    }

    /// The reference geometry.
    #[must_use]
    pub const fn reference(&self) -> &ReferenceGeometry {
        &self.reference
    }

    /// Attach control-surface definitions (set by [`VlmModel::from_surfaces`]).
    #[must_use]
    pub fn with_controls(mut self, controls: Vec<ControlRegion>) -> Self {
        self.controls = controls;
        self
    }

    /// The control surfaces of this model, in the order deflections are given
    /// to [`VlmModel::control_twists`].
    #[must_use]
    pub fn controls(&self) -> &[ControlRegion] {
        &self.controls
    }

    /// Index of the control called `name`.
    #[must_use]
    pub fn control_index(&self, name: &str) -> Option<usize> {
        self.controls.iter().position(|c| c.control.name == name)
    }

    /// The per-element twist slice for a set of control deflections (rad, one
    /// per control in [`VlmModel::controls`] order), to pass as the `twists`
    /// argument of [`VlmModel::solve`].
    ///
    /// # Errors
    /// [`VlmError::DimensionMismatch`] if `deflections.len()` differs from
    /// the number of controls.
    pub fn control_twists(&self, deflections: &[f64]) -> Result<Vec<f64>, VlmError> {
        if deflections.len() != self.controls.len() {
            return Err(VlmError::DimensionMismatch);
        }
        let tol = 1.0e-9;
        Ok(self
            .lattice
            .elements()
            .iter()
            .map(|e| {
                self.controls
                    .iter()
                    .zip(deflections)
                    .filter(|(r, _)| {
                        let c = &r.control;
                        e.surface == r.surface
                            && e.chord_frac >= c.hinge
                            && e.eta.abs() >= c.y_in - tol
                            && e.eta.abs() <= c.y_out + tol
                    })
                    .map(|(r, d)| match r.control.symmetry {
                        ControlSymmetry::Symmetric => *d,
                        ControlSymmetry::Antisymmetric => {
                            if e.eta >= 0.0 {
                                *d
                            } else {
                                -*d
                            }
                        }
                    })
                    .sum()
            })
            .collect())
    }

    /// Override the reference geometry (e.g. to match a wind-tunnel model).
    ///
    /// The reference is given in physical units. On a compressible model
    /// (`β < 1`) the area and span are stretched by `β` internally, the same
    /// as the default reference, so that `S` in the coefficients is the true
    /// wing area and is not counted twice.
    #[must_use]
    pub fn with_reference(mut self, reference: ReferenceGeometry) -> Self {
        let b = self.beta;
        self.reference = ReferenceGeometry {
            area: reference.area * b,
            span: reference.span * b,
            ..reference
        };
        self
    }

    /// Assemble the AIC system in `R`, solve for the horseshoe strengths, and
    /// reconstruct near-field forces — the differentiable path.
    ///
    /// `alpha_rad` is nose-up positive, `v_inf` the freestream speed (m/s),
    /// `rho` the density (kg/m³), and `twists` the per-element geometric
    /// twist (rad). An empty `twists` slice means zero twist everywhere; for
    /// gradient work pass exactly one entry per lattice element, seeded as
    /// `Dual` variables.
    ///
    /// # Errors
    /// [`VlmError::SingularAic`] when the AIC matrix is numerically singular;
    /// [`VlmError::DimensionMismatch`] on inconsistent input sizes.
    pub fn solve_generic<R: Real>(
        &self,
        alpha_rad: R,
        v_inf: R,
        rho: R,
        twists: &[R],
    ) -> Result<VlmSolutionFields<R>, VlmError> {
        self.solve_core(alpha_rad, v_inf, rho, twists)
            .map(|core| core.fields)
    }

    /// The shared solve: AIC assembly, linear solve, and a single pass of
    /// local velocities that force, moment and strip loads all reuse.
    fn solve_core<R: Real>(
        &self,
        alpha_rad: R,
        v_inf: R,
        rho: R,
        twists: &[R],
    ) -> Result<SolveCore<R>, VlmError> {
        // Same flight-condition contract as `solve`, checked on the primal
        // value so the generic (dual) path cannot silently return a zero-speed
        // or non-finite solution.
        let (a_val, v_val, rho_val) = (alpha_rad.to_f64(), v_inf.to_f64(), rho.to_f64());
        if !(a_val.is_finite()
            && v_val.is_finite()
            && v_val > 0.0
            && rho_val.is_finite()
            && rho_val > 0.0)
        {
            return Err(VlmError::InvalidFlightCondition);
        }
        let elements = self.lattice.elements();
        let n = elements.len();
        if !twists.is_empty() && twists.len() != n {
            return Err(VlmError::DimensionMismatch);
        }
        self.check_ground()?;
        let zero = R::zero();
        // Effective twist = the planform's own geometric twist plus the
        // per-solve perturbation (the differentiable design variable).
        let twist_at =
            |i: usize| R::from_f64(elements[i].twist) + twists.get(i).copied().unwrap_or(zero);

        // Freestream (body frame) and trailing-leg direction. The default body
        // +x wake keeps the AIC matrix independent of the flow state, so the
        // model is exactly linear in V∞ (freestream reversal ⇒ exact force
        // reversal) and avoids upstream-leg singularities at large |α|. The
        // opt-in freestream wake aligns the legs with the flow instead.
        let (s, c) = alpha_rad.sin_cos();
        let v = Vec3::new(v_inf * c, zero, v_inf * s);
        let d = match self.wake {
            WakeModel::BodyAxis => Vec3::new(R::one(), zero, zero),
            WakeModel::Freestream => Vec3::new(c, zero, s),
        };

        // AIC: a_ij = n_i · v_j(cp_i) [1/m]; rhs_i = −V∞·n_i [m/s].
        let mut a: Vec<Vec<R>> = vec![vec![zero; n]; n];
        let mut rhs = vec![zero; n];
        for (i, ei) in elements.iter().enumerate() {
            let ni = ei.normal_for_twist(twist_at(i));
            let cp = ei.cp::<R>();
            rhs[i] = v.dot(&ni).neg();
            for (j, ej) in elements.iter().enumerate() {
                let (p0, p1) = ej.bounds::<R>();
                let vel = horseshoe_velocity(&p0, &p1, &d, &cp, R::one())
                    + self.image_horseshoe(&p0, &p1, &d, &cp, R::one());
                a[i][j] = ni.dot(&vel);
            }
        }

        // Native solve: f64 uses the tpt-math dense partial-pivot LU; other
        // Real types (dual numbers) use the from-scratch generic LU.
        let gamma: Vec<R> = if core::any::TypeId::of::<R>() == core::any::TypeId::of::<f64>() {
            let af: Vec<Vec<f64>> = a
                .iter()
                .map(|row| row.iter().map(|x| x.to_f64()).collect())
                .collect();
            let bf: Vec<f64> = rhs.iter().map(|x| x.to_f64()).collect();
            lu_solve(&af, &bf)?.into_iter().map(R::from_f64).collect()
        } else {
            lu_solve_generic::<R>(&a, &rhs)?
        };

        let vlocs = self.local_velocities::<R>(&gamma, &v, &d);
        let force = self.total_force::<R>(&gamma, &vlocs, rho);
        let moment = self.total_moment::<R>(&gamma, &vlocs, rho);
        let (cl, cdi, cm, cl_roll, cn) =
            self.coefficients::<R>(&force, &moment, &alpha_rad, v_inf, rho);
        Ok(SolveCore {
            fields: VlmSolutionFields {
                cl,
                cdi,
                cm,
                cl_roll,
                cn,
                gamma,
            },
            vlocs,
            force,
        })
    }

    /// Do the geometry-only work once for repeated zero-twist solves; see
    /// [`PreparedVlm`].
    ///
    /// # Errors
    /// [`VlmError::SingularAic`] when the AIC matrix is numerically singular.
    pub fn prepare(&self) -> Result<PreparedVlm<'_>, VlmError> {
        self.prepare_with_twists(&[])
    }

    /// [`VlmModel::prepare`] with per-element twist perturbations added to the
    /// geometric twist (rad, one per element, as from
    /// [`VlmModel::control_twists`]). Empty `twists` means zero perturbation.
    /// Only the surface normals and the AIC factorization depend on twist, so
    /// this is the entry point for control-deflected stability loads.
    ///
    /// # Errors
    /// [`VlmError::FreestreamWakeNotPrepared`], [`VlmError::InvalidGroundHeight`],
    /// [`VlmError::DimensionMismatch`] for a wrong-length twist slice, or
    /// [`VlmError::SingularAic`].
    pub fn prepare_with_twists(&self, twists: &[f64]) -> Result<PreparedVlm<'_>, VlmError> {
        let kernel = self.influence_kernel()?;
        self.prepare_with_kernel(&kernel, twists)
    }

    /// The twist-independent geometric influence of every horseshoe, computed
    /// once per model. Twist only changes the surface normals, so a caller that
    /// re-prepares many times (a six-DOF run changing control deflections, say)
    /// builds this once and passes it to [`VlmModel::prepare_with_kernel`].
    ///
    /// # Errors
    /// [`VlmError::FreestreamWakeNotPrepared`] for the freestream wake, and
    /// [`VlmError::InvalidGroundHeight`] / ground-clearance errors as in
    /// [`VlmModel::prepare`].
    pub fn influence_kernel(&self) -> Result<InfluenceKernel, VlmError> {
        if self.wake == WakeModel::Freestream {
            return Err(VlmError::FreestreamWakeNotPrepared);
        }
        self.check_ground()?;
        let elements = self.lattice.elements();
        let n = elements.len();
        let d = Vec3::new(1.0, 0.0, 0.0);
        // Velocity at collocation point i from element j, per unit Γ (row-major).
        let mut aic = Vec::with_capacity(n * n);
        for ei in elements {
            let cp = ei.cp::<f64>();
            for ej in elements {
                let (p0, p1) = ej.bounds::<f64>();
                aic.push(
                    horseshoe_velocity(&p0, &p1, &d, &cp, 1.0)
                        + self.image_horseshoe(&p0, &p1, &d, &cp, 1.0),
                );
            }
        }
        // Velocity at element j's bound midpoint from element k (row-major).
        let mut load = Vec::with_capacity(n * n);
        for (j, ej) in elements.iter().enumerate() {
            let (p0, p1) = ej.bounds::<f64>();
            let mid = (p0 + p1).scaled(0.5);
            for (k, ek) in elements.iter().enumerate() {
                let (ka, kb) = ek.bounds::<f64>();
                // Same load model as `local_velocity`: an element's own bound
                // segment is excluded, its trailing legs are kept.
                let direct = if k == j {
                    semi_infinite_velocity(&kb, &d, &mid, 1.0)
                        + semi_infinite_velocity(&ka, &d, &mid, -1.0)
                } else {
                    horseshoe_velocity(&ka, &kb, &d, &mid, 1.0)
                };
                load.push(direct + self.image_horseshoe(&ka, &kb, &d, &mid, 1.0));
            }
        }
        Ok(InfluenceKernel { n, aic, load })
    }

    /// Prepare the model for `twists` from a kernel built by
    /// [`VlmModel::influence_kernel`] on this same model. Only the normals and
    /// one LU factorisation depend on the twist; the horseshoe geometry is reused.
    ///
    /// # Errors
    /// [`VlmError::DimensionMismatch`] for a wrong-length twist slice or a
    /// kernel from a different lattice, or [`VlmError::SingularAic`].
    pub fn prepare_with_kernel(
        &self,
        kernel: &InfluenceKernel,
        twists: &[f64],
    ) -> Result<PreparedVlm<'_>, VlmError> {
        let elements = self.lattice.elements();
        let n = elements.len();
        if kernel.n != n || (!twists.is_empty() && twists.len() != n) {
            return Err(VlmError::DimensionMismatch);
        }
        let normals: Vec<Vec3<f64>> = elements
            .iter()
            .enumerate()
            .map(|(i, e)| e.normal_for_twist(e.twist + twists.get(i).copied().unwrap_or(0.0)))
            .collect();
        let mut a = vec![vec![0.0; n]; n];
        for (i, row) in a.iter_mut().enumerate() {
            for (j, aij) in row.iter_mut().enumerate() {
                *aij = normals[i].dot(&kernel.aic[i * n + j]);
            }
        }
        let lu = LuFactors::factor(&a)?;
        Ok(PreparedVlm {
            model: self,
            lu,
            normals,
            influence: kernel.load.clone(),
        })
    }

    /// Solve a whole angle-of-attack sweep (zero twist) with one
    /// [`VlmModel::prepare`]; typically an order of magnitude or more faster
    /// than calling [`VlmModel::solve`] per angle.
    ///
    /// # Errors
    /// As [`VlmModel::prepare`] and [`PreparedVlm::solve`].
    pub fn solve_sweep(
        &self,
        alphas_rad: &[f64],
        v_inf: f64,
        density: f64,
    ) -> Result<Vec<VlmSolution>, VlmError> {
        if self.wake == WakeModel::Freestream {
            // The wake depends on α, so there is nothing to cache across angles.
            return alphas_rad
                .iter()
                .map(|&a| self.solve(a, v_inf, density, &[]))
                .collect();
        }
        let prepared = self.prepare()?;
        alphas_rad
            .iter()
            .map(|&a| prepared.solve(a, v_inf, density))
            .collect()
    }

    /// Convenience primal solve: `f64` in, [`VlmSolution`] out.
    ///
    /// `alpha_rad` is nose-up positive, `v_inf` in m/s, `density` in kg/m³,
    /// `twists` per-element (empty ⇒ zero twist).
    ///
    /// # Errors
    /// See [`VlmModel::solve_generic`].
    pub fn solve(
        &self,
        alpha_rad: f64,
        v_inf: f64,
        density: f64,
        twists: &[f64],
    ) -> Result<VlmSolution, VlmError> {
        if !(alpha_rad.is_finite()
            && v_inf.is_finite()
            && v_inf > 0.0
            && density.is_finite()
            && density > 0.0)
        {
            return Err(VlmError::InvalidFlightCondition);
        }
        self.check_ground()?;
        let n = self.lattice.len();
        if !twists.is_empty() && twists.len() != n {
            return Err(VlmError::DimensionMismatch);
        }
        let SolveCore {
            fields,
            vlocs,
            force,
        } = self.solve_core::<f64>(alpha_rad, v_inf, density, twists)?;
        let strip_cl = self.strip_coefficients(&fields.gamma, &vlocs, v_inf, density, alpha_rad);
        Ok(VlmSolution {
            cl: fields.cl,
            cdi: fields.cdi,
            cm: fields.cm,
            cl_roll: fields.cl_roll,
            cn: fields.cn,
            gamma: fields.gamma,
            strip_cl,
            force,
        })
    }

    /// Velocity at the bound-vortex midpoint of element `j`: freestream plus
    /// every horseshoe influence except the panel's *own bound segment* (its
    /// own trailing legs are kept — they represent the wake). This is the
    /// standard VLM load model: the own-bound singularity is the only singular
    /// contribution at the midpoint, and dropping it leaves the physical
    /// induction. The 2D limit (single wide strip) then yields exactly zero
    /// induced drag, as it must.
    fn local_velocity<R: Real>(&self, j: usize, gamma: &[R], v: &Vec3<R>, d: &Vec3<R>) -> Vec3<R> {
        let elements = self.lattice.elements();
        let (p0, p1) = elements[j].bounds::<R>();
        let mid = (p0 + p1).scaled(R::half());
        let mut vloc = *v;
        for (k, ek) in elements.iter().enumerate() {
            let (ka, kb) = ek.bounds::<R>();
            if k == j {
                vloc = vloc
                    + semi_infinite_velocity(&kb, d, &mid, gamma[k])
                    + semi_infinite_velocity(&ka, d, &mid, gamma[k].neg());
            } else {
                vloc = vloc + horseshoe_velocity(&ka, &kb, d, &mid, gamma[k]);
            }
            vloc = vloc + self.image_horseshoe(&ka, &kb, d, &mid, gamma[k]);
        }
        vloc
    }

    /// [`Self::local_velocity`] for every element, computed once per solve.
    fn local_velocities<R: Real>(&self, gamma: &[R], v: &Vec3<R>, d: &Vec3<R>) -> Vec<Vec3<R>> {
        (0..self.lattice.len())
            .map(|j| self.local_velocity::<R>(j, gamma, v, d))
            .collect()
    }

    /// Total dimensional body-frame force (N): F = Σ ρ Γ (V_loc × (b − a)).
    fn total_force<R: Real>(&self, gamma: &[R], vlocs: &[Vec3<R>], rho: R) -> Vec3<R> {
        let elements = self.lattice.elements();
        let mut total = Vec3::<R>::zero();
        for (j, ej) in elements.iter().enumerate() {
            let (p0, p1) = ej.bounds::<R>();
            let vloc = vlocs[j];
            // b − a points +y for a lift-up Γ.
            let ell = p1 - p0;
            total = total + vloc.cross(&ell).scaled(rho * gamma[j]);
        }
        total
    }

    /// Moment vector (N·m) about the reference point; each panel force acts
    /// at its bound-vortex midpoint.
    fn total_moment<R: Real>(&self, gamma: &[R], vlocs: &[Vec3<R>], rho: R) -> Vec3<R> {
        let elements = self.lattice.elements();
        let r0 = Vec3::new(
            R::from_f64(self.reference.moment_center.x),
            R::from_f64(self.reference.moment_center.y),
            R::from_f64(self.reference.moment_center.z),
        );
        let mut total = Vec3::<R>::zero();
        for (j, ej) in elements.iter().enumerate() {
            let (p0, p1) = ej.bounds::<R>();
            let mid = (p0 + p1).scaled(R::half());
            let vloc = vlocs[j];
            let ell = p1 - p0;
            let f = vloc.cross(&ell).scaled(rho * gamma[j]);
            let r = mid - r0;
            total = total + r.cross(&f);
        }
        total
    }

    /// Decompose force/moment into (cl, cdi, cm, cl_roll, cn).
    #[allow(clippy::type_complexity)]
    fn coefficients<R: Real>(
        &self,
        force: &Vec3<R>,
        moment: &Vec3<R>,
        alpha_rad: &R,
        v_inf: R,
        rho: R,
    ) -> (R, R, R, R, R) {
        let (s, c) = alpha_rad.sin_cos();
        let lift_dir = Vec3::new(s.neg(), R::zero(), c);
        let drag_dir = Vec3::new(c, R::zero(), s);
        let lift = force.dot(&lift_dir);
        let drag = force.dot(&drag_dir);
        let q_s = R::half() * rho * v_inf * v_inf * R::from_f64(self.reference.area * self.beta);
        let q_sc = q_s * R::from_f64(self.reference.chord);
        let q_sb = q_s * R::from_f64(self.reference.span);
        (
            lift / q_s,
            drag / q_s,
            moment.y / q_sc,
            moment.x / q_sb,
            moment.z / q_sb,
        )
    }

    /// Section lift coefficients per spanwise strip (midpoint load model).
    fn strip_coefficients(
        &self,
        gamma: &[f64],
        vlocs: &[Vec3<f64>],
        v_inf: f64,
        density: f64,
        alpha_rad: f64,
    ) -> Vec<f64> {
        let elements = self.lattice.elements();
        let (s, c) = alpha_rad.sin_cos();
        let lift_dir = Vec3::new(-s, 0.0, c);
        let q = 0.5 * density * v_inf * v_inf * self.beta;
        let mut force_per_strip: Vec<f64> = vec![0.0; self.lattice.strip_count()];
        let mut ref_len_per_strip: Vec<f64> = vec![0.0; self.lattice.strip_count()];
        for (j, ej) in elements.iter().enumerate() {
            let vloc = vlocs[j];
            let (p0, p1) = ej.bounds::<f64>();
            let ell = p1 - p0;
            let f = vloc.cross(&ell).scaled(density * gamma[j]);
            force_per_strip[ej.strip] += f.dot(&lift_dir);
        }
        // Reference length per strip: q . chord . span (counted once).
        for e in elements.iter() {
            ref_len_per_strip[e.strip] = q * e.strip_chord * e.strip_dy;
        }
        force_per_strip
            .iter()
            .zip(ref_len_per_strip.iter())
            .map(|(f, r)| if r.abs() < 1.0e-12 { 0.0 } else { f / r })
            .collect()
    }
}

/// Whether a flight condition lies inside the regime the linear VLM is
/// validated for. Each flag is independent, so a result can fail one limit and
/// pass the others.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ValidityFlags {
    /// `|α| ≤ 15°`: attached flow, where the linear lift slope holds.
    pub alpha_ok: bool,
    /// `|β| ≤ 15°`: small sideslip, where the linearised side loads hold.
    pub beta_ok: bool,
    /// `M < 0.7`: the Göthert compressibility correction is used below this.
    pub mach_ok: bool,
    /// `Re ≥ 1e5`: below this, viscous and laminar effects dominate the
    /// inviscid lattice.
    pub reynolds_ok: bool,
}

impl ValidityFlags {
    /// True when every flag is inside its limit.
    #[must_use]
    pub const fn all_ok(&self) -> bool {
        self.alpha_ok && self.beta_ok && self.mach_ok && self.reynolds_ok
    }
}

/// Validity limits for [`validity_flags`], in radians, Mach and Reynolds number.
pub const VALIDITY_ALPHA_LIMIT: f64 = 15.0 * core::f64::consts::PI / 180.0;
/// See [`ValidityFlags::mach_ok`].
pub const VALIDITY_MACH_LIMIT: f64 = 0.7;
/// See [`ValidityFlags::reynolds_ok`].
pub const VALIDITY_REYNOLDS_MIN: f64 = 1.0e5;

/// Check a flight condition against the linear VLM's validity limits.
///
/// Non-finite inputs fail every flag they touch, so a NaN never reads as valid.
#[must_use]
pub fn validity_flags(alpha_rad: f64, beta_rad: f64, mach: f64, reynolds: f64) -> ValidityFlags {
    ValidityFlags {
        alpha_ok: alpha_rad.abs() <= VALIDITY_ALPHA_LIMIT,
        beta_ok: beta_rad.abs() <= VALIDITY_ALPHA_LIMIT,
        mach_ok: mach.is_finite() && mach < VALIDITY_MACH_LIMIT,
        reynolds_ok: reynolds >= VALIDITY_REYNOLDS_MIN,
    }
}
