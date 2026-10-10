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

use crate::lattice::Lattice;
use crate::linalg::{lu_solve, lu_solve_generic, SolveError};
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
        }
    }
}

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
    /// Section lift coefficients per spanwise strip (root-side first).
    pub strip_cl: Vec<f64>,
    /// Total dimensional force (N) in body axes.
    pub force: Vec3<f64>,
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

/// A horseshoe-vortex VLM model: lattice + reference geometry.
#[derive(Debug, Clone)]
pub struct VlmModel {
    lattice: Lattice,
    reference: ReferenceGeometry,
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
        let lattice =
            Lattice::build(planform, n_span, n_chord, camber).map_err(VlmError::Lattice)?;
        Ok(Self {
            lattice,
            reference: ReferenceGeometry::from_planform(planform),
        })
    }

    /// Build from a pre-made lattice; reference geometry is estimated from
    /// the lattice panels (use [`VlmModel::with_reference`] to override).
    #[must_use]
    pub fn from_lattice(lattice: Lattice) -> Self {
        let area: f64 = lattice.elements().iter().map(|e| e.area).sum();
        let span: f64 = lattice
            .elements()
            .iter()
            .filter(|e| e.strip == 0)
            .map(|e| e.strip_dy)
            .sum();
        let chord = if area.abs() < 1.0e-12 {
            1.0
        } else {
            area / span.max(1.0e-12)
        };
        Self {
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

    /// Override the reference geometry (e.g. to match a wind-tunnel model).
    #[must_use]
    pub fn with_reference(mut self, reference: ReferenceGeometry) -> Self {
        self.reference = reference;
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
        let elements = self.lattice.elements();
        let n = elements.len();
        if !twists.is_empty() && twists.len() != n {
            return Err(VlmError::DimensionMismatch);
        }
        let zero = R::zero();
        let twist_at = |i: usize| twists.get(i).copied().unwrap_or(zero);

        // Freestream (body frame). The horseshoe trailing legs convect along
        // the fixed body +x axis: this keeps the AIC matrix independent of
        // the flow state, so the model is exactly linear in V∞ (freestream
        // reversal ⇒ exact force reversal) and avoids upstream-leg
        // singularities at large |α|. For the small incidences VLM is valid
        // for, this is indistinguishable from a freestream-aligned wake.
        let (s, c) = alpha_rad.sin_cos();
        let v = Vec3::new(v_inf * c, zero, v_inf * s);
        let d = Vec3::new(R::one(), zero, zero);

        // AIC: a_ij = n_i · v_j(cp_i) [1/m]; rhs_i = −V∞·n_i [m/s].
        let mut a: Vec<Vec<R>> = vec![vec![zero; n]; n];
        let mut rhs = vec![zero; n];
        for (i, ei) in elements.iter().enumerate() {
            let ni = ei.normal_for_twist(twist_at(i));
            let cp = ei.cp::<R>();
            rhs[i] = v.dot(&ni).neg();
            for (j, ej) in elements.iter().enumerate() {
                let (ba, bb) = ej.bounds::<R>();
                let vel = horseshoe_velocity(&ba, &bb, &d, &cp, R::one());
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

        let force = self.total_force::<R>(&gamma, &v, &d, rho);
        let moment = self.total_moment::<R>(&gamma, &v, &d, rho);
        let (cl, cdi, cm, cl_roll, cn) =
            self.coefficients::<R>(&force, &moment, &alpha_rad, v_inf, rho);
        Ok(VlmSolutionFields {
            cl,
            cdi,
            cm,
            cl_roll,
            cn,
            gamma,
        })
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
        let n = self.lattice.len();
        if !twists.is_empty() && twists.len() != n {
            return Err(VlmError::DimensionMismatch);
        }
        let fields = self.solve_generic::<f64>(alpha_rad, v_inf, density, twists)?;
        let (s, c) = alpha_rad.sin_cos();
        let v = Vec3::new(v_inf * c, 0.0, v_inf * s);
        let d = Vec3::new(1.0, 0.0, 0.0);
        let force = self.total_force::<f64>(&fields.gamma, &v, &d, density);
        let strip_cl = self.strip_coefficients(&fields.gamma, &v, &d, density, alpha_rad);
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
        let (ba, bb) = elements[j].bounds::<R>();
        let mid = (ba + bb).scaled(R::half());
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
        }
        vloc
    }

    /// Total dimensional body-frame force (N): F = Σ ρ Γ (V_loc × (b − a)).
    fn total_force<R: Real>(&self, gamma: &[R], v: &Vec3<R>, d: &Vec3<R>, rho: R) -> Vec3<R> {
        let elements = self.lattice.elements();
        let mut total = Vec3::<R>::zero();
        for (j, ej) in elements.iter().enumerate() {
            let (ba, bb) = ej.bounds::<R>();
            let vloc = self.local_velocity::<R>(j, gamma, v, d);
            // b − a points +y for a lift-up Γ.
            let ell = bb - ba;
            total = total + vloc.cross(&ell).scaled(rho * gamma[j]);
        }
        total
    }

    /// Moment vector (N·m) about the reference point; each panel force acts
    /// at its bound-vortex midpoint.
    fn total_moment<R: Real>(&self, gamma: &[R], v: &Vec3<R>, d: &Vec3<R>, rho: R) -> Vec3<R> {
        let elements = self.lattice.elements();
        let r0 = Vec3::new(
            R::from_f64(self.reference.moment_center.x),
            R::from_f64(self.reference.moment_center.y),
            R::from_f64(self.reference.moment_center.z),
        );
        let mut total = Vec3::<R>::zero();
        for (j, ej) in elements.iter().enumerate() {
            let (ba, bb) = ej.bounds::<R>();
            let mid = (ba + bb).scaled(R::half());
            let vloc = self.local_velocity::<R>(j, gamma, v, d);
            let ell = bb - ba;
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
        let q_s = R::half() * rho * v_inf * v_inf * R::from_f64(self.reference.area);
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
        v: &Vec3<f64>,
        d: &Vec3<f64>,
        density: f64,
        alpha_rad: f64,
    ) -> Vec<f64> {
        let elements = self.lattice.elements();
        let (s, c) = alpha_rad.sin_cos();
        let lift_dir = Vec3::new(-s, 0.0, c);
        let q = 0.5 * density * v.norm().powi(2);
        let mut force_per_strip: Vec<f64> = vec![0.0; self.lattice.strip_count()];
        let mut ref_len_per_strip: Vec<f64> = vec![0.0; self.lattice.strip_count()];
        for (j, ej) in elements.iter().enumerate() {
            let vloc = self.local_velocity::<f64>(j, gamma, v, d);
            let (ba, bb) = ej.bounds::<f64>();
            let ell = bb - ba;
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
