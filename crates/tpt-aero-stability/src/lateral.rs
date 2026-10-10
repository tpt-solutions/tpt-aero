//! Lateral-directional stability derivatives from the VLM.
//!
//! Extracted by central finite differences of [`PreparedVlm::solve_state`]
//! (sideslip and body rates enter the lattice as kinematic velocity changes,
//! see the VLM crate), so the whole set costs one AIC factorization. All
//! signs follow aircraft convention: `β` positive for wind from the right,
//! `p` right-wing-down, `r` nose-right, side force positive right, roll
//! positive right-wing-down, yaw positive nose-right.
//!
//! Rate derivatives are per non-dimensional rate: `p̂ = p·b/(2V)`,
//! `r̂ = r·b/(2V)`; `β` derivatives are per radian.

use alloc::vec::Vec;

use tpt_aero_vlm::{FlightState, PreparedVlm, StateLoads, VlmError, VlmModel};

use crate::controls::ControlDerivatives;
use crate::modes::{eigenvalues4, Complex};
use crate::sim::StateSpace;

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// Lateral-directional derivative set (per rad, or per unit non-dimensional rate).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LateralDerivatives {
    /// Side force due to sideslip, `C_{Yβ}`.
    pub cy_beta: f64,
    /// Rolling moment due to sideslip (dihedral effect), `C_{lβ}`.
    pub cl_beta: f64,
    /// Yawing moment due to sideslip (weathercock stability), `C_{nβ}`.
    pub cn_beta: f64,
    /// Side force due to roll rate, `C_{Yp}`.
    pub cy_p: f64,
    /// Roll damping, `C_{lp}`.
    pub cl_p: f64,
    /// Yawing moment due to roll rate, `C_{np}`.
    pub cn_p: f64,
    /// Side force due to yaw rate, `C_{Yr}`.
    pub cy_r: f64,
    /// Rolling moment due to yaw rate, `C_{lr}`.
    pub cl_r: f64,
    /// Yaw damping, `C_{nr}`.
    pub cn_r: f64,
}

/// Step for angles (rad) and normalized rates.
const STEP: f64 = 1.0e-3;

fn central(
    prepared: &PreparedVlm<'_>,
    state_at: impl Fn(f64) -> FlightState,
    v_inf: f64,
    density: f64,
) -> Result<[f64; 3], VlmError> {
    let plus: StateLoads = prepared.solve_state(&state_at(STEP), v_inf, density)?;
    let minus: StateLoads = prepared.solve_state(&state_at(-STEP), v_inf, density)?;
    let d = |a: f64, b: f64| (a - b) / (2.0 * STEP);
    Ok([
        d(plus.c_side, minus.c_side),
        d(plus.c_roll, minus.c_roll),
        d(plus.c_yaw, minus.c_yaw),
    ])
}

/// Extract the lateral-directional derivatives of `model` at trim angle of
/// attack `alpha_rad`, speed `v_inf` (m/s) and density (kg/m³).
///
/// The moment reference point is the model's own; wing-alone results have no
/// vertical-tail or fuselage contribution (`C_{nβ}` in particular is only the
/// small wing part).
///
/// # Errors
/// Any [`VlmError`] from preparing or solving the model, including
/// [`VlmError::InvalidFlightCondition`].
pub fn lateral_fd(
    model: &VlmModel,
    alpha_rad: f64,
    v_inf: f64,
    density: f64,
) -> Result<LateralDerivatives, VlmError> {
    let prepared = model.prepare()?;
    let span = model.reference().span;
    // p = p̂·2V/b, r = r̂·2V/b.
    let rate = 2.0 * v_inf / span;
    let base = FlightState {
        alpha: alpha_rad,
        ..FlightState::default()
    };
    let [cy_beta, cl_beta, cn_beta] = central(
        &prepared,
        |b| FlightState { beta: b, ..base },
        v_inf,
        density,
    )?;
    let [cy_p, cl_p, cn_p] = central(
        &prepared,
        |ph| FlightState {
            p: ph * rate,
            ..base
        },
        v_inf,
        density,
    )?;
    let [cy_r, cl_r, cn_r] = central(
        &prepared,
        |rh| FlightState {
            r: rh * rate,
            ..base
        },
        v_inf,
        density,
    )?;
    Ok(LateralDerivatives {
        cy_beta,
        cl_beta,
        cn_beta,
        cy_p,
        cl_p,
        cn_p,
        cy_r,
        cl_r,
        cn_r,
    })
}

/// Mass properties and flight condition for dimensionalizing the lateral plant.
#[derive(Debug, Clone, Copy)]
pub struct LateralParameters {
    /// True airspeed (m/s).
    pub v_inf: f64,
    /// Reference span `b` (m).
    pub span: f64,
    /// Reference area `S` (m^2).
    pub area: f64,
    /// Air density (kg/m^3).
    pub density: f64,
    /// Mass (kg).
    pub mass: f64,
    /// Roll inertia `I_x` (kg m^2), stability axes.
    pub inertia_x: f64,
    /// Yaw inertia `I_z` (kg m^2), stability axes.
    pub inertia_z: f64,
    /// Product of inertia `I_xz` (kg m^2), stability axes (0 if unknown).
    pub inertia_xz: f64,
}

/// Linearized lateral-directional plant `x' = A x` for level flight with
/// state `[beta, p, r, phi]` (rad, rad/s, rad/s, rad), in stability axes:
///
/// ```text
/// beta' = (qS/mV)(CY_b beta + CY_p p b/2V + CY_r r b/2V) - r + (g/V) phi
/// [Ix  -Ixz][p']   [qSb (Cl_b beta + Cl_p p^ + Cl_r r^)]
/// [-Ixz  Iz][r'] = [qSb (Cn_b beta + Cn_p p^ + Cn_r r^)]
/// phi' = p
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LateralPlant {
    /// Dynamics matrix (4x4, row-major).
    pub a: [[f64; 4]; 4],
}

/// Standard gravity (m/s^2).
const GRAVITY: f64 = 9.80665;

impl LateralPlant {
    /// Assemble the plant from derivatives and mass properties. Returns
    /// `None` for non-positive speed, span, area, density, mass or inertias,
    /// or an inertia matrix that is not positive definite.
    #[must_use]
    pub fn from_derivatives(d: &LateralDerivatives, p: &LateralParameters) -> Option<Self> {
        let positive = |x: f64| x.is_finite() && x > 0.0;
        let det = p.inertia_x * p.inertia_z - p.inertia_xz * p.inertia_xz;
        if !(positive(p.v_inf)
            && positive(p.span)
            && positive(p.area)
            && positive(p.density)
            && positive(p.mass)
            && positive(p.inertia_x)
            && positive(p.inertia_z)
            && p.inertia_xz.is_finite()
            && det > 0.0)
        {
            return None;
        }
        let q_s = 0.5 * p.density * p.v_inf * p.v_inf * p.area;
        let rate = p.span / (2.0 * p.v_inf);
        let side = q_s / (p.mass * p.v_inf);
        let moment = q_s * p.span;

        // Moments per state: rows (L, N), columns (beta, p, r).
        let l = [
            d.cl_beta * moment,
            d.cl_p * moment * rate,
            d.cl_r * moment * rate,
        ];
        let n = [
            d.cn_beta * moment,
            d.cn_p * moment * rate,
            d.cn_r * moment * rate,
        ];
        // [p'; r'] = M^-1 [L; N], M = [[Ix, -Ixz], [-Ixz, Iz]].
        let (ix, iz, ixz) = (p.inertia_x, p.inertia_z, p.inertia_xz);
        let inv = |lk: f64, nk: f64| [(iz * lk + ixz * nk) / det, (ixz * lk + ix * nk) / det];
        let (col_b, col_p, col_r) = (inv(l[0], n[0]), inv(l[1], n[1]), inv(l[2], n[2]));
        let a = [
            [
                side * d.cy_beta,
                side * d.cy_p * rate,
                side * d.cy_r * rate - 1.0,
                GRAVITY / p.v_inf,
            ],
            [col_b[0], col_p[0], col_r[0], 0.0],
            [col_b[1], col_p[1], col_r[1], 0.0],
            [0.0, 1.0, 0.0, 0.0],
        ];
        Some(Self { a })
    }

    /// Eigenvalues of the plant.
    #[must_use]
    pub fn eigenvalues(&self) -> [Complex; 4] {
        eigenvalues4(&self.a)
    }

    /// The plant as a [`StateSpace`] with inputs `[aileron, rudder, v_gust]`:
    /// the two control commands (rad, VLM sign convention: positive aileron
    /// rolls left, positive rudder yaws right) with their effectiveness
    /// `aileron` and `rudder`, and the sidewash `v_gust` (m/s, positive for
    /// air arriving from the right), which adds `v_gust / V` to the
    /// aerodynamic sideslip.
    #[must_use]
    pub fn state_space(
        &self,
        p: &LateralParameters,
        aileron: &ControlDerivatives,
        rudder: &ControlDerivatives,
    ) -> StateSpace<4, 3> {
        let q_s = 0.5 * p.density * p.v_inf * p.v_inf * p.area;
        let side = q_s / (p.mass * p.v_inf);
        let moment = q_s * p.span;
        let det = p.inertia_x * p.inertia_z - p.inertia_xz * p.inertia_xz;
        // [p'; r'] = M^-1 [L; N].
        let apply = |l: f64, n: f64| {
            [
                (p.inertia_z * l + p.inertia_xz * n) / det,
                (p.inertia_xz * l + p.inertia_x * n) / det,
            ]
        };
        let ail = apply(moment * aileron.c_roll, moment * aileron.c_yaw);
        let rud = apply(moment * rudder.c_roll, moment * rudder.c_yaw);
        let mut b = [[0.0; 3]; 4];
        b[0][0] = side * aileron.cy;
        b[0][1] = side * rudder.cy;
        b[1][0] = ail[0];
        b[1][1] = rud[0];
        b[2][0] = ail[1];
        b[2][1] = rud[1];
        // The beta column of the dynamics is purely aerodynamic.
        for i in 0..3 {
            b[i][2] = self.a[i][0] / p.v_inf;
        }
        StateSpace { a: self.a, b }
    }

    /// Classify the eigenvalues into the three classical lateral modes.
    ///
    /// Returns `None` when the spectrum does not have the usual shape (one
    /// complex pair plus two real roots), e.g. when the Dutch roll has
    /// degenerated into two real roots.
    #[must_use]
    pub fn modes(&self) -> Option<LateralModes> {
        let eig = self.eigenvalues();
        let pair = eig.iter().find(|z| z.im > 0.0)?;
        let mut reals: [f64; 2] = {
            let r: Vec<f64> = eig.iter().filter(|z| z.im == 0.0).map(|z| z.re).collect();
            if r.len() != 2 {
                return None;
            }
            [r[0], r[1]]
        };
        // Roll subsidence is the fast real root, spiral the slow one.
        reals.sort_by(|x, y| {
            x.abs()
                .partial_cmp(&y.abs())
                .unwrap_or(core::cmp::Ordering::Equal)
        });
        let wn = pair.abs();
        Some(LateralModes {
            dutch_roll_wn: wn,
            dutch_roll_zeta: -pair.re / wn,
            roll_root: reals[1],
            spiral_root: reals[0],
        })
    }
}

/// The classical lateral-directional modes of a [`LateralPlant`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LateralModes {
    /// Dutch-roll natural frequency (rad/s).
    pub dutch_roll_wn: f64,
    /// Dutch-roll damping ratio (negative = divergent).
    pub dutch_roll_zeta: f64,
    /// Roll-subsidence eigenvalue (1/s; negative = stable, time constant
    /// `-1/root`).
    pub roll_root: f64,
    /// Spiral eigenvalue (1/s; negative = stable, positive = divergent with
    /// time to double `ln 2 / root`).
    pub spiral_root: f64,
}

impl LateralModes {
    /// Roll-subsidence time constant (s), `None` if the mode is unstable.
    #[must_use]
    pub fn roll_time_constant(&self) -> Option<f64> {
        (self.roll_root < 0.0).then(|| -1.0 / self.roll_root)
    }

    /// Time for the spiral mode to double in amplitude (s), `None` if stable.
    #[must_use]
    pub fn spiral_time_to_double(&self) -> Option<f64> {
        (self.spiral_root > 0.0).then(|| core::f64::consts::LN_2 / self.spiral_root)
    }
}
