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

use tpt_aero_vlm::{FlightState, PreparedVlm, StateLoads, VlmError, VlmModel};

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
