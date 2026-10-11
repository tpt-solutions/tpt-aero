//! Steady turning and steady sideslip trim, solved nonlinearly on the VLM.
//!
//! A steady flight is defined by the bank angle `φ`, the turn rate `ψ̇`, and
//! the airspeed `V`. The pitch angle follows from zero flight-path rate,
//! `tan θ = (v sin φ + w cos φ) / u`. The body rates are then fixed by the
//! steady-turn kinematics (aircraft axes, x forward, y right, z down):
//!
//! * `p = −ψ̇ sin θ`, `q = ψ̇ sin φ cos θ`, `r = ψ̇ cos φ cos θ`
//!
//! Steady force balance `F = m (ω × V − g)` then gives the required aircraft
//! force, which is projected onto the wind-axis lift and side directions.
//! Thrust acts along the airspeed vector, so it drops out of those two
//! projections. The unknowns are `α`, `β` and the three control deflections.
//! The equations are these two force balances plus the pitch, roll and yaw
//! moments, all zero.
//!
//! A coordinated turn has `ψ̇ = g tan φ / V` and `β ≈ 0`. A steady sideslip
//! has `ψ̇ = 0`. The solve is Newton iteration with a finite-difference
//! Jacobian, continued from level flight, and it reuses
//! [`tpt_aero_vlm::PreparedVlm::solve_state`] with control-deflected twists.

use alloc::vec;
use alloc::vec::Vec;

use tpt_aero_core::linalg::lu_solve;
use tpt_aero_vlm::{FlightState, VlmModel};

use crate::trim::TrimError;

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// Standard gravity (m/s²).
/// Standard gravity (m/s²), re-exported from the crate root.
pub const G: f64 = crate::STANDARD_GRAVITY;

/// Turn rate `ψ̇` (rad/s) of a coordinated level turn at `bank_rad` and speed
/// `v_inf`: `g tan φ / V`. Positive for a right turn.
#[must_use]
pub fn coordinated_turn_rate(v_inf: f64, bank_rad: f64) -> f64 {
    G * bank_rad.tan() / v_inf
}

/// What to trim for a steady turn or sideslip.
#[derive(Debug, Clone, Copy)]
pub struct TurnRequest<'a> {
    /// Weight (N).
    pub weight_n: f64,
    /// True airspeed (m/s).
    pub v_inf: f64,
    /// Air density (kg/m³).
    pub density: f64,
    /// Bank angle (rad), positive right wing down.
    pub bank_rad: f64,
    /// Turn rate `ψ̇` (rad/s), positive for a right turn. Use
    /// [`coordinated_turn_rate`] for a coordinated turn, `0.0` for a steady
    /// sideslip.
    pub turn_rate: f64,
    /// Symmetric pitch control (e.g. `"elevator"`).
    pub elevator: &'a str,
    /// Antisymmetric roll control (e.g. `"aileron"`).
    pub aileron: &'a str,
    /// Yaw control (e.g. `"rudder"`).
    pub rudder: &'a str,
}

/// A converged steady turn or sideslip.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TurnTrim {
    /// Angle of attack (rad).
    pub alpha_rad: f64,
    /// Sideslip angle (rad), positive for wind from the right.
    pub beta_rad: f64,
    /// Elevator command (rad).
    pub elevator_rad: f64,
    /// Aileron command (rad).
    pub aileron_rad: f64,
    /// Rudder command (rad).
    pub rudder_rad: f64,
    /// Load factor, the wind-axis lift over weight.
    pub load_factor: f64,
    /// Newton iterations used, summed over the continuation steps.
    pub iterations: usize,
}

/// Number of unknowns and residuals: `[α, β, δe, δa, δr]`.
const N: usize = 5;

/// Number of continuation steps from level flight to the requested turn.
/// Each step warm-starts Newton from the previous root.
const CONTINUATION_STEPS: usize = 8;

/// Problem data shared by the residual and the Newton solve.
struct Setup {
    weight_n: f64,
    v_inf: f64,
    density: f64,
    q_s: f64,
    ie: usize,
    ia: usize,
    ir: usize,
}

/// Body rates `(p, q, r)` (rad/s, aircraft axes) of a steady turn.
fn steady_rates(bank: f64, pitch: f64, turn_rate: f64) -> [f64; 3] {
    let (sphi, cphi) = bank.sin_cos();
    let (sth, cth) = pitch.sin_cos();
    [
        -turn_rate * sth,
        turn_rate * sphi * cth,
        turn_rate * cphi * cth,
    ]
}

/// Stability-axis rates `(p_s, q_s, r_s)` that `solve_state` takes, from
/// aircraft-axis body rates at angle of attack `alpha`.
fn stability_rates(alpha: f64, body: [f64; 3]) -> [f64; 3] {
    // Model axes: ω_m = (−p, q, −r); stability axes project onto forward/down.
    let (sa, ca) = alpha.sin_cos();
    let [p, q, r] = body;
    [ca * p + sa * r, q, -sa * p + ca * r]
}

/// Residuals `[lift, side, pitch, roll, yaw]` at unknowns `x = [α, β, δe, δa, δr]`
/// for bank `bank` (rad) and turn rate `turn_rate` (rad/s).
fn residual(
    model: &VlmModel,
    s: &Setup,
    x: &[f64; N],
    bank: f64,
    turn_rate: f64,
) -> Result<[f64; N], TrimError> {
    let (alpha, beta) = (x[0], x[1]);
    let mut deflections = vec![0.0; model.controls().len()];
    deflections[s.ie] = x[2];
    deflections[s.ia] = x[3];
    deflections[s.ir] = x[4];
    let twists = model.control_twists(&deflections)?;
    let prepared = model.prepare_with_twists(&twists)?;

    // Steady kinematics and the required aircraft-axis force.
    let (sa, ca) = alpha.sin_cos();
    let (sb, cb) = beta.sin_cos();
    let (u, v, w) = (s.v_inf * ca * cb, s.v_inf * sb, s.v_inf * sa * cb);
    let (sphi, cphi) = bank.sin_cos();
    // atan2 form of atan(·/u), valid for u > 0.
    let pitch = (v * sphi + w * cphi).atan2(u);
    let (sth, cth) = pitch.sin_cos();
    let [p, q, r] = steady_rates(bank, pitch, turn_rate);
    let gravity = [-G * sth, G * cth * sphi, G * cth * cphi];
    let mass = s.weight_n / G;
    // F = m (ω × V − g), aircraft axes.
    let f_aircraft = [
        mass * (q * w - r * v - gravity[0]),
        mass * (r * u - p * w - gravity[1]),
        mass * (p * v - q * u - gravity[2]),
    ];
    // Aircraft → model axes: R = diag(−1, 1, −1).
    let f = [-f_aircraft[0], f_aircraft[1], -f_aircraft[2]];

    let [p_s, q_s, r_s] = stability_rates(alpha, [p, q, r]);
    let loads = prepared.solve_state(
        &FlightState {
            alpha,
            beta,
            p: p_s,
            q: q_s,
            r: r_s,
        },
        s.v_inf,
        s.density,
    )?;

    // Wind-axis lift and side directions, as in `PreparedVlm::solve_state`.
    let drag_dir = [ca * cb, -sb, sa * cb];
    let lift_raw = [-drag_dir[2], 0.0, drag_dir[0]];
    let lift_norm = (lift_raw[0] * lift_raw[0] + lift_raw[2] * lift_raw[2]).sqrt();
    let lift_dir = [lift_raw[0] / lift_norm, 0.0, lift_raw[2] / lift_norm];
    let side_raw = [
        -lift_dir[2] * drag_dir[1],
        lift_dir[2] * drag_dir[0] - lift_dir[0] * drag_dir[2],
        lift_dir[0] * drag_dir[1],
    ];
    let side_norm =
        (side_raw[0] * side_raw[0] + side_raw[1] * side_raw[1] + side_raw[2] * side_raw[2]).sqrt();
    let side_dir = [
        side_raw[0] / side_norm,
        side_raw[1] / side_norm,
        side_raw[2] / side_norm,
    ];
    let dot = |a: &[f64; 3], b: &[f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let req_lift = dot(&f, &lift_dir) / s.q_s;
    let req_side = dot(&f, &side_dir) / s.q_s;
    Ok([
        loads.c_lift - req_lift,
        loads.c_side - req_side,
        loads.c_pitch,
        loads.c_roll,
        loads.c_yaw,
    ])
}

/// Newton iteration at fixed bank and turn rate, from `x0`. Returns the root
/// and the iteration count.
fn newton(
    model: &VlmModel,
    s: &Setup,
    x0: [f64; N],
    bank: f64,
    turn_rate: f64,
) -> Result<([f64; N], usize), TrimError> {
    let mut x = x0;
    let h = 1.0e-6;
    for iteration in 0..=25 {
        let f = residual(model, s, &x, bank, turn_rate)?;
        if f.iter().all(|fi| fi.abs() < 1.0e-9) {
            return Ok((x, iteration));
        }
        if iteration == 25 {
            break;
        }
        // Forward-difference Jacobian, columns = ∂f/∂x_k.
        let mut jac: Vec<Vec<f64>> = vec![vec![0.0; N]; N];
        for k in 0..N {
            let mut xp = x;
            xp[k] += h;
            let fp = residual(model, s, &xp, bank, turn_rate)?;
            for i in 0..N {
                jac[i][k] = (fp[i] - f[i]) / h;
            }
        }
        let neg_f: Vec<f64> = f.iter().map(|fi| -fi).collect();
        let step = lu_solve(&jac, &neg_f).map_err(|_| TrimError::Singular)?;
        // Damped step: shrink it until the residual norm falls, so an overshoot
        // does not jump to a different root or leave the valid envelope.
        let norm = |r: &[f64; N]| r.iter().map(|v| v * v).sum::<f64>().sqrt();
        let current = norm(&f);
        let mut lambda = 1.0_f64;
        let mut accepted = x;
        for _ in 0..20 {
            let mut trial = x;
            for k in 0..N {
                trial[k] += lambda * step[k];
            }
            if let Ok(rt) = residual(model, s, &trial, bank, turn_rate) {
                if norm(&rt) < current || lambda < 1.0e-3 {
                    accepted = trial;
                    break;
                }
            }
            lambda *= 0.5;
            accepted = trial;
        }
        x = accepted;
        if !x.iter().all(|v| v.is_finite()) || x[0].abs() > 1.0 || x[1].abs() > 1.0 {
            return Err(TrimError::NoConvergence);
        }
    }
    Err(TrimError::NoConvergence)
}

/// Trim `model` for a steady turn or sideslip described by `req`.
///
/// The solve continues from level flight: bank and turn rate are ramped to
/// their requested values in a fixed number of steps, each warm-started
/// from the previous root.
///
/// # Errors
/// [`TrimError::InvalidRequest`] for non-finite or non-positive inputs or a
/// bank angle at or beyond ±90°; [`TrimError::NoSuchControl`] if any of the
/// three named controls is missing; [`TrimError::Singular`] if the Jacobian
/// is singular (e.g. a control with no authority); [`TrimError::NoConvergence`]
/// otherwise; [`TrimError::Solve`] for VLM errors.
pub fn trim_steady_turn(model: &VlmModel, req: &TurnRequest<'_>) -> Result<TurnTrim, TrimError> {
    let positive = |x: f64| x.is_finite() && x > 0.0;
    if !(positive(req.weight_n)
        && positive(req.v_inf)
        && positive(req.density)
        && req.bank_rad.is_finite()
        && req.bank_rad.abs() < core::f64::consts::FRAC_PI_2
        && req.turn_rate.is_finite())
    {
        return Err(TrimError::InvalidRequest);
    }
    let setup = Setup {
        weight_n: req.weight_n,
        v_inf: req.v_inf,
        density: req.density,
        q_s: 0.5 * req.density * req.v_inf * req.v_inf * model.reference().area,
        ie: model
            .control_index(req.elevator)
            .ok_or(TrimError::NoSuchControl)?,
        ia: model
            .control_index(req.aileron)
            .ok_or(TrimError::NoSuchControl)?,
        ir: model
            .control_index(req.rudder)
            .ok_or(TrimError::NoSuchControl)?,
    };

    let mut x = [0.05, 0.0, 0.0, 0.0, 0.0];
    let mut iterations = 0;
    // Adaptive continuation: ramp the bank and turn rate from level flight. A
    // ramp step that fails is halved and retried, down to a floor, so a hard
    // case gets finer steps instead of a fixed schedule that gives up.
    let max_step = 1.0 / CONTINUATION_STEPS as f64;
    let min_step = max_step / 128.0;
    let mut done = 0.0_f64;
    let mut increment = max_step;
    while done < 1.0 - 1.0e-12 {
        let target = (done + increment).min(1.0);
        match newton(
            model,
            &setup,
            x,
            req.bank_rad * target,
            req.turn_rate * target,
        ) {
            Ok((root, its)) => {
                x = root;
                iterations += its;
                done = target;
                increment = (increment * 1.5).min(max_step);
            }
            Err(TrimError::Singular) => return Err(TrimError::Singular),
            Err(e) => {
                if increment <= min_step {
                    return Err(e);
                }
                increment *= 0.5;
            }
        }
    }

    // Load factor: wind-axis lift over weight, from the converged state.
    let (sa, ca) = x[0].sin_cos();
    let (sb, cb) = x[1].sin_cos();
    let (u, v, w) = (req.v_inf * ca * cb, req.v_inf * sb, req.v_inf * sa * cb);
    let (sphi, cphi) = req.bank_rad.sin_cos();
    // atan2 form of atan(·/u), valid for u > 0.
    let pitch = (v * sphi + w * cphi).atan2(u);
    let (sth, cth) = pitch.sin_cos();
    let [p, q, r] = steady_rates(req.bank_rad, pitch, req.turn_rate);
    let gravity = [-G * sth, G * cth * sphi, G * cth * cphi];
    let mass = req.weight_n / G;
    // Wind-axis lift of the required aircraft-axis force: in model axes the
    // lift direction is (−sin α, 0, cos α) (after normalisation, β drops out),
    // so L = F_xa sin α − F_za cos α.
    let f_x = mass * (q * w - r * v - gravity[0]);
    let f_z = mass * (p * v - q * u - gravity[2]);
    let lift = f_x * sa - f_z * ca;
    let load_factor = lift / req.weight_n;

    Ok(TurnTrim {
        alpha_rad: x[0],
        beta_rad: x[1],
        elevator_rad: x[2],
        aileron_rad: x[3],
        rudder_rad: x[4],
        load_factor,
        iterations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coordinated_turn_rate_matches_textbook() {
        // 20° bank at 60 m/s: ψ̇ = g tan φ / V = 9.80665 × 0.36397 / 60 ≈ 0.05949 rad/s.
        let r = coordinated_turn_rate(60.0, 20.0_f64.to_radians());
        assert!((r - 0.059_49).abs() < 1e-4, "ψ̇ = {r}");
    }

    #[test]
    fn steady_rates_reduce_to_yaw_in_level_flight() {
        let [p, q, r] = steady_rates(0.0, 0.0, 0.05);
        assert!(p.abs() < 1e-15 && q.abs() < 1e-15 && (r - 0.05).abs() < 1e-15);
    }
}
