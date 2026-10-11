//! Level-flight longitudinal trim: angle of attack and elevator.
//!
//! Solves `C_L(α, δe) = W/(q S)` and `C_m(α, δe) = 0` about the model's
//! reference point (the centre of gravity), by Newton iteration with a
//! finite-difference Jacobian on the VLM. Thrust is assumed to act through
//! the CG, so it only has to balance drag.

use alloc::vec;

use tpt_aero_vlm::{VlmError, VlmModel};

/// What to trim.
#[derive(Debug, Clone, Copy)]
pub struct TrimRequest<'a> {
    /// Weight (N) = required lift in level flight.
    pub weight_n: f64,
    /// True airspeed (m/s).
    pub v_inf: f64,
    /// Air density (kg/m³).
    pub density: f64,
    /// Name of the symmetric pitch control to trim with (e.g. `"elevator"`).
    pub elevator: &'a str,
    /// Zero-lift (profile and parasite) drag coefficient to add to the
    /// induced drag when reporting the thrust required (0 if unknown).
    pub cd0: f64,
}

/// A converged trim point.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Trim {
    /// Angle of attack (rad).
    pub alpha_rad: f64,
    /// Elevator command (rad; positive = trailing edge down, see
    /// [`tpt_aero_vlm::Control`] for the sign convention).
    pub elevator_rad: f64,
    /// Trimmed lift coefficient (= `W/(qS)`).
    pub cl: f64,
    /// Induced drag coefficient from the VLM.
    pub cdi: f64,
    /// Thrust required for level flight (N): `(C_Di + cd0)·q·S`.
    pub thrust_n: f64,
    /// Newton iterations used.
    pub iterations: usize,
}

/// Why a trim failed.
#[derive(Debug, Clone, PartialEq)]
pub enum TrimError {
    /// Non-finite or non-positive weight, speed, density or `cd0 < 0`.
    InvalidRequest,
    /// The model has no control with the requested name.
    NoSuchControl,
    /// The VLM solve failed.
    Solve(VlmError),
    /// The Jacobian is singular: the controls cannot change `C_L` and `C_m`
    /// independently (a tail-less wing, or a control with no authority).
    Singular,
    /// Newton iteration did not converge.
    NoConvergence,
    /// The trimmed angle of attack is beyond ±20°, outside the range where the
    /// linear, attached-flow VLM is meaningful.
    OutsideEnvelope,
}

/// Largest angle of attack (rad) a trim may return: 20°.
pub const TRIM_ALPHA_LIMIT: f64 = 20.0 * core::f64::consts::PI / 180.0;

impl core::fmt::Display for TrimError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidRequest => write!(f, "trim request must be finite and positive"),
            Self::NoSuchControl => write!(f, "model has no control with that name"),
            Self::Solve(e) => write!(f, "VLM solve failed: {e}"),
            Self::Singular => write!(f, "angle of attack and elevator cannot trim independently"),
            Self::NoConvergence => write!(f, "trim iteration did not converge"),
            Self::OutsideEnvelope => write!(f, "trimmed angle of attack is beyond ±20 deg"),
        }
    }
}

impl core::error::Error for TrimError {}

impl From<VlmError> for TrimError {
    fn from(e: VlmError) -> Self {
        Self::Solve(e)
    }
}

/// Trim `model` for steady level flight.
///
/// # Errors
/// See [`TrimError`].
pub fn trim_level_flight(model: &VlmModel, req: &TrimRequest<'_>) -> Result<Trim, TrimError> {
    let ok = |x: f64| x.is_finite() && x > 0.0;
    if !(ok(req.weight_n)
        && ok(req.v_inf)
        && ok(req.density)
        && req.cd0.is_finite()
        && req.cd0 >= 0.0)
    {
        return Err(TrimError::InvalidRequest);
    }
    let index = model
        .control_index(req.elevator)
        .ok_or(TrimError::NoSuchControl)?;
    let q_s = 0.5 * req.density * req.v_inf * req.v_inf * model.reference().area;
    let cl_required = req.weight_n / q_s;

    // Residuals (CL - CL_required, Cm) at (alpha, elevator), and the induced drag.
    let evaluate = |alpha: f64, elevator: f64| -> Result<(f64, f64, f64), TrimError> {
        let mut deflections = vec![0.0; model.controls().len()];
        deflections[index] = elevator;
        let twists = model.control_twists(&deflections)?;
        let s = model.solve(alpha, req.v_inf, req.density, &twists)?;
        Ok((s.cl - cl_required, s.cm, s.cdi))
    };

    let (mut alpha, mut elevator) = (0.03_f64, 0.0_f64);
    let h = 1.0e-4;
    // Relative tolerance on the lift residual; the moment is dimensionless.
    let tol_cl = 1.0e-10 * cl_required.abs().max(1.0);
    for iteration in 1..=30 {
        let (f_cl, f_cm, cdi) = evaluate(alpha, elevator)?;
        if f_cl.abs() < tol_cl && f_cm.abs() < 1.0e-10 {
            if alpha.abs() > TRIM_ALPHA_LIMIT {
                return Err(TrimError::OutsideEnvelope);
            }
            return Ok(Trim {
                alpha_rad: alpha,
                elevator_rad: elevator,
                cl: cl_required + f_cl,
                cdi,
                thrust_n: (cdi + req.cd0) * q_s,
                iterations: iteration - 1,
            });
        }
        let (a_p, a_n) = (
            evaluate(alpha + h, elevator)?,
            evaluate(alpha - h, elevator)?,
        );
        let (e_p, e_n) = (
            evaluate(alpha, elevator + h)?,
            evaluate(alpha, elevator - h)?,
        );
        let j = [
            [(a_p.0 - a_n.0) / (2.0 * h), (e_p.0 - e_n.0) / (2.0 * h)],
            [(a_p.1 - a_n.1) / (2.0 * h), (e_p.1 - e_n.1) / (2.0 * h)],
        ];
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        let scale = j[0][0].abs() * j[1][1].abs() + j[0][1].abs() * j[1][0].abs();
        if det.is_nan() || det.abs() <= 1.0e-9 * scale.max(1.0e-30) {
            return Err(TrimError::Singular);
        }
        let d_alpha = (j[1][1] * f_cl - j[0][1] * f_cm) / det;
        let d_elev = (-j[1][0] * f_cl + j[0][0] * f_cm) / det;
        // Damped step: halve until the residual norm decreases, so a Newton
        // step that overshoots does not throw the iterate away.
        // Squared norm: only compared, never reported, and `hypot` is not in `no_std`.
        let norm = |r: (f64, f64)| {
            let a = r.0 / cl_required.abs().max(1.0);
            a * a + r.1 * r.1
        };
        let current = norm((f_cl, f_cm));
        let mut lambda = 1.0_f64;
        let mut accepted = (alpha - lambda * d_alpha, elevator - lambda * d_elev);
        for _ in 0..20 {
            let (a, e) = accepted;
            if let Ok((c, m, _)) = evaluate(a, e) {
                if norm((c, m)) < current || lambda < 1.0e-3 {
                    break;
                }
            }
            lambda *= 0.5;
            accepted = (alpha - lambda * d_alpha, elevator - lambda * d_elev);
        }
        alpha = accepted.0;
        elevator = accepted.1;
        if !(alpha.is_finite() && elevator.is_finite() && alpha.abs() < 1.0) {
            return Err(TrimError::NoConvergence);
        }
    }
    Err(TrimError::NoConvergence)
}
