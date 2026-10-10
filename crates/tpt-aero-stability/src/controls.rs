//! Control-effectiveness derivatives for any named control surface.

use alloc::vec;

use tpt_aero_vlm::{VlmError, VlmModel};

/// Sensitivity of the force and moment coefficients to a control command
/// (per radian), about the stability axes, in aircraft-standard signs.
///
/// The command convention is the VLM's (positive = trailing edge down on the
/// `+y` half): a positive aileron command rolls left and a positive rudder
/// command yaws right. The aircraft-standard `C_{lδa}` and `C_{nδr}` are the
/// negatives of the fields below for those two controls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ControlDerivatives {
    /// Lift, `∂C_L/∂δ`.
    pub cl: f64,
    /// Induced drag, `∂C_{Di}/∂δ`.
    pub cd: f64,
    /// Side force, `∂C_Y/∂δ`.
    pub cy: f64,
    /// Pitching moment, `∂C_m/∂δ`.
    pub cm: f64,
    /// Rolling moment, `∂C_l/∂δ`.
    pub c_roll: f64,
    /// Yawing moment, `∂C_n/∂δ`.
    pub c_yaw: f64,
}

/// Why control derivatives could not be extracted.
#[derive(Debug, Clone, PartialEq)]
pub enum ControlError {
    /// The model has no control with that name.
    NoSuchControl,
    /// The VLM solve failed (including invalid flight conditions).
    Solve(VlmError),
}

impl core::fmt::Display for ControlError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoSuchControl => write!(f, "model has no control with that name"),
            Self::Solve(e) => write!(f, "VLM solve failed: {e}"),
        }
    }
}

impl core::error::Error for ControlError {}

impl From<VlmError> for ControlError {
    fn from(e: VlmError) -> Self {
        Self::Solve(e)
    }
}

/// Extract the derivatives of the control called `name` at angle of attack
/// `alpha_rad`, speed `v_inf` (m/s) and density (kg/m³), by central
/// difference of the VLM about the undeflected state.
///
/// # Errors
/// [`ControlError::NoSuchControl`] for an unknown name; otherwise any VLM
/// error, e.g. [`VlmError::InvalidFlightCondition`].
pub fn control_derivatives(
    model: &VlmModel,
    name: &str,
    alpha_rad: f64,
    v_inf: f64,
    density: f64,
) -> Result<ControlDerivatives, ControlError> {
    let index = model
        .control_index(name)
        .ok_or(ControlError::NoSuchControl)?;
    let q_s = 0.5 * density * v_inf * v_inf * model.reference().area;
    let h = 1.0e-3;
    let at = |delta: f64| -> Result<[f64; 6], ControlError> {
        let mut command = vec![0.0; model.controls().len()];
        command[index] = delta;
        let twists = model.control_twists(&command)?;
        let s = model.solve(alpha_rad, v_inf, density, &twists)?;
        let (roll, yaw) = s.stability_roll_yaw(alpha_rad);
        Ok([s.cl, s.cdi, s.force.y / q_s, s.cm, roll, yaw])
    };
    let (hi, lo) = (at(h)?, at(-h)?);
    let d = |i: usize| (hi[i] - lo[i]) / (2.0 * h);
    Ok(ControlDerivatives {
        cl: d(0),
        cd: d(1),
        cy: d(2),
        cm: d(3),
        c_roll: d(4),
        c_yaw: d(5),
    })
}
