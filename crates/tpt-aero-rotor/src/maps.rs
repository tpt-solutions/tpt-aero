//! Propeller / rotor performance mapping over the advance-ratio axis.

use alloc::vec::Vec;

use crate::bemt::{BemtSolver, OperatingCondition};
use crate::blade::RotorError;

/// One point of the performance map.
#[derive(Debug, Clone, Copy)]
pub struct MapPoint {
    /// Advance ratio `J = V/(nD)` (0 at hover).
    pub j: f64,
    /// Thrust coefficient.
    pub ct: f64,
    /// Torque coefficient.
    pub cq: f64,
    /// Power coefficient `C_P = C_Q + λ_c C_T`.
    pub cp: f64,
    /// Propeller efficiency `η = J C_T / (2π C_Q)`.
    pub efficiency: f64,
    /// Figure of merit (meaningful near hover).
    pub figure_of_merit: f64,
}

/// Sweep a rotor/propeller over a range of axial speeds at fixed shaft speed.
///
/// # Errors
/// Propagates [`RotorError`] from any individual BEMT solve.
pub fn performance_map(
    solver: &BemtSolver,
    omega: f64,
    density: f64,
    axial_speeds: &[f64],
    collective: f64,
) -> Result<Vec<MapPoint>, RotorError> {
    axial_speeds
        .iter()
        .map(|&v| {
            let op = OperatingCondition {
                axial_speed: v,
                omega,
                density,
                collective,
            };
            let sol = solver.solve(&op)?;
            Ok(MapPoint {
                j: core::f64::consts::PI * v / (omega * solver.geometry().tip_radius),
                ct: sol.loads.ct,
                cq: sol.loads.cq,
                cp: sol.loads.cp,
                efficiency: sol.loads.propeller_efficiency,
                figure_of_merit: sol.loads.figure_of_merit,
            })
        })
        .collect()
}
