//! Linear plant-matrix construction for flight-control synthesis.
//!
//! The longitudinal short-period/heave-pitch plant in dimensional
//! linearized form, assembled from the extracted derivatives. This is the
//! handoff structure for a future `tpt-control` crate: state
//! `[Δu, Δα, q, Δθ]` (m/s, rad, rad/s, rad) with the standard small-disturbance
//! rigid-body dynamics decoupled from the lateral-directional set.

use crate::longitudinal::LongitudinalDerivatives;

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// Linearized longitudinal plant `Δẋ = A·Δx + B·Δδ` with state
/// `[Δu, Δα, q, Δθ]` and single control `Δδ` (flap/elevator).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LongitudinalPlant {
    /// Dynamics matrix (4×4, row-major).
    pub a: [[f64; 4]; 4],
    /// Control matrix (4×1).
    pub b: [[f64; 1]; 4],
}

/// Parameters used to dimensionalize the derivative-based plant.
#[derive(Debug, Clone, Copy)]
pub struct PlantParameters {
    /// Flight speed (m/s).
    pub v_inf: f64,
    /// Reference chord (m).
    pub chord: f64,
    /// Wing reference area (m²).
    pub area: f64,
    /// Mass (kg).
    pub mass: f64,
    /// Pitch inertia about the reference point (kg·m²).
    pub inertia_y: f64,
    /// Air density (kg/m³).
    pub density: f64,
}

impl LongitudinalPlant {
    /// Assemble the classical derivative-based longitudinal plant.
    ///
    /// The heave/pitch (short-period-like) block uses the standard
    /// quasi-steady relations with the derivatives' dimensionalization
    /// `q̂ = q c/(2V)`; the speed-damping row uses the simplified drag
    /// slope `C_{Du}` supplied by the caller (0 if unknown).
    ///
    /// `pitch_damping` multiplies q̇ by `C_{mq}` + `C_{Lα}`-driven pitch
    /// coupling — pass the raw extracted values via
    /// [`LongitudinalDerivatives`].
    #[must_use]
    pub fn from_derivatives(
        d: &LongitudinalDerivatives,
        p: &PlantParameters,
        drag_slope_u: f64,
    ) -> Self {
        let q_bar_s = 0.5 * p.density * p.v_inf * p.v_inf * p.area;
        let m_dot = p.mass * p.v_inf;
        // Δα row: m(V α̇ + u α_trim…) — linearized heave:
        // m V α̇ = −q S (C_{Lα} V α/V + C_{Lq} c q/(2V)) − m g sinθ ≈
        // folded into the standard matrix form below.
        let a_11 = -drag_slope_u * q_bar_s / m_dot;
        let a_12 = -d.cl_alpha * q_bar_s / m_dot;
        let a_21 = 0.0; // angle-of-attack equation insensitive to Δu at this order
        let a_22 = -d.cl_alpha * q_bar_s / m_dot;
        let q_damp = d.cm_q * q_bar_s * p.chord / (2.0 * p.v_inf) / p.inertia_y * p.chord;
        let alpha_to_q = d.cm_alpha * q_bar_s * p.chord / p.inertia_y;
        let a: [[f64; 4]; 4] = [
            [a_11, a_12, 0.0, 0.0],
            [a_21, a_22, 1.0, 0.0],
            [0.0, alpha_to_q, q_damp, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ];
        // Control: flap effectiveness drives both lift and pitch.
        let b_lift = d.cl_delta * q_bar_s / m_dot;
        let b_pitch = d.cm_delta * q_bar_s * p.area * p.chord / p.inertia_y;
        let b: [[f64; 1]; 4] = [[b_lift], [0.0], [b_pitch], [0.0]];
        Self { a, b }
    }

    /// Short-period characteristics: natural frequency and damping ratio of
    /// the (α, q) subsystem — the primary handling metrics.
    #[must_use]
    pub fn short_period(&self) -> (f64, f64) {
        let (a22, a23, a32, a33) = (self.a[1][1], self.a[1][2], self.a[2][1], self.a[2][2]);
        // (α, q) block: λ² − bλ + c = 0 with b = −(a22+a33), c = a22a33 − a23a32.
        let b = -(a22 + a33);
        let c = a22 * a33 - a23 * a32;
        // Standard second-order mapping (also valid overdamped, ζ > 1).
        let wn = c.max(0.0).sqrt();
        let zeta = if wn > 1.0e-9 {
            b / (2.0 * wn)
        } else {
            f64::NAN
        };
        (wn, zeta)
    }
}

/// Hook for a future `tpt-control` crate: exposes the plant in the standard
/// state-space form with a doc-stable interface.
impl LongitudinalPlant {
    /// The dynamics matrix as a flat row-major slice (for control-law crates
    /// that want raw buffers).
    #[must_use]
    pub fn a_flat(&self) -> [f64; 16] {
        let mut out = [0.0; 16];
        for (i, row) in self.a.iter().enumerate() {
            out[i * 4..i * 4 + 4].copy_from_slice(row);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::longitudinal::longitudinal_fd;
    use tpt_aero_core::WingPlanform;
    use tpt_aero_vlm::{CamberModel, VlmModel};

    #[test]
    fn short_period_is_damped() {
        // Reference the moment at the wing LEADING edge (ahead of the
        // aerodynamic center) so the airframe is statically stable — a
        // wing-alone plant referenced at its AC is neutrally stable by
        // construction.
        let wing = WingPlanform::rectangular(10.0, 1.0);
        let model = VlmModel::from_planform(&wing, 16, 4, CamberModel::None)
            .expect("valid")
            .with_reference(tpt_aero_vlm::ReferenceGeometry {
                area: wing.area(),
                chord: wing.mac(),
                span: wing.span(),
                moment_center: tpt_aero_core::Vec3::new(0.0, 0.0, 0.0),
            });
        let d = longitudinal_fd(&model, 2.0_f64.to_radians(), 50.0, 1.225, 1e-4);
        assert!(
            d.cm_alpha < 0.0,
            "C_mα = {} must be negative ahead of the AC",
            d.cm_alpha
        );
        let p = PlantParameters {
            v_inf: 50.0,
            chord: 1.0,
            area: 10.0,
            mass: 500.0,
            inertia_y: 400.0,
            density: 1.225,
        };
        let plant = LongitudinalPlant::from_derivatives(&d, &p, 0.05);
        let (wn, zeta) = plant.short_period();
        assert!(wn > 0.0, "ω_n = {wn}");
        assert!(zeta > 0.0, "damping ζ = {zeta} must be positive (stable)");
    }
}
