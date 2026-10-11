//! Linear plant-matrix construction for flight-control synthesis.
//!
//! The longitudinal short-period/heave-pitch plant in dimensional
//! linearized form, assembled from the extracted derivatives. This is the
//! handoff structure for a future `tpt-control` crate: state
//! `[Δu, Δα, q, Δθ]` (m/s, rad, rad/s, rad) with the standard small-disturbance
//! rigid-body dynamics decoupled from the lateral-directional set.

use crate::longitudinal::LongitudinalDerivatives;

/// Standard gravity (m/s²).
const STANDARD_GRAVITY: f64 = 9.80665;

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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
    /// Small-disturbance rigid-body dynamics about level, steady flight
    /// (`γ₀ = 0`, `θ₀ = α₀`), in wind axes with state `[Δu, Δα, q, Δθ]`:
    ///
    /// ```text
    /// Δu̇ = −(C_Du q̄S/(mV)) Δu − g Δθ
    /// Δα̇ = q − (q̄S/(mV)) (C_Lα Δα + C_Lq (c/2V) q + C_Lδ Δδ)
    /// q̇  = (q̄S c/I_y) (C_mα Δα + C_mq (c/2V) q + C_mδ Δδ)
    /// Δθ̇ = q
    /// ```
    ///
    /// `drag_slope_u` is the caller-supplied `C_{Du}` (0 if unknown). Speed
    /// effects on lift and pitching moment and the `C_{Dα}` coupling are
    /// neglected — the usual short-period/phugoid first-order model.
    #[must_use]
    pub fn from_derivatives(
        d: &LongitudinalDerivatives,
        p: &PlantParameters,
        drag_slope_u: f64,
    ) -> Self {
        let q_bar_s = 0.5 * p.density * p.v_inf * p.v_inf * p.area;
        let mv = p.mass * p.v_inf;
        let half_c_over_v = p.chord / (2.0 * p.v_inf);
        let lift_per_alpha = q_bar_s / mv;

        let x_u = -drag_slope_u * lift_per_alpha;
        let z_alpha = d.cl_alpha * lift_per_alpha;
        let z_q = d.cl_q * lift_per_alpha * half_c_over_v;
        let m_scale = q_bar_s * p.chord / p.inertia_y;
        let m_alpha = d.cm_alpha * m_scale;
        let m_q = d.cm_q * m_scale * half_c_over_v;

        let a: [[f64; 4]; 4] = [
            [x_u, 0.0, 0.0, -STANDARD_GRAVITY],
            [0.0, -z_alpha, 1.0 - z_q, 0.0],
            [0.0, m_alpha, m_q, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ];
        // Flap/elevator effectiveness: positive δ raises lift (so α̇ falls)
        // and moves the pitching moment by C_mδ.
        let b_alpha = -d.cl_delta * lift_per_alpha;
        let b_pitch = d.cm_delta * m_scale;
        let b: [[f64; 1]; 4] = [[0.0], [b_alpha], [b_pitch], [0.0]];
        Self { a, b }
    }

    /// Short-period characteristics `(ω_n, ζ)` of the `(α, q)` subsystem —
    /// the primary handling metrics.
    ///
    /// Returns `None` when the subsystem has no oscillatory/stable-pair
    /// mapping, i.e. `det ≤ 0`: the airframe is statically unstable (a real
    /// positive root), so `ω_n` and `ζ` are undefined.
    #[must_use]
    pub fn short_period(&self) -> Option<(f64, f64)> {
        let (a22, a23, a32, a33) = (self.a[1][1], self.a[1][2], self.a[2][1], self.a[2][2]);
        // (α, q) block: λ² − bλ + c = 0 with b = −(a22+a33), c = a22a33 − a23a32.
        let b = -(a22 + a33);
        let c = a22 * a33 - a23 * a32;
        if c.is_nan() || c <= 0.0 {
            return None;
        }
        let wn = c.sqrt();
        Some((wn, b / (2.0 * wn)))
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
        let d = longitudinal_fd(&model, 2.0_f64.to_radians(), 50.0, 1.225, 1e-4).expect("solve");
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
        let (wn, zeta) = plant.short_period().expect("statically stable");
        assert!(wn > 0.0, "ω_n = {wn}");
        assert!(zeta > 0.0, "damping ζ = {zeta} must be positive (stable)");
    }

    #[test]
    fn plant_matrix_is_dimensionally_consistent() {
        let d = LongitudinalDerivatives {
            cl_alpha: 5.0,
            cm_alpha: -1.0,
            cl_q: 3.0,
            cm_q: -10.0,
            cl_delta: 2.0,
            cm_delta: -1.5,
        };
        let p = PlantParameters {
            v_inf: 50.0,
            chord: 1.0,
            area: 10.0,
            mass: 500.0,
            inertia_y: 400.0,
            density: 1.225,
        };
        let plant = LongitudinalPlant::from_derivatives(&d, &p, 0.1);
        let q_bar_s = 0.5 * 1.225 * 50.0 * 50.0 * 10.0;
        // Control effectiveness: S enters once, via q̄S.
        let b_pitch = -1.5 * q_bar_s * 1.0 / 400.0;
        assert!((plant.b[2][0] - b_pitch).abs() < 1e-12);
        // Lift control acts through the α row, with opposite sign.
        assert_eq!(plant.b[0][0], 0.0);
        assert!((plant.b[1][0] + 2.0 * q_bar_s / (500.0 * 50.0)).abs() < 1e-12);
        // Gravity couples pitch attitude into speed.
        assert!((plant.a[0][3] + 9.80665).abs() < 1e-12);
        // C_Lq enters the α row's q coupling.
        let z_q = 3.0 * q_bar_s / (500.0 * 50.0) * 0.01;
        assert!((plant.a[1][2] - (1.0 - z_q)).abs() < 1e-12);
    }

    #[test]
    fn statically_unstable_has_no_short_period() {
        let d = LongitudinalDerivatives {
            cl_alpha: 5.0,
            cm_alpha: 0.5,
            cl_q: 0.0,
            cm_q: -10.0,
            cl_delta: 0.0,
            cm_delta: 0.0,
        };
        let p = PlantParameters {
            v_inf: 50.0,
            chord: 1.0,
            area: 10.0,
            mass: 500.0,
            inertia_y: 400.0,
            density: 1.225,
        };
        let plant = LongitudinalPlant::from_derivatives(&d, &p, 0.0);
        assert!(plant.short_period().is_none());
    }
}
