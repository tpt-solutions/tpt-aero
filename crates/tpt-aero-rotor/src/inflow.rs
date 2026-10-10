//! Inflow models beyond steady axial BEMT: skewed wake modulation and
//! first-order dynamic inflow (Pitt–Peters axial mode).

// Brings the `no_std` float methods (`atan`, `powf`, ...) into scope.
#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// Skewed-wake inflow modulation for oblique inflow (forward flight,
/// helicopter/VTOL transitions).
///
/// The first-harmonic Glauert series: with shaft inflow ratio `λ₀` and wake
/// skew angle `χ` measured from the shaft axis, the local inflow at azimuth
/// ψ and span fraction x is modulated as
/// `λ(ψ) ≈ λ₀ · (1 + E · x · cos ψ)` with `E = tan χ` (first term of the
/// Glauert/Pitt–Peters expansion). Engineering accuracy for |χ| ≤ 45°.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SkewedWake {
    /// Wake skew angle χ (rad) from the shaft axis; 0 = axial.
    pub skew_angle: f64,
}

impl SkewedWake {
    /// Axial wake (no modulation).
    #[must_use]
    pub const fn axial() -> Self {
        Self { skew_angle: 0.0 }
    }

    /// The modulation factor for azimuth `ψ` (rad; 0 = retreating side
    /// convention of the standard rotor diagram) and span fraction `x`.
    #[must_use]
    pub fn modulation(&self, psi: f64, x: f64) -> f64 {
        1.0 + self.skew_angle.tan() * x * psi.cos()
    }

    /// Modulated inflow ratio at `(ψ, x)` for the axial value `λ₀`.
    #[must_use]
    pub fn inflow(&self, lambda0: f64, psi: f64, x: f64) -> f64 {
        lambda0 * self.modulation(psi, x)
    }
}

/// First-order dynamic inflow (axial Pitt–Peters mode).
///
/// Models the rotor-induced inflow as a first-order lag toward the steady
/// BEMT value with the actuator-disk time constant
/// `τ = 1.13 R / (V + 2 v_ss)` (Leishman ch. 4 approximation of the
/// Peters–He finite-state model's axial pole).
///
/// clamps nothing away from physicality: the state is clamped to `[0, 1]`
/// like the steady induction factor it tracks.
#[derive(Debug, Clone, Copy)]
pub struct DynamicInflow {
    /// Rotor radius (m).
    radius: f64,
    /// Current induced inflow ratio state.
    lambda: f64,
}

impl DynamicInflow {
    /// Initialize the lag state (commonly the steady value at trim).
    #[must_use]
    pub const fn new(radius: f64, initial_lambda: f64) -> Self {
        Self {
            radius,
            lambda: initial_lambda,
        }
    }

    /// Current induced inflow ratio.
    #[must_use]
    pub const fn lambda(&self) -> f64 {
        self.lambda
    }

    /// The actuator-disk lag time constant (s) for axial speed `v` (m/s) and
    /// steady induced velocity `v_ss` (m/s).
    #[must_use]
    pub fn time_constant(&self, axial_speed: f64, steady_induced: f64) -> f64 {
        let denom = (axial_speed + 2.0 * steady_induced).abs().max(1.0e-6);
        1.13 * self.radius / denom
    }

    /// Advance the lag one time step `dt` (s) toward a new steady state.
    ///
    /// Uses the exact ZOH solution of `τ λ̇ = λ_ss − λ`, i.e. unconditionally
    /// stable for any `dt > 0`.
    ///
    /// `steady_induced_ratio` is the dimensionless `λ_ss = v_ss/(ΩR)`;
    /// `tip_speed` (`ΩR`, m/s) converts it to the m/s the time constant needs.
    pub fn step(&mut self, axial_speed: f64, steady_induced_ratio: f64, tip_speed: f64, dt: f64) {
        let tau = self.time_constant(axial_speed, steady_induced_ratio * tip_speed);
        let alpha = 1.0 - (-dt.max(0.0) / tau).exp();
        let target = steady_induced_ratio.clamp(0.0, 1.0);
        self.lambda += (target - self.lambda) * alpha;
        self.lambda = self.lambda.clamp(0.0, 1.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn axial_wake_has_no_modulation() {
        let sw = SkewedWake::axial();
        assert!((sw.modulation(1.234, 0.7) - 1.0).abs() < 1e-12);
        assert!((sw.inflow(0.08, 2.0, 0.5) - 0.08).abs() < 1e-12);
    }

    #[test]
    fn skewed_wake_first_harmonic() {
        let sw = SkewedWake { skew_angle: 0.4 };
        let m = sw.modulation(0.0, 0.7);
        let expected = 1.0 + 0.4_f64.tan() * 0.7;
        assert!((m - expected).abs() < 1e-12);
        // Advancing side (ψ=π) modulates opposite to retreating (ψ=0).
        assert!(sw.modulation(0.0, 0.7) > 1.0);
        assert!(sw.modulation(core::f64::consts::PI, 0.7) < 1.0);
    }

    #[test]
    fn dynamic_inflow_converges_to_steady_state() {
        let mut di = DynamicInflow::new(1.0, 0.0);
        let dt = 0.01;
        for _ in 0..20000 {
            di.step(0.0, 0.07, 100.0, dt);
        }
        assert!((di.lambda() - 0.07).abs() < 1e-4, "λ = {}", di.lambda());
        // Stability: a huge step must not overshoot the target.
        let mut di = DynamicInflow::new(1.0, 0.0);
        di.step(0.0, 0.07, 100.0, 1.0e6);
        assert!((di.lambda() - 0.07).abs() < 1e-9);
        // Bound: never leaves [0, 1].
        let mut di = DynamicInflow::new(1.0, 0.999);
        di.step(0.0, 2.0, 100.0, 0.1);
        assert!(di.lambda() <= 1.0);
    }

    #[test]
    fn time_constant_scales_physically() {
        let di = DynamicInflow::new(2.0, 0.0);
        // Hover: τ = 1.13 R / (2 v_i).
        let tau_hover = di.time_constant(0.0, 5.65);
        assert!((tau_hover - 1.13 * 2.0 / 11.3).abs() < 1e-9);
        // Faster axial flow → shorter lag.
        assert!(di.time_constant(50.0, 1.0) < tau_hover);
    }
}
