//! Continuous atmospheric turbulence: Dryden and von Kármán models
//! (MIL-HDBK-1797 / MIL-F-8785C spectral forms) with discretized state-space
//! filter realizations for flight-dynamics simulation.

use alloc::vec::Vec;

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// Feet per metre: multiply a length in metres by this to get feet.
const FEET_PER_METRE: f64 = 3.280_84;

/// Turbulence spectral model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TurbulenceModel {
    /// Dryden spectra (rational transfer functions — exact realization).
    Dryden,
    /// von Kármán spectra (rational approximations per MIL-HDBK-1797).
    VonKarman,
}

/// One axis of the turbulence field, driven by unit white noise.
///
/// Two realizations, both with the stationary variance set to `σ²` exactly:
///
/// * **First order** (Dryden longitudinal, and von Kármán as an approximation):
///   `dx = −x/τ dt + σ√(2/τ) dW`, integrated with the exact discrete update
///   `x ← e^{−dt/τ} x + σ√(1 − e^{−2dt/τ}) ξ`, which is stable for any step.
/// * **Second order** (Dryden lateral and vertical, `vertical = true`): the
///   MIL-HDBK-1797 form `H(s) = σ√(3/τ)·(1 + τs/√3)/(1 + τs)²`, realized as two
///   cascaded first-order states, `v₁ = (1 + τs)⁻¹ w` and `v₂ = (1 + τs)⁻¹ v₁`,
///   with output `σ√(3/τ)·(v₁/√3 + (1 − 1/√3) v₂)`. The stationary variance is
///   `σ²` by the Lyapunov solution of this two-state system.
///
/// The von Kármán spectrum is approximated by the first-order filter with the
/// 1.339 pole factor. That is a scale-length correction, not the exact von
/// Kármán spectrum.
#[derive(Debug, Clone)]
pub struct TurbulenceAxis {
    /// Filter time constant τ (s).
    tau: f64,
    /// Output standard deviation (m/s).
    sigma: f64,
    /// Second-order (Dryden v/w) realization?
    second_order: bool,
    /// First-order state, or the first stage `v₁` of the second-order filter.
    x1: f64,
    /// Second stage `v₂` (second-order only).
    x2: f64,
    /// Current output (m/s).
    output: f64,
}

impl TurbulenceAxis {
    /// Intensity σ (m/s) and integral scale L (m) at true airspeed `v` (m/s).
    ///
    /// `vertical` selects the second-order Dryden form for the Dryden model
    /// (the lateral and vertical axes). A negative or non-finite `sigma` is
    /// treated as zero, and `length` and `speed` are floored at 10⁻⁶.
    #[must_use]
    pub fn new(
        model: TurbulenceModel,
        sigma: f64,
        length: f64,
        speed: f64,
        vertical: bool,
    ) -> Self {
        let sigma = if sigma.is_finite() && sigma > 0.0 {
            sigma
        } else {
            0.0
        };
        let v = speed.max(1.0e-6);
        let l = length.max(1.0e-6);
        let second_order = vertical && model == TurbulenceModel::Dryden;
        // von Kármán spectra carry a 1.339 pole-factor; Dryden does not.
        let pole_factor = match model {
            TurbulenceModel::Dryden => 1.0,
            TurbulenceModel::VonKarman => 1.339,
        };
        let tau = pole_factor * l / v;
        Self {
            tau,
            sigma,
            second_order,
            x1: 0.0,
            x2: 0.0,
            output: 0.0,
        }
    }

    /// Advance the filter one time step with a standard normal sample `w`.
    /// A non-positive `dt` leaves the state unchanged.
    pub fn step(&mut self, w: f64, dt: f64) {
        if !(dt.is_finite() && dt > 0.0) {
            return;
        }
        let tau = self.tau;
        if self.second_order {
            // Euler–Maruyama on the cascade: dv₁ = −v₁/τ dt + dW, dv₂ = (v₁−v₂)/τ dt.
            // Stable for dt < 2τ; the scale lengths keep τ well above typical steps.
            let dw = w * dt.sqrt();
            let v1_old = self.x1;
            self.x1 += -v1_old / tau * dt + dw;
            self.x2 += (v1_old - self.x2) / tau * dt;
            let k = self.sigma * (3.0 / tau).sqrt();
            let a = 1.0 / 3.0_f64.sqrt();
            self.output = k * (a * self.x1 + (1.0 - a) * self.x2);
        } else {
            // Exact update of the Ornstein–Uhlenbeck process with stationary σ².
            let phi = (-dt / tau).exp();
            self.x1 = phi * self.x1 + self.sigma * (1.0 - phi * phi).sqrt() * w;
            self.output = self.x1;
        }
    }

    /// Current turbulence velocity (m/s).
    #[must_use]
    pub const fn value(&self) -> f64 {
        self.output
    }

    /// The design output standard deviation (m/s).
    #[must_use]
    pub const fn sigma(&self) -> f64 {
        self.sigma
    }
}

/// A three-axis turbulence field.
#[derive(Debug, Clone)]
pub struct TurbulenceField {
    /// Longitudinal (along-track) component.
    pub u: TurbulenceAxis,
    /// Lateral component.
    pub v: TurbulenceAxis,
    /// Vertical component.
    pub w: TurbulenceAxis,
}

impl TurbulenceField {
    /// Build a field for the flight condition: true airspeed (m/s),
    /// altitude (m, for low-altitude scale lengths), exceedance intensity
    /// `sigma` (m/s, "light" ≈ 0.5, "moderate" ≈ 1.5, "severe" ≈ 3.0).
    ///
    /// Low-altitude scale lengths (MIL-HDBK-1797 §3.7.2.1):
    /// `L_w = h`, `L_u = L_v = h/(0.177 + 0.000823 h)^1.1` (capped at 150 m
    /// above 760 m where the high-altitude values take over).
    #[must_use]
    pub fn low_altitude(model: TurbulenceModel, speed: f64, altitude: f64, sigma: f64) -> Self {
        // MIL-HDBK-1797 low-altitude scales, in feet:
        // L_w = h, L_u = L_v = h / (0.177 + 0.000823 h)^1.2 (h in ft).
        // Valid to 2000 ft (610 m); above that the medium-altitude scale applies.
        let h_m = altitude.clamp(3.0, 610.0);
        let h_ft = h_m * FEET_PER_METRE;
        let l_u_ft = h_ft / (0.177 + 0.000_823 * h_ft).powf(1.2);
        let l_u = l_u_ft / FEET_PER_METRE;
        let l_w = h_m;
        Self {
            u: TurbulenceAxis::new(model, sigma, l_u, speed, false),
            v: TurbulenceAxis::new(model, sigma, l_u, speed, true),
            w: TurbulenceAxis::new(model, sigma, l_w, speed, true),
        }
    }

    /// Medium- and high-altitude isotropic scale (MIL-HDBK-1797): L = 1750 ft
    /// (533.4 m) on all three axes.
    #[must_use]
    pub fn medium_altitude(model: TurbulenceModel, speed: f64, sigma: f64) -> Self {
        const L: f64 = 533.4;
        Self {
            u: TurbulenceAxis::new(model, sigma, L, speed, false),
            v: TurbulenceAxis::new(model, sigma, L, speed, true),
            w: TurbulenceAxis::new(model, sigma, L, speed, true),
        }
    }

    /// Advance all three axes one step with independent white-noise samples.
    pub fn step(&mut self, wu: f64, wv: f64, ww: f64, dt: f64) {
        self.u.step(wu, dt);
        self.v.step(wv, dt);
        self.w.step(ww, dt);
    }

    /// Current `(u, v, w)` turbulence velocities (m/s).
    #[must_use]
    pub fn values(&self) -> (f64, f64, f64) {
        (self.u.value(), self.v.value(), self.w.value())
    }
}

/// Simple deterministic pseudo-random white-noise generator (xorshift64*),
/// so turbulence simulations are reproducible without external crates and
/// work in `no_std`.
#[derive(Debug, Clone)]
pub struct WhiteNoise {
    state: u64,
}

impl WhiteNoise {
    /// Seed the generator (a zero seed becomes 1).
    #[must_use]
    pub const fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 { 1 } else { seed },
        }
    }

    /// Next sample in `[-1, 1)` (uniform).
    pub fn next_uniform(&mut self) -> f64 {
        // xorshift64*
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        let v = x.wrapping_mul(0x2545_F491_4F6C_DD1D);
        (v >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
    }

    /// Next standard-normal sample (Box–Muller with cached second value).
    pub fn next_gaussian(&mut self) -> f64 {
        let u1 = (self.next_uniform() + 1.0) * 0.5; // (0, 1]
        let u2 = (self.next_uniform() + 1.0) * 0.5;
        (-2.0 * u1.ln()).sqrt() * (2.0 * core::f64::consts::PI * u2).cos()
    }
}

/// Simulate a turbulence axis for `n` steps; returns the sample series.
/// Statistical helper for tests and examples (drives the filter with unit
/// Gaussian noise; the axis scales internally).
#[must_use]
pub fn simulate_axis(
    axis: &mut TurbulenceAxis,
    noise: &mut WhiteNoise,
    dt: f64,
    n: usize,
) -> Vec<f64> {
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        axis.step(noise.next_gaussian(), dt);
        out.push(axis.value());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dryden_filter_is_stable_and_bounded() {
        let mut axis = TurbulenceAxis::new(TurbulenceModel::Dryden, 1.5, 100.0, 50.0, false);
        let mut noise = WhiteNoise::new(42);
        let series = simulate_axis(&mut axis, &mut noise, 0.02, 100_000);
        // Correlation time τ = L/V = 2 s = 100 samples ⇒ ~1000 independent
        // samples: the standard error of the mean is ~0.05.
        let mean: f64 = series.iter().sum::<f64>() / series.len() as f64;
        assert!(mean.abs() < 0.3, "mean {mean}");
        let var: f64 = series.iter().map(|x| x * x).sum::<f64>() / series.len() as f64;
        let std = var.sqrt();
        assert!(
            (1.5_f64 * 0.75..1.5 * 1.35).contains(&std),
            "std {std} should match the designed intensity (σ=1.5)"
        );
    }

    #[test]
    fn von_karman_filter_runs() {
        let mut axis = TurbulenceAxis::new(TurbulenceModel::VonKarman, 1.0, 200.0, 80.0, true);
        let mut noise = WhiteNoise::new(7);
        let series = simulate_axis(&mut axis, &mut noise, 0.01, 10_000);
        assert!(series.iter().all(|x| x.is_finite()));
    }

    #[test]
    fn white_noise_statistics() {
        let mut noise = WhiteNoise::new(123);
        let n = 50_000;
        let sum: f64 = (0..n).map(|_| noise.next_gaussian()).sum();
        assert!((sum / n as f64).abs() < 0.05, "mean drift");
        let var: f64 = (0..n)
            .map(|_| {
                let g = noise.next_gaussian();
                g * g
            })
            .sum::<f64>()
            / n as f64;
        assert!((0.8..1.2).contains(&var), "variance {var}");
    }

    #[test]
    fn field_low_altitude_scales() {
        // At 100 m altitude: L_w = 100 m, L_u ≈ 263 m.
        let field = TurbulenceField::low_altitude(TurbulenceModel::Dryden, 60.0, 100.0, 1.0);
        assert!(field.u.sigma() > 0.0);
        let mut field = field;
        field.step(0.1, 0.2, -0.1, 0.02);
        let (u, v, w) = field.values();
        assert!(u.is_finite() && v.is_finite() && w.is_finite());
    }
}

#[cfg(test)]
mod realization_tests {
    use super::*;

    fn std_of(axis: &mut TurbulenceAxis, dt: f64, n: usize, seed: u64) -> f64 {
        let mut noise = WhiteNoise::new(seed);
        let series = simulate_axis(axis, &mut noise, dt, n);
        let burn = n / 10;
        let tail = &series[burn..];
        let mean = tail.iter().sum::<f64>() / tail.len() as f64;
        (tail.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / tail.len() as f64).sqrt()
    }

    #[test]
    fn second_order_dryden_vertical_has_the_designed_intensity() {
        let mut axis = TurbulenceAxis::new(TurbulenceModel::Dryden, 1.5, 100.0, 50.0, true);
        let s = std_of(&mut axis, 0.02, 400_000, 11);
        assert!((s - 1.5).abs() < 0.1, "σ_w realized {s}");
    }

    #[test]
    fn first_order_is_stable_for_steps_much_larger_than_tau() {
        // τ = 2 s; a step of 50 s used to blow up under Euler integration.
        let mut axis = TurbulenceAxis::new(TurbulenceModel::Dryden, 1.0, 100.0, 50.0, false);
        let mut noise = WhiteNoise::new(3);
        for _ in 0..2_000 {
            axis.step(noise.next_gaussian(), 50.0);
            assert!(axis.value().is_finite() && axis.value().abs() < 20.0);
        }
    }

    #[test]
    fn negative_intensity_is_treated_as_calm_air() {
        let mut axis = TurbulenceAxis::new(TurbulenceModel::Dryden, -2.0, 100.0, 50.0, false);
        axis.step(1.0, 0.1);
        assert_eq!(axis.value(), 0.0);
    }

    #[test]
    fn low_altitude_longitudinal_scale_uses_feet() {
        // At 100 m (328 ft): L_u = 328 / (0.177 + 0.000823·328)^1.2 ft ≈ 862 ft ≈ 263 m.
        let f = TurbulenceField::low_altitude(TurbulenceModel::Dryden, 60.0, 100.0, 1.0);
        let l_u = f.u.tau * 60.0;
        assert!((l_u - 263.0).abs() < 3.0, "L_u = {l_u} m");
    }
}
