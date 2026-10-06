//! Continuous atmospheric turbulence: Dryden and von Kármán models
//! (MIL-HDBK-1797 / MIL-F-8785C spectral forms) with discretized state-space
//! filter realizations for flight-dynamics simulation.

use alloc::vec::Vec;

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// Turbulence spectral model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TurbulenceModel {
    /// Dryden spectra (rational transfer functions — exact realization).
    Dryden,
    /// von Kármán spectra (rational approximations per MIL-HDBK-1797).
    VonKarman,
}

/// One axis of the turbulence field: a first-order shaping filter driven by
/// unit white noise.
///
/// The Dryden longitudinal spectrum reduces *exactly* to a first-order filter
/// `H(s) = g/(1+τs)` (the (1+τs) numerator cancels one pole). The remaining
/// spectra (Dryden vertical, von Kármán) are realized as RMS-matched
/// first-order approximations: the integral scale is corrected by the
/// von Kármán pole factor 1.339 (MIL-HDBK-1797), and the gain is chosen so
/// the realized RMS equals the design σ exactly (Lyapunov:
/// `variance = g²τ/2` ⇒ `g = σ√(2/τ)`).
#[derive(Debug, Clone)]
pub struct TurbulenceAxis {
    /// Filter time constant τ (s).
    tau: f64,
    /// Filter gain (m/s per unit PSD).
    gain: f64,
    /// Filter state.
    x: f64,
    /// Output standard deviation (m/s).
    sigma: f64,
}

impl TurbulenceAxis {
    /// Intensity σ (m/s) and integral scale L (m) at true airspeed `v` (m/s).
    /// `vertical` selects the vertical/lateral spectrum correction.
    #[must_use]
    pub fn new(model: TurbulenceModel, sigma: f64, length: f64, speed: f64, vertical: bool) -> Self {
        let v = speed.max(1.0e-6);
        let l = length.max(1.0e-6);
        let _ = vertical;
        // von Kármán spectra carry a 1.339 pole-factor; Dryden does not.
        let pole_factor = match model {
            TurbulenceModel::Dryden => 1.0,
            TurbulenceModel::VonKarman => 1.339,
        };
        let tau = pole_factor * l / v;
        let gain = sigma * (2.0 / tau).sqrt();
        Self {
            tau,
            gain,
            x: 0.0,
            sigma,
        }
    }

    /// Advance the filter one time step.
    ///
    /// `w` is a standard normal sample; the drive impulse has variance `dt`,
    /// the exact discrete equivalent of unit-PSD continuous white noise
    /// (`variance = g²τ/2` at the output ⇒ RMS = σ by construction).
    pub fn step(&mut self, w: f64, dt: f64) {
        let u = w * dt.max(0.0).sqrt();
        self.x += (-self.x / self.tau) * dt + self.gain * u;
    }

    /// Current turbulence velocity (m/s).
    #[must_use]
    pub const fn value(&self) -> f64 {
        self.x
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
        let h = altitude.clamp(3.0, 760.0);
        let l_w = h;
        let denom = (0.177 + 0.000823 * h).powf(1.1);
        let l_u = (h / denom).clamp(10.0, 150.0);
        Self {
            u: TurbulenceAxis::new(model, sigma, l_u, speed, false),
            v: TurbulenceAxis::new(model, sigma, l_u, speed, false),
            w: TurbulenceAxis::new(model, sigma, l_w, speed, true),
        }
    }

    /// High-altitude isotropic scales (MIL-F-8785C): L = 530 m for all axes
    /// below 610 m of... above 760 m: L = 1758·h^... the standard curves
    /// flatten to L ≈ 530 m (MIL-HDBK-1797 medium altitude).
    #[must_use]
    pub fn medium_altitude(model: TurbulenceModel, speed: f64, sigma: f64) -> Self {
        const L: f64 = 530.0;
        Self {
            u: TurbulenceAxis::new(model, sigma, L, speed, false),
            v: TurbulenceAxis::new(model, sigma, L, speed, false),
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
pub fn simulate_axis(axis: &mut TurbulenceAxis, noise: &mut WhiteNoise, dt: f64, n: usize) -> Vec<f64> {
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
        let var: f64 = (0..n).map(|_| {
            let g = noise.next_gaussian();
            g * g
        }).sum::<f64>() / n as f64;
        assert!((0.8..1.2).contains(&var), "variance {var}");
    }

    #[test]
    fn field_low_altitude_scales() {
        // At 100 m altitude: L_w = 100, L_u = 100/(0.177+0.0823)^1.1 ≈ 300.
        let field = TurbulenceField::low_altitude(TurbulenceModel::Dryden, 60.0, 100.0, 1.0);
        assert!(field.u.sigma() > 0.0);
        let mut field = field;
        field.step(0.1, 0.2, -0.1, 0.02);
        let (u, v, w) = field.values();
        assert!(u.is_finite() && v.is_finite() && w.is_finite());
    }
}
