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
    /// von Kármán spectra (MIL-HDBK-1797): rational shaping filters fitted to
    /// the handbook spectra, third order for the lateral and vertical axes and
    /// second order for the longitudinal axis. Accurate to a few per cent for
    /// `L·ω/V ≤ 3`; the handbook says the true spectrum needs higher-order
    /// filters above that.
    VonKarman,
}

/// One axis of the turbulence field, driven by unit white noise.
///
/// Two realizations, both with the stationary variance set to `σ²` exactly:
///
/// * **First order** (Dryden longitudinal):
///   `dx = −x/τ dt + σ√(2/τ) dW`, integrated with the exact discrete update
///   `x ← e^{−dt/τ} x + σ√(1 − e^{−2dt/τ}) ξ`, which is stable for any step.
/// * **Second order** (Dryden lateral and vertical, `vertical = true`): the
///   MIL-HDBK-1797 form `H(s) = σ√(3/τ)·(1 + τs/√3)/(1 + τs)²`, realized as two
///   cascaded first-order states, `v₁ = (1 + τs)⁻¹ w` and `v₂ = (1 + τs)⁻¹ v₁`,
///   with output `σ√(3/τ)·(v₁/√3 + (1 − 1/√3) v₂)`. The stationary variance is
///   `σ²` by the Lyapunov solution of this two-state system.
///
/// The von Kármán model uses internal rational shaping filters (`ShapingFilter`),
/// whose spectra follow the MIL-HDBK-1797 von Kármán forms.
#[derive(Debug, Clone)]
pub struct TurbulenceAxis {
    /// Filter time constant τ = L/V (s).
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
    /// The von Kármán shaping filter (`None` for Dryden).
    shaping: Option<ShapingFilter>,
}

impl TurbulenceAxis {
    /// Intensity σ (m/s) and scale length `length` (m) at true airspeed `v`
    /// (m/s).
    ///
    /// `vertical` selects the lateral and vertical axes: the second-order
    /// Dryden form, or the third-order von Kármán filter. Otherwise the axis is
    /// longitudinal. For von Kármán, `length` is the filter scale: `L_u` on the
    /// longitudinal axis and `2L_v`, `2L_w` on the lateral and vertical axes,
    /// the MIL-HDBK-1797 convention in which `L_u = 2L_v = 2L_w`. A negative or
    /// non-finite `sigma` is treated as zero, and `length` and `speed` are
    /// floored at 10⁻⁶.
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
        let tau = l / v;
        let shaping = match model {
            TurbulenceModel::Dryden => None,
            TurbulenceModel::VonKarman => Some(ShapingFilter::von_karman(vertical, sigma, tau)),
        };
        Self {
            tau,
            sigma,
            second_order,
            x1: 0.0,
            x2: 0.0,
            output: 0.0,
            shaping,
        }
    }

    /// Advance the filter one time step with a standard normal sample `w`.
    /// A non-positive `dt` leaves the state unchanged.
    ///
    /// The von Kármán filter is integrated by Euler–Maruyama, so `dt` should
    /// stay below about `0.1·τ` (the fastest lateral pole is at `11/τ`).
    pub fn step(&mut self, w: f64, dt: f64) {
        if !(dt.is_finite() && dt > 0.0) {
            return;
        }
        if let Some(filter) = &mut self.shaping {
            self.output = filter.step(w, dt);
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

/// Von Kármán lateral and vertical shaping numerator, `N(s)` with `s = (L/V)·p`
/// and constant term 1. Least-squares fit to the MIL-HDBK-1797 spectrum
/// `(1 + (8/3)(1.339x)²)/(1 + (1.339x)²)^(11/6)`, `x = Lω/V`.
const VK_LATERAL_NUM: [f64; 3] = [1.0, 2.7478, 0.3398];
/// Von Kármán lateral and vertical shaping denominator `D(s)` (same fit).
const VK_LATERAL_DEN: [f64; 4] = [1.0, 2.9958, 1.9754, 0.1539];
/// Von Kármán longitudinal shaping numerator (fit to `(1 + (1.339x)²)^(−5/6)`).
const VK_LONGITUDINAL_NUM: [f64; 2] = [1.0, 0.25];
/// Von Kármán longitudinal shaping denominator `D(s)` (same fit).
const VK_LONGITUDINAL_DEN: [f64; 3] = [1.0, 1.357, 0.1987];

/// `x^k` for a small non-negative integer power.
fn pow_int(x: f64, k: usize) -> f64 {
    (0..k).fold(1.0, |acc, _| acc * x)
}

/// `|P(jx)|²` for a real polynomial `P(s) = Σ c_k s^k`.
fn poly_mag2(c: &[f64], x: f64) -> f64 {
    let (mut re, mut im) = (0.0, 0.0);
    for (k, &ck) in c.iter().enumerate() {
        let term = ck * pow_int(x, k);
        // s^k = j^k x^k: the powers cycle through 1, j, −1, −j.
        match k % 4 {
            0 => re += term,
            1 => im += term,
            2 => re -= term,
            _ => im -= term,
        }
    }
    re * re + im * im
}

/// `∫₀^∞ |N(jx)/D(jx)|² dx`, by the midpoint rule in `θ` with `x = tan θ`.
fn shaping_integral(num: &[f64], den: &[f64]) -> f64 {
    const STEPS: usize = 4096;
    let h = core::f64::consts::FRAC_PI_2 / STEPS as f64;
    (0..STEPS)
        .map(|i| {
            let theta = (i as f64 + 0.5) * h;
            let x = theta.tan();
            // dx = sec²θ dθ.
            let ratio = poly_mag2(num, x) / poly_mag2(den, x);
            ratio * h / (theta.cos() * theta.cos())
        })
        .sum()
}

/// A rational shaping filter `H(p) = g·N(τp)/D(τp)` driven by unit white noise,
/// in controllable companion form.
///
/// The states satisfy `ẋᵢ = xᵢ₊₁` and `ẋₙ₋₁ = −Σ aₖxₖ + w`, and the output is
/// `Σ bₖxₖ`. The gain `g` sits in the `b` coefficients. The state update is
/// Euler–Maruyama, the same scheme as the second-order Dryden branch.
#[derive(Debug, Clone)]
struct ShapingFilter {
    /// Coefficients `a₀..aₙ₋₁` of the monic denominator `pⁿ + Σ aₖpᵏ`.
    a: Vec<f64>,
    /// Coefficients `b₀..bₙ₋₁` of the numerator, gain included.
    b: Vec<f64>,
    /// States `x₀..xₙ₋₁`.
    x: Vec<f64>,
}

impl ShapingFilter {
    /// The von Kármán filter for the lateral and vertical axes (`vertical`) or
    /// the longitudinal axis, with stationary standard deviation `sigma` and
    /// time constant `tau = L/V`.
    fn von_karman(vertical: bool, sigma: f64, tau: f64) -> Self {
        let (num, den): (&[f64], &[f64]) = if vertical {
            (&VK_LATERAL_NUM, &VK_LATERAL_DEN)
        } else {
            (&VK_LONGITUDINAL_NUM, &VK_LONGITUDINAL_DEN)
        };
        // With unit white noise, the variance of g·N/D is (g²/π)(V/L)∫|N/D|²dx.
        // Setting that to σ² gives g = σ√(πτ/I), with I the integral of the
        // rational filter itself. This makes the realised intensity exact.
        let gain = sigma * (core::f64::consts::PI * tau / shaping_integral(num, den)).sqrt();
        let n = den.len() - 1;
        // Divide through by the leading coefficient of D(τp) to make it monic.
        let lead = den[n] * pow_int(tau, n);
        let a = (0..n).map(|k| den[k] * pow_int(tau, k) / lead).collect();
        let b = (0..n)
            .map(|k| {
                let nk = num.get(k).copied().unwrap_or(0.0);
                gain * nk * pow_int(tau, k) / lead
            })
            .collect();
        Self {
            a,
            b,
            x: alloc::vec![0.0; n],
        }
    }

    /// Advance one step with white-noise sample `w` and return the output.
    fn step(&mut self, w: f64, dt: f64) -> f64 {
        let n = self.x.len();
        // Feedback from the old states, before any of them is updated.
        let feedback: f64 = self.a.iter().zip(&self.x).map(|(a, x)| -a * x).sum();
        for i in 0..n - 1 {
            self.x[i] += self.x[i + 1] * dt;
        }
        self.x[n - 1] += feedback * dt + w * dt.sqrt();
        self.b.iter().zip(&self.x).map(|(b, x)| b * x).sum()
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
    /// `L_w = h`, `L_u = L_v = h/(0.177 + 0.000823 h)^1.2` (h in feet). The
    /// altitude is clamped to 3 to 610 m; above 610 m the medium-altitude
    /// preset applies instead.
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

    /// Medium- and high-altitude scales (MIL-HDBK-1797): 1750 ft (533.4 m) on
    /// all three axes for Dryden, and 2500 ft (762 m) for von Kármán, where
    /// `L_u = 2L_v = 2L_w = 2500 ft` puts the filter length at 2500 ft on every
    /// axis.
    #[must_use]
    pub fn medium_altitude(model: TurbulenceModel, speed: f64, sigma: f64) -> Self {
        let l = match model {
            TurbulenceModel::Dryden => 533.4,
            TurbulenceModel::VonKarman => 762.0,
        };
        Self {
            u: TurbulenceAxis::new(model, sigma, l, speed, false),
            v: TurbulenceAxis::new(model, sigma, l, speed, true),
            w: TurbulenceAxis::new(model, sigma, l, speed, true),
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

    #[test]
    fn von_karman_shaping_matches_handbook_spectrum_in_band() {
        // MIL-HDBK-1797 shapes, x = Lω/V; the fitted filters are checked over
        // 0 < x ≤ 3 against the spectrum they approximate.
        let lateral = |x: f64| {
            (1.0 + 8.0 / 3.0 * (1.339 * x).powi(2)) / (1.0 + (1.339 * x).powi(2)).powf(11.0 / 6.0)
        };
        let longitudinal = |x: f64| (1.0 + (1.339 * x).powi(2)).powf(-5.0 / 6.0);
        let mut worst_lat: f64 = 0.0;
        let mut worst_lon: f64 = 0.0;
        for i in 1..=300 {
            let x = f64::from(i) * 0.01;
            let r_lat = poly_mag2(&VK_LATERAL_NUM, x) / poly_mag2(&VK_LATERAL_DEN, x);
            let r_lon = poly_mag2(&VK_LONGITUDINAL_NUM, x) / poly_mag2(&VK_LONGITUDINAL_DEN, x);
            worst_lat = worst_lat.max((r_lat / lateral(x) - 1.0).abs());
            worst_lon = worst_lon.max((r_lon / longitudinal(x) - 1.0).abs());
        }
        assert!(worst_lat < 0.02, "lateral error {worst_lat}");
        assert!(worst_lon < 0.04, "longitudinal error {worst_lon}");
    }

    #[test]
    fn von_karman_lateral_has_the_designed_intensity() {
        // τ = 762 m / 50 m/s ≈ 15 s; the slowest pole is at 0.48/τ.
        let mut axis = TurbulenceAxis::new(TurbulenceModel::VonKarman, 1.5, 762.0, 50.0, true);
        let s = std_of(&mut axis, 0.1, 400_000, 5);
        assert!((s - 1.5).abs() < 0.12, "σ_v realized {s}");
    }

    #[test]
    fn von_karman_longitudinal_has_the_designed_intensity() {
        let mut axis = TurbulenceAxis::new(TurbulenceModel::VonKarman, 1.0, 762.0, 50.0, false);
        let s = std_of(&mut axis, 0.1, 400_000, 17);
        assert!((s - 1.0).abs() < 0.08, "σ_u realized {s}");
    }

    #[test]
    fn von_karman_filter_is_stable_at_a_tenth_of_its_time_constant() {
        // τ = 100 m / 50 m/s = 2 s; dt = 0.2 s. The fastest lateral pole gives
        // Euler amplification factor 1 − 11.14·0.1 ≈ −0.11, which is stable.
        let mut axis = TurbulenceAxis::new(TurbulenceModel::VonKarman, 1.0, 100.0, 50.0, true);
        let mut noise = WhiteNoise::new(9);
        let series = simulate_axis(&mut axis, &mut noise, 0.2, 50_000);
        assert!(series.iter().all(|x| x.is_finite() && x.abs() < 20.0));
    }

    #[test]
    fn von_karman_medium_altitude_uses_2500_ft_on_every_axis() {
        let f = TurbulenceField::medium_altitude(TurbulenceModel::VonKarman, 50.0, 1.0);
        for axis in [&f.u, &f.v, &f.w] {
            assert!(
                (axis.tau * 50.0 - 762.0).abs() < 1e-9,
                "L = {} m",
                axis.tau * 50.0
            );
        }
    }
}
