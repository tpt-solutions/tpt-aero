//! Time-domain simulation of the linear plants.
//!
//! A fixed-step RK4 integrator for `x' = A x + B u(t)`. Inputs are closures of
//! time, so control histories, discrete gusts and turbulence samples all plug
//! in the same way, and this crate needs no dependency on the atmosphere
//! crate.

use alloc::vec::Vec;

/// A linear time-invariant system `x' = A x + B u` with `N` states and `M`
/// inputs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StateSpace<const N: usize, const M: usize> {
    /// Dynamics matrix.
    pub a: [[f64; N]; N],
    /// Input matrix.
    pub b: [[f64; M]; N],
}

impl<const N: usize, const M: usize> StateSpace<N, M> {
    /// `A x + B u`.
    #[must_use]
    pub fn derivative(&self, x: &[f64; N], u: &[f64; M]) -> [f64; N] {
        core::array::from_fn(|i| {
            let ax: f64 = (0..N).map(|j| self.a[i][j] * x[j]).sum();
            let bu: f64 = (0..M).map(|j| self.b[i][j] * u[j]).sum();
            ax + bu
        })
    }

    /// Advance one RK4 step of size `dt` from time `t`; the input is sampled
    /// at `t`, `t + dt/2` and `t + dt`.
    #[must_use]
    pub fn step(
        &self,
        x: &[f64; N],
        t: f64,
        dt: f64,
        input: &impl Fn(f64) -> [f64; M],
    ) -> [f64; N] {
        let shifted = |base: &[f64; N], k: &[f64; N], h: f64| -> [f64; N] {
            core::array::from_fn(|i| base[i] + h * k[i])
        };
        let k1 = self.derivative(x, &input(t));
        let k2 = self.derivative(&shifted(x, &k1, 0.5 * dt), &input(t + 0.5 * dt));
        let k3 = self.derivative(&shifted(x, &k2, 0.5 * dt), &input(t + 0.5 * dt));
        let k4 = self.derivative(&shifted(x, &k3, dt), &input(t + dt));
        core::array::from_fn(|i| x[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]))
    }
}

/// One recorded sample of a simulation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sample<const N: usize> {
    /// Time (s).
    pub t: f64,
    /// State at that time.
    pub x: [f64; N],
}

/// Simulate `steps` RK4 steps of size `dt` from `x0` at `t = 0`, returning
/// `steps + 1` samples (including the initial state).
///
/// Returns an empty history for a non-finite or non-positive `dt`, or when
/// the state stops being finite (a diverging run is cut off at the first
/// non-finite sample, which is not recorded).
#[must_use]
pub fn simulate<const N: usize, const M: usize>(
    system: &StateSpace<N, M>,
    x0: [f64; N],
    dt: f64,
    steps: usize,
    input: impl Fn(f64) -> [f64; M],
) -> Vec<Sample<N>> {
    if !(dt.is_finite() && dt > 0.0) {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(steps + 1);
    let mut x = x0;
    out.push(Sample { t: 0.0, x });
    for k in 0..steps {
        let t = k as f64 * dt;
        x = system.step(&x, t, dt, &input);
        if !x.iter().all(|v| v.is_finite()) {
            break;
        }
        out.push(Sample {
            t: (k + 1) as f64 * dt,
            x,
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rk4_reproduces_an_exponential_decay() {
        let sys = StateSpace::<1, 1> {
            a: [[-1.0]],
            b: [[0.0]],
        };
        let hist = simulate(&sys, [1.0], 0.01, 300, |_| [0.0]);
        let last = hist.last().unwrap();
        assert!((last.t - 3.0).abs() < 1e-12);
        assert!((last.x[0] - (-3.0_f64).exp()).abs() < 1e-9);
    }

    #[test]
    fn rk4_preserves_a_harmonic_oscillator() {
        // x'' = -4x: period pi.
        let sys = StateSpace::<2, 1> {
            a: [[0.0, 1.0], [-4.0, 0.0]],
            b: [[0.0], [0.0]],
        };
        let hist = simulate(&sys, [1.0, 0.0], 0.005, 2000, |_| [0.0]);
        let energy = |s: &Sample<2>| 0.5 * s.x[1] * s.x[1] + 2.0 * s.x[0] * s.x[0];
        let e0 = energy(&hist[0]);
        assert!(hist.iter().all(|s| (energy(s) - e0).abs() < 1e-6 * e0));
        // After exactly two periods (t = 2 pi) the state returns.
        let k = (2.0 * core::f64::consts::PI / 0.005).round() as usize;
        assert!((hist[k].x[0] - 1.0).abs() < 1e-4);
    }

    #[test]
    fn constant_input_drives_the_expected_steady_state_and_time_varying_input_is_sampled() {
        // x' = -2x + 3u: steady state 1.5 u.
        let sys = StateSpace::<1, 1> {
            a: [[-2.0]],
            b: [[3.0]],
        };
        let hist = simulate(&sys, [0.0], 0.01, 1000, |_| [2.0]);
        assert!((hist.last().unwrap().x[0] - 3.0).abs() < 1e-6);
        // Ramp input u = t, from rest: x(t) = 1.5 t - 0.75 + 0.75 e^{-2t}.
        let hist = simulate(&sys, [0.0], 0.01, 200, |t| [t]);
        let t = hist.last().unwrap().t;
        let exact = 1.5 * t - 0.75 + 0.75 * (-2.0 * t).exp();
        assert!((hist.last().unwrap().x[0] - exact).abs() < 1e-8);
    }

    #[test]
    fn bad_step_and_divergence_are_handled() {
        let sys = StateSpace::<1, 1> {
            a: [[1.0]],
            b: [[0.0]],
        };
        assert!(simulate(&sys, [1.0], 0.0, 10, |_| [0.0]).is_empty());
        assert!(simulate(&sys, [1.0], f64::NAN, 10, |_| [0.0]).is_empty());
        // Violent growth overflows to infinity and the run is cut short.
        let hist = simulate(&sys, [1.0], 10.0, 10_000, |_| [0.0]);
        assert!(hist.len() < 10_001);
        assert!(hist.iter().all(|s| s.x[0].is_finite()));
    }
}
