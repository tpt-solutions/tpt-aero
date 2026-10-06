//! Numeric abstraction layer: one trait implemented by `f64` and by the
//! forward-mode dual numbers from `tpt-math-autodiff-fwd`.
//!
//! Solver kernels in tpt-aero are generic over [`Real`], so the same VLM or
//! BEMT code path runs on plain `f64` for analysis and on
//! `Dual<f64, N>` for gradient-based shape optimization (e.g. induced-drag
//! minimization via twist distribution). This is the "differentiable AIC
//! assembly" requirement from `spec.txt` §4 (tpt-aero-vlm).

use core::ops::{Add, Div, Mul, Sub};

use tpt_math_autodiff_fwd::Dual;

/// A floating-point-like scalar usable in solver kernels.
///
/// Implemented for [`f64`] directly and for `Dual<f64, N>` (any number of
/// derivative directions) via the `tpt-math` forward-mode AD type. Comparisons
/// on the primal value go through [`Real::lt_f`] etc.; methods that are
/// non-smooth in exact arithmetic (`abs`, `max`, `min`) follow the standard AD
/// convention of taking the derivative branch of the selected operand.
pub trait Real:
    Copy + 'static + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Div<Output = Self>
{
    /// Additive identity.
    fn zero() -> Self;
    /// Multiplicative identity.
    fn one() -> Self;
    /// Two.
    fn two() -> Self;
    /// One half.
    fn half() -> Self;
    /// π.
    fn pi() -> Self;

    /// Lift a primitive constant into the scalar type.
    fn from_f64(v: f64) -> Self;

    /// The primal value as `f64` (derivative information is dropped).
    fn to_f64(self) -> f64;

    /// Negation.
    fn neg(self) -> Self;

    /// Square root.
    fn sqrt(self) -> Self;
    /// Absolute value.
    fn abs(self) -> Self;
    /// Sine.
    fn sin(self) -> Self;
    /// Cosine.
    fn cos(self) -> Self;
    /// Sine and cosine together.
    fn sin_cos(self) -> (Self, Self);
    /// Tangent.
    fn tan(self) -> Self;
    /// Arc sine.
    fn asin(self) -> Self;
    /// Arc cosine.
    fn acos(self) -> Self;
    /// Four-quadrant arc tangent of `self / other`.
    fn atan2(self, other: Self) -> Self;
    /// Raise to an integer power.
    fn powi(self, n: i32) -> Self;
    /// Raise to a real power.
    fn powf(self, n: Self) -> Self;
    /// Natural logarithm.
    fn ln(self) -> Self;
    /// Natural exponential.
    fn exp(self) -> Self;
    /// Larger of the two values (primal comparison).
    fn max(self, other: Self) -> Self;
    /// Smaller of the two values (primal comparison).
    fn min(self, other: Self) -> Self;
    /// Multiplicative inverse.
    fn recip(self) -> Self {
        Self::one() / self
    }
    /// True if the primal value (and any derivatives) are finite.
    fn is_finite(self) -> bool;

    /// Primal comparison: `self.to_f64() < other.to_f64()`.
    fn lt_f(self, other: Self) -> bool {
        self.to_f64() < other.to_f64()
    }

    /// Primal comparison: `self.to_f64() <= other.to_f64()`.
    fn le_f(self, other: Self) -> bool {
        self.to_f64() <= other.to_f64()
    }
}

impl Real for f64 {
    #[inline]
    fn zero() -> Self {
        0.0
    }

    #[inline]
    fn one() -> Self {
        1.0
    }

    #[inline]
    fn two() -> Self {
        2.0
    }

    #[inline]
    fn half() -> Self {
        0.5
    }

    #[inline]
    fn pi() -> Self {
        core::f64::consts::PI
    }

    #[inline]
    fn from_f64(v: f64) -> Self {
        v
    }

    #[inline]
    fn to_f64(self) -> f64 {
        self
    }

    #[inline]
    fn neg(self) -> Self {
        -self
    }

    #[inline]
    fn sqrt(self) -> Self {
        tpt_math_numeric::Float::sqrt(self)
    }

    #[inline]
    fn abs(self) -> Self {
        tpt_math_numeric::Float::abs(self)
    }

    #[inline]
    fn sin(self) -> Self {
        tpt_math_numeric::Float::sin(self)
    }

    #[inline]
    fn cos(self) -> Self {
        tpt_math_numeric::Float::cos(self)
    }

    #[inline]
    fn sin_cos(self) -> (Self, Self) {
        (
            tpt_math_numeric::Float::sin(self),
            tpt_math_numeric::Float::cos(self),
        )
    }

    #[inline]
    fn tan(self) -> Self {
        tpt_math_numeric::Float::tan(self)
    }

    #[inline]
    fn asin(self) -> Self {
        tpt_math_numeric::Float::asin(self)
    }

    #[inline]
    fn acos(self) -> Self {
        tpt_math_numeric::Float::acos(self)
    }

    #[inline]
    fn atan2(self, other: Self) -> Self {
        tpt_math_numeric::Float::atan2(self, other)
    }

    #[inline]
    fn powi(self, n: i32) -> Self {
        tpt_math_numeric::Float::powi(self, n)
    }

    #[inline]
    fn powf(self, n: Self) -> Self {
        tpt_math_numeric::Float::powf(self, n)
    }

    #[inline]
    fn ln(self) -> Self {
        tpt_math_numeric::Float::ln(self)
    }

    #[inline]
    fn exp(self) -> Self {
        tpt_math_numeric::Float::exp(self)
    }

    #[inline]
    fn max(self, other: Self) -> Self {
        tpt_math_numeric::Float::max(self, other)
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        tpt_math_numeric::Float::min(self, other)
    }

    #[inline]
    fn is_finite(self) -> bool {
        f64::is_finite(self)
    }
}

impl<const N: usize> Real for Dual<f64, N> {
    #[inline]
    fn zero() -> Self {
        Dual::constant(0.0)
    }

    #[inline]
    fn one() -> Self {
        Dual::constant(1.0)
    }

    #[inline]
    fn two() -> Self {
        Dual::constant(2.0)
    }

    #[inline]
    fn half() -> Self {
        Dual::constant(0.5)
    }

    #[inline]
    fn pi() -> Self {
        Dual::constant(core::f64::consts::PI)
    }

    #[inline]
    fn from_f64(v: f64) -> Self {
        Dual::constant(v)
    }

    #[inline]
    fn to_f64(self) -> f64 {
        self.re()
    }

    #[inline]
    fn neg(self) -> Self {
        self.negate()
    }

    #[inline]
    fn sqrt(self) -> Self {
        let r = self.re().sqrt();
        // d(sqrt(x)) = du / (2 sqrt(x)); at x = 0 the derivative limit is +inf.
        let du = core::array::from_fn(|i| self.du(i) * 0.5 / r);
        Dual::new(r, du)
    }

    #[inline]
    fn abs(self) -> Self {
        if self.re() < 0.0 {
            self.negate()
        } else {
            self
        }
    }

    #[inline]
    fn sin(self) -> Self {
        Dual::sin(self)
    }

    #[inline]
    fn cos(self) -> Self {
        Dual::cos(self)
    }

    #[inline]
    fn sin_cos(self) -> (Self, Self) {
        (self.sin(), self.cos())
    }

    #[inline]
    fn tan(self) -> Self {
        self.sin() * self.cos().recip()
    }

    #[inline]
    fn asin(self) -> Self {
        let s = self.re();
        let c = (1.0 - s * s).sqrt();
        let du = core::array::from_fn(|i| self.du(i) / c);
        Dual::new(s.asin(), du)
    }

    #[inline]
    fn acos(self) -> Self {
        let s = self.re();
        let c = (1.0 - s * s).sqrt();
        let du = core::array::from_fn(|i| -self.du(i) / c);
        Dual::new(s.acos(), du)
    }

    #[inline]
    fn atan2(self, other: Self) -> Self {
        let (a, b) = (self.re(), other.re());
        let denom = a * a + b * b;
        let du = core::array::from_fn(|i| (b * self.du(i) - a * other.du(i)) / denom);
        Dual::new(a.atan2(b), du)
    }

    #[inline]
    fn powi(self, n: i32) -> Self {
        let p = self.re().powi(n);
        // d(x^n) = n x^(n-1) du
        let du = core::array::from_fn(|i| self.du(i) * (n as f64) * self.re().powi(n - 1));
        Dual::new(p, du)
    }

    #[inline]
    fn powf(self, n: Self) -> Self {
        let (x, a) = (self.re(), n.re());
        // d(x^a) = x^a (a du_x / x + ln(x) du_a)
        let p = x.powf(a);
        let du = core::array::from_fn(|i| p * (a * self.du(i) / x + x.ln() * n.du(i)));
        Dual::new(p, du)
    }

    #[inline]
    fn ln(self) -> Self {
        Dual::ln(self)
    }

    #[inline]
    fn exp(self) -> Self {
        Dual::exp(self)
    }

    #[inline]
    fn max(self, other: Self) -> Self {
        if self.re() >= other.re() {
            self
        } else {
            other
        }
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        if self.re() <= other.re() {
            self
        } else {
            other
        }
    }

    #[inline]
    fn is_finite(self) -> bool {
        self.re().is_finite() && self.deriv().iter().all(|d| d.is_finite())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f64_circle() {
        assert!((f64::pi() - 0.0_f64.atan2(-1.0)).abs() < 1e-12);
    }

    #[test]
    fn dual_sqrt_derivative() {
        let x = Dual::<f64>::variable(4.0, 0);
        let y: Dual<f64> = Real::sqrt(x);
        assert!((y.re() - 2.0).abs() < 1e-12);
        assert!((y.du(0) - 0.25).abs() < 1e-12);
    }

    #[test]
    fn dual_atan2_gradient() {
        // atan2(y, x) partials: x/(x^2+y^2), -y/(x^2+y^2)
        let y = Dual::<f64>::variable(3.0, 0);
        let x = Dual::<f64, 1>::constant(4.0);
        let a = y.atan2(x);
        assert!((a.re() - 3.0_f64.atan2(4.0)).abs() < 1e-12);
        assert!((a.du(0) - 4.0 / 25.0).abs() < 1e-12);
    }

    #[test]
    fn dual_powi_gradient() {
        let x = Dual::<f64>::variable(2.0, 0);
        let y = x.powi(3);
        assert!((y.re() - 8.0).abs() < 1e-12);
        assert!((y.du(0) - 12.0).abs() < 1e-12);
    }

    #[test]
    fn dual_abs_max_min() {
        let x = Dual::<f64>::variable(-2.0, 0);
        assert_eq!(x.abs().re(), 2.0);
        assert_eq!(x.abs().du(0), -1.0);
        let a = Dual::<f64, 1>::constant(1.0);
        let b = Dual::<f64, 1>::constant(3.0);
        assert_eq!(a.max(b).re(), 3.0);
        assert_eq!(a.min(b).re(), 1.0);
    }

    #[test]
    fn dual_vec3_kernel_runs() {
        // Vec3<R> works on duals: a smoke test for the generic kernel path.
        let a = crate::geometry::Vec3::<Dual<f64, 1>>::new(
            Dual::variable(1.0, 0),
            Dual::constant(0.0),
            Dual::constant(0.0),
        );
        let b = crate::geometry::Vec3::<Dual<f64, 1>>::new(
            Dual::constant(0.0),
            Dual::variable(1.0, 0),
            Dual::constant(0.0),
        );
        let c = a.cross(&b);
        assert!((c.z.re() - 1.0).abs() < 1e-12);
        // c_z = x_a*y_b - y_a*x_b; both x_a and y_b are seeded in direction
        // 0, so d(c_z)/d0 = y_b + x_a = 2.
        assert!((c.z.du(0) - 2.0).abs() < 1e-12);
    }
}
