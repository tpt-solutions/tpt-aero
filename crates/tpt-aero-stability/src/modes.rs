//! Eigenvalues of small plant matrices, for flight-dynamics modes.
//!
//! A 4×4 real matrix is reduced to its characteristic polynomial
//! (Faddeev–LeVerrier) whose roots come from the Durand–Kerner iteration:
//! no external linear-algebra dependency, accurate to roughly 1e-10 relative
//! for the well-scaled plant matrices used here.

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// A complex number (real, imaginary).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Complex {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}

impl Complex {
    /// Construct from parts.
    #[must_use]
    pub const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    /// Modulus.
    #[must_use]
    pub fn abs(self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }

    fn mul(self, o: Self) -> Self {
        Self::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }

    fn sub(self, o: Self) -> Self {
        Self::new(self.re - o.re, self.im - o.im)
    }

    fn div(self, o: Self) -> Self {
        let d = o.re * o.re + o.im * o.im;
        Self::new(
            (self.re * o.re + self.im * o.im) / d,
            (self.im * o.re - self.re * o.im) / d,
        )
    }
}

/// Coefficients `[c4, c3, c2, c1, c0]` (highest power first, `c4 = 1`) of
/// `det(λI − A)` for a 4×4 matrix, by the Faddeev–LeVerrier recursion.
#[must_use]
pub fn characteristic_polynomial(a: &[[f64; 4]; 4]) -> [f64; 5] {
    let mul = |x: &[[f64; 4]; 4], y: &[[f64; 4]; 4]| {
        let mut out = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                out[i][j] = (0..4).map(|k| x[i][k] * y[k][j]).sum();
            }
        }
        out
    };
    let mut coeffs = [0.0; 5];
    coeffs[0] = 1.0;
    let mut m = [[0.0; 4]; 4];
    for k in 1..=4 {
        // M_k = A M_{k-1} + c_{k-1} I
        let am = mul(a, &m);
        m = am;
        for (i, row) in m.iter_mut().enumerate() {
            row[i] += coeffs[k - 1];
        }
        let am = mul(a, &m);
        let trace: f64 = (0..4).map(|i| am[i][i]).sum();
        coeffs[k] = -trace / k as f64;
    }
    coeffs
}

/// Roots of the monic quartic with coefficients `[1, c3, c2, c1, c0]`.
fn quartic_roots(c: &[f64; 5]) -> [Complex; 4] {
    // Cauchy bound on the root magnitudes, used to scale the starting circle.
    let bound = 1.0 + c[1..].iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    let eval = |z: Complex| {
        let mut acc = Complex::new(1.0, 0.0);
        for coeff in &c[1..] {
            acc = acc.mul(z);
            acc.re += coeff;
        }
        acc
    };
    let seed = Complex::new(0.4, 0.9);
    let mut roots = [Complex::default(); 4];
    let mut power = Complex::new(1.0, 0.0);
    for r in &mut roots {
        *r = Complex::new(power.re * 0.5 * bound, power.im * 0.5 * bound);
        power = power.mul(seed);
    }
    for _ in 0..2000 {
        let mut largest_step = 0.0_f64;
        for i in 0..4 {
            let mut denom = Complex::new(1.0, 0.0);
            for j in 0..4 {
                if i != j {
                    denom = denom.mul(roots[i].sub(roots[j]));
                }
            }
            if denom.abs() == 0.0 {
                roots[i].re += 1.0e-8;
                continue;
            }
            let step = eval(roots[i]).div(denom);
            roots[i] = roots[i].sub(step);
            largest_step = largest_step.max(step.abs());
        }
        if largest_step < 1.0e-14 * bound {
            break;
        }
    }
    roots
}

/// Eigenvalues of a 4×4 real matrix, sorted by descending real part, with
/// conjugate pairs adjacent (positive imaginary part first). Results within
/// `1e-9 · scale` of the real axis are snapped onto it.
#[must_use]
pub fn eigenvalues4(a: &[[f64; 4]; 4]) -> [Complex; 4] {
    let poly = characteristic_polynomial(a);
    let mut roots = quartic_roots(&poly);
    let scale = roots.iter().fold(1.0_f64, |m, r| m.max(r.abs()));
    for r in &mut roots {
        if r.im.abs() < 1.0e-9 * scale {
            r.im = 0.0;
        }
    }
    roots.sort_by(|x, y| {
        y.re.partial_cmp(&x.re)
            .unwrap_or(core::cmp::Ordering::Equal)
            .then(
                y.im.partial_cmp(&x.im)
                    .unwrap_or(core::cmp::Ordering::Equal),
            )
    });
    roots
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec::Vec;

    fn sorted_re(roots: &[Complex; 4]) -> [f64; 4] {
        let mut v = roots.map(|r| r.re);
        v.sort_by(|a, b| b.partial_cmp(a).unwrap());
        v
    }

    #[test]
    fn diagonal_matrix_returns_its_diagonal() {
        let a = [
            [-1.0, 0.0, 0.0, 0.0],
            [0.0, -2.0, 0.0, 0.0],
            [0.0, 0.0, -30.0, 0.0],
            [0.0, 0.0, 0.0, 0.05],
        ];
        let r = eigenvalues4(&a);
        let want = [0.05, -1.0, -2.0, -30.0];
        for (got, want) in sorted_re(&r).iter().zip(want) {
            assert!((got - want).abs() < 1e-9, "{got} vs {want}");
        }
        assert!(r.iter().all(|z| z.im == 0.0));
    }

    #[test]
    fn complex_pair_from_a_companion_matrix() {
        // (s^2 + 0.4 s + 4)(s + 1)(s + 20): roots -0.2 +- i 1.99, -1, -20.
        let poly = {
            let a = [1.0, 0.4, 4.0];
            let b = [1.0, 21.0, 20.0];
            let mut p = [0.0; 5];
            for (i, x) in a.iter().enumerate() {
                for (j, y) in b.iter().enumerate() {
                    p[i + j] += x * y;
                }
            }
            p
        };
        // Companion matrix of the monic quartic.
        let a = [
            [-poly[1], -poly[2], -poly[3], -poly[4]],
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ];
        let r = eigenvalues4(&a);
        let pair: Vec<_> = r.iter().filter(|z| z.im != 0.0).collect();
        assert_eq!(pair.len(), 2);
        let expected_im = (4.0_f64 - 0.04).sqrt();
        for z in pair {
            assert!((z.re + 0.2).abs() < 1e-9, "re {}", z.re);
            assert!((z.im.abs() - expected_im).abs() < 1e-9, "im {}", z.im);
        }
        let reals = sorted_re(&r);
        assert!((reals[2] + 1.0).abs() < 1e-9 && (reals[3] + 20.0).abs() < 1e-9);
    }

    #[test]
    fn characteristic_polynomial_matches_trace_and_determinant() {
        let a = [
            [1.0, 2.0, 0.0, -1.0],
            [0.5, -3.0, 1.0, 0.0],
            [0.0, 4.0, 2.0, 1.0],
            [2.0, 0.0, -1.0, 0.5],
        ];
        let p = characteristic_polynomial(&a);
        let trace: f64 = (0..4).map(|i| a[i][i]).sum();
        assert!((p[1] + trace).abs() < 1e-12, "c3 = -trace");
        // det(A) = c0 for a quartic (even degree).
        let det = {
            // Cofactor expansion along the first row.
            let minor = |skip: usize| {
                let cols: Vec<usize> = (0..4).filter(|&c| c != skip).collect();
                let m = |r: usize, c: usize| a[r + 1][cols[c]];
                m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1))
                    - m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0))
                    + m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0))
            };
            a[0][0] * minor(0) - a[0][1] * minor(1) + a[0][2] * minor(2) - a[0][3] * minor(3)
        };
        assert!((p[4] - det).abs() < 1e-9, "c0 {} vs det {det}", p[4]);
    }
}
