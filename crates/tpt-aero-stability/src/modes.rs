//! Eigenvalues of small plant matrices, for flight-dynamics modes.
//!
//! A 4×4 real matrix is reduced to upper Hessenberg form and solved with the
//! Francis double-shift QR iteration: no external linear-algebra dependency.
//! Repeated and defective eigenvalues converge like any others, and a failed
//! iteration is reported as NaN rather than as a plausible-looking result.

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

/// Eigenvalues of a 4×4 real matrix by reduction to upper Hessenberg form
/// (Gaussian elimination with pivoting) followed by the Francis double-shift
/// QR iteration. Returns `None` for non-finite input or if the iteration
/// fails to deflate within its budget.
///
/// The arrays are 1-indexed (index 0 unused) to follow the textbook
/// formulation, which keeps the index arithmetic auditable.
fn francis_eigenvalues(a_in: &[[f64; 4]; 4]) -> Option<[Complex; 4]> {
    const N: usize = 4;
    if !a_in.iter().flatten().all(|x| x.is_finite()) {
        return None;
    }
    let mut a = [[0.0_f64; N + 1]; N + 1];
    for i in 0..N {
        for j in 0..N {
            a[i + 1][j + 1] = a_in[i][j];
        }
    }

    // Reduction to Hessenberg form by stabilised elimination.
    for m in 2..N {
        let mut x = 0.0;
        let mut piv = m;
        for j in m..=N {
            if a[j][m - 1].abs() > x.abs() {
                x = a[j][m - 1];
                piv = j;
            }
        }
        if piv != m {
            for j in (m - 1)..=N {
                let t = a[piv][j];
                a[piv][j] = a[m][j];
                a[m][j] = t;
            }
            for row in a.iter_mut().skip(1) {
                row.swap(piv, m);
            }
        }
        if x != 0.0 {
            for i in (m + 1)..=N {
                let mut y = a[i][m - 1];
                if y != 0.0 {
                    y /= x;
                    a[i][m - 1] = y;
                    for j in m..=N {
                        a[i][j] -= y * a[m][j];
                    }
                    for j in 1..=N {
                        a[j][m] += y * a[j][i];
                    }
                }
            }
        }
    }
    // The elimination leaves multipliers below the subdiagonal; clear them.
    for i in 3..=N {
        for j in 1..(i - 1) {
            a[i][j] = 0.0;
        }
    }

    let mut wr = [0.0_f64; N + 1];
    let mut wi = [0.0_f64; N + 1];
    let mut anorm = 0.0_f64;
    for i in 1..=N {
        for j in i.saturating_sub(1).max(1)..=N {
            anorm += a[i][j].abs();
        }
    }
    let sign = |a: f64, b: f64| if b >= 0.0 { a.abs() } else { -a.abs() };

    let mut nn = N;
    let mut t = 0.0_f64;
    while nn >= 1 {
        let mut its = 0;
        loop {
            // Look for a single small subdiagonal element.
            let mut l = nn;
            while l >= 2 {
                let mut s = a[l - 1][l - 1].abs() + a[l][l].abs();
                if s == 0.0 {
                    s = anorm;
                }
                if a[l][l - 1].abs() + s == s {
                    a[l][l - 1] = 0.0;
                    break;
                }
                l -= 1;
            }
            let x = a[nn][nn];
            if l == nn {
                // One root found.
                wr[nn] = x + t;
                wi[nn] = 0.0;
                nn -= 1;
                break;
            }
            let y = a[nn - 1][nn - 1];
            let w = a[nn][nn - 1] * a[nn - 1][nn];
            if l == nn - 1 {
                // Two roots found.
                let p = 0.5 * (y - x);
                let q = p * p + w;
                let z = q.abs().sqrt();
                let x = x + t;
                if q >= 0.0 {
                    let z = p + sign(z, p);
                    wr[nn - 1] = x + z;
                    wr[nn] = x + z;
                    if z != 0.0 {
                        wr[nn] = x - w / z;
                    }
                    wi[nn - 1] = 0.0;
                    wi[nn] = 0.0;
                } else {
                    wr[nn - 1] = x + p;
                    wr[nn] = x + p;
                    wi[nn - 1] = -z;
                    wi[nn] = z;
                }
                nn = nn.saturating_sub(2);
                break;
            }
            if its == 60 {
                return None;
            }
            // Exceptional shifts to break cycles.
            let (mut x, mut y, mut w) = (x, y, w);
            if its == 10 || its == 20 {
                t += x;
                for i in 1..=nn {
                    a[i][i] -= x;
                }
                let s = a[nn][nn - 1].abs() + a[nn - 1][nn - 2].abs();
                x = 0.75 * s;
                y = x;
                w = -0.4375 * s * s;
            }
            its += 1;

            // Look for two consecutive small subdiagonal elements.
            let mut m = nn - 2;
            let (mut p, mut q, mut r);
            loop {
                let z = a[m][m];
                r = x - z;
                let s = y - z;
                p = (r * s - w) / a[m + 1][m] + a[m][m + 1];
                q = a[m + 1][m + 1] - z - r - s;
                r = a[m + 2][m + 1];
                let s = p.abs() + q.abs() + r.abs();
                p /= s;
                q /= s;
                r /= s;
                if m == l {
                    break;
                }
                let u = a[m][m - 1].abs() * (q.abs() + r.abs());
                let v = p.abs() * (a[m - 1][m - 1].abs() + z.abs() + a[m + 1][m + 1].abs());
                if u + v == v {
                    break;
                }
                m -= 1;
            }
            for i in (m + 2)..=nn {
                a[i][i - 2] = 0.0;
                if i != m + 2 {
                    a[i][i - 3] = 0.0;
                }
            }

            // Double QR step on rows l..nn and columns m..nn.
            let mut k = m;
            while k < nn {
                if k != m {
                    p = a[k][k - 1];
                    q = a[k + 1][k - 1];
                    r = 0.0;
                    if k != nn - 1 {
                        r = a[k + 2][k - 1];
                    }
                    let xx = p.abs() + q.abs() + r.abs();
                    x = xx;
                    if xx != 0.0 {
                        p /= xx;
                        q /= xx;
                        r /= xx;
                    }
                }
                let s = sign((p * p + q * q + r * r).sqrt(), p);
                if s != 0.0 {
                    if k == m {
                        if l != m {
                            a[k][k - 1] = -a[k][k - 1];
                        }
                    } else {
                        a[k][k - 1] = -s * x;
                    }
                    p += s;
                    x = p / s;
                    y = q / s;
                    let z = r / s;
                    q /= p;
                    r /= p;
                    for j in k..=nn {
                        let mut pp = a[k][j] + q * a[k + 1][j];
                        if k != nn - 1 {
                            pp += r * a[k + 2][j];
                            a[k + 2][j] -= pp * z;
                        }
                        a[k + 1][j] -= pp * y;
                        a[k][j] -= pp * x;
                    }
                    let mmin = if nn < k + 3 { nn } else { k + 3 };
                    for i in l..=mmin {
                        let mut pp = x * a[i][k] + y * a[i][k + 1];
                        if k != nn - 1 {
                            pp += z * a[i][k + 2];
                            a[i][k + 2] -= pp * r;
                        }
                        a[i][k + 1] -= pp * q;
                        a[i][k] -= pp;
                    }
                }
                k += 1;
            }
        }
    }

    let mut out = [Complex::default(); N];
    for i in 0..N {
        if !(wr[i + 1].is_finite() && wi[i + 1].is_finite()) {
            return None;
        }
        out[i] = Complex::new(wr[i + 1], wi[i + 1]);
    }
    Some(out)
}

/// Eigenvalues of a 4×4 real matrix, sorted by descending real part, with
/// conjugate pairs adjacent (positive imaginary part first). Results within
/// `1e-9 · scale` of the real axis are snapped onto it.
///
/// Returns four NaN eigenvalues if the input is non-finite or the QR iteration
/// does not converge, so a failed computation cannot pass for a valid mode.
#[must_use]
pub fn eigenvalues4(a: &[[f64; 4]; 4]) -> [Complex; 4] {
    let Some(mut roots) = francis_eigenvalues(a) else {
        return [Complex::new(f64::NAN, f64::NAN); 4];
    };
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

    #[test]
    fn defective_repeated_eigenvalue_is_found() {
        // Jordan block: all four eigenvalues equal -2. Defective, so the
        // error scales like eps^(1/4); the bound is loose on purpose.
        let a = [
            [-2.0, 1.0, 0.0, 0.0],
            [0.0, -2.0, 1.0, 0.0],
            [0.0, 0.0, -2.0, 1.0],
            [0.0, 0.0, 0.0, -2.0],
        ];
        for z in eigenvalues4(&a) {
            assert!((z.re + 2.0).abs() < 1e-3 && z.im.abs() < 1e-3, "{z:?}");
        }
    }

    #[test]
    fn eigenvalues_match_trace_and_determinant() {
        let a = [
            [1.0, 2.0, 0.0, -1.0],
            [0.5, -3.0, 1.0, 0.0],
            [0.0, 4.0, 2.0, 1.0],
            [2.0, 0.0, -1.0, 0.5],
        ];
        let r = eigenvalues4(&a);
        let trace: f64 = (0..4).map(|i| a[i][i]).sum();
        let sum_re: f64 = r.iter().map(|z| z.re).sum();
        assert!(
            (sum_re - trace).abs() < 1e-9,
            "sum {sum_re} vs trace {trace}"
        );
        // Product of eigenvalues (complex) must equal det(A), which is real.
        let mut prod = Complex::new(1.0, 0.0);
        for z in r {
            let (re, im) = (
                prod.re * z.re - prod.im * z.im,
                prod.re * z.im + prod.im * z.re,
            );
            prod = Complex::new(re, im);
        }
        let p = characteristic_polynomial(&a);
        assert!(
            (prod.re - p[4]).abs() < 1e-8 && prod.im.abs() < 1e-8,
            "{prod:?} vs {}",
            p[4]
        );
    }

    #[test]
    fn non_finite_input_gives_nan_not_a_mode() {
        let mut a = [[0.0; 4]; 4];
        a[1][2] = f64::NAN;
        assert!(eigenvalues4(&a).iter().all(|z| z.re.is_nan()));
    }
}
