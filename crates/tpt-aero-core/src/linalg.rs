//! Native linear solves.
//!
//! * [`lu_solve`] — dense partial-pivot LU via the in-house
//!   `tpt-math-linalg-dense` backend (no external solver crates).
//! * [`lu_solve_generic`] — from-scratch partial-pivot LU generic over
//!   [`Real`], used by the differentiable assembly path (pivoting compares
//!   primal values only, which is exactly right for AD).
//! * [`cholesky_spd`] — from-scratch Cholesky–Banachiewicz for symmetric
//!   positive-definite systems (e.g. the Trefftz-plane induced-drag Gram
//!   matrix used in twist optimization).

use alloc::{vec, vec::Vec};

use crate::num::Real;

/// Errors from linear solves.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SolveError {
    /// Matrix is not square or dimensions mismatch.
    DimensionMismatch,
    /// A pivot fell below tolerance — the matrix is numerically singular.
    Singular,
    /// Cholesky applied to a matrix that is not symmetric positive definite.
    NotPositiveDefinite,
}

impl core::fmt::Display for SolveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            SolveError::DimensionMismatch => write!(f, "dimension mismatch"),
            SolveError::Singular => write!(f, "numerically singular matrix"),
            SolveError::NotPositiveDefinite => {
                write!(f, "matrix is not symmetric positive definite")
            }
        }
    }
}

impl core::error::Error for SolveError {}

/// Solve `A x = b` on `f64` via the native partial-pivot LU in
/// `tpt-math-linalg-dense`.
///
/// # Errors
/// [`SolveError::Singular`] if the matrix is numerically singular;
/// [`SolveError::DimensionMismatch`] on ragged input.
pub fn lu_solve(a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SolveError> {
    let n = b.len();
    if a.len() != n || a.iter().any(|row| row.len() != n) {
        return Err(SolveError::DimensionMismatch);
    }
    // One pivot test and one NaN rule for every solve: the same factorization
    // the prepared (cached) path uses, so the two cannot disagree.
    if !a.iter().flatten().chain(b).all(|v| v.is_finite()) {
        return Err(SolveError::Singular);
    }
    LuFactors::factor(a)?.solve(b)
}

/// From-scratch partial-pivot LU solve generic over [`Real`].
///
/// # Errors
/// [`SolveError::Singular`] when a pivot magnitude drops below `eps * max|a|`.
pub fn lu_solve_generic<R: Real>(a: &[Vec<R>], b: &[R]) -> Result<Vec<R>, SolveError> {
    let n = b.len();
    if n == 0 || a.len() != n || a.iter().any(|row| row.len() != n) {
        return Err(SolveError::DimensionMismatch);
    }

    // Partial-pivot LU in place: rows of `m` and entries of `x` are swapped
    // together, so after elimination + back substitution `x` is the solution
    // in the original variable order (PA = LU with Pb applied implicitly).
    let mut m: Vec<Vec<R>> = a.to_vec();
    let mut x: Vec<R> = b.to_vec();

    // A NaN entry would otherwise pass every comparison below and propagate.
    if !a.iter().flatten().chain(b).all(|v| v.to_f64().is_finite()) {
        return Err(SolveError::Singular);
    }

    // Global magnitude scale for the singularity tolerance.
    let mut scale = R::zero();
    for row in a {
        for &v in row {
            scale = scale.max(v.abs());
        }
    }
    let eps = R::from_f64(1.0e-12) * scale.max(R::one());

    for k in 0..n {
        // Pivot: largest |m[i][k]| for i ≥ k (primal comparison only).
        let mut p = k;
        let mut best = m[k][k].abs();
        for i in (k + 1)..n {
            let c = m[i][k].abs();
            if best.lt_f(c) {
                best = c;
                p = i;
            }
        }
        // `!(… > …)` rather than `<=` so a NaN pivot is also rejected.
        let (best_f, eps_f) = (best.to_f64(), eps.to_f64());
        if best_f.is_nan() || best_f <= eps_f {
            return Err(SolveError::Singular);
        }
        if p != k {
            m.swap(k, p);
            x.swap(k, p);
        }
        for i in (k + 1)..n {
            let f = m[i][k] / m[k][k];
            m[i][k] = R::zero();
            for j in (k + 1)..n {
                m[i][j] = m[i][j] - f * m[k][j];
            }
            x[i] = x[i] - f * x[k];
        }
    }
    // Back substitution.
    for i in (0..n).rev() {
        let mut s = x[i];
        for j in (i + 1)..n {
            s = s - m[i][j] * x[j];
        }
        x[i] = s / m[i][i];
    }
    // Undo the row permutation.
    Ok(x)
}

/// Solve `A x = b` for symmetric positive-definite `A` via the
/// Cholesky–Banachiewicz factorization `A = L Lᵀ`, entirely from scratch.
///
/// # Errors
/// [`SolveError::NotPositiveDefinite`] if a diagonal term stays ≤ 0;
/// [`SolveError::DimensionMismatch`] on ragged input.
pub fn cholesky_spd(a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, SolveError> {
    let n = b.len();
    if n == 0 || a.len() != n || a.iter().any(|row| row.len() != n) {
        return Err(SolveError::DimensionMismatch);
    }
    // NaN or infinite input would otherwise slip through the pivot test below.
    if !a.iter().flatten().chain(b).all(|v| v.is_finite()) {
        return Err(SolveError::NotPositiveDefinite);
    }
    let mut l = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut s = a[i][j];
            for k in 0..j {
                s -= l[i][k] * l[j][k];
            }
            if i == j {
                // NaN fails `s > 0` but is not caught by `s <= 0`, so test it too.
                if s.is_nan() || s <= 0.0 {
                    return Err(SolveError::NotPositiveDefinite);
                }
                l[i][j] = s.sqrt();
            } else {
                l[i][j] = s / l[j][j];
            }
        }
    }
    // Forward solve L y = b.
    let mut y = vec![0.0_f64; n];
    for i in 0..n {
        let mut s = b[i];
        for k in 0..i {
            s -= l[i][k] * y[k];
        }
        y[i] = s / l[i][i];
    }
    // Backward solve Lᵀ x = y.
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut s = y[i];
        for k in (i + 1)..n {
            s -= l[k][i] * x[k];
        }
        x[i] = s / l[i][i];
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dense(n: usize, f: impl Fn(usize, usize) -> f64) -> Vec<Vec<f64>> {
        (0..n).map(|i| (0..n).map(|j| f(i, j)).collect()).collect()
    }

    #[test]
    fn lu_matches_dense_backend() {
        // Diagonally dominant, well-conditioned.
        let a = dense(3, |i, j| if i == j { 6.0 } else { 1.0 + (i + j) as f64 });
        let b = [1.0, -2.0, 3.5];
        let x = lu_solve(&a, &b).expect("solvable");
        // Verify A x = b.
        for i in 0..3 {
            let r: f64 = (0..3).map(|j| a[i][j] * x[j]).sum();
            assert!((r - b[i]).abs() < 1e-9);
        }
    }

    #[test]
    fn generic_lu_matches_f64() {
        let a = dense(4, |i, j| {
            if i == j {
                8.0
            } else {
                0.5 * (1.0 + (i + j) as f64)
            }
        });
        let b = [2.0, -1.0, 0.5, 4.0];
        let x1 = lu_solve(&a, &b).expect("f64");
        let x2 = lu_solve_generic::<f64>(&a, &b).expect("generic");
        for (p, q) in x1.iter().zip(x2.iter()) {
            assert!((p - q).abs() < 1e-10);
        }
    }

    #[test]
    fn generic_lu_handles_row_pivoting() {
        // Zero diagonal forces heavy pivoting; the permutation undo must
        // reconstruct the true solution order.
        let a = vec![
            vec![0.0, 1.0, 2.0],
            vec![3.0, 0.0, 4.0],
            vec![5.0, 6.0, 0.0],
        ];
        let b = [5.0, 11.0, 17.0];
        // Solve by hand: x = [1, 1, 2] (0+1+4=5, 3+0+8=11, 5+6+0=11? no)
        // Verify by residual instead of a hand-picked vector.
        let x = lu_solve_generic::<f64>(&a, &b).expect("solvable");
        for i in 0..3 {
            let r: f64 = (0..3).map(|j| a[i][j] * x[j]).sum();
            assert!((r - b[i]).abs() < 1e-9, "row {i}: {r} != {}", b[i]);
        }
    }

    #[test]
    fn generic_lu_detects_singular() {
        // Row 1 == Row 2 → singular.
        let a: Vec<Vec<f64>> = vec![
            vec![1.0, 0.0, 0.0],
            vec![1.0, 2.0, 3.0],
            vec![1.0, 2.0, 3.0],
        ];
        assert_eq!(lu_solve(&a, &[1.0, 1.0, 1.0]), Err(SolveError::Singular));
        assert_eq!(
            lu_solve_generic::<f64>(&a, &[1.0, 1.0, 1.0]),
            Err(SolveError::Singular)
        );
    }

    #[test]
    fn cached_and_generic_lu_share_singularity_threshold() {
        // The prepared (cached) path and the generic path must reject the same
        // matrices: pivot 1e-13 is below 1e-12 * max(|a|, 1), pivot 1e-11 is not.
        for (pivot, expect_ok) in [(1.0e-13, false), (1.0e-11, true)] {
            let a = vec![vec![pivot, 0.0], vec![0.0, 1.0]];
            let b = [1.0, 1.0];
            let cached = lu_solve(&a, &b);
            let generic = lu_solve_generic::<f64>(&a, &b);
            assert_eq!(cached.is_ok(), expect_ok, "cached, pivot {pivot}");
            assert_eq!(generic.is_ok(), expect_ok, "generic, pivot {pivot}");
        }
    }

    #[test]
    fn cholesky_solves_spd_system() {
        // A = L Lᵀ for a known L → guaranteed SPD.
        let l = [
            vec![2.0, 0.0, 0.0],
            vec![0.5, 1.5, 0.0],
            vec![1.0, -0.5, 1.0],
        ];
        let mut a = vec![vec![0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                a[i][j] = (0..3).map(|k| l[i][k] * l[j][k]).sum();
            }
        }
        let b = [4.0, -1.0, 2.0];
        let x = cholesky_spd(&a, &b).expect("spd");
        for i in 0..3 {
            let r: f64 = (0..3).map(|j| a[i][j] * x[j]).sum();
            assert!((r - b[i]).abs() < 1e-9);
        }
    }

    #[test]
    fn cholesky_rejects_indefinite() {
        let a = vec![vec![1.0, 0.0], vec![0.0, -1.0]];
        assert_eq!(
            cholesky_spd(&a, &[1.0, 1.0]),
            Err(SolveError::NotPositiveDefinite)
        );
    }
}

/// A reusable partial-pivot LU factorization of a dense `f64` matrix.
///
/// Factor once with [`LuFactors::factor`] (O(n^3)), then solve any number of
/// right-hand sides with [`LuFactors::solve`] (O(n^2) each): the right tool
/// for sweeps where only the right-hand side changes.
#[derive(Debug, Clone)]
pub struct LuFactors {
    n: usize,
    /// Row-major combined factors: unit-lower L below the diagonal, U above.
    lu: Vec<f64>,
    /// `piv[k]` is the row swapped with row `k` at elimination step `k`.
    piv: Vec<usize>,
}

impl LuFactors {
    /// Factor the square matrix `a`.
    ///
    /// # Errors
    /// [`SolveError::DimensionMismatch`] if `a` is empty or not square;
    /// [`SolveError::Singular`] when a pivot magnitude drops below
    /// `1e-12 * max(max|a|, 1)`.
    pub fn factor(a: &[Vec<f64>]) -> Result<Self, SolveError> {
        let n = a.len();
        if n == 0 || a.iter().any(|row| row.len() != n) {
            return Err(SolveError::DimensionMismatch);
        }
        let mut lu = Vec::with_capacity(n * n);
        let mut scale = 0.0_f64;
        for row in a {
            for &v in row {
                scale = scale.max(v.abs());
            }
            lu.extend_from_slice(row);
        }
        let eps = 1.0e-12 * scale.max(1.0);
        let mut piv = vec![0; n];
        for k in 0..n {
            let mut p = k;
            let mut best = lu[k * n + k].abs();
            for i in (k + 1)..n {
                let c = lu[i * n + k].abs();
                if c > best {
                    best = c;
                    p = i;
                }
            }
            if best.is_nan() || best <= eps {
                return Err(SolveError::Singular);
            }
            piv[k] = p;
            if p != k {
                for j in 0..n {
                    lu.swap(k * n + j, p * n + j);
                }
            }
            let pivot = lu[k * n + k];
            for i in (k + 1)..n {
                let factor = lu[i * n + k] / pivot;
                lu[i * n + k] = factor;
                for j in (k + 1)..n {
                    lu[i * n + j] -= factor * lu[k * n + j];
                }
            }
        }
        Ok(Self { n, lu, piv })
    }

    /// Solve `A x = b` using the stored factors.
    ///
    /// # Errors
    /// [`SolveError::DimensionMismatch`] if `b.len()` differs from the matrix size.
    pub fn solve(&self, b: &[f64]) -> Result<Vec<f64>, SolveError> {
        let n = self.n;
        if b.len() != n {
            return Err(SolveError::DimensionMismatch);
        }
        let mut x = b.to_vec();
        for k in 0..n {
            x.swap(k, self.piv[k]);
        }
        for i in 1..n {
            let mut acc = x[i];
            for j in 0..i {
                acc -= self.lu[i * n + j] * x[j];
            }
            x[i] = acc;
        }
        for i in (0..n).rev() {
            let mut acc = x[i];
            for j in (i + 1)..n {
                acc -= self.lu[i * n + j] * x[j];
            }
            x[i] = acc / self.lu[i * n + i];
        }
        Ok(x)
    }
}

#[cfg(test)]
mod lu_factor_tests {
    use super::*;

    #[test]
    fn matches_one_shot_solver_for_many_rhs() {
        let a = vec![
            vec![0.0, 2.0, 1.0],
            vec![1.0, -1.0, 3.0],
            vec![4.0, 0.5, -2.0],
        ];
        let lu = LuFactors::factor(&a).expect("nonsingular");
        for b in [[1.0, 2.0, 3.0], [-4.0, 0.0, 7.5], [0.0, 0.0, 1.0]] {
            let x = lu.solve(&b).expect("sized");
            let reference = lu_solve(&a, &b).expect("solvable");
            for (u, v) in x.iter().zip(&reference) {
                assert!((u - v).abs() < 1e-12, "{u} vs {v}");
            }
        }
    }

    #[test]
    fn rejects_singular_and_ragged() {
        let singular = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        assert_eq!(
            LuFactors::factor(&singular).unwrap_err(),
            SolveError::Singular
        );
        let ragged = vec![vec![1.0, 2.0], vec![3.0]];
        assert_eq!(
            LuFactors::factor(&ragged).unwrap_err(),
            SolveError::DimensionMismatch
        );
    }
}

#[cfg(test)]
mod nan_tests {
    use super::*;

    #[test]
    fn cholesky_rejects_nan_input() {
        let a = vec![vec![4.0, f64::NAN], vec![f64::NAN, 2.0]];
        assert_eq!(
            cholesky_spd(&a, &[1.0, 1.0]),
            Err(SolveError::NotPositiveDefinite)
        );
        let good = vec![vec![4.0, 1.0], vec![1.0, 3.0]];
        assert!(cholesky_spd(&good, &[1.0, 2.0]).is_ok());
    }
}

#[cfg(test)]
mod generic_nan_tests {
    use super::*;

    #[test]
    fn generic_lu_rejects_nan_entries() {
        let a = vec![vec![2.0, f64::NAN], vec![1.0, 3.0]];
        assert_eq!(
            lu_solve_generic::<f64>(&a, &[1.0, 1.0]),
            Err(SolveError::Singular)
        );
        let good = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        assert!(lu_solve_generic::<f64>(&good, &[1.0, 1.0]).is_ok());
    }
}

#[cfg(test)]
mod single_criterion_tests {
    use super::*;

    #[test]
    fn one_shot_and_cached_solves_agree_on_singular_and_nan_input() {
        let singular = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        assert_eq!(lu_solve(&singular, &[1.0, 1.0]), Err(SolveError::Singular));
        assert!(LuFactors::factor(&singular).is_err());
        let nan = vec![vec![1.0, f64::NAN], vec![0.0, 1.0]];
        assert_eq!(lu_solve(&nan, &[1.0, 1.0]), Err(SolveError::Singular));
        let good = vec![vec![3.0, 1.0], vec![1.0, 2.0]];
        let x = lu_solve(&good, &[9.0, 8.0]).expect("solves");
        assert!((x[0] - 2.0).abs() < 1e-12 && (x[1] - 3.0).abs() < 1e-12);
    }
}
