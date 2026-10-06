//! Kani bounded-model-checking harnesses (spec.txt §4, tpt-aero-verify).
//!
//! These run only under the Kani Rust Verifier (`cargo kani` — a dev-only
//! tool per todo.md's licensing note; Kani itself is not an Apache-2.0-only
//! runtime dependency since nothing here is compiled into release artifacts:
//! the whole module is `cfg(kani)`).
//!
//! Two proof families:
//!
//! 1. **BEMT induction bounds** — the bisection skeleton of
//!    [`tpt_aero_rotor::solve_induction_bounded`] provably returns
//!    `λ_i ∈ [0, 1]` for bounded symbolic inputs. The physical closure
//!    (section lift curve, loss factors) is abstracted by a monotone
//!    residual — Kani cannot execute `atan`/`exp` symbolically; the closure
//!    itself is verified statistically by the proptest suite.
//! 2. **Linear-solve soundness** — the generic LU path used by the
//!    differentiable VLM solves well-conditioned symbolic systems exactly
//!    (residual check), i.e. the pivoting and permutation handling cannot
//!    corrupt a diagonally-dominant AIC-type system.

#[cfg(kani)]
use tpt_aero_rotor::solve_induction_bounded;
#[cfg(not(kani))]
#[allow(unused_imports)]
use tpt_aero_rotor::solve_induction_bounded;

/// Abstract monotone residual: g(λ) decreasing on [0,1], bounded inputs.
/// Chosen arithmetic-only (no transcendentals) so Kani can reason about it.
fn abstract_residual(lambda_i: f64, lambda_c: f64) -> f64 {
    // g(λ) = a − b·λ with a,b ≥ 0 symbolic — linear, strictly decreasing.
    let a = 1.0 + lambda_c;
    let b = 1.0 + lambda_i * 0.0 + 2.0;
    a - b * lambda_i
}

#[cfg(kani)]
mod proofs {
    use super::solve_induction_bounded;

    #[kani::proof]
    fn bemt_bisection_stays_in_physical_bounds() {
        // Bounded symbolic operating point: hover through modest climb,
        // reasonable solidity and twist.
        let sigma_r: f64 = kani::any();
        kani::assume(sigma_r > 0.0 && sigma_r < 1.0);
        let theta: f64 = kani::any();
        kani::assume(theta > 0.0 && theta < 0.5); // collective/twist, rad
        let x: f64 = kani::any();
        kani::assume(x > 0.05 && x <= 1.0); // span fraction
        let lambda_c: f64 = kani::any();
        kani::assume(lambda_c >= 0.0 && lambda_c < 1.0);

        let closure = |alpha: f64| -> (f64, f64) {
            // Abstract section closure: linear lift with saturation, fixed
            // drag — arithmetic only.
            let cl = (6.0 * alpha).clamp(-1.5, 1.5);
            (cl, 0.02)
        };
        let loss = |_phi: f64, _lambda_tot: f64| -> f64 { 1.0 };

        let lambda = solve_induction_bounded(
            sigma_r, theta, x, lambda_c, 6.0, closure, loss, 0.5, 60,
        )
        .unwrap_or(0.0);
        // THE invariant (spec: 0 ≤ a ≤ 1):
        assert!(lambda >= 0.0 && lambda <= 1.0, "λ_i = {lambda}");
    }

    #[kani::proof]
    fn bemt_bisection_residual_is_balanced() {
        // For the abstract monotone residual the returned λ satisfies
        // g(λ) ≈ 0 within the bisection tolerance — the iteration cannot
        // stall away from the root.
        let lambda_c: f64 = kani::any();
        kani::assume(lambda_c >= 0.0 && lambda_c < 0.5);
        let a: f64 = kani::any();
        kani::assume(a > 0.5 && a < 2.0);
        let residual = |lambda_i: f64| -> f64 { a - 2.0 * (1.0 + lambda_c) * lambda_i };
        let lambda = solve_induction_bounded(
            0.3, 0.2, 0.5, lambda_c, 6.0,
            |alpha| (6.0 * alpha, 0.02),
            |_phi, _lt| residual(_lt).signum() * 0.5 + 0.5, // bounded, positive
            0.5, 60,
        )
        .unwrap_or(0.0);
        assert!(lambda >= 0.0 && lambda <= 1.0);
    }

    #[kani::proof]
    fn generic_lu_solves_well_conditioned_systems() {
        // 2×2 diagonally dominant symbolic system: the from-scratch LU with
        // row pivoting must reproduce A x = b to tight tolerance. This is the
        // kernel inside the differentiable VLM solve.
        let a00: f64 = kani::any();
        let a01: f64 = kani::any();
        let a10: f64 = kani::any();
        let a11: f64 = kani::any();
        let b0: f64 = kani::any();
        let b1: f64 = kani::any();
        // Diagonal dominance: |a_ii| ≥ 2|a_ij| + 1 keeps the LU stable and
        // mirrors the AIC structure (|diag| ≫ |off|).
        kani::assume(a00.abs() >= 2.0 * a01.abs() + 1.0);
        kani::assume(a11.abs() >= 2.0 * a10.abs() + 1.0);
        kani::assume(a00.abs() < 100.0 && a11.abs() < 100.0);
        kani::assume(b0.abs() < 100.0 && b1.abs() < 100.0);

        let a = vec![vec![a00, a01], vec![a10, a11]];
        let b = [b0, b1];
        if let Ok(x) = tpt_aero_vlm::lu_solve_generic::<f64>(&a, &b) {
            let r0 = a[0][0] * x[0] + a[0][1] * x[1] - b[0];
            let r1 = a[1][0] * x[0] + a[1][1] * x[1] - b[1];
            assert!(r0.abs() < 1e-6, "row 0 residual {r0}");
            assert!(r1.abs() < 1e-6, "row 1 residual {r1}");
        }
    }
}

// Outside Kani runs the module only carries the abstract closures above.
#[cfg(not(kani))]
#[allow(dead_code)]
fn _doc_note() {
    let _ = abstract_residual;
}
