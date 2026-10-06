//! Subsonic compressibility corrections for pressure coefficients.

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// Prandtl–Glauert correction: `Cp_c = Cp_i / β` with `β = √(1 − M²)`.
///
/// First-order subsonic scaling; valid while the flow stays subsonic
/// everywhere (M < ~0.7 for typical bodies).
#[must_use]
pub fn prandtl_glauert(cp_incompressible: f64, mach: f64) -> f64 {
    let beta = beta(mach);
    cp_incompressible / beta
}

/// The Prandtl–Glauert factor β = √(1 − M²) (1 at M = 0).
///
/// # Panics
/// For Mach ≥ 1 (no real β) — the caller is expected to gate on
/// [`is_applicable`].
#[must_use]
pub fn beta(mach: f64) -> f64 {
    assert!((0.0..1.0).contains(&mach), "β requires 0 ≤ M < 1, got {mach}");
    (1.0 - mach * mach).sqrt()
}

/// True when the correction is applicable (0 ≤ M < 1).
#[must_use]
pub const fn is_applicable(mach: f64) -> bool {
    mach >= 0.0 && mach < 1.0
}

/// von Kármán–Tsien correction (second-order subsonic):
///
/// `Cp_c = Cp_i / (β + (M²/(2(1+β)))·Cp_i)`
///
/// Accounts for the finite thickness effect on the disturbance level; it
/// predicts weaker suction amplification than Prandtl–Glauert for negative
/// Cp and stronger compression for positive Cp.
#[must_use]
pub fn karman_tsien(cp_incompressible: f64, mach: f64) -> f64 {
    let b = beta(mach);
    cp_incompressible / (b + mach * mach / (2.0 * (1.0 + b)) * cp_incompressible)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pg_identity_at_zero_mach() {
        assert!((prandtl_glauert(-1.2, 0.0) + 1.2).abs() < 1e-12);
        assert!((beta(0.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn pg_amplifies_suction() {
        let cp = prandtl_glauert(-1.0, 0.6);
        let expected = -1.0 / (1.0_f64 - 0.36).sqrt();
        assert!((cp - expected).abs() < 1e-12);
        assert!(cp < -1.0);
    }

    #[test]
    fn kt_amplifies_suction_beyond_pg() {
        // For negative Cp the Kármán–Tsien denominator shrinks ⇒ stronger
        // suction than Prandtl–Glauert; for positive Cp it is weaker.
        let m = 0.6;
        let pg = prandtl_glauert(-1.0, m);
        let kt = karman_tsien(-1.0, m);
        assert!(kt < pg, "kt {kt} should be stronger suction than pg {pg}");
        let pg_p = prandtl_glauert(0.8, m);
        let kt_p = karman_tsien(0.8, m);
        assert!(kt_p < pg_p, "kt {kt_p} should be weaker than pg {pg_p}");
    }

    #[test]
    fn kt_reduces_to_cp_at_zero_mach() {
        assert!((karman_tsien(-0.7, 0.0) + 0.7).abs() < 1e-12);
    }

    #[test]
    #[should_panic]
    fn beta_rejects_supersonic() {
        let _ = beta(1.1);
    }
}
