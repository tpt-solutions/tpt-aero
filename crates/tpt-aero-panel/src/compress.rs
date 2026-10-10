//! Subsonic compressibility corrections for pressure coefficients.

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// Failure of a subsonic compressibility correction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompressibilityError {
    /// Mach number outside the subsonic range `0 ≤ M < 1` (or not finite).
    MachOutOfRange,
    /// The correction's denominator vanished or changed sign (Cp too far
    /// negative for the model; the flow is beyond the correction's validity).
    Singular,
}

impl core::fmt::Display for CompressibilityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MachOutOfRange => write!(f, "Mach number must satisfy 0 <= M < 1"),
            Self::Singular => write!(f, "compressibility correction is singular for this Cp"),
        }
    }
}

impl core::error::Error for CompressibilityError {}

/// Prandtl–Glauert correction: `Cp_c = Cp_i / β` with `β = √(1 − M²)`.
///
/// First-order subsonic scaling; valid while the flow stays subsonic
/// everywhere (M < ~0.7 for typical bodies).
///
/// # Errors
/// [`CompressibilityError::MachOutOfRange`] unless `0 ≤ M < 1`.
pub fn prandtl_glauert(cp_incompressible: f64, mach: f64) -> Result<f64, CompressibilityError> {
    Ok(cp_incompressible / beta(mach)?)
}

/// The Prandtl–Glauert factor β = √(1 − M²) (1 at M = 0).
///
/// # Errors
/// [`CompressibilityError::MachOutOfRange`] for M ≥ 1, M < 0 or non-finite M
/// (no real β); [`is_applicable`] is the equivalent gate.
pub fn beta(mach: f64) -> Result<f64, CompressibilityError> {
    if is_applicable(mach) {
        Ok((1.0 - mach * mach).sqrt())
    } else {
        Err(CompressibilityError::MachOutOfRange)
    }
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
///
/// # Errors
/// [`CompressibilityError::MachOutOfRange`] unless `0 ≤ M < 1`;
/// [`CompressibilityError::Singular`] if the denominator is not positive.
pub fn karman_tsien(cp_incompressible: f64, mach: f64) -> Result<f64, CompressibilityError> {
    let b = beta(mach)?;
    let denom = b + mach * mach / (2.0 * (1.0 + b)) * cp_incompressible;
    if denom.is_nan() || denom <= 0.0 {
        return Err(CompressibilityError::Singular);
    }
    Ok(cp_incompressible / denom)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pg_identity_at_zero_mach() {
        assert!((prandtl_glauert(-1.2, 0.0).unwrap() + 1.2).abs() < 1e-12);
        assert!((beta(0.0).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn pg_amplifies_suction() {
        let cp = prandtl_glauert(-1.0, 0.6).unwrap();
        let expected = -1.0 / (1.0_f64 - 0.36).sqrt();
        assert!((cp - expected).abs() < 1e-12);
        assert!(cp < -1.0);
    }

    #[test]
    fn kt_amplifies_suction_beyond_pg() {
        // For negative Cp the Kármán–Tsien denominator shrinks ⇒ stronger
        // suction than Prandtl–Glauert; for positive Cp it is weaker.
        let m = 0.6;
        let pg = prandtl_glauert(-1.0, m).unwrap();
        let kt = karman_tsien(-1.0, m).unwrap();
        assert!(kt < pg, "kt {kt} should be stronger suction than pg {pg}");
        let pg_p = prandtl_glauert(0.8, m).unwrap();
        let kt_p = karman_tsien(0.8, m).unwrap();
        assert!(kt_p < pg_p, "kt {kt_p} should be weaker than pg {pg_p}");
    }

    #[test]
    fn kt_reduces_to_cp_at_zero_mach() {
        assert!((karman_tsien(-0.7, 0.0).unwrap() + 0.7).abs() < 1e-12);
    }

    #[test]
    fn rejects_out_of_range_mach() {
        for m in [1.0, 1.1, -0.1, f64::NAN, f64::INFINITY] {
            assert_eq!(beta(m), Err(CompressibilityError::MachOutOfRange));
            assert!(prandtl_glauert(-1.0, m).is_err());
            assert!(karman_tsien(-1.0, m).is_err());
        }
    }

    #[test]
    fn karman_tsien_flags_singular_denominator() {
        // Very large suction drives the denominator through zero.
        assert_eq!(
            karman_tsien(-50.0, 0.9),
            Err(CompressibilityError::Singular)
        );
    }
}
