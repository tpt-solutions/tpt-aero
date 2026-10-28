//! Viscous (parasitic) drag: flat-plate skin friction and a component
//! drag build-up for the zero-lift drag `C_D0`.
//!
//! The VLM and panel solvers are inviscid, so their drag is induced drag
//! only. Adding [`cd0_buildup`] to the induced drag gives a total drag
//! estimate for performance and trim work.

#[allow(unused_imports)]
use crate::num::Real;

/// Reynolds number below which the boundary layer is treated as laminar when
/// fixing transition on a flat plate.
pub const TRANSITION_RE: f64 = 5.0e5;

/// Laminar flat-plate skin-friction coefficient (Blasius), `C_f = 1.328/√Re`.
///
/// Returns NaN for a Reynolds number that is not finite and positive.
#[must_use]
pub fn laminar_flat_plate_cf(re: f64) -> f64 {
    if !(re.is_finite() && re > 0.0) {
        return f64::NAN;
    }
    1.328 / re.sqrt()
}

/// Flat-plate skin-friction coefficient with laminar leading edge and
/// turbulent trailing flow (Schlichting, transition at `Re = 5·10⁵`):
/// `C_f = 0.074 Re^{-1/5} − 1700/Re` for `Re ≥ 5·10⁵`, Blasius below.
///
/// `re` is the Reynolds number based on the plate length (or chord) in the
/// flow direction. Returns NaN unless `re` is finite and positive.
#[must_use]
pub fn flat_plate_cf(re: f64) -> f64 {
    if !(re.is_finite() && re > 0.0) {
        return f64::NAN;
    }
    if re < TRANSITION_RE {
        laminar_flat_plate_cf(re)
    } else {
        0.074 * re.powf(-0.2) - 1700.0 / re
    }
}

/// Van Driest II compressible skin friction on a flat plate, from the
/// incompressible Coles–Fernholz relation.
///
/// The compressible coefficient is `C_f = C_f,i / F_C`, with the incompressible
/// coefficient at the transformed Reynolds number `F_θ·Re_θ`:
/// `√(2/C_f,i) = ln(F_θ·Re_θ)/κ + C` with `κ = 0.384` and `C = 4.127`. The
/// transformation factors are (Zhu et al., arXiv:2609.02650, eqs. 2.2, A.8 and
/// the Walz relation):
///
/// * `F_C = (T_r/T_∞ − 1) / (arcsin α₀ + arcsin β₀)²`, with
///   `α₀ = (2A₀² − B₀)/√(B₀² + 4A₀²)`, `β₀ = B₀/√(B₀² + 4A₀²)`,
///   `A₀² = (T_r − T_∞)/T_w`, `B₀ = T_r/T_w − 1`;
/// * `F_θ = μ_∞/μ_w`, with Sutherland viscosity (`S = 110.4 K`, `T_∞ = 288.15 K`);
/// * `T_r = T_∞ [1 + r(γ−1)M²/2]`, `r = Pr^{1/3}`, `Pr = 0.72`, `γ = 1.4`.
///
/// Arguments: `re_theta` is the momentum-thickness Reynolds number of the
/// compressible boundary layer, `mach` the freestream Mach number, and
/// `wall_temperature_ratio` the wall temperature over the freestream
/// temperature (`T_w/T_∞`). For an adiabatic wall use `T_r/T_∞`.
///
/// The paper reports a mean error of 3.67 % against its DNS database for this
/// form (Mach 0.30 to 13.6, `T_w/T_r` from 0.18 to 1.89). That check is not
/// reproduced here. The incompressible relation is stated for
/// `Re_θ,i ≳ 5000`, and the function does not check that.
///
/// Returns NaN unless every argument is finite and positive. `mach = 0` is
/// excluded: both the numerator and denominator of `F_C` vanish there.
#[must_use]
pub fn van_driest_ii_cf(re_theta: f64, mach: f64, wall_temperature_ratio: f64) -> f64 {
    const PRANDTL: f64 = 0.72;
    const GAMMA: f64 = 1.4;
    const SUTHERLAND_OVER_TINF: f64 = 110.4 / 288.15;
    const KAPPA: f64 = 0.384;
    const C_LOG: f64 = 4.127;

    let valid = |x: f64| x.is_finite() && x > 0.0;
    if !(valid(re_theta) && valid(mach) && valid(wall_temperature_ratio)) {
        return f64::NAN;
    }
    let recovery = PRANDTL.powf(1.0 / 3.0);
    // Recovery and wall temperatures as ratios to the freestream temperature.
    let t_r = 1.0 + recovery * (GAMMA - 1.0) * mach * mach / 2.0;
    let t_w = wall_temperature_ratio;
    let a2 = (t_r - 1.0) / t_w;
    let b = t_r / t_w - 1.0;
    let root = (b * b + 4.0 * a2).sqrt();
    let alpha = (2.0 * a2 - b) / root;
    let beta = b / root;
    let f_c = (t_r - 1.0) / (alpha.asin() + beta.asin()).powi(2);
    // Sutherland: μ_w/μ_∞ = (T_w/T_∞)^{3/2} (1 + S')/(T_w/T_∞ + S'), S' = S/T_∞.
    let mu_w = t_w.powf(1.5) * (1.0 + SUTHERLAND_OVER_TINF) / (t_w + SUTHERLAND_OVER_TINF);
    let f_theta = 1.0 / mu_w;
    let root_cf = (f_theta * re_theta).ln() / KAPPA + C_LOG;
    let cf_incompressible = 2.0 / (root_cf * root_cf);
    cf_incompressible / f_c
}

/// Hoerner's form factor for a streamlined body or lifting surface of
/// thickness ratio `t/c`: `FF = 1 + 2(t/c) + 60(t/c)⁴`.
///
/// Multiplies the flat-plate friction to account for pressure drag due to
/// thickness. Returns NaN for a thickness ratio that is negative or not
/// finite.
#[must_use]
pub fn form_factor(thickness_ratio: f64) -> f64 {
    let t = thickness_ratio;
    if !(t.is_finite() && t >= 0.0) {
        return f64::NAN;
    }
    1.0 + 2.0 * t + 60.0 * t.powf(4.0)
}

/// One component of a drag build-up.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DragComponent {
    /// Reynolds number based on the component's reference length.
    pub reynolds: f64,
    /// Wetted area divided by the reference area of the aircraft (`S_wet/S`).
    pub wetted_ratio: f64,
    /// Thickness ratio used for the form factor (`t/c`; `0.0` for a body
    /// whose form factor is supplied directly via `form_override`).
    pub thickness_ratio: f64,
    /// Overrides the Hoerner form factor when set (e.g. for a fuselage).
    pub form_override: Option<f64>,
    /// Interference factor (≥ 1; 1.0 for an isolated component).
    pub interference: f64,
}

impl DragComponent {
    /// A lifting surface (wing or tail) with its own thickness ratio.
    #[must_use]
    pub const fn lifting_surface(reynolds: f64, wetted_ratio: f64, thickness_ratio: f64) -> Self {
        Self {
            reynolds,
            wetted_ratio,
            thickness_ratio,
            form_override: None,
            interference: 1.0,
        }
    }

    /// Contribution to `C_D0` referenced to the aircraft reference area.
    #[must_use]
    pub fn cd0(&self) -> f64 {
        let ff = self
            .form_override
            .unwrap_or_else(|| form_factor(self.thickness_ratio));
        flat_plate_cf(self.reynolds) * ff * self.interference * self.wetted_ratio
    }
}

/// Zero-lift drag coefficient `C_D0` as the sum of component contributions,
/// each `C_f · FF · Q · S_wet/S` (`Q` = interference factor).
#[must_use]
pub fn cd0_buildup(components: &[DragComponent]) -> f64 {
    components.iter().map(DragComponent::cd0).sum()
}

/// Total drag coefficient: `C_D = C_D0 + C_Di`, where `C_Di` comes from the
/// inviscid solver (e.g. the VLM induced drag).
#[must_use]
pub fn total_drag(cd0: f64, cd_induced: f64) -> f64 {
    cd0 + cd_induced
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn laminar_blasius_reference() {
        // Re = 10⁶: Cf = 1.328/1000 = 0.001328.
        assert!((laminar_flat_plate_cf(1.0e6) - 0.001_328).abs() < 1e-9);
    }

    #[test]
    fn transition_is_continuous_enough_and_turbulent_is_larger() {
        let below = flat_plate_cf(TRANSITION_RE * 0.999);
        let above = flat_plate_cf(TRANSITION_RE * 1.001);
        // Schlichting's mixed law jumps at the switch; the jump is bounded.
        assert!((above - below).abs() < 6e-4, "jump {}", above - below);
        // At Re = 10⁷ the turbulent value is well above the laminar one.
        assert!(flat_plate_cf(1.0e7) > laminar_flat_plate_cf(1.0e7));
        // Turbulent value at Re = 10⁶ is ≈ 0.00297.
        assert!((flat_plate_cf(1.0e6) - 0.002_97).abs() < 5e-5);
    }

    #[test]
    fn form_factor_reference_values() {
        assert!((form_factor(0.0) - 1.0).abs() < 1e-12);
        // t/c = 0.12: 1 + 0.24 + 60·0.12⁴ ≈ 1.2524.
        assert!((form_factor(0.12) - 1.252_4).abs() < 1e-4);
    }

    #[test]
    fn buildup_sums_components() {
        let wing = DragComponent::lifting_surface(3.4e6, 2.0, 0.12);
        let tail = DragComponent {
            wetted_ratio: 0.4,
            ..DragComponent::lifting_surface(1.0e6, 0.0, 0.09)
        };
        let total = cd0_buildup(&[wing, tail]);
        assert!((total - (wing.cd0() + tail.cd0())).abs() < 1e-12);
        // A reasonable light-aircraft-scale C_D0 (order 0.02–0.03).
        assert!(total > 0.005 && total < 0.05, "C_D0 = {total}");
    }

    #[test]
    fn form_override_replaces_hoerner() {
        let body = DragComponent {
            form_override: Some(1.5),
            ..DragComponent::lifting_surface(5.0e6, 1.0, 0.0)
        };
        assert!((body.cd0() - flat_plate_cf(5.0e6) * 1.5).abs() < 1e-12);
    }
}

#[cfg(test)]
mod validation_tests {
    use super::*;

    #[test]
    fn invalid_inputs_give_nan_not_plausible_values() {
        assert!(flat_plate_cf(0.0).is_nan());
        assert!(flat_plate_cf(-1.0e6).is_nan());
        assert!(flat_plate_cf(f64::NAN).is_nan());
        assert!(laminar_flat_plate_cf(0.0).is_nan());
        assert!(form_factor(-0.1).is_nan());
        assert!(form_factor(f64::INFINITY).is_nan());
        assert!(flat_plate_cf(1.0e6).is_finite());
    }

    fn incompressible_cf(re_theta: f64) -> f64 {
        let root = re_theta.ln() / 0.384 + 4.127;
        2.0 / (root * root)
    }

    #[test]
    fn van_driest_ii_reduces_to_the_incompressible_relation_at_low_mach() {
        // M → 0 with an isothermal wall: F_C = F_θ = 1.
        let cf = van_driest_ii_cf(1.0e4, 1.0e-3, 1.0);
        assert!((cf / incompressible_cf(1.0e4) - 1.0).abs() < 1e-4, "{cf}");
    }

    #[test]
    fn van_driest_ii_adiabatic_factor_matches_the_closed_form() {
        // Adiabatic wall (B₀ = 0): α₀ = A₀, β₀ = 0, so
        // F_C = (T_r/T_∞ − 1)/arcsin²(√(1 − T_∞/T_r)).
        let mach: f64 = 3.0;
        let t_r = 1.0 + 0.72_f64.powf(1.0 / 3.0) * 0.4 * mach * mach / 2.0;
        let f_c = (t_r - 1.0) / (1.0 - 1.0 / t_r).sqrt().asin().powi(2);
        let cf = van_driest_ii_cf(1.0e4, mach, t_r);
        let f_theta_adiabatic =
            1.0 / (t_r.powf(1.5) * (1.0 + 110.4 / 288.15) / (t_r + 110.4 / 288.15));
        let expected = incompressible_cf(f_theta_adiabatic * 1.0e4) / f_c;
        assert!((cf / expected - 1.0).abs() < 1e-9, "{cf} vs {expected}");
    }

    #[test]
    fn van_driest_ii_mach_3_adiabatic_regression() {
        // Value from the closed form above: C_f/C_f,i(Re_θ = 10⁴) ≈ 0.506 at M = 3.
        let t_r = 1.0 + 0.72_f64.powf(1.0 / 3.0) * 0.4 * 9.0 / 2.0;
        let cf = van_driest_ii_cf(1.0e4, 3.0, t_r);
        assert!((cf - 0.00146).abs() < 2.0e-5, "{cf}");
    }

    #[test]
    fn van_driest_ii_cooled_wall_has_more_friction_than_adiabatic() {
        // At M = 3 friction falls as the wall warms toward the recovery temperature.
        let adiabatic = 1.0 + 0.72_f64.powf(1.0 / 3.0) * 0.4 * 9.0 / 2.0;
        let cooled = van_driest_ii_cf(1.0e4, 3.0, 0.5);
        let mid = van_driest_ii_cf(1.0e4, 3.0, 1.0);
        let hot = van_driest_ii_cf(1.0e4, 3.0, adiabatic);
        assert!(cooled > mid && mid > hot, "{cooled} {mid} {hot}");
    }

    #[test]
    fn van_driest_ii_rejects_invalid_inputs() {
        assert!(van_driest_ii_cf(1.0e4, 0.0, 1.0).is_nan());
        assert!(van_driest_ii_cf(0.0, 3.0, 1.0).is_nan());
        assert!(van_driest_ii_cf(1.0e4, 3.0, -1.0).is_nan());
        assert!(van_driest_ii_cf(f64::NAN, 3.0, 1.0).is_nan());
    }
}
