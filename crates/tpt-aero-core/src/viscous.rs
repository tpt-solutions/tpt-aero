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
}
