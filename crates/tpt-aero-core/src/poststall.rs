//! Post-stall lift and drag: the empirical AERODAS model (Spera, NASA/CR-2008-215434,
//! equations 8 to 12). It is fitted to airfoil test data and to full-scale
//! rotor power data.
//!
//! The functions cover the post-stall branches only. Below the lift stall angle
//! `ACL1` the caller uses the attached-flow curve, and below the drag stall
//! angle `ACD1` it uses the pre-stall drag curve. The two branches are not
//! blended here.
//!
//! The report's own fit to 231 post-stall lift points (NACA 0012, 4415 and
//! 4418) has a mean deviation of 0.012 and a standard deviation of 0.038. Its
//! S809 parameter table is reproduced in the tests of this module.
//!
//! Angles are in radians at the interface. The report's constants are in
//! degrees, and the conversion happens inside.

#[allow(unused_imports)]
use crate::num::Real;

/// Radians to degrees.
const DEG_PER_RAD: f64 = 180.0 / core::f64::consts::PI;

/// Post-stall lift: the maximum lift coefficient `CL2max` at a thickness ratio
/// `t/c` (a fraction, `0 ≤ t/c ≤ 1`) and aspect ratio `AR`.
///
/// `CL2max = F1·F2`, with `F1 = 1.190·(1 − (t/c)²)` and
/// `F2 = 0.65 + 0.35·exp(−(9/AR)^2.3)` (equations 8a, 9a, 9b). `F2` tends to 1
/// as `AR` grows, and the report's constants are fitted to the test data
/// above.
///
/// Returns NaN for a thickness ratio outside `[0, 1]` or a non-positive
/// `aspect_ratio`.
#[must_use]
pub fn cl_max(thickness_ratio: f64, aspect_ratio: f64) -> f64 {
    if !((0.0..=1.0).contains(&thickness_ratio) && aspect_ratio.is_finite() && aspect_ratio > 0.0) {
        return f64::NAN;
    }
    let f1 = 1.190 * (1.0 - thickness_ratio * thickness_ratio);
    let f2 = 0.65 + 0.35 * (-(9.0 / aspect_ratio).powf(2.3)).exp();
    f1 * f2
}

/// Post-stall drag: the maximum drag coefficient `CD2max` at 90° for a
/// thickness ratio `t/c` (a fraction) and aspect ratio `AR`.
///
/// `CD2max = G1·G2`, with `G1 = 2.300·exp(−(0.65·t/c)^0.9)` and
/// `G2 = 0.52 + 0.48·exp(−(6.5/AR)^1.1)` (equations 8b, 10b, 10c). `G1` is a
/// default. The report's S809 table gives `G1` as an input that differs from
/// this default by about 1 %, so a measured `G1` can replace it by using
/// `G1·G2` directly.
///
/// Returns NaN for a thickness ratio outside `[0, 1]` or a non-positive
/// `aspect_ratio`.
#[must_use]
pub fn cd_max(thickness_ratio: f64, aspect_ratio: f64) -> f64 {
    if !((0.0..=1.0).contains(&thickness_ratio) && aspect_ratio.is_finite() && aspect_ratio > 0.0) {
        return f64::NAN;
    }
    let g1 = 2.300 * (-(0.65 * thickness_ratio).powf(0.9)).exp();
    let g2 = 0.52 + 0.48 * (-(6.5 / aspect_ratio).powf(1.1)).exp();
    g1 * g2
}

/// Post-stall lift coefficient at `alpha` (rad), with `cl2_max` from
/// [`cl_max`] (equations 11a to 11e).
///
/// Below `alpha_stall` the post-stall lift is zero, as in equation 11a. From
/// the stall angle the lift follows a power law that falls from `cl2_max` at
/// 41° to zero at 92°, and turns negative beyond that. The report gives the
/// value at the stall angle as about 0.49 for the S809 parameters, which is
/// well below its pre-stall peak of 1.07. So the post-stall curve does not
/// join the attached curve continuously at `alpha_stall`.
///
/// Returns NaN for a non-finite angle, or for `cl2_max` outside `(0, 1.632)`:
/// the report's constant 1.632 bounds the exponent's denominator.
#[must_use]
pub fn post_stall_lift(alpha: f64, alpha_stall: f64, cl2_max: f64) -> f64 {
    let r = 1.632 - cl2_max;
    if !(alpha.is_finite()
        && alpha_stall.is_finite()
        && cl2_max.is_finite()
        && cl2_max > 0.0
        && r > 0.0)
    {
        return f64::NAN;
    }
    let a = alpha * DEG_PER_RAD;
    if alpha < alpha_stall {
        return 0.0;
    }
    let n = 1.0 + cl2_max / r;
    let linear = -0.032 * (a - 92.0);
    if a <= 92.0 {
        linear - r * ((92.0 - a) / 51.0).powf(n)
    } else {
        linear + r * ((a - 92.0) / 51.0).powf(n)
    }
}

/// Post-stall drag coefficient at `alpha` (rad) for `alpha_drag_stall ≤ alpha ≤
/// π/2`, from the pre-stall drag `cd1_max` at the drag stall angle and the
/// maximum `cd2_max` at 90° (equation 12b).
///
/// The curve is a sine that starts at `cd1_max` at `alpha_drag_stall` and peaks
/// at `cd2_max` at 90°, so it meets the pre-stall curve there.
///
/// Returns NaN outside `[alpha_drag_stall, π/2]` or for non-finite inputs.
#[must_use]
pub fn post_stall_drag(alpha: f64, alpha_drag_stall: f64, cd1_max: f64, cd2_max: f64) -> f64 {
    let lo = alpha_drag_stall * DEG_PER_RAD;
    let a = alpha * DEG_PER_RAD;
    if !(a.is_finite()
        && lo.is_finite()
        && lo < 90.0
        && cd1_max.is_finite()
        && cd2_max.is_finite()
        && a >= lo
        && a <= 90.0)
    {
        return f64::NAN;
    }
    // Both angles are in degrees, so their ratio is dimensionless.
    let phase = core::f64::consts::FRAC_PI_2 * (a - lo) / (90.0 - lo);
    cd1_max + (cd2_max - cd1_max) * phase.sin()
}

#[cfg(test)]
mod tests {
    use super::*;

    const RAD: f64 = core::f64::consts::PI / 180.0;

    #[test]
    fn s809_cl_max_reproduces_the_report_table() {
        // S809: t/c = 0.21, infinite aspect ratio, table CL2max = 1.138.
        let cl2 = cl_max(0.21, 1.0e4);
        assert!((cl2 - 1.138).abs() < 0.002, "CL2max {cl2}");
    }

    #[test]
    fn blade_cl_max_reproduces_the_report_table() {
        // Blade of the S809 rotor: AR 15.28, table CL2max = 1.036.
        let cl2 = cl_max(0.21, 15.28);
        assert!((cl2 - 1.036).abs() < 0.002, "CL2max {cl2}");
    }

    #[test]
    fn s809_lift_hits_its_anchors() {
        let cl2 = cl_max(0.21, 1.0e4);
        // The S809 table has RCL2 = 0.494 and N2 = 3.30.
        assert!((1.632 - cl2 - 0.494).abs() < 0.002);
        // Equation 11b at 41° is CL2max exactly, and at 92° it is zero.
        assert!((post_stall_lift(41.0 * RAD, 14.0 * RAD, cl2) - cl2).abs() < 1e-9);
        assert!(post_stall_lift(92.0 * RAD, 14.0 * RAD, cl2).abs() < 1e-12);
        // The report's value at the stall angle, about 0.49, sits well below
        // the pre-stall peak of 1.07.
        let at_stall = post_stall_lift(14.0 * RAD, 14.0 * RAD, cl2);
        assert!((at_stall - 0.49).abs() < 0.01, "CL2 at 14° = {at_stall}");
    }

    #[test]
    fn post_stall_lift_is_zero_before_the_stall_angle() {
        assert_eq!(post_stall_lift(10.0 * RAD, 14.0 * RAD, 1.1), 0.0);
    }

    #[test]
    fn post_stall_lift_falls_after_its_peak() {
        let cl2 = 1.138;
        let at_41 = post_stall_lift(41.0 * RAD, 14.0 * RAD, cl2);
        let at_70 = post_stall_lift(70.0 * RAD, 14.0 * RAD, cl2);
        let at_90 = post_stall_lift(90.0 * RAD, 14.0 * RAD, cl2);
        assert!(at_41 > at_70 && at_70 > at_90 && at_90 > 0.0);
    }

    #[test]
    fn cd_max_reproduces_the_blade_table_with_the_table_g1() {
        // Table G1 = 1.922 for the S809 blade (AR 15.28): CD2max = 1.624.
        let g2 = 0.52 + 0.48 * (-(6.5_f64 / 15.28).powf(1.1)).exp();
        assert!((1.922 * g2 - 1.624).abs() < 0.002, "{}", 1.922 * g2);
        // The default G1 lies within 2 % of that input.
        let default = cd_max(0.21, 15.28);
        assert!((default / (1.922 * g2) - 1.0).abs() < 0.02, "{default}");
    }

    #[test]
    fn post_stall_drag_meets_the_pre_stall_curve_and_peaks_at_ninety() {
        let acd1 = 20.1 * RAD;
        assert!((post_stall_drag(acd1, acd1, 0.2, 1.921) - 0.2).abs() < 1e-12);
        assert!((post_stall_drag(90.0 * RAD, acd1, 0.2, 1.921) - 1.921).abs() < 1e-12);
        let mut last = 0.0;
        for d in 21..=90 {
            let cd = post_stall_drag(f64::from(d) * RAD, acd1, 0.2, 1.921);
            assert!(cd >= last - 1e-12, "drag fell at {d}°");
            last = cd;
        }
    }

    #[test]
    fn invalid_inputs_are_nan() {
        assert!(cl_max(1.5, 10.0).is_nan());
        assert!(cl_max(0.12, 0.0).is_nan());
        assert!(cd_max(-0.1, 10.0).is_nan());
        assert!(post_stall_lift(0.5, 0.2, 2.0).is_nan());
        assert!(post_stall_drag(0.1, 0.5, 0.2, 1.9).is_nan());
    }
}
