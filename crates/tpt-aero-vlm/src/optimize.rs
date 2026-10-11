//! Induced-drag optimisation over the spanwise twist distribution.
//!
//! The induced drag of a planar wing at fixed lift is least for the elliptic
//! loading, which needs a particular twist. This module searches a one-parameter
//! family: a linear washout from the root (0 twist) to the tips (`τ`), with the
//! tip twist as the design variable. For each `τ` the angle of attack is solved
//! so the wing carries the target lift coefficient, and the induced drag is
//! minimised over `τ` by golden-section search.
//!
//! One parameter is a deliberately narrow search: it finds the best linear
//! washout, not the globally optimal twist, which is a general shape. The
//! result says how far the best washout is from the elliptic bound.

use alloc::vec::Vec;

use tpt_aero_core::{WingPlanform, WingSection};

use crate::lattice::CamberModel;
use crate::solver::{VlmError, VlmModel};

/// The result of [`optimize_linear_twist`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TwistOptimum {
    /// Best tip twist (deg): nose-up at the tips relative to the root (negative is washout).
    pub tip_twist_deg: f64,
    /// Angle of attack that carries the target lift at that twist (deg).
    pub alpha_deg: f64,
    /// Induced drag coefficient at the optimum.
    pub cdi: f64,
    /// Induced drag of the untwisted wing at the same lift, for comparison.
    pub cdi_untwisted: f64,
    /// Elliptic-loading induced drag `C_L²/(π·AR)`. A reference, not a bound: the
    /// VLM's near-field result can sit a few per cent below it (Oswald e above 1).
    pub cdi_elliptic: f64,
}

/// The same planform with a linear washout: each section keeps its geometry and
/// gets the twist `tip_twist · |y| / half-span`, so the root (at `y = 0`) is
/// untwisted and the tips carry `tip_twist`.
fn washout(planform: &WingPlanform, tip_twist: f64) -> Result<WingPlanform, VlmError> {
    let secs = planform.sections();
    let half = 0.5 * (secs[secs.len() - 1].le.y - secs[0].le.y).abs();
    if half.is_nan() || half <= 0.0 {
        return Err(VlmError::Lattice(
            crate::lattice::LatticeError::InvalidPlanform,
        ));
    }
    let sections: Vec<WingSection> = secs
        .iter()
        .map(|s| WingSection {
            twist: tip_twist * (s.le.y.abs() / half),
            ..*s
        })
        .collect();
    WingPlanform::from_sections(sections)
        .map_err(|_| VlmError::Lattice(crate::lattice::LatticeError::InvalidPlanform))
}

/// `(C_Di, α)` at the target lift for one twist, with α found by secant on `C_L`.
fn induced_at(
    planform: &WingPlanform,
    n_span: usize,
    n_chord: usize,
    tip_twist: f64,
    cl_target: f64,
    speed: f64,
    density: f64,
) -> Result<(f64, f64), VlmError> {
    let p = washout(planform, tip_twist)?;
    let model = VlmModel::from_planform(&p, n_span, n_chord, CamberModel::None)?
        .with_reference(crate::solver::ReferenceGeometry::from_planform(planform));
    let cl = |a: f64| -> Result<f64, VlmError> { Ok(model.solve(a, speed, density, &[])?.cl) };
    let (mut a0, mut a1) = (0.0_f64, 0.05_f64);
    let (mut f0, mut f1) = (cl(a0)? - cl_target, cl(a1)? - cl_target);
    for _ in 0..30 {
        if (f1 - f0).abs() < 1.0e-14 || f1.abs() < 1.0e-10 {
            break;
        }
        let a2 = a1 - f1 * (a1 - a0) / (f1 - f0);
        a0 = a1;
        f0 = f1;
        a1 = a2;
        f1 = cl(a1)? - cl_target;
    }
    let sol = model.solve(a1, speed, density, &[])?;
    Ok((sol.cdi, a1))
}

/// Minimise the induced drag over the tip twist of a linear washout, at the
/// lift coefficient `cl_target`.
///
/// `planform` is a tapered wing (two end sections, symmetric about `y = 0`).
/// The search covers tip twists from −10° to +5° (negative is washout).
///
/// # Errors
/// [`VlmError::InvalidFlightCondition`] for a non-finite target or a non-positive
/// speed or density; the VLM's own errors otherwise.
pub fn optimize_linear_twist(
    planform: &WingPlanform,
    n_span: usize,
    n_chord: usize,
    cl_target: f64,
    speed: f64,
    density: f64,
) -> Result<TwistOptimum, VlmError> {
    if !(cl_target.is_finite()
        && speed.is_finite()
        && speed > 0.0
        && density.is_finite()
        && density > 0.0)
    {
        return Err(VlmError::InvalidFlightCondition);
    }
    let eval = |deg: f64| {
        induced_at(
            planform,
            n_span,
            n_chord,
            deg.to_radians(),
            cl_target,
            speed,
            density,
        )
    };

    // Golden-section search on the tip twist, in degrees.
    let (mut lo, mut hi) = (-10.0_f64, 5.0_f64);
    // Golden ratio conjugate, (√5 − 1)/2, as a literal (`sqrt` is not in `no_std`).
    let phi = 0.618_033_988_749_894_9_f64;
    let mut x1 = hi - phi * (hi - lo);
    let mut x2 = lo + phi * (hi - lo);
    let mut f1 = eval(x1)?.0;
    let mut f2 = eval(x2)?.0;
    for _ in 0..24 {
        if f1 < f2 {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = hi - phi * (hi - lo);
            f1 = eval(x1)?.0;
        } else {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = lo + phi * (hi - lo);
            f2 = eval(x2)?.0;
        }
        if hi - lo < 1.0e-4 {
            break;
        }
    }
    let best = 0.5 * (lo + hi);
    let (cdi, alpha) = eval(best)?;
    let (cdi_untwisted, _) = eval(0.0)?;
    let span = planform.span();
    let aspect = span * span / planform.area();
    Ok(TwistOptimum {
        tip_twist_deg: best,
        alpha_deg: alpha.to_degrees(),
        cdi,
        cdi_untwisted,
        cdi_elliptic: cl_target * cl_target / (core::f64::consts::PI * aspect),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn washout_keeps_the_planform_area_and_span() {
        let wing = WingPlanform::tapered(10.0, 1.2, 0.8, 0.0, 0.0);
        let w = washout(&wing, 3.0_f64.to_radians()).expect("valid");
        assert!((w.area() - wing.area()).abs() < 1e-9);
        assert!((w.span() - wing.span()).abs() < 1e-9);
    }

    #[test]
    fn twist_optimum_beats_the_untwisted_wing_and_respects_the_bound() {
        let wing = WingPlanform::tapered(10.0, 1.2, 0.6, 0.0, 0.0);
        let opt = optimize_linear_twist(&wing, 16, 3, 0.5, 50.0, 1.225).expect("optimises");
        assert!(opt.cdi <= opt.cdi_untwisted + 1e-9, "{opt:?}");
        // The VLM stays within 10 % of the elliptic reference. It can be a little
        // below it (Oswald e above 1, a known near-field effect), so this is a band.
        assert!((opt.cdi / opt.cdi_elliptic - 1.0).abs() < 0.10, "{opt:?}");
        // A tapered wing wants washout at the tips.
        assert!(opt.tip_twist_deg < 0.0, "tip twist {}", opt.tip_twist_deg);
    }
}
