//! Induced-drag optimisation over the spanwise twist distribution.
//!
//! The induced drag of a planar wing at fixed lift is least for the elliptic
//! loading, which needs a particular twist. This module searches a one-parameter
//! family: a linear washout from the root (0 twist) to the tips (`τ`), with the
//! tip twist as the design variable. For each `τ` the angle of attack is solved
//! so the wing carries the target lift coefficient, and the induced drag is
//! minimised over `τ` by golden-section search.
//!
//! [`optimize_linear_twist`] is the narrow search: the best linear washout, with
//! the tip twist as the only design variable. [`optimize_twist_shape`] widens it
//! to a cubic spanwise twist law. It is still a local search, not a global one,
//! but it contains the washout family, so it can only do as well or better. Both
//! results say how far the best twist is from the elliptic bound.

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

/// The same planform with a spanwise twist law: each section keeps its geometry
/// and gets the twist `shape(η)`, where `η = |y| / half-span` runs from 0 at the
/// root to 1 at the tips. Twist is in radians.
fn with_twist_shape(
    planform: &WingPlanform,
    shape: impl Fn(f64) -> f64,
) -> Result<WingPlanform, VlmError> {
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
            twist: shape(s.le.y.abs() / half),
            ..*s
        })
        .collect();
    WingPlanform::from_sections(sections)
        .map_err(|_| VlmError::Lattice(crate::lattice::LatticeError::InvalidPlanform))
}

/// `(C_Di, α)` at the target lift for one twist law, with α found by secant on `C_L`.
fn induced_at(
    planform: &WingPlanform,
    n_span: usize,
    n_chord: usize,
    shape: impl Fn(f64) -> f64,
    cl_target: f64,
    speed: f64,
    density: f64,
) -> Result<(f64, f64), VlmError> {
    let p = with_twist_shape(planform, shape)?;
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
        let tip = deg.to_radians();
        induced_at(
            planform,
            n_span,
            n_chord,
            |eta| tip * eta,
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

/// The result of [`optimize_twist_shape`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TwistShapeOptimum {
    /// Coefficients `c₁, c₂, c₃` (deg) of the twist law
    /// `twist(η) = c₁ η + c₂ η² + c₃ η³`, with `η = |y| / half-span`.
    pub coefficients_deg: [f64; 3],
    /// Angle of attack that carries the target lift at that twist law (deg).
    pub alpha_deg: f64,
    /// Induced drag coefficient at the optimum.
    pub cdi: f64,
    /// Induced drag of the untwisted wing at the same lift, for comparison.
    pub cdi_untwisted: f64,
    /// Elliptic-loading induced drag `C_L²/(π·AR)`. A reference, not a bound.
    pub cdi_elliptic: f64,
    /// VLM solves used by the search (each is one planform build plus a secant
    /// on the angle of attack).
    pub evaluations: usize,
}

/// Minimise the induced drag over a cubic spanwise twist law `c₁ η + c₂ η² +
/// c₃ η³` at the lift coefficient `cl_target`, with `η = |y| / half-span`.
///
/// The root stays untwisted (the law has no constant term). The search is a
/// compass pattern search from the untwisted wing: it steps each coefficient by
/// ±`step` (starting at 2°), keeps any improvement, and halves the step when no
/// direction improves, down to 0.01°. It is derivative-free and deterministic, but it finds
/// a local minimum. The linear washout of [`optimize_linear_twist`] is the slice
/// `c₂ = c₃ = 0`, so this search can only do as well or better on the same
/// planform, up to the search tolerance.
///
/// # Errors
/// [`VlmError::InvalidFlightCondition`] for a non-finite target or a non-positive
/// speed or density; the VLM's own errors otherwise.
pub fn optimize_twist_shape(
    planform: &WingPlanform,
    n_span: usize,
    n_chord: usize,
    cl_target: f64,
    speed: f64,
    density: f64,
) -> Result<TwistShapeOptimum, VlmError> {
    if !(cl_target.is_finite()
        && speed.is_finite()
        && speed > 0.0
        && density.is_finite()
        && density > 0.0)
    {
        return Err(VlmError::InvalidFlightCondition);
    }
    // Coefficients in degrees; the shape converts to radians for the planform.
    let eval = |c: [f64; 3]| {
        let c = c.map(f64::to_radians);
        induced_at(
            planform,
            n_span,
            n_chord,
            |eta| c[0] * eta + c[1] * eta * eta + c[2] * eta * eta * eta,
            cl_target,
            speed,
            density,
        )
    };

    const MAX_EVALUATIONS: usize = 400;
    const MIN_STEP_DEG: f64 = 0.01;
    let mut evaluations = 0_usize;
    let mut best = [0.0_f64; 3];
    let (mut best_cdi, mut best_alpha) = eval(best)?;
    evaluations += 1;
    let mut step = 2.0_f64;
    while step > MIN_STEP_DEG && evaluations < MAX_EVALUATIONS {
        let mut improved = false;
        for i in 0..3 {
            for sign in [1.0_f64, -1.0] {
                let mut trial = best;
                trial[i] += sign * step;
                let (cdi, alpha) = eval(trial)?;
                evaluations += 1;
                if cdi < best_cdi {
                    best = trial;
                    best_cdi = cdi;
                    best_alpha = alpha;
                    improved = true;
                }
            }
        }
        if !improved {
            step *= 0.5;
        }
    }
    let (cdi_untwisted, _) = induced_at(
        planform,
        n_span,
        n_chord,
        |_| 0.0,
        cl_target,
        speed,
        density,
    )?;
    let span = planform.span();
    let aspect = span * span / planform.area();
    Ok(TwistShapeOptimum {
        coefficients_deg: best,
        alpha_deg: best_alpha.to_degrees(),
        cdi: best_cdi,
        cdi_untwisted,
        cdi_elliptic: cl_target * cl_target / (core::f64::consts::PI * aspect),
        evaluations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn washout_keeps_the_planform_area_and_span() {
        let wing = WingPlanform::tapered(10.0, 1.2, 0.8, 0.0, 0.0);
        let tip = 3.0_f64.to_radians();
        let w = with_twist_shape(&wing, |eta| tip * eta).expect("valid");
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

    #[test]
    fn twist_shape_search_matches_or_beats_the_linear_washout() {
        let wing = WingPlanform::tapered(10.0, 1.2, 0.6, 0.0, 0.0);
        let linear = optimize_linear_twist(&wing, 16, 3, 0.5, 50.0, 1.225).expect("linear");
        let shaped = optimize_twist_shape(&wing, 16, 3, 0.5, 50.0, 1.225).expect("shaped");
        assert!(shaped.cdi <= shaped.cdi_untwisted + 1e-9, "{shaped:?}");
        // The linear washout is inside the cubic family, so the search should not
        // lose to it by more than a small fraction of the induced drag.
        assert!(
            shaped.cdi <= linear.cdi * 1.01,
            "shaped {shaped:?} linear {linear:?}"
        );
        // The root stays untwisted by construction: no constant term in the law.
        assert!(shaped.evaluations > 1 && shaped.evaluations <= 400);
    }
}
