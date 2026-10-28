//! Lift-limited (post-stall plateau) VLM solve: invariants only. There is no
//! reference lift-limited solution here, so the checks are that the limit is
//! met, that it never adds lift, and that it reduces to the linear solve
//! below the limit.

use tpt_aero_core::{AirfoilPolar, WingPlanform};
use tpt_aero_vlm::{CamberModel, VlmError, VlmModel};

const V_INF: f64 = 50.0;
const RHO: f64 = 1.225;
const DEG: f64 = core::f64::consts::PI / 180.0;

fn wing() -> VlmModel {
    // Aspect ratio 8, flat plate, 16 spanwise strips of 4 chordwise panels.
    let planform = WingPlanform::rectangular(8.0, 1.0);
    VlmModel::from_planform(&planform, 16, 4, CamberModel::None).unwrap()
}

#[test]
fn below_the_limit_the_result_is_the_linear_solution() {
    let model = wing();
    let linear = model.solve(4.0 * DEG, V_INF, RHO, &[]).unwrap();
    let limited = model
        .solve_stall_limited(4.0 * DEG, V_INF, RHO, &[], 10.0)
        .unwrap();
    assert!(
        (limited.cl - linear.cl).abs() < 1e-9,
        "{} vs {}",
        limited.cl,
        linear.cl
    );
    for (a, b) in limited.gamma.iter().zip(&linear.gamma) {
        assert!((a - b).abs() < 1e-9);
    }
}

#[test]
fn strips_above_the_limit_sit_on_the_limit() {
    let model = wing();
    let cl_max = 0.9;
    let linear = model.solve(14.0 * DEG, V_INF, RHO, &[]).unwrap();
    // The linear solution must actually exceed the limit, or the test is vacuous.
    assert!(linear.strip_cl.iter().any(|c| c.abs() > cl_max));
    let limited = model
        .solve_stall_limited(14.0 * DEG, V_INF, RHO, &[], cl_max)
        .unwrap();
    let peak = limited.strip_cl.iter().fold(0.0_f64, |m, c| m.max(c.abs()));
    assert!((peak - cl_max).abs() < 1e-6, "peak strip cl {peak}");
    for c in &limited.strip_cl {
        assert!(c.abs() <= cl_max + 1e-6, "strip above the limit: {c}");
    }
}

#[test]
fn the_limit_never_adds_lift() {
    let model = wing();
    for deg in [6.0, 10.0, 14.0, 18.0] {
        let linear = model.solve(deg * DEG, V_INF, RHO, &[]).unwrap();
        let limited = model
            .solve_stall_limited(deg * DEG, V_INF, RHO, &[], 0.9)
            .unwrap();
        assert!(
            limited.cl <= linear.cl + 1e-9,
            "α={deg}: {} > {}",
            limited.cl,
            linear.cl
        );
        for (l, lin) in limited.strip_cl.iter().zip(&linear.strip_cl) {
            // Each strip is capped on its lift-producing side only.
            assert!(l.abs() <= lin.abs() + 1e-6, "α={deg}: strip gained lift");
        }
    }
}

#[test]
fn total_lift_is_monotone_on_the_plateau() {
    let model = wing();
    let mut previous = f64::NEG_INFINITY;
    for deg in (8..=22).step_by(2) {
        let cl = model
            .solve_stall_limited(f64::from(deg) * DEG, V_INF, RHO, &[], 0.9)
            .unwrap()
            .cl;
        assert!(
            cl >= previous - 1e-6,
            "CL fell at α={deg}: {cl} < {previous}"
        );
        previous = cl;
    }
}

#[test]
fn rejects_a_non_positive_limit() {
    let model = wing();
    for bad in [0.0, -0.5, f64::NAN] {
        assert_eq!(
            model
                .solve_stall_limited(5.0 * DEG, V_INF, RHO, &[], bad)
                .unwrap_err(),
            VlmError::InvalidStallLimit
        );
    }
}

#[test]
fn per_strip_uniform_caps_match_the_scalar_solve() {
    let model = wing();
    let scalar = model
        .solve_stall_limited(14.0 * DEG, V_INF, RHO, &[], 0.9)
        .unwrap();
    let caps = vec![0.9; model.lattice().strip_count()];
    let per_strip = model
        .solve_stall_limited_per_strip(14.0 * DEG, V_INF, RHO, &[], &caps)
        .unwrap();
    assert!((per_strip.cl - scalar.cl).abs() < 1e-10);
}

#[test]
fn each_strip_respects_its_own_cap_from_a_polar() {
    // Flat-plate-like polar: linear to a peak of 1.0 at 10°, then falling.
    let alphas: Vec<f64> = (-20..=30).map(|d| f64::from(d).to_radians()).collect();
    let cl: Vec<f64> = alphas
        .iter()
        .map(|a| {
            let lin = 2.0 * core::f64::consts::PI * a;
            if lin > 1.0 {
                1.0 - (lin - 1.0)
            } else {
                lin.max(-1.0)
            }
        })
        .collect();
    let polar = AirfoilPolar::new(alphas, cl.clone(), cl.clone(), cl, 0.2, 1.0e6).unwrap();
    let peak = polar.cl_max();
    let model = wing();
    let n = model.lattice().strip_count();
    // Taper the cap across the span, from 60 % of the peak at one tip to the
    // full peak at the other.
    let caps: Vec<f64> = (0..n)
        .map(|s| peak * (0.6 + 0.4 * s as f64 / (n - 1) as f64))
        .collect();
    let linear = model.solve(16.0 * DEG, V_INF, RHO, &[]).unwrap();
    let limited = model
        .solve_stall_limited_per_strip(16.0 * DEG, V_INF, RHO, &[], &caps)
        .unwrap();
    // No strip may exceed its own cap.
    for s in 0..n {
        assert!(
            limited.strip_cl[s].abs() <= caps[s] + 1e-6,
            "strip {s} above its cap"
        );
    }
    // Strips that the linear solution over-loads are held on their caps, unless
    // their neighbours shed enough lift to unload them: then they sit below the
    // cap with no incidence change. Some strip must be held on its cap.
    let held = (0..n)
        .filter(|&s| linear.strip_cl[s].abs() > caps[s])
        .filter(|&s| (limited.strip_cl[s].abs() - caps[s]).abs() < 1e-6)
        .count();
    let over = (0..n)
        .filter(|&s| linear.strip_cl[s].abs() > caps[s])
        .count();
    assert!(
        over > 0,
        "the test must exercise at least one over-loaded strip"
    );
    assert!(held > 0, "no over-loaded strip was held on its cap");
}

#[test]
fn per_strip_rejects_wrong_length_and_bad_entries() {
    let model = wing();
    let n = model.lattice().strip_count();
    let short = vec![0.9; n - 1];
    assert_eq!(
        model
            .solve_stall_limited_per_strip(5.0 * DEG, V_INF, RHO, &[], &short)
            .unwrap_err(),
        VlmError::DimensionMismatch
    );
    let mut bad = vec![0.9; n];
    bad[3] = -1.0;
    assert_eq!(
        model
            .solve_stall_limited_per_strip(5.0 * DEG, V_INF, RHO, &[], &bad)
            .unwrap_err(),
        VlmError::InvalidStallLimit
    );
}
