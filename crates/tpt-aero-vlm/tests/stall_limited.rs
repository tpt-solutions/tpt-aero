//! Lift-limited (post-stall plateau) VLM solve: invariants only. There is no
//! reference lift-limited solution here, so the checks are that the limit is
//! met, that it never adds lift, and that it reduces to the linear solve
//! below the limit.

use tpt_aero_core::WingPlanform;
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
