//! Physical-invariant property tests (spec.txt §2): these encode the laws the
//! library must never violate, over randomized valid inputs.

use proptest::prelude::*;
use tpt_aero_verify::{angle_of_attack, freestream_speed, lattice_resolution, symmetric_wing};
use tpt_aero_vlm::{CamberModel, VlmModel};

const RHO: f64 = 1.225;
const V_INF: f64 = 50.0;

/// Generate a symmetric closed contour (ellipse-like) with random axes.
fn symmetric_contour(a: f64, b: f64, n: usize) -> Vec<(f64, f64)> {
    (0..n)
        .map(|i| {
            let t = 2.0 * core::f64::consts::PI * i as f64 / n as f64;
            (a * t.cos(), b * t.sin())
        })
        .collect()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// spec.txt §2: symmetric geometries at zero AoA produce zero rolling and
    /// yawing moments — over randomized symmetric planforms and lattices.
    #[test]
    fn symmetric_wing_zero_alpha_has_zero_roll_yaw(
        params in symmetric_wing(),
        (n_span, n_chord) in lattice_resolution(),
    ) {
        let wing = params.planform();
        let model = VlmModel::from_planform(&wing, n_span, n_chord, CamberModel::None)
            .expect("valid lattice");
        let sol = model.solve(0.0, V_INF, RHO, &[]).expect("solvable");
        prop_assert!(sol.cl.abs() < 1e-9, "cl = {}", sol.cl);
        prop_assert!(sol.cl_roll.abs() < 1e-9, "cl_roll = {}", sol.cl_roll);
        prop_assert!(sol.cn.abs() < 1e-9, "cn = {}", sol.cn);
        prop_assert!(sol.cdi.abs() < 1e-12, "cdi = {}", sol.cdi);
    }

    /// spec.txt §2: reversing the freestream velocity vector reverses the
    /// lift vector (exact for the linear, fixed-wake VLM model). The moment
    /// vector is invariant under reversal; only the coefficient projections
    /// (lift, drag, roll/yaw via their reference directions) transform.
    #[test]
    fn freestream_reversal_reverses_lift(
        params in symmetric_wing(),
        alpha in angle_of_attack(),
        (n_span, n_chord) in lattice_resolution(),
    ) {
        let wing = params.planform();
        let model = VlmModel::from_planform(&wing, n_span, n_chord, CamberModel::None)
            .expect("valid lattice");
        let fwd = model.solve(alpha, V_INF, RHO, &[]).expect("solvable");
        let rev = model
            .solve(alpha + core::f64::consts::PI, V_INF, RHO, &[])
            .expect("solvable");
        let scale = fwd.cl.abs().max(1.0);
        prop_assert!(
            (rev.cl + fwd.cl).abs() < 1e-9 * scale,
            "cl fwd {} rev {}",
            fwd.cl,
            rev.cl
        );
        // Drag coefficient flips (it is the force component along the
        // reversed freestream direction).
        prop_assert!(
            (rev.cdi + fwd.cdi).abs() < 1e-9 * fwd.cdi.abs().max(1e-6),
            "cdi fwd {} rev {}",
            fwd.cdi,
            rev.cdi
        );
        // The dimensional force vector is invariant under reversal.
        prop_assert!(
            (rev.force.x - fwd.force.x).abs() < 1e-6 * V_INF,
            "force x fwd {} rev {}",
            fwd.force.x,
            rev.force.x
        );
    }

    /// Lift must scale with dynamic pressure at fixed geometry and α.
    #[test]
    fn lift_scales_with_dynamic_pressure(
        params in symmetric_wing(),
        alpha in angle_of_attack(),
        speed in freestream_speed(),
        (n_span, n_chord) in lattice_resolution(),
    ) {
        let wing = params.planform();
        let model = VlmModel::from_planform(&wing, n_span, n_chord, CamberModel::None)
            .expect("valid lattice");
        let fast = model.solve(alpha, speed, RHO, &[]).expect("solvable");
        let slow = model.solve(alpha, speed * 0.5, RHO, &[]).expect("solvable");
        prop_assert!(
            (fast.cl - slow.cl).abs() < 1e-9 * fast.cl.abs().max(1.0),
            "coefficients must be speed-invariant: {} vs {}",
            fast.cl,
            slow.cl
        );
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]

    /// spec.txt §4 invariant (panel methods): a symmetric body at zero
    /// incidence produces zero net vertical force — the pressure-integrated
    /// lift of the source-panel solution vanishes over randomized ellipses.
    #[test]
    fn panel_symmetric_body_zero_lift(
        a in 0.5_f64..3.0,
        b in 0.3_f64..1.5,
        n in 16_usize..48,
    ) {
        let pts = symmetric_contour(a, b, n);
        let sol = tpt_aero_panel::source::solve_source_panels(&pts, 0.0, 1.0)
            .expect("solvable");
        let mut lift = 0.0_f64;
        for (i, &cp) in sol.cp.iter().enumerate() {
            let (x1, y1) = pts[i];
            let (x2, y2) = pts[(i + 1) % n];
            let dx = x2 - x1;
            let dy = y2 - y1;
            let l = (dx * dx + dy * dy).sqrt();
            let (nx, _ny) = (dy / l, -dx / l);
            let _ = nx;
            // Vertical force = −Cp·n_y·l (y is the lift direction).
            lift += -cp * (-dx / l) * l;
        }
        prop_assert!(lift.abs() < 1e-6, "lift = {lift}");
    }
}
