//! Physical-invariant property tests (spec.txt §2): these encode the laws the
//! library must never violate, over randomized valid inputs.

use proptest::prelude::*;
use tpt_aero_verify::{angle_of_attack, freestream_speed, lattice_resolution, symmetric_wing};
use tpt_aero_vlm::{CamberModel, VlmModel};

const RHO: f64 = 1.225;
const V_INF: f64 = 50.0;

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
