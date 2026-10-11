//! Using the verify strategies outside a test suite: draw random wings and
//! angles from the proptest strategies, solve each one, and check the physical
//! invariants that the property tests enforce.
//!
//! Run with: `cargo run -p tpt-aero-verify --example property_sweep`

use proptest::strategy::{Strategy, ValueTree};
use proptest::test_runner::TestRunner;
use tpt_aero_verify::{angle_of_attack, symmetric_wing};
use tpt_aero_vlm::{CamberModel, VlmModel};

fn main() {
    let mut runner = TestRunner::deterministic();
    let wings = symmetric_wing();
    let alphas = angle_of_attack();

    let (mut checked, mut worst_roll, mut worst_yaw) = (0, 0.0_f64, 0.0_f64);
    for _ in 0..25 {
        let params = wings.new_tree(&mut runner).expect("sample").current();
        let alpha = alphas.new_tree(&mut runner).expect("sample").current();

        // Symmetric wing, no sideslip: roll and yaw moments must vanish.
        let model = VlmModel::from_planform(&params.planform(), 12, 3, CamberModel::None)
            .expect("strategy produces valid wings");
        let sol = model.solve(alpha, 40.0, 1.225, &[]).expect("solvable");
        assert!(sol.cl.is_finite(), "non-finite lift for {params:?}");

        worst_roll = worst_roll.max(sol.cl_roll.abs());
        worst_yaw = worst_yaw.max(sol.cn.abs());
        checked += 1;
    }

    println!("checked {checked} random symmetric wings at random angles of attack");
    println!("largest |roll moment| {worst_roll:.2e}, largest |yaw moment| {worst_yaw:.2e}");
    println!("(symmetric flight: both should be at round-off level)");
}
