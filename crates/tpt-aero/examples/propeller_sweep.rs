//! Sample propeller BEMT sweep: a small UAV propeller swept across its
//! advance-ratio range, printing the performance map and comparing hover
//! induced velocity against ideal actuator-disk momentum theory.

use tpt_aero::rotor::bemt::{BemtSolver, OperatingCondition};
use tpt_aero::rotor::blade::{RotorGeometry, SectionModel};
use tpt_aero::rotor::loss::{HubLoss, TipLoss};
use tpt_aero::rotor::performance_map;

fn main() {
    // 10-inch (0.254 m) UAV propeller: 2 blades, chord 20 mm, 25 deg root
    // twist tapering to 5 deg at the tip.
    let geom = RotorGeometry::rectangular(0.127, 0.0127, 2, 0.020, 30)
        .with_linear_twist(25.0_f64.to_radians(), 5.0_f64.to_radians());
    let solver = BemtSolver::new(
        geom,
        SectionModel::Linear {
            cl_slope: 2.0 * core::f64::consts::PI,
            cl_max: 1.2,
            cd0: 0.012,
            cd_cl2: 0.02,
        },
        TipLoss::Prandtl,
        HubLoss::Prandtl { x_hub: 0.1 },
    )
    .expect("valid geometry");

    let omega = 800.0_f64; // rad/s (~7600 rpm)
    let rho = 1.225_f64;

    println!("=== Propeller BEMT sweep ===");
    println!(
        "R = 0.127 m, 2 blades, c = 20 mm, {} rpm\n",
        omega * 60.0 / core::f64::consts::TAU
    );
    println!(
        "{:>6} {:>8} {:>8} {:>8} {:>8} {:>8}",
        "J", "C_T", "C_Q", "C_P", "eta", "FOM"
    );
    let speeds: Vec<f64> = (0..=16).map(|i| i as f64 * 5.0).collect();
    let map =
        performance_map(&solver, omega, rho, &speeds, 20.0_f64.to_radians()).expect("map solved");
    for p in &map {
        println!(
            "{:>6.2} {:>8.4} {:>8.5} {:>8.5} {:>8.3} {:>8.3}",
            p.j, p.ct, p.cq, p.cp, p.efficiency, p.figure_of_merit
        );
    }

    // Momentum theory cross-check at hover.
    let op = OperatingCondition::hover(omega, rho);
    let sol = solver.solve(&op).expect("hover solved");
    let vi_ideal = (sol.loads.ct / 2.0).sqrt() * omega * 0.127;
    let vi_mid = sol.sections[sol.sections.len() / 2].lambda_induced * omega * 0.127;
    println!("\nHover cross-check:");
    println!("  ideal actuator disk v_i  = {vi_ideal:.3} m/s");
    println!("  BEMT mid-span v_i        = {vi_mid:.3} m/s");
    println!(
        "  figure of merit          = {:.3}",
        sol.loads.figure_of_merit
    );

    // Dynamic inflow lag after a collective step.
    let mut inflow = tpt_aero::rotor::DynamicInflow::new(0.127, 0.0);
    let steady = vi_ideal / (omega * 0.127);
    let dt = 0.002;
    for step in 0..=500 {
        inflow.step(0.0, steady, omega * 0.127, dt);
        if step % 100 == 0 {
            println!(
                "  t = {:5.2} s: lambda_i = {:.4}",
                step as f64 * dt,
                inflow.lambda()
            );
        }
    }
}
