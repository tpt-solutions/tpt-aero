//! Blade Element Momentum Theory for a small UAV propeller: one hover point,
//! then an advance-ratio performance map.
//!
//! Run with: `cargo run -p tpt-aero-rotor --example hover_bemt`

use tpt_aero_rotor::{
    performance_map, BemtSolver, HubLoss, OperatingCondition, RotorGeometry, SectionModel, TipLoss,
};

fn main() {
    // 10 inch propeller: 2 blades, 20 mm chord, root twist 25 deg to 5 deg at the tip.
    let geometry = RotorGeometry::rectangular(0.127, 0.0127, 2, 0.020, 30)
        .with_linear_twist(25.0_f64.to_radians(), 5.0_f64.to_radians());
    let solver = BemtSolver::new(
        geometry,
        SectionModel::Linear {
            cl_slope: 2.0 * std::f64::consts::PI,
            cl_max: 1.2,
            cd0: 0.012,
            cd_cl2: 0.02,
        },
        TipLoss::Prandtl,
        HubLoss::Prandtl { x_hub: 0.1 },
    )
    .expect("valid geometry");

    // Hover at 800 rad/s (about 7600 rpm).
    let omega = 800.0;
    let hover = OperatingCondition {
        axial_speed: 0.0,
        omega,
        density: 1.225,
        collective: 0.0,
    };
    let sol = solver.solve(&hover).expect("converges");
    println!(
        "hover: C_T {:.4}, C_P {:.5}, figure of merit {:.3}, {} iterations",
        sol.loads.ct, sol.loads.cp, sol.loads.figure_of_merit, sol.iterations
    );

    // Thrust in newtons: T = C_T rho n^2 D^4, with n in rev/s and D the diameter.
    let n = omega / std::f64::consts::TAU;
    let thrust = sol.loads.ct * 1.225 * n * n * 0.254_f64.powi(4);
    println!("hover thrust: {:.2} N", thrust);

    // Advance-ratio map: J = V / (n D).
    let speeds: Vec<f64> = (0..=8).map(|i| f64::from(i) * 5.0).collect();
    let map = performance_map(&solver, omega, 1.225, &speeds, 0.0).expect("map");
    println!("\n{:>6} {:>8} {:>8} {:>8}", "J", "C_T", "C_P", "eta");
    for p in map {
        println!(
            "{:>6.3} {:>8.4} {:>8.5} {:>8.3}",
            p.j, p.ct, p.cp, p.efficiency
        );
    }
}
