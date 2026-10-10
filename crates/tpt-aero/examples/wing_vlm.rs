//! Sample wing VLM analysis: rectangular AR = 10 wing swept through an
//! angle-of-attack sweep, printing the force breakdown and spanwise loading.

use tpt_aero::{CamberModel, VlmModel, WingPlanform};

fn main() {
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let model = VlmModel::from_planform(&wing, 32, 4, CamberModel::Naca4 { m: 0.02, p: 0.4 })
        .expect("valid planform");
    let model = model.with_reference(tpt_aero::vlm::ReferenceGeometry::from_planform(&wing));

    println!("=== Wing VLM analysis ===");
    println!("AR = 10 rectangular, NACA 24-series camber, 32x4 lattice\n");
    println!(
        "{:>6} {:>8} {:>8} {:>8} {:>8}",
        "a[deg]", "C_L", "C_Di", "C_m", "e"
    );
    for deg in [-4.0_f64, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0] {
        let sol = model
            .solve(deg.to_radians(), 50.0, 1.225, &[])
            .expect("solvable");
        // The Oswald factor divides C_L² by a small C_Di, so it is only
        // meaningful away from zero lift; near C_L = 0 it is reported as n/a.
        let e = if sol.cl.abs() >= 0.1 {
            format!(
                "{:.3}",
                sol.cl * sol.cl / (core::f64::consts::PI * 10.0 * sol.cdi)
            )
        } else {
            "n/a".to_string()
        };
        println!(
            "{deg:>6.1} {:>8.4} {:>8.5} {:>8.4} {:>8}",
            sol.cl, sol.cdi, sol.cm, e
        );
    }

    let sol = model
        .solve(5.0_f64.to_radians(), 50.0, 1.225, &[])
        .expect("solvable");
    println!("\nSpanwise loading at alpha = 5 deg (32 strips):");
    let n = sol.strip_cl.len();
    for (i, &cl) in sol.strip_cl.iter().enumerate() {
        let eta = (i as f64 + 0.5) / n as f64 * 2.0 - 1.0;
        let bar = "#".repeat((cl.abs() * 40.0) as usize);
        println!("  y/b={eta:+.3}  c_l={cl:.4}  {bar}");
    }
}
