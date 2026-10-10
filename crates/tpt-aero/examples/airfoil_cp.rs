//! Pressure distribution over a NACA 4-digit airfoil with source panels,
//! plus a Prandtl–Glauert compressibility correction.
//!
//! Run: `cargo run --example airfoil_cp`

use tpt_aero::prelude::*;

fn main() {
    let code = "0012";
    // Closed trailing edge: a sharp TE with no gap panel.
    let airfoil = naca4_closed(code, 60).expect("valid NACA code");
    let points: Vec<(f64, f64)> = airfoil
        .x
        .iter()
        .zip(&airfoil.y)
        .map(|(&x, &y)| (x, y))
        .collect();

    // Source panels give the thickness (symmetric) pressure field; at zero
    // incidence a symmetric section carries no lift.
    let sol = solve_source_panels(&points, 0.0, 1.0).expect("solvable contour");

    println!("=== NACA {code}, source panels, alpha = 0 ===\n");
    println!("{:>7} {:>9} {:>14}", "x/c", "Cp (M=0)", "Cp (M=0.6, PG)");
    let n = points.len();
    // Upper surface runs TE -> LE, i.e. the first half of the contour.
    for i in (0..n / 2).step_by(6) {
        let x_mid = 0.5 * (points[i].0 + points[(i + 1) % n].0);
        let cp = sol.cp[i];
        let cp_c = compress::prandtl_glauert(cp, 0.6).expect("subsonic");
        println!("{x_mid:>7.3} {cp:>9.4} {cp_c:>14.4}");
    }
    let cp_min = sol.cp.iter().cloned().fold(f64::INFINITY, f64::min);
    println!("\nPeak suction Cp_min = {cp_min:.3} (incompressible; literature ~ -0.43)");
}
