//! Two-dimensional panel methods on a NACA 0012, plus a Prandtl-Glauert
//! correction to Mach 0.6.
//!
//! Two results are shown, each the one the crate validates:
//!
//! * the lift from the doublet (lifting) solver at 4 degrees, which agrees
//!   with thick-airfoil theory;
//! * the pressure distribution from the source (non-lifting) solver at zero
//!   incidence, whose peak suction matches the published inviscid value.
//!
//! The lifting solver's surface pressure at incidence is also checked against
//! thin-airfoil theory in the crate's tests; see the crate README.
//!
//! Run with: `cargo run -p tpt-aero-panel --example airfoil_pressure`

use tpt_aero_core::airfoil::naca4_closed;
use tpt_aero_panel::{compress, solve_doublet_panels, solve_source_panels};

fn main() {
    // Closed trailing edge, counter-clockwise ordering (required by the solver).
    let foil = naca4_closed("0012", 100).expect("valid NACA code");
    let points: Vec<(f64, f64)> = foil.x.iter().zip(&foil.y).map(|(&x, &y)| (x, y)).collect();

    // Lift at 4 degrees, 50 m/s, 1 m chord.
    let alpha = 4.0_f64.to_radians();
    let lift = solve_doublet_panels(&points, alpha, 50.0, 1.0).expect("solvable");
    let thick_theory = 2.0 * std::f64::consts::PI * (1.0 + 0.77 * 0.12) * alpha;
    println!(
        "NACA 0012 at 4 deg: cl = {:.3} (thick-airfoil theory {:.3})",
        lift.cl, thick_theory
    );

    // Pressure at zero incidence, from source panels.
    let source = solve_source_panels(&points, 0.0, 1.0).expect("solvable");
    let (i_min, cp_min) =
        source.cp.iter().enumerate().fold(
            (0, f64::MAX),
            |best, (i, &c)| if c < best.1 { (i, c) } else { best },
        );
    // Each panel is between two consecutive points; take its start point's x.
    println!(
        "\nzero incidence: peak suction Cp {:.3} at x/c {:.3}",
        cp_min, points[i_min].0
    );
    println!("(published inviscid reference: about -0.43 near x/c 0.1)");

    // Compressibility: correct that peak suction to Mach 0.6.
    let cp_mach = compress::prandtl_glauert(cp_min, 0.6).expect("subsonic");
    println!("Prandtl-Glauert at M 0.6: {:.3}", cp_mach);
}
