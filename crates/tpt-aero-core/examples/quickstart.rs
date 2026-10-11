//! Tour of `tpt-aero-core`: geometry, airfoils, unit-checked quantities,
//! a drag build-up and the native linear solver.
//!
//! Run with: `cargo run -p tpt-aero-core --example quickstart`

use tpt_aero_core::airfoil::naca4;
use tpt_aero_core::viscous::{cd0_buildup, flat_plate_cf, DragComponent};
use tpt_aero_core::{lu_solve, MachNumber, WingPlanform};

fn main() {
    // Geometry: a 10 m span, 1 m chord rectangular wing.
    let wing = WingPlanform::rectangular(10.0, 1.0);
    println!(
        "wing: span {:.1} m, area {:.1} m^2, MAC {:.2} m",
        wing.span(),
        wing.area(),
        wing.mac()
    );

    // Airfoil: NACA 2412 with 61 stations per surface. `x` is chordwise, `y` the ordinate.
    let foil = naca4("2412", 61).expect("valid NACA code");
    let y_max = foil.y.iter().cloned().fold(f64::MIN, f64::max);
    println!(
        "NACA 2412: {} points, max ordinate {:.4} (about 0.12 of the chord)",
        foil.len(),
        y_max
    );

    // Unit-checked quantities reject nonsense on construction.
    match MachNumber::try_new(0.78) {
        Some(m) => println!("Mach {:.2} accepted", m.value()),
        None => println!("Mach rejected"),
    }
    println!(
        "Mach -0.1 accepted? {}",
        MachNumber::try_new(-0.1).is_some()
    );

    // Skin friction: laminar-to-turbulent flat plate at Re = 5e6.
    println!("flat-plate Cf at Re 5e6: {:.5}", flat_plate_cf(5.0e6));

    // Zero-lift drag build-up: wing and fuselage, each with its own form factor.
    let components = [
        DragComponent {
            reynolds: 5.0e6,
            wetted_ratio: 2.0,
            thickness_ratio: 0.12,
            form_override: None,
            interference: 1.1,
        },
        DragComponent {
            reynolds: 2.0e7,
            wetted_ratio: 3.5,
            thickness_ratio: 0.0,
            form_override: Some(1.3),
            interference: 1.0,
        },
    ];
    println!("C_D0 build-up: {:.4}", cd0_buildup(&components));

    // Native dense LU solve of  [2 1; 1 3] x = [3; 5]  ->  x = [0.8, 1.4].
    let a = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
    let x = lu_solve(&a, &[3.0, 5.0]).expect("non-singular");
    println!(
        "lu_solve: x = [{:.3}, {:.3}] (expected [0.800, 1.400])",
        x[0], x[1]
    );
}
