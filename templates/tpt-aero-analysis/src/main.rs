//! A starting point for a tpt-aero analysis: a tapered wing swept through
//! angle of attack, reporting lift, induced drag and the pitching moment.
//!
//! Run with `cargo run --release`.

use tpt_aero::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Planform: 10 m span, 1.4 m root chord tapering to 0.8 m, 0.5 m aft sweep.
    let wing = WingPlanform::tapered(10.0, 1.4, 0.8, 0.5, 0.0);

    // Discretize into 24 spanwise by 4 chordwise horseshoe vortices, with a
    // NACA-style camber line (m = 2 %, p = 40 %).
    let model = VlmModel::from_planform(
        &wing,
        24,
        4,
        CamberModel::Naca4 { m: 0.02, p: 0.4 },
    )?;

    let speed = 50.0; // m/s
    let density = 1.225; // kg/m³ (sea-level ISA)

    println!("{:>8} {:>10} {:>10} {:>10}", "alpha", "CL", "CDi", "Cm");
    for deg in (-4..=10).step_by(2) {
        let sol = model.solve(f64::from(deg).to_radians(), speed, density, &[])?;
        println!("{deg:>8} {:>10.4} {:>10.5} {:>10.4}", sol.cl, sol.cdi, sol.cm);
    }
    Ok(())
}
