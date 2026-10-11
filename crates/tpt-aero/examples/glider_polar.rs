//! Glider polar: lift and drag over an angle-of-attack sweep, with the best
//! lift-to-drag ratio and the minimum-sink condition.
//!
//! Scenario: a 15 m, taper-free glider at 20 m/s. The parasite drag is a fixed
//! `C_D0`, the induced drag comes from the VLM, and the polar is `C_D = C_D0 +
//! C_Di`. The best glide ratio is `max C_L/C_D`, and the minimum sink rate is
//! `max C_L^{3/2}/C_D`, which sets the sink rate `w = √(2W/ρS) · (C_D/C_L^{3/2})`.
//!
//! Run with `cargo run --release -p tpt-aero --example glider_polar`.

use tpt_aero::vlm::{CamberModel, ReferenceGeometry, VlmModel};
use tpt_aero::WingPlanform;

const SPAN: f64 = 15.0;
const CHORD_ROOT: f64 = 1.0;
const CHORD_TIP: f64 = 0.5;
const SPEED: f64 = 20.0;
const RHO: f64 = 1.225;
const MASS_KG: f64 = 250.0;
const G: f64 = 9.806_65;
const CD0: f64 = 0.008;

fn main() {
    let wing = WingPlanform::tapered(SPAN, CHORD_ROOT, CHORD_TIP, 0.0, 3.0_f64.to_radians());
    let area = wing.area();
    let model = VlmModel::from_planform(&wing, 30, 4, CamberModel::Naca4 { m: 0.015, p: 0.4 })
        .expect("valid planform")
        .with_reference(ReferenceGeometry::from_planform(&wing));
    let weight = MASS_KG * G;
    let aspect = SPAN * SPAN / area;

    println!("=== Glider polar ===");
    println!(
        "span {SPAN} m, S = {area:.2} m², AR = {aspect:.1}, mass {MASS_KG} kg, V = {SPEED} m/s\n"
    );
    println!(
        "{:>8} {:>8} {:>9} {:>9} {:>8} {:>10}",
        "alpha[deg]", "C_L", "C_Di", "C_D", "L/D", "sink[m/s]"
    );

    let mut best_ld = (0.0_f64, 0.0_f64);
    let mut min_sink = (f64::INFINITY, 0.0_f64);
    let mut alpha = -2.0_f64;
    while alpha <= 12.0 + 1.0e-9 {
        let sol = model
            .solve(alpha.to_radians(), SPEED, RHO, &[])
            .expect("solvable");
        let cd = sol.cdi + CD0;
        if sol.cl > 0.05 {
            let ld = sol.cl / cd;
            // Sink rate at this lift: w = √(2W/(ρS)) · C_D / C_L^{3/2}.
            let sink = (2.0 * weight / (RHO * area)).sqrt() * cd / sol.cl.powf(1.5);
            // Sink is the trim sink at the speed that carries the weight at this C_L,
            // so the table uses the same W and S for every row.
            if ld > best_ld.0 {
                best_ld = (ld, alpha);
            }
            if sink < min_sink.0 {
                min_sink = (sink, alpha);
            }
            println!(
                "{alpha:>8.1} {:>8.4} {:>9.5} {cd:>9.5} {ld:>8.2} {sink:>10.3}",
                sol.cl, sol.cdi
            );
        } else {
            println!(
                "{alpha:>8.1} {:>8.4} {:>9.5} {cd:>9.5} {:>8} {:>10}",
                sol.cl, sol.cdi, "n/a", "n/a"
            );
        }
        alpha += 1.0;
    }
    println!(
        "\nBest glide ratio {:.1} at alpha {:.1}° (C_L is about {:.2} there).",
        best_ld.0,
        best_ld.1,
        model
            .solve(best_ld.1.to_radians(), SPEED, RHO, &[])
            .expect("solvable")
            .cl
    );
    println!(
        "Minimum sink {:.2} m/s at alpha {:.1}°. The sink figure assumes the flight speed",
        min_sink.0, min_sink.1
    );
    println!("implied by the lift coefficient; read it as a trend, not a glider's sink polar.");
}
