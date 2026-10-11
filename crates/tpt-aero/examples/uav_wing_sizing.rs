//! UAV wing sizing: pick the wing area for a cruise lift coefficient, then
//! trade span against aspect ratio and report the induced drag and lift-to-drag
//! ratio at each option.
//!
//! Scenario: a 2.5 kg fixed-wing UAV cruising at 15 m/s at sea level. The wing
//! area is `S = W / (½ ρ V² C_L)` for the target `C_L`. For each span the chord
//! is `c = S / b`, and the angle of attack that gives the target `C_L` comes
//! from a secant solve on the VLM. The parasite drag is a fixed `C_D0`, so the
//! table shows how aspect ratio trades induced drag against structure.
//!
//! Run with `cargo run --release -p tpt-aero --example uav_wing_sizing`.

use tpt_aero::vlm::{CamberModel, ReferenceGeometry, VlmModel};
use tpt_aero::WingPlanform;

const MASS_KG: f64 = 2.5;
const SPEED: f64 = 15.0;
const RHO: f64 = 1.225;
const G: f64 = 9.806_65;
const CL_TARGET: f64 = 0.6;
const CD0: f64 = 0.02;

/// Lift coefficient of a rectangular wing at `alpha` (rad).
fn cl_at(span: f64, chord: f64, alpha: f64) -> (f64, f64) {
    let wing = WingPlanform::rectangular(span, chord);
    let model = VlmModel::from_planform(&wing, 16, 3, CamberModel::Naca4 { m: 0.02, p: 0.4 })
        .expect("valid planform")
        .with_reference(ReferenceGeometry::from_planform(&wing));
    let sol = model.solve(alpha, SPEED, RHO, &[]).expect("solvable");
    (sol.cl, sol.cdi)
}

/// Secant solve for the angle of attack (rad) that gives `CL_TARGET`.
fn alpha_for_target(span: f64, chord: f64) -> (f64, f64) {
    let f = |a: f64| cl_at(span, chord, a).0 - CL_TARGET;
    let (mut a0, mut a1) = (0.0_f64, 0.08_f64);
    let (mut f0, mut f1) = (f(a0), f(a1));
    for _ in 0..12 {
        if (f1 - f0).abs() < 1.0e-14 {
            break;
        }
        let a2 = a1 - f1 * (a1 - a0) / (f1 - f0);
        a0 = a1;
        f0 = f1;
        a1 = a2;
        f1 = f(a1);
        if f1.abs() < 1.0e-10 {
            break;
        }
    }
    (a1, cl_at(span, chord, a1).1)
}

fn main() {
    let weight = MASS_KG * G;
    let area = weight / (0.5 * RHO * SPEED * SPEED * CL_TARGET);
    println!("=== UAV wing sizing ===");
    println!("mass {MASS_KG} kg, cruise {SPEED} m/s, C_L target {CL_TARGET}");
    println!("required wing area S = {area:.3} m²\n");
    println!(
        "{:>7} {:>7} {:>6} {:>9} {:>9} {:>9} {:>8}",
        "span[m]", "chord", "AR", "alpha[deg]", "C_Di", "C_D", "L/D"
    );
    for span in [1.6_f64, 1.8, 2.0, 2.2, 2.4] {
        let chord = area / span;
        let aspect = span * span / area;
        let (alpha, cdi) = alpha_for_target(span, chord);
        let cd = cdi + CD0;
        println!(
            "{span:>7.2} {chord:>7.3} {aspect:>6.2} {:>9.2} {cdi:>9.5} {cd:>9.5} {:>8.2}",
            alpha.to_degrees(),
            CL_TARGET / cd
        );
    }
    println!("\nA higher aspect ratio lowers C_Di for the same lift. Structure and the");
    println!("span limit of the launch and transport set the upper end.");
}
