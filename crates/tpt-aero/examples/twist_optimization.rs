//! Minimum induced drag over a linear washout, at several lift coefficients.
//!
//! A tapered wing is not elliptically loaded, so a linear washout (less twist at
//! the root than at the tips) can move its loading toward elliptic and cut the
//! induced drag. This example finds the tip twist that minimises the induced
//! drag at each lift coefficient, and compares it with the elliptic reference
//! `C_L²/(π·AR)`. It searches one parameter, the linear washout; a full twist
//! shape can do better.
//!
//! Run with `cargo run --release -p tpt-aero --example twist_optimization`.

use tpt_aero::vlm::optimize_linear_twist;
use tpt_aero::WingPlanform;

fn main() {
    // A tapered wing: 10 m span, root 1.2 m, tip 0.6 m, no sweep or dihedral.
    let wing = WingPlanform::tapered(10.0, 1.2, 0.6, 0.0, 0.0);
    println!("=== Induced-drag twist optimisation ===");
    println!(
        "tapered wing, span {} m, root {} m, tip {} m, AR {:.2}\n",
        10.0,
        1.2,
        0.6,
        100.0 / wing.area()
    );
    println!(
        "{:>6} {:>12} {:>10} {:>11} {:>11} {:>11} {:>10}",
        "C_L", "tip twist°", "alpha°", "C_Di", "C_Di(0°)", "C_Di ellip", "vs ellip"
    );
    for cl in [0.3_f64, 0.5, 0.7, 0.9] {
        let opt = optimize_linear_twist(&wing, 20, 4, cl, 50.0, 1.225).expect("optimises");
        println!(
            "{cl:>6.2} {:>12.2} {:>10.2} {:>11.6} {:>11.6} {:>11.6} {:>9.1}%",
            opt.tip_twist_deg,
            opt.alpha_deg,
            opt.cdi,
            opt.cdi_untwisted,
            opt.cdi_elliptic,
            100.0 * (opt.cdi / opt.cdi_elliptic - 1.0)
        );
    }
    println!("\n'vs ellip' is the optimum against the elliptic reference. A negative value means");
    println!("the VLM's near-field induced drag is below the textbook formula (Oswald e above 1,");
    println!("see todo.md). The washout gain is the C_Di(0°) column against the optimum: a");
    println!("fraction of a per cent for this planform.");
}
