//! A wing polar from a VLM alpha sweep, checked against lifting-line theory.
//!
//! For an elliptic-limit wing the lift slope is `2 pi AR / (AR + 2)` per
//! radian, so the printed slope should land close to that value.
//!
//! Run with: `cargo run -p tpt-aero-vlm --example wing_polar`

use tpt_aero_core::WingPlanform;
use tpt_aero_vlm::{Surface, VlmModel};

fn main() {
    // Rectangular wing, aspect ratio 10, at 50 m/s, sea-level density.
    let (span, chord) = (10.0, 1.0);
    let aspect_ratio = span * span / (span * chord);
    let model = VlmModel::from_surfaces(&[Surface::horizontal(
        WingPlanform::rectangular(span, chord),
        24,
        4,
    )])
    .expect("valid configuration");

    // The sweep factors the influence matrix once, then solves each angle.
    let degrees: Vec<i32> = (-2..=10).step_by(2).collect();
    let alphas: Vec<f64> = degrees.iter().map(|&d| f64::from(d).to_radians()).collect();
    let sweep = model.solve_sweep(&alphas, 50.0, 1.225).expect("solvable");

    println!("{:>6} {:>8} {:>10} {:>9}", "alpha", "CL", "CDi", "Cm");
    for (d, s) in degrees.iter().zip(&sweep) {
        println!("{:>6} {:>8.4} {:>10.5} {:>9.4}", d, s.cl, s.cdi, s.cm);
    }

    // Least-squares lift slope against alpha (per radian).
    let n = alphas.len() as f64;
    let (sx, sy) = (
        alphas.iter().sum::<f64>(),
        sweep.iter().map(|s| s.cl).sum::<f64>(),
    );
    let sxy: f64 = alphas.iter().zip(&sweep).map(|(a, s)| a * s.cl).sum();
    let sxx: f64 = alphas.iter().map(|a| a * a).sum();
    let slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    let theory = 2.0 * std::f64::consts::PI * aspect_ratio / (aspect_ratio + 2.0);
    println!();
    println!("lift slope (VLM) {:.3} /rad", slope);
    println!(
        "elliptic theory  {:.3} /rad  ({:+.1} %)",
        theory,
        100.0 * (slope / theory - 1.0)
    );
}
