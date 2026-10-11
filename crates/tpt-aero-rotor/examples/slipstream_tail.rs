//! Propeller slipstream over a horizontal tail: how much the tail's dynamic
//! pressure rises, and what that does to the pitch stiffness.
//!
//! The disk sits at the nose (x = 0). A horizontal tail 4 m behind it spans
//! ±2 m on the propeller axis height (z = 0), so it sits in the slipstream
//! tube. Run with `cargo run -p tpt-aero-rotor --example slipstream_tail`.

use tpt_aero_rotor::Slipstream;

fn main() {
    let rho = 1.225_f64;
    let v = 25.0_f64; // m/s
    let radius = 0.9_f64; // m
    let x_tail = 4.0_f64; // m behind the disk
    let tail_area = 4.0 * 0.6; // m²
    let wing_area = 12.0_f64; // m²
    let wing_mac = 1.2_f64; // m
    let tail_arm = 4.6_f64; // m, tail quarter chord from the CG
    let tail_cla = 4.0_f64; // per rad, tail lift-curve slope (assumed)
    let q = 0.5 * rho * v * v;

    // Tail volume coefficient and the classical pitch-stiffness contribution.
    let v_h = tail_area * tail_arm / (wing_area * wing_mac);
    println!("propeller radius {radius} m, freestream {v} m/s, tail at x = {x_tail} m");
    println!("tail volume V_H = {v_h:.3}, tail lift slope {tail_cla:.2} /rad\n");
    println!("thrust (N)  v_d (m/s)  far wake (m/s)  tube R_w at tail (m)  mean dq/q over tail  dCm_alpha (/rad)");

    for &thrust in &[0.0, 50.0, 100.0, 150.0] {
        let s = Slipstream::new(radius, thrust, v, rho).expect("valid propeller");
        let span_points = 9;
        let mut dq_sum = 0.0;
        for i in 0..span_points {
            let y = -2.0 + 4.0 * i as f64 / (span_points - 1) as f64;
            dq_sum += s.dynamic_pressure_increment(x_tail, y);
        }
        let dq_over_q = dq_sum / span_points as f64 / q;
        // Tail lift rises by Δq/q, so C_mα changes by −(Δq/q)·V_H·a_t.
        let dcm_alpha = -dq_over_q * v_h * tail_cla;
        println!(
            "{thrust:>9.0}  {:>9.3}  {:>14.3}  {:>20.3}  {:>19.3}  {:>16.3}",
            s.disk_induced_velocity(),
            v + 2.0 * s.disk_induced_velocity(),
            s.tube_radius(x_tail),
            dq_over_q,
            dcm_alpha
        );
    }
}
