//! The atmosphere, a wind-shear profile and a discrete 1-cos gust.
//!
//! Run with: `cargo run -p tpt-aero-atmo --example atmosphere_gusts`

use tpt_aero_atmo::isa;
use tpt_aero_atmo::shear::{power_law, DiscreteGust, GustShape};

fn main() {
    // Standard atmosphere, with a hot day for comparison.
    println!(
        "{:>8} {:>10} {:>10} {:>10}",
        "alt (m)", "T (K)", "p (Pa)", "rho"
    );
    for h in [0.0, 1000.0, 5000.0, 11000.0] {
        let s = isa::at(h).expect("within range");
        println!(
            "{:>8.0} {:>10.2} {:>10.0} {:>10.4}",
            h, s.temperature, s.pressure, s.density
        );
    }
    let hot = isa::at_isa_offset(0.0, 15.0).expect("within range");
    println!(
        "sea level, ISA+15 K: rho {:.4} (standard 1.2250)",
        hot.density
    );

    // Power-law wind profile: 10 m/s at 10 m, exponent 1/7.
    println!("\nwind shear, 10 m/s at 10 m:");
    for z in [10.0, 50.0, 100.0, 200.0] {
        println!(
            "  {:>5.0} m  {:>6.2} m/s",
            z,
            power_law(10.0, 10.0, z, 1.0 / 7.0)
        );
    }

    // 1-cos gust: 3 m/s peak, gradient distance 15 m.
    let gust = DiscreteGust::new(3.0, 15.0, GustShape::OneMinusCosine).expect("valid gust");
    println!("\n1-cos gust (3 m/s peak, H = 15 m), vertical speed vs penetration:");
    for k in 0..=8 {
        let x = f64::from(k) * 7.5; // 0 .. 60 m, spanning the 2H gust length
        println!("  x = {:>5.1} m  w = {:>6.3} m/s", x, gust.speed_at(x));
    }
}
