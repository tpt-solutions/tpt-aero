//! International Standard Atmosphere table with Mach-to-airspeed conversion.
//!
//! Run: `cargo run --example isa_table`

use tpt_aero::prelude::*;

fn main() {
    println!("=== ISA (US Standard Atmosphere 1976) ===\n");
    println!(
        "{:>8} {:>8} {:>10} {:>9} {:>8} {:>12}",
        "h [m]", "T [K]", "p [Pa]", "rho", "a [m/s]", "V @ M0.3 [m/s]"
    );
    for h in [
        0.0, 1_000.0, 3_000.0, 5_000.0, 8_000.0, 11_000.0, 15_000.0, 20_000.0,
    ] {
        let s = isa::at(h).expect("altitude in range");
        println!(
            "{h:>8.0} {:>8.2} {:>10.1} {:>9.4} {:>8.2} {:>12.2}",
            s.temperature,
            s.pressure,
            s.density,
            s.speed_of_sound,
            0.3 * s.speed_of_sound
        );
    }
}
