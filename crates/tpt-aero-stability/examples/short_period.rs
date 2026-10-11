//! From a VLM model to a short-period mode: stability derivatives, a
//! linearised longitudinal plant, and the natural frequency and damping of
//! the short-period oscillation.
//!
//! Run with: `cargo run -p tpt-aero-stability --example short_period`

use tpt_aero_core::{Vec3, WingPlanform};
use tpt_aero_stability::longitudinal::longitudinal_fd;
use tpt_aero_stability::plant::{LongitudinalPlant, PlantParameters};
use tpt_aero_vlm::{CamberModel, ReferenceGeometry, VlmModel};

fn main() {
    // A 10 m x 1 m wing. Moments are taken about the leading edge, which sits
    // ahead of the neutral point, so the wing alone is statically stable.
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let model = VlmModel::from_planform(&wing, 16, 4, CamberModel::None)
        .expect("valid")
        .with_reference(ReferenceGeometry {
            area: wing.area(),
            chord: wing.mac(),
            span: wing.span(),
            moment_center: Vec3::new(0.0, 0.0, 0.0),
        });

    // Derivatives at 2 degrees, 50 m/s, sea level.
    let (v, rho) = (50.0, 1.225);
    let d = longitudinal_fd(&model, 2.0_f64.to_radians(), v, rho, 1e-4).expect("solves");
    println!("C_L,alpha = {:.3} /rad", d.cl_alpha);
    println!(
        "C_m,alpha = {:.3} /rad  (negative: statically stable)",
        d.cm_alpha
    );
    println!("C_m,q     = {:.3}       (pitch damping)", d.cm_q);

    // Plant: a 500 kg aircraft, 10 m^2 reference area, 1 m chord.
    let params = PlantParameters {
        v_inf: v,
        chord: 1.0,
        area: 10.0,
        mass: 500.0,
        inertia_y: 400.0,
        density: rho,
    };
    let plant = LongitudinalPlant::from_derivatives(&d, &params, 0.05);
    match plant.short_period() {
        Some((wn, zeta)) => {
            println!(
                "\nshort period: omega_n = {:.2} rad/s (period {:.2} s), zeta = {:.3}",
                wn,
                std::f64::consts::TAU / wn,
                zeta
            );
        }
        None => println!("\nno short-period mode: the plant is statically unstable"),
    }
}
