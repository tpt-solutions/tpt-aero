//! Longitudinal stability derivatives of a wing and its short-period mode.
//!
//! The pitching moment is referenced ahead of the aerodynamic centre so the
//! airframe is statically stable; slide `x_ref` aft to watch C_m_alpha go
//! positive and the short-period mode vanish.
//!
//! Run: `cargo run --example stability_derivatives`

use tpt_aero::prelude::*;

fn main() {
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let (v, rho) = (50.0, 1.225);

    println!("=== Longitudinal derivatives (rectangular AR 10) ===\n");
    println!(
        "{:>8} {:>8} {:>9} {:>9} {:>9} {:>9}",
        "x_ref/c", "CLa", "Cma", "Cmq", "wn[rad/s]", "zeta"
    );
    for x_ref in [0.0, 0.1, 0.2, 0.25, 0.3] {
        let model = VlmModel::from_planform(&wing, 16, 4, CamberModel::None)
            .expect("valid planform")
            .with_reference(ReferenceGeometry {
                area: wing.area(),
                chord: wing.mac(),
                span: wing.span(),
                moment_center: Vec3::new(x_ref * wing.mac(), 0.0, 0.0),
            });
        let d = longitudinal_fd(&model, 2.0_f64.to_radians(), v, rho, 1e-4).expect("solve");
        let params = PlantParameters {
            v_inf: v,
            chord: wing.mac(),
            area: wing.area(),
            mass: 500.0,
            inertia_y: 400.0,
            density: rho,
        };
        let plant = LongitudinalPlant::from_derivatives(&d, &params, 0.05);
        match plant.short_period() {
            Some((wn, zeta)) => println!(
                "{x_ref:>8.2} {:>8.3} {:>9.3} {:>9.3} {wn:>9.3} {zeta:>9.3}",
                d.cl_alpha, d.cm_alpha, d.cm_q
            ),
            None => println!(
                "{x_ref:>8.2} {:>8.3} {:>9.3} {:>9.3} {:>9} {:>9}",
                d.cl_alpha, d.cm_alpha, d.cm_q, "unstable", "-"
            ),
        }
    }
}
