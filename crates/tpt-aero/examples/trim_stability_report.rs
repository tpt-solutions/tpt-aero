//! Trim and stability report for a small conventional aircraft: a wing, a
//! tailplane with an elevator, and a fin with a rudder. Trims the aircraft in
//! level flight, then reports the static margin, the short-period and Dutch-roll
//! modes, and the main derivatives.
//!
//! The inertias are estimated from the span (a uniform-mass rule of thumb), so
//! the mode frequencies are illustrative. Replace them with measured values for
//! a real analysis. The aerodynamics and the trim are the library's.
//!
//! Run with `cargo run --release -p tpt-aero --example trim_stability_report`.

use tpt_aero::aero_core::{Vec3, WingPlanform};
use tpt_aero::stability::lateral::{lateral_fd, LateralParameters, LateralPlant};
use tpt_aero::stability::longitudinal::longitudinal_fd;
use tpt_aero::stability::plant::{LongitudinalPlant, PlantParameters};
use tpt_aero::stability::trim::{trim_level_flight, TrimRequest};
use tpt_aero::vlm::{ControlSymmetry, ReferenceGeometry, Surface, VlmModel};

const MASS: f64 = 500.0;
const G: f64 = 9.806_65;
const SPEED: f64 = 45.0;
const RHO: f64 = 1.225;
const CD0: f64 = 0.03;
/// Centre of gravity, streamwise from the wing root leading edge (m).
const CG_X: f64 = 0.4;

fn build() -> VlmModel {
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let tail = WingPlanform::rectangular(4.0, 0.8);
    let fin = WingPlanform::rectangular(1.2, 0.8);
    let model = VlmModel::from_surfaces(&[
        Surface::horizontal(wing.clone(), 16, 4),
        Surface::horizontal(tail, 8, 3)
            .at(4.5, 0.0, 0.0)
            .with_control("elevator", 0.6, 0.0, 2.0, ControlSymmetry::Symmetric),
        Surface::vertical(fin, 6, 3).at(4.5, 0.0, 0.0).with_control(
            "rudder",
            0.6,
            0.0,
            0.6,
            ControlSymmetry::Symmetric,
        ),
    ])
    .expect("valid aircraft");
    // The CG is the moment reference: it sets C_m and the static margin.
    let r = *model.reference();
    model.with_reference(ReferenceGeometry {
        area: r.area,
        chord: r.chord,
        span: r.span,
        moment_center: Vec3::new(CG_X, 0.0, 0.0),
    })
}

fn main() {
    let model = build();
    let weight = MASS * G;
    let span = model.reference().span;
    let chord = model.reference().chord;

    println!("=== Trim and stability report ===");
    println!("mass {MASS} kg, level flight at {SPEED} m/s, sea level, CG at x = {CG_X} m\n");

    let trim = trim_level_flight(
        &model,
        &TrimRequest {
            weight_n: weight,
            v_inf: SPEED,
            density: RHO,
            elevator: "elevator",
            cd0: CD0,
        },
    )
    .expect("trims");
    let q_s = 0.5 * RHO * SPEED * SPEED * model.reference().area;
    println!("-- trim --");
    println!("alpha            {:>8.3} deg", trim.alpha_rad.to_degrees());
    println!(
        "elevator         {:>8.3} deg",
        trim.elevator_rad.to_degrees()
    );
    println!("C_L              {:>8.4}", weight / q_s);
    println!("C_D (T/qS)       {:>8.4}", trim.thrust_n / q_s);
    println!("L/D              {:>8.2}", weight / trim.thrust_n);

    let d = longitudinal_fd(&model, trim.alpha_rad, SPEED, RHO, 1.0e-4).expect("derivatives");
    let static_margin = -d.cm_alpha / d.cl_alpha;
    println!("\n-- longitudinal --");
    println!("C_L,alpha        {:>8.3} /rad", d.cl_alpha);
    println!("C_m,alpha        {:>8.3} /rad", d.cm_alpha);
    println!(
        "static margin    {:>8.3} of the reference chord ({:.2} m)",
        static_margin,
        static_margin * chord
    );
    println!("C_m,q            {:>8.3}", d.cm_q);
    println!("C_m,delta_e      {:>8.3} /rad", d.cm_delta);

    // Illustrative inertias from a uniform-mass rule of thumb (see the header).
    let inertia_y = MASS * (span / 3.0).powi(2);
    let plant = LongitudinalPlant::from_derivatives(
        &d,
        &PlantParameters {
            v_inf: SPEED,
            chord,
            area: model.reference().area,
            mass: MASS,
            inertia_y,
            density: RHO,
        },
        0.0,
    );
    println!("\n-- short period (illustrative inertia I_y = {inertia_y:.0} kg·m²) --");
    match plant.short_period() {
        Some((wn, zeta)) => println!("ω_n {wn:.3} rad/s, ζ {zeta:.3}"),
        None => println!("not oscillatory: the airframe is statically unstable"),
    }

    let lat = lateral_fd(&model, trim.alpha_rad, SPEED, RHO).expect("lateral derivatives");
    println!("\n-- lateral-directional --");
    println!("C_Y,beta         {:>8.3} /rad", lat.cy_beta);
    println!(
        "C_l,beta         {:>8.3} /rad   (dihedral effect)",
        lat.cl_beta
    );
    println!("C_n,beta         {:>8.3} /rad   (weathercock)", lat.cn_beta);

    let inertia_x = MASS * (span / 4.0).powi(2);
    let inertia_z = MASS * (span / 2.5).powi(2);
    let lateral = LateralPlant::from_derivatives(
        &lat,
        &LateralParameters {
            v_inf: SPEED,
            span,
            area: model.reference().area,
            density: RHO,
            mass: MASS,
            inertia_x,
            inertia_z,
            inertia_xz: 0.0,
        },
    );
    match lateral.and_then(|p| p.modes()) {
        Some(m) => {
            println!(
                "Dutch roll       ω_n {:.3} rad/s, ζ {:.3}",
                m.dutch_roll_wn, m.dutch_roll_zeta
            );
            println!("roll subsidence  {:.3} 1/s", m.roll_root);
            println!(
                "spiral           {:.4} 1/s (negative is stable)",
                m.spiral_root
            );
        }
        None => println!("lateral modes not available for these inputs"),
    }
    println!("\nInertias are estimates; the trim, derivatives and aerodynamics are computed.");
}
