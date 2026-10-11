//! Two classic stability-augmentation designs on a light aircraft: a yaw
//! damper (rudder from yaw rate) that fixes a lightly damped Dutch roll, and
//! a pitch-attitude hold (elevator from pitch angle and rate).
//!
//! Everything is analysed through eigenvalues and checked in the time domain
//! with the RK4 simulator. Sign convention: a positive rudder command yaws
//! the nose right, a positive elevator command pitches it down.
//!
//! Run: `cargo run --example control_design`

use tpt_aero::core::WingSection;
use tpt_aero::prelude::*;
use tpt_aero::stability::modes::Complex;
use tpt_aero::stability::sim::{simulate, StateSpace};

const MASS: f64 = 600.0;
const G: f64 = 9.80665;

fn build_aircraft(dihedral_deg: f64) -> VlmModel {
    let section = |y: f64, chord: f64| WingSection {
        le: Vec3::new(0.0, y, 0.0),
        chord,
        twist: 0.0,
    };
    let fin = WingPlanform::from_sections(vec![section(0.0, 0.9), section(1.5, 0.5)])
        .expect("two ordered sections");
    let wing = WingPlanform::tapered(10.0, 1.2, 1.2, 0.0, dihedral_deg.to_radians());
    VlmModel::from_surfaces(&[
        Surface::horizontal(wing, 20, 3).with_control(
            "aileron",
            0.75,
            2.5,
            4.8,
            ControlSymmetry::Antisymmetric,
        ),
        Surface::horizontal(WingPlanform::rectangular(4.0, 0.6), 10, 3)
            .at(4.0, 0.0, 0.0)
            .with_control("elevator", 0.6, 0.0, 2.0, ControlSymmetry::Symmetric),
        Surface::vertical(fin, 6, 3).at(3.8, 0.0, 0.0).with_control(
            "rudder",
            0.6,
            0.0,
            1.5,
            ControlSymmetry::Symmetric,
        ),
    ])
    .expect("valid configuration")
}

fn describe(eig: &[Complex; 4]) -> String {
    let mut parts = Vec::new();
    for z in eig {
        if z.im > 0.0 {
            parts.push(format!(
                "oscillation wn {:.2} rad/s zeta {:.2}",
                z.abs(),
                -z.re / z.abs()
            ));
        } else if z.im == 0.0 {
            parts.push(format!("real {:+.3} /s", z.re));
        }
    }
    parts.join("; ")
}

fn main() {
    let (v, rho) = (50.0, 1.225);
    let alpha = 3.0_f64.to_radians();

    // ---------------------------------------------------------------- yaw damper
    let model = build_aircraft(3.0);
    let d = lateral_fd(&model, alpha, v, rho).expect("solvable");
    let params = LateralParameters {
        v_inf: v,
        span: model.reference().span,
        area: model.reference().area,
        density: rho,
        mass: MASS,
        inertia_x: 1000.0,
        inertia_z: 1800.0,
        inertia_xz: 0.0,
    };
    let plant = LateralPlant::from_derivatives(&d, &params).expect("valid plant");
    let aileron = control_derivatives(&model, "aileron", alpha, v, rho).expect("control exists");
    let rudder = control_derivatives(&model, "rudder", alpha, v, rho).expect("control exists");
    let open = plant.state_space(&params, &aileron, &rudder);

    println!("=== Yaw damper (rudder = -k * yaw rate) ===");
    println!("open loop:   {}", describe(&open.eigenvalues()));
    let gain = 0.5;
    let closed = open.with_feedback(1, &[0.0, 0.0, gain, 0.0]);
    println!("k = {gain}:     {}", describe(&closed.eigenvalues()));

    let kick = [0.05, 0.0, 0.0, 0.0]; // 2.9 deg of sideslip
    let peak_yaw_rate = |sys: &StateSpace<4, 3>, t0: f64, t1: f64| {
        simulate(sys, kick, 0.005, 2400, |_| [0.0; 3])
            .iter()
            .filter(|s| s.t >= t0 && s.t <= t1)
            .map(|s| s.x[2].abs())
            .fold(0.0, f64::max)
    };
    println!(
        "yaw-rate peak after a sideslip kick, 3-6 s: open {:.4} rad/s, closed {:.5} rad/s",
        peak_yaw_rate(&open, 3.0, 6.0),
        peak_yaw_rate(&closed, 3.0, 6.0)
    );

    // ------------------------------------------------------- pitch attitude hold
    let model = build_aircraft(0.0);
    let trim = trim_level_flight(
        &model,
        &TrimRequest {
            weight_n: MASS * G,
            v_inf: v,
            density: rho,
            elevator: "elevator",
            cd0: 0.03,
        },
    )
    .expect("trims");
    let r = model.reference();
    let long_params = FullLongitudinalParameters {
        v_inf: v,
        chord: r.chord,
        area: r.area,
        density: rho,
        mass: MASS,
        inertia_y: 900.0,
        cl: trim.cl,
        cd: trim.cdi + 0.03,
        cd_alpha: drag_alpha_fd(&model, trim.alpha_rad, v, rho).expect("solvable"),
    };
    let d = longitudinal_fd(&model, trim.alpha_rad, v, rho, 1e-4).expect("solve");
    let elevator =
        control_derivatives(&model, "elevator", trim.alpha_rad, v, rho).expect("control exists");
    let plant = FullLongitudinalPlant::from_derivatives(&d, &long_params).expect("valid plant");
    let open = plant.state_space(&long_params, &elevator);

    println!("\n=== Pitch attitude hold (elevator = k_theta*(theta - cmd) + k_q*q) ===");
    println!("open loop:   {}", describe(&open.eigenvalues()));
    let (k_theta, k_q) = (2.0, 0.3);
    let closed = open.with_feedback(0, &[0.0, 0.0, -k_q, -k_theta]);
    println!(
        "k_theta {k_theta}, k_q {k_q}: {}",
        describe(&closed.eigenvalues())
    );

    let theta_cmd = 5.0_f64.to_radians();
    let command = -k_theta * theta_cmd; // so that u = k_theta (theta - cmd) + k_q q
    let dt = 0.005;
    let hist = simulate(&closed, [0.0; 4], dt, 3000, |_| [command, 0.0]);
    println!("5 deg pitch-attitude step:");
    println!(
        "{:>7} {:>11} {:>14}",
        "t [s]", "theta [deg]", "elevator [deg]"
    );
    for t in [0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 15.0] {
        let s = &hist[(t / dt) as usize];
        let u = command + k_q * s.x[2] + k_theta * s.x[3];
        println!(
            "{t:>7.2} {:>11.2} {:>14.2}",
            s.x[3].to_degrees(),
            u.to_degrees()
        );
    }
    let overshoot = hist.iter().map(|s| s.x[3]).fold(0.0, f64::max) / theta_cmd - 1.0;
    println!(
        "overshoot {:.1} %; a proportional hold leaves a steady-state error (add an integrator or autotrim to remove it)",
        100.0 * overshoot.max(0.0)
    );
}
