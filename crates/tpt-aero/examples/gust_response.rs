//! Gust and turbulence response of a complete light aircraft: a 1-cos
//! discrete gust (peak load factor against the quasi-static estimate) and a
//! Dryden turbulence record (load-factor statistics).
//!
//! Linear short-period/phugoid model of wing + tailplane + fin; the numbers
//! show the workflow, not a certified gust-load analysis.
//!
//! Run: `cargo run --example gust_response`

use tpt_aero::atmo::{DiscreteGust, GustShape, TurbulenceField, TurbulenceModel, WhiteNoise};
use tpt_aero::core::WingSection;
use tpt_aero::prelude::*;
use tpt_aero::stability::sim::simulate;

const MASS: f64 = 600.0;
const G: f64 = 9.80665;

fn build_aircraft() -> VlmModel {
    let section = |y: f64, chord: f64| WingSection {
        le: Vec3::new(0.0, y, 0.0),
        chord,
        twist: 0.0,
    };
    let fin = WingPlanform::from_sections(vec![section(0.0, 0.9), section(1.5, 0.5)])
        .expect("two ordered sections");
    VlmModel::from_surfaces(&[
        Surface::horizontal(WingPlanform::rectangular(10.0, 1.2), 20, 3),
        Surface::horizontal(WingPlanform::rectangular(4.0, 0.6), 10, 3)
            .at(4.0, 0.0, 0.0)
            .with_control("elevator", 0.6, 0.0, 2.0, ControlSymmetry::Symmetric),
        Surface::vertical(fin, 6, 3).at(3.8, 0.0, 0.0),
    ])
    .expect("valid configuration")
}

fn main() {
    let (v, rho) = (50.0, 1.225);
    let model = build_aircraft();
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
    let derivs = longitudinal_fd(&model, trim.alpha_rad, v, rho, 1e-4);
    let elevator =
        control_derivatives(&model, "elevator", trim.alpha_rad, v, rho).expect("control exists");
    let r = model.reference();
    let params = FullLongitudinalParameters {
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
    let plant = FullLongitudinalPlant::from_derivatives(&derivs, &params).expect("valid plant");
    let ss = plant.state_space(&params, &elevator);

    // Load-factor increment from the state: dn = V (q - alpha') / g.
    let load_factor = |x: &[f64; 4], w_gust: f64| {
        let xdot = ss.derivative(x, &[0.0, w_gust]);
        v * (x[2] - xdot[1]) / G
    };

    println!(
        "=== Gust response: {MASS} kg light aircraft at {v} m/s (CL {:.2}, CL_alpha {:.2}) ===\n",
        trim.cl, derivs.cl_alpha
    );

    // 1-cos discrete gust (speed as a function of penetration x = V t), swept
    // over the gradient distance H. The design condition is H = 12.5 chords,
    // where FAR 25 gives the alleviation factor Kg = 0.88 mu / (5.3 + mu).
    let u_ds = 7.6;
    let wing_loading = MASS * G / r.area;
    let quasi_static = rho * v * derivs.cl_alpha * u_ds / (2.0 * wing_loading);
    let mu = 2.0 * wing_loading / (rho * r.chord * derivs.cl_alpha * G);
    println!(
        "1-cos gust, U_ds {u_ds} m/s. Quasi-static sharp-gust estimate {quasi_static:.2} g; \
         mass ratio mu {mu:.1}, FAR 25 Kg {:.2}",
        0.88 * mu / (5.3 + mu)
    );
    println!(
        "{:>9} {:>8} {:>8} {:>12} {:>10}",
        "H [m]", "chords", "2H/V [s]", "peak [g]", "ratio"
    );
    let dt = 0.002;
    for h in [6.0, 9.0, 15.0, 30.0, 60.0, 107.0] {
        let gust = DiscreteGust::new(u_ds, h, GustShape::OneMinusCosine).expect("valid gust");
        let duration = 2.0 * h / v + 6.0;
        let hist = simulate(&ss, [0.0; 4], dt, (duration / dt) as usize, |t| {
            [0.0, gust.speed_at(v * t)]
        });
        let peak = hist
            .iter()
            .map(|s| load_factor(&s.x, gust.speed_at(v * s.t)).abs())
            .fold(0.0_f64, f64::max);
        println!(
            "{h:>9.1} {:>8.1} {:>8.2} {peak:>12.2} {:>10.2}",
            h / r.chord,
            2.0 * h / v,
            peak / quasi_static
        );
    }

    let dt = 0.005;

    // Dryden vertical turbulence, moderate intensity at 500 m.
    let sigma_w = 1.5;
    let mut field = TurbulenceField::low_altitude(TurbulenceModel::Dryden, v, 500.0, sigma_w);
    let mut noise = WhiteNoise::new(2024);
    let t_end = 120.0;
    let n_steps = (t_end / dt) as usize;
    let mut w_series = Vec::with_capacity(n_steps + 2);
    for _ in 0..=n_steps + 1 {
        field.step(0.0, 0.0, noise.next_gaussian(), dt);
        w_series.push(field.values().2);
    }
    let w_at = |t: f64| w_series[((t / dt) as usize).min(w_series.len() - 1)];
    let hist = simulate(&ss, [0.0; 4], dt, n_steps, |t| [0.0, w_at(t)]);
    let loads: Vec<f64> = hist.iter().map(|s| load_factor(&s.x, w_at(s.t))).collect();
    let rms = (loads.iter().map(|n| n * n).sum::<f64>() / loads.len() as f64).sqrt();
    let peak = loads.iter().fold(0.0_f64, |m, n| m.max(n.abs()));
    let w_rms = (w_series.iter().map(|w| w * w).sum::<f64>() / w_series.len() as f64).sqrt();
    println!("\nDryden turbulence, sigma_w {sigma_w} m/s at 500 m, {t_end} s record:");
    println!("  realized w rms {w_rms:.2} m/s");
    println!("  load factor rms {rms:.3} g, peak {peak:.2} g");
}
