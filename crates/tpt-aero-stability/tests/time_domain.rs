//! Time-domain responses of the linear plants, checked against the modes
//! and steady states the analysis predicts.

use tpt_aero_core::{Vec3, WingPlanform, WingSection};
use tpt_aero_stability::controls::{control_derivatives, ControlDerivatives};
use tpt_aero_stability::lateral::{lateral_fd, LateralParameters, LateralPlant};
use tpt_aero_stability::longitudinal::longitudinal_fd;
use tpt_aero_stability::modes::Complex;
use tpt_aero_stability::phugoid::{
    drag_alpha_fd, FullLongitudinalParameters, FullLongitudinalPlant,
};
use tpt_aero_stability::sim::{simulate, Sample};
use tpt_aero_stability::trim::{trim_level_flight, TrimRequest};
use tpt_aero_vlm::{ControlSymmetry, Surface, VlmModel};

const V: f64 = 50.0;
const RHO: f64 = 1.225;
const MASS: f64 = 600.0;
const G: f64 = 9.80665;

fn aircraft(dihedral_deg: f64) -> VlmModel {
    let section = |y: f64, chord: f64| WingSection {
        le: Vec3::new(0.0, y, 0.0),
        chord,
        twist: 0.0,
    };
    let fin = WingPlanform::from_sections(vec![section(0.0, 0.9), section(1.5, 0.5)]).unwrap();
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
    .unwrap()
}

struct Longitudinal {
    plant: FullLongitudinalPlant,
    params: FullLongitudinalParameters,
    elevator: ControlDerivatives,
}

fn longitudinal() -> Longitudinal {
    let m = aircraft(0.0);
    let t = trim_level_flight(
        &m,
        &TrimRequest {
            weight_n: MASS * G,
            v_inf: V,
            density: RHO,
            elevator: "elevator",
            cd0: 0.03,
        },
    )
    .unwrap();
    let r = m.reference();
    let params = FullLongitudinalParameters {
        v_inf: V,
        chord: r.chord,
        area: r.area,
        density: RHO,
        mass: MASS,
        inertia_y: 900.0,
        cl: t.cl,
        cd: t.cdi + 0.03,
        cd_alpha: drag_alpha_fd(&m, t.alpha_rad, V, RHO).unwrap(),
    };
    let d = longitudinal_fd(&m, t.alpha_rad, V, RHO, 1e-4);
    Longitudinal {
        plant: FullLongitudinalPlant::from_derivatives(&d, &params).unwrap(),
        params,
        elevator: control_derivatives(&m, "elevator", t.alpha_rad, V, RHO).unwrap(),
    }
}

/// Times of upward zero crossings of `value`, linearly interpolated.
fn upward_zero_crossings<const N: usize>(
    hist: &[Sample<N>],
    value: impl Fn(&Sample<N>) -> f64,
) -> Vec<f64> {
    hist.windows(2)
        .filter_map(|w| {
            let (a, b) = (value(&w[0]), value(&w[1]));
            (a < 0.0 && b >= 0.0).then(|| w[0].t + (w[1].t - w[0].t) * (-a) / (b - a))
        })
        .collect()
}

fn mean_spacing(times: &[f64]) -> f64 {
    (times[times.len() - 1] - times[0]) / (times.len() - 1) as f64
}

/// The complex pair of `eig` with the smallest (or largest) modulus.
fn pair(eig: &[Complex; 4], slowest: bool) -> Complex {
    let mut pairs: Vec<Complex> = eig.iter().copied().filter(|z| z.im > 0.0).collect();
    pairs.sort_by(|a, b| a.abs().partial_cmp(&b.abs()).unwrap());
    if slowest {
        pairs[0]
    } else {
        pairs[pairs.len() - 1]
    }
}

#[test]
fn elevator_step_starts_with_the_input_matrix_response() {
    let l = longitudinal();
    let ss = l.plant.state_space(&l.params, &l.elevator);
    let delta = 0.02; // 1.15 deg trailing edge down
    let hist = simulate(&ss, [0.0; 4], 1e-6, 5, |_| [delta, 0.0]);
    let q_s = 0.5 * RHO * V * V * l.params.area;
    let q_dot = q_s * l.params.chord / l.params.inertia_y * l.elevator.cm * delta;
    let alpha_dot = -q_s / (MASS * V) * l.elevator.cl * delta;
    let t = hist.last().unwrap().t;
    assert!((hist.last().unwrap().x[2] - q_dot * t).abs() < 1e-3 * (q_dot * t).abs());
    assert!((hist.last().unwrap().x[1] - alpha_dot * t).abs() < 1e-2 * (alpha_dot * t).abs());
    // Trailing-edge-down elevator pitches the nose down and lifts the tail.
    assert!(q_dot < 0.0 && alpha_dot < 0.0);
}

#[test]
fn simulated_phugoid_has_the_period_and_decay_of_its_eigenvalue() {
    let l = longitudinal();
    let ss = l.plant.state_space(&l.params, &l.elevator);
    let phugoid = pair(&l.plant.eigenvalues(), true);
    let period = 2.0 * core::f64::consts::PI / phugoid.im;

    // 3 m/s speed perturbation, no input, 150 s.
    let hist = simulate(&ss, [3.0, 0.0, 0.0, 0.0], 0.01, 15_000, |_| [0.0, 0.0]);
    let crossings = upward_zero_crossings(&hist, |s| s.x[0]);
    assert!(crossings.len() >= 4, "{} crossings", crossings.len());
    let measured = mean_spacing(&crossings);
    assert!(
        (measured / period - 1.0).abs() < 0.02,
        "period {measured:.2} s vs eigenvalue {period:.2} s"
    );

    // Amplitude decay per cycle follows exp(re * period).
    let peak = |lo: f64, hi: f64| {
        hist.iter()
            .filter(|s| s.t >= lo && s.t < hi)
            .map(|s| s.x[0].abs())
            .fold(0.0, f64::max)
    };
    let (first, later) = (peak(0.0, period), peak(3.0 * period, 4.0 * period));
    let expected = (phugoid.re * 3.0 * period).exp();
    assert!(
        ((later / first) / expected - 1.0).abs() < 0.12,
        "decay {} vs {expected}",
        later / first
    );
}

#[test]
fn a_sustained_updraft_leaves_the_airplane_aligned_with_the_relative_wind() {
    let l = longitudinal();
    let ss = l.plant.state_space(&l.params, &l.elevator);
    let w_gust = 2.0;
    // Long enough for the lightly damped phugoid to die out (~55 s time constant).
    let hist = simulate(&ss, [0.0; 4], 0.02, 30_000, |_| [0.0, w_gust]);
    let end = hist.last().unwrap().x;
    // Equilibrium: no aerodynamic angle of attack, so alpha = -w/V; no
    // residual speed change, pitch rate or pitch attitude.
    assert!(
        (end[1] + w_gust / V).abs() < 0.02 * w_gust / V,
        "alpha {}",
        end[1]
    );
    assert!(end[0].abs() < 0.05, "u {}", end[0]);
    assert!(end[2].abs() < 1e-4 && end[3].abs() < 1e-3, "{end:?}");
}

#[test]
fn integration_converges_with_the_step_size() {
    let l = longitudinal();
    let ss = l.plant.state_space(&l.params, &l.elevator);
    let input = |t: f64| [0.01 * (3.0 * t).sin(), 1.0 * (0.5 * t).cos()];
    let coarse = simulate(&ss, [0.0; 4], 0.02, 500, input);
    let fine = simulate(&ss, [0.0; 4], 0.01, 1000, input);
    let (a, b) = (coarse.last().unwrap(), fine.last().unwrap());
    assert!((a.t - b.t).abs() < 1e-12);
    for i in 0..4 {
        assert!(
            (a.x[i] - b.x[i]).abs() < 1e-6,
            "state {i}: {} vs {}",
            a.x[i],
            b.x[i]
        );
    }
}

struct Lateral {
    plant: LateralPlant,
    params: LateralParameters,
    aileron: ControlDerivatives,
    rudder: ControlDerivatives,
    cl_p: f64,
}

fn lateral() -> Lateral {
    let m = aircraft(3.0);
    let alpha = 3.0_f64.to_radians();
    let d = lateral_fd(&m, alpha, V, RHO).unwrap();
    let params = LateralParameters {
        v_inf: V,
        span: m.reference().span,
        area: m.reference().area,
        density: RHO,
        mass: MASS,
        inertia_x: 1000.0,
        inertia_z: 1800.0,
        inertia_xz: 0.0,
    };
    Lateral {
        plant: LateralPlant::from_derivatives(&d, &params).unwrap(),
        params,
        aileron: control_derivatives(&m, "aileron", alpha, V, RHO).unwrap(),
        rudder: control_derivatives(&m, "rudder", alpha, V, RHO).unwrap(),
        cl_p: d.cl_p,
    }
}

#[test]
fn simulated_dutch_roll_matches_its_eigenvalue() {
    let l = lateral();
    let ss = l.plant.state_space(&l.params, &l.aileron, &l.rudder);
    let dutch_roll = pair(&l.plant.eigenvalues(), false);
    let period = 2.0 * core::f64::consts::PI / dutch_roll.im;

    let hist = simulate(&ss, [0.05, 0.0, 0.0, 0.0], 0.005, 4000, |_| [0.0; 3]);
    let crossings = upward_zero_crossings(&hist, |s| s.x[2]);
    assert!(crossings.len() >= 3, "{} crossings", crossings.len());
    let measured = mean_spacing(&crossings);
    assert!(
        (measured / period - 1.0).abs() < 0.03,
        "period {measured:.3} vs {period:.3}"
    );
}

#[test]
fn aileron_step_settles_to_the_roll_rate_set_by_roll_damping() {
    let l = lateral();
    let ss = l.plant.state_space(&l.params, &l.aileron, &l.rudder);
    let delta = 0.02;
    let hist = simulate(&ss, [0.0; 4], 0.002, 500, |_| [delta, 0.0, 0.0]);
    let p_at_1s = hist.last().unwrap().x[1];
    // Positive aileron command rolls left (negative p). Steady roll rate from
    // the roll-damping balance: p = -(Cl_delta delta) / (Cl_p b / 2V).
    let p_ss = -(l.aileron.c_roll * delta) / (l.cl_p * l.params.span / (2.0 * V));
    assert!(p_at_1s < 0.0 && p_ss < 0.0);
    assert!(
        (p_at_1s / p_ss - 1.0).abs() < 0.15,
        "p(1 s) {p_at_1s:.4} vs balance {p_ss:.4}"
    );
}

#[test]
fn rudder_step_yaws_the_nose_right_and_sideslips_into_the_turn() {
    let l = lateral();
    let ss = l.plant.state_space(&l.params, &l.aileron, &l.rudder);
    let hist = simulate(&ss, [0.0; 4], 0.005, 400, |_| [0.0, 0.02, 0.0]);
    let early = &hist[40]; // 0.2 s
    assert!(
        early.x[2] > 0.0,
        "positive rudder command yaws right: r = {}",
        early.x[2]
    );
    // Sideslip builds in response: the weathercock restores it, so it stays
    // small relative to the yaw excursion.
    assert!(hist.iter().all(|s| s.x.iter().all(|v| v.is_finite())));
}

/// Peak load-factor increment for a 1-cos gust of gradient distance `h` (m).
fn peak_gust_load(l: &Longitudinal, u_ds: f64, h: f64) -> f64 {
    let ss = l.plant.state_space(&l.params, &l.elevator);
    let gust = |t: f64| {
        let x = V * t;
        if x <= 0.0 || x >= 2.0 * h {
            0.0
        } else {
            0.5 * u_ds * (1.0 - (core::f64::consts::PI * x / h).cos())
        }
    };
    let dt = 0.002;
    let hist = simulate(
        &ss,
        [0.0; 4],
        dt,
        ((2.0 * h / V + 6.0) / dt) as usize,
        |t| [0.0, gust(t)],
    );
    hist.iter()
        .map(|s| {
            let xdot = ss.derivative(&s.x, &[0.0, gust(s.t)]);
            (V * (s.x[2] - xdot[1]) / G).abs()
        })
        .fold(0.0, f64::max)
}

#[test]
fn gust_alleviation_matches_the_far25_factor_at_the_design_gradient() {
    let l = longitudinal();
    let wing_loading = MASS * G / l.params.area;
    let cl_alpha = {
        // Lift-curve slope entering the quasi-static estimate.
        let m = aircraft(0.0);
        longitudinal_fd(&m, 3.0_f64.to_radians(), V, RHO, 1e-4).cl_alpha
    };
    let u_ds = 7.6;
    let quasi_static = RHO * V * cl_alpha * u_ds / (2.0 * wing_loading);
    let mu = 2.0 * wing_loading / (RHO * l.params.chord * cl_alpha * G);
    let kg = 0.88 * mu / (5.3 + mu);

    // Design gradient: 12.5 chords.
    let ratio = peak_gust_load(&l, u_ds, 12.5 * l.params.chord) / quasi_static;
    assert!(
        (ratio / kg - 1.0).abs() < 0.15,
        "alleviation {ratio:.3} vs FAR 25 Kg {kg:.3}"
    );

    // Longer gradients are alleviated more: the aircraft heaves and
    // weathercocks with the air.
    let loads: Vec<f64> = [6.0, 15.0, 30.0, 60.0, 107.0]
        .iter()
        .map(|h| peak_gust_load(&l, u_ds, *h))
        .collect();
    assert!(loads.windows(2).all(|w| w[1] < w[0]), "{loads:?}");
    assert!(loads[4] < 0.1 * quasi_static);
}
