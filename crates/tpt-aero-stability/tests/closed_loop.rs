//! Closed-loop control on the linear plants: a yaw damper and a pitch
//! attitude hold, analysed through eigenvalues and simulated step responses.

use tpt_aero_core::{Vec3, WingPlanform, WingSection};
use tpt_aero_stability::controls::{control_derivatives, ControlDerivatives};
use tpt_aero_stability::lateral::{lateral_fd, LateralParameters, LateralPlant};
use tpt_aero_stability::longitudinal::longitudinal_fd;
use tpt_aero_stability::modes::Complex;
use tpt_aero_stability::phugoid::{
    drag_alpha_fd, FullLongitudinalParameters, FullLongitudinalPlant,
};
use tpt_aero_stability::sim::{simulate, StateSpace};
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

fn lateral_system() -> StateSpace<4, 3> {
    let m = aircraft(3.0);
    let alpha = 3.0_f64.to_radians();
    let d = lateral_fd(&m, alpha, V, RHO).unwrap();
    let p = LateralParameters {
        v_inf: V,
        span: m.reference().span,
        area: m.reference().area,
        density: RHO,
        mass: MASS,
        inertia_x: 1000.0,
        inertia_z: 1800.0,
        inertia_xz: 0.0,
    };
    let plant = LateralPlant::from_derivatives(&d, &p).unwrap();
    let aileron: ControlDerivatives = control_derivatives(&m, "aileron", alpha, V, RHO).unwrap();
    let rudder = control_derivatives(&m, "rudder", alpha, V, RHO).unwrap();
    plant.state_space(&p, &aileron, &rudder)
}

fn longitudinal_system() -> StateSpace<4, 2> {
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
    let p = FullLongitudinalParameters {
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
    let plant = FullLongitudinalPlant::from_derivatives(&d, &p).unwrap();
    plant.state_space(
        &p,
        &control_derivatives(&m, "elevator", t.alpha_rad, V, RHO).unwrap(),
    )
}

fn max_real(eig: &[Complex; 4]) -> f64 {
    eig.iter().map(|z| z.re).fold(f64::MIN, f64::max)
}

fn peak_abs(hist: &[tpt_aero_stability::sim::Sample<4>], state: usize, t0: f64, t1: f64) -> f64 {
    hist.iter()
        .filter(|s| s.t >= t0 && s.t <= t1)
        .map(|s| s.x[state].abs())
        .fold(0.0, f64::max)
}

#[test]
fn yaw_damper_stabilizes_the_dutch_roll_and_the_spiral() {
    let open = lateral_system();
    // Open loop: lightly damped Dutch roll and a (slightly) divergent spiral.
    assert!(max_real(&open.eigenvalues()) > 0.0);

    // Rudder command (positive yaws right) opposing the yaw rate.
    let closed = open.with_feedback(1, &[0.0, 0.0, 0.5, 0.0]);
    assert!(
        max_real(&closed.eigenvalues()) < 0.0,
        "{:?}",
        closed.eigenvalues()
    );
    // The Dutch-roll oscillation is gone or at least heavily damped.
    let ev = closed.eigenvalues();
    for z in ev.iter().filter(|z| z.im > 0.0) {
        let zeta = -z.re / z.abs();
        assert!(zeta > 0.7, "remaining oscillation has zeta {zeta}");
    }
}

#[test]
fn yaw_damper_suppresses_the_yaw_oscillation_in_time() {
    let open = lateral_system();
    let closed = open.with_feedback(1, &[0.0, 0.0, 0.5, 0.0]);
    let kick = [0.05, 0.0, 0.0, 0.0]; // 2.9 deg of sideslip
    let run = |sys: &StateSpace<4, 3>| simulate(sys, kick, 0.005, 2400, |_| [0.0; 3]);
    let (h_open, h_closed) = (run(&open), run(&closed));
    // After 3 s the open-loop yaw rate still rings; the damped one is quiet.
    let r_open = peak_abs(&h_open, 2, 3.0, 6.0);
    let r_closed = peak_abs(&h_closed, 2, 3.0, 6.0);
    assert!(
        r_closed < 0.15 * r_open,
        "closed {r_closed} vs open {r_open}"
    );
}

#[test]
fn pitch_attitude_hold_damps_the_short_period_and_replaces_the_phugoid() {
    let open = longitudinal_system();
    let ev = open.eigenvalues();
    let open_zeta_sp = ev
        .iter()
        .filter(|z| z.im > 0.0)
        .map(|z| (z.abs(), -z.re / z.abs()))
        .fold((0.0, 0.0), |a, b| if b.0 > a.0 { b } else { a })
        .1;

    // Elevator command (positive = nose down): u = kth theta + kq q.
    let closed = open.with_feedback(0, &[0.0, 0.0, -0.3, -2.0]);
    let cev = closed.eigenvalues();
    assert!(max_real(&cev) < 0.0, "{cev:?}");
    let closed_zeta_sp = cev
        .iter()
        .filter(|z| z.im > 0.0)
        .map(|z| (z.abs(), -z.re / z.abs()))
        .fold((0.0, 0.0), |a, b| if b.0 > a.0 { b } else { a })
        .1;
    assert!(closed_zeta_sp > open_zeta_sp);
    assert!(closed_zeta_sp > 0.8, "short period zeta {closed_zeta_sp}");
    // The lightly damped phugoid pair no longer exists as an oscillation at
    // its old frequency: any remaining slow oscillation is well damped.
    for z in cev.iter().filter(|z| z.im > 0.0 && z.abs() < 2.0) {
        assert!(-z.re / z.abs() > 0.5, "{z:?}");
    }
}

#[test]
fn pitch_attitude_hold_tracks_a_step_command() {
    let open = longitudinal_system();
    let (kth, kq) = (2.0, 0.3);
    let closed = open.with_feedback(0, &[0.0, 0.0, -kq, -kth]);
    let theta_cmd = 5.0_f64.to_radians();
    // u = r - K x with K = [0, 0, -kq, -kth]: choose r = -kth theta_cmd so that
    // u = kth (theta - theta_cmd) + kq q.
    let r = -kth * theta_cmd;
    let hist = simulate(&closed, [0.0; 4], 0.005, 2400, |_| [r, 0.0]);

    let theta_at = |t: f64| hist[(t / 0.005) as usize].x[3];
    let overshoot = hist.iter().map(|s| s.x[3]).fold(0.0, f64::max) / theta_cmd - 1.0;
    // Gain-limited proportional hold: the loop closes most of the error
    // quickly and the response is well damped.
    assert!(overshoot < 0.1, "overshoot {overshoot}");
    // Within 3% of the command after 2 s, 8% after 8 s: constant thrust lets
    // the airspeed bleed off slowly, so a pure proportional hold drifts.
    let early = theta_at(2.0);
    assert!(
        (early / theta_cmd - 1.0).abs() < 0.03,
        "theta(2 s) = {:.2} deg",
        early.to_degrees()
    );
    let settled = theta_at(8.0);
    assert!(
        (settled / theta_cmd - 1.0).abs() < 0.08,
        "theta(8 s) = {:.2} deg",
        settled.to_degrees()
    );
    // Monotone-ish approach: no sustained oscillation after 3 s.
    let tail: Vec<f64> = hist.iter().filter(|s| s.t > 3.0).map(|s| s.x[3]).collect();
    let (lo, hi) = tail
        .iter()
        .fold((f64::MAX, f64::MIN), |(a, b), x| (a.min(*x), b.max(*x)));
    assert!(hi - lo < 0.5 * theta_cmd, "tail range {}", hi - lo);
}
