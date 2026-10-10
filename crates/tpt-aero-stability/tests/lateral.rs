//! Lateral-directional derivative checks against physical expectations and
//! the independent longitudinal extraction.

use tpt_aero_core::WingPlanform;
use tpt_aero_stability::lateral::lateral_fd;
use tpt_aero_stability::longitudinal::longitudinal_fd;
use tpt_aero_vlm::{CamberModel, FlightState, VlmError, VlmModel};

const V: f64 = 50.0;
const RHO: f64 = 1.225;
const ALPHA: f64 = 5.0 * core::f64::consts::PI / 180.0;

fn model(w: &WingPlanform) -> VlmModel {
    VlmModel::from_planform(w, 16, 4, CamberModel::None).expect("valid planform")
}

fn dihedral_wing(deg: f64) -> WingPlanform {
    WingPlanform::tapered(10.0, 1.0, 1.0, 0.0, deg.to_radians())
}

#[test]
fn flat_symmetric_wing_has_no_sideslip_roll_or_yaw_response() {
    let d = lateral_fd(&model(&WingPlanform::rectangular(10.0, 1.0)), ALPHA, V, RHO).unwrap();
    // No dihedral or sweep: a straight bound-vortex lattice cannot make roll
    // or yaw moment from sideslip, nor side force from rates.
    assert!(d.cl_beta.abs() < 1e-9, "Clb = {}", d.cl_beta);
    assert!(d.cn_beta.abs() < 1e-9, "Cnb = {}", d.cn_beta);
    assert!(d.cy_p.abs() < 1e-9, "Cyp = {}", d.cy_p);
    assert!(d.cy_r.abs() < 1e-9, "Cyr = {}", d.cy_r);
}

#[test]
fn dihedral_effect_is_negative_and_linear_in_dihedral() {
    let d6 = lateral_fd(&model(&dihedral_wing(6.0)), ALPHA, V, RHO).unwrap();
    let d12 = lateral_fd(&model(&dihedral_wing(12.0)), ALPHA, V, RHO).unwrap();
    assert!(d6.cl_beta < -0.05, "Clb(6 deg) = {}", d6.cl_beta);
    let ratio = d12.cl_beta / d6.cl_beta;
    assert!((1.8..=2.2).contains(&ratio), "doubling dihedral: {ratio}");
    // Wind from the right pushes the airplane to the left: Cy_beta < 0.
    assert!(d6.cy_beta < 0.0, "Cyb = {}", d6.cy_beta);
}

#[test]
fn sweepback_adds_dihedral_effect() {
    let straight = lateral_fd(&model(&dihedral_wing(0.0)), ALPHA, V, RHO).unwrap();
    let swept = WingPlanform::tapered(10.0, 1.0, 1.0, 5.0 * 30f64.to_radians().tan(), 0.0);
    let swept = lateral_fd(&model(&swept), ALPHA, V, RHO).unwrap();
    assert!(
        swept.cl_beta < straight.cl_beta - 0.03,
        "swept {} vs straight {}",
        swept.cl_beta,
        straight.cl_beta
    );
}

#[test]
fn roll_damping_is_negative_with_textbook_magnitude() {
    let rect = lateral_fd(&model(&WingPlanform::rectangular(10.0, 1.0)), ALPHA, V, RHO).unwrap();
    // Rectangular AR 10: published roll damping is about -0.5 to -0.65.
    assert!((-0.75..=-0.45).contains(&rect.cl_p), "Clp = {}", rect.cl_p);
    // Taper moves lift inboard and reduces damping.
    let tapered = WingPlanform::tapered(10.0, 1.33, 0.67, 0.0, 0.0);
    let tapered = lateral_fd(&model(&tapered), ALPHA, V, RHO).unwrap();
    assert!(tapered.cl_p < 0.0 && tapered.cl_p > rect.cl_p);
    // A lifting wing yaws against a roll (adverse) and rolls with yaw rate.
    assert!(rect.cn_p < 0.0, "Cnp = {}", rect.cn_p);
    assert!(rect.cl_r > 0.0, "Clr = {}", rect.cl_r);
}

#[test]
fn pitch_rate_derivatives_agree_with_incidence_method() {
    // Two independent routes to q-derivatives: the longitudinal extraction
    // perturbs per-element incidence, `solve_state` applies -omega x r.
    let m = model(&WingPlanform::rectangular(10.0, 1.0));
    let alpha = 2f64.to_radians();
    let old = longitudinal_fd(&m, alpha, V, RHO, 1e-4);
    let prepared = m.prepare().unwrap();
    let c = m.reference().chord;
    let qhat = 1e-3;
    let q = qhat * 2.0 * V / c;
    let at = |q: f64| {
        prepared
            .solve_state(
                &FlightState {
                    alpha,
                    q,
                    ..FlightState::default()
                },
                V,
                RHO,
            )
            .unwrap()
    };
    let (hi, lo) = (at(q), at(-q));
    let cm_q = (hi.c_pitch - lo.c_pitch) / (2.0 * qhat);
    let cl_q = (hi.c_lift - lo.c_lift) / (2.0 * qhat);
    assert!(
        (cm_q / old.cm_q - 1.0).abs() < 0.01,
        "Cmq {cm_q} vs {}",
        old.cm_q
    );
    assert!(
        (cl_q / old.cl_q - 1.0).abs() < 0.01,
        "CLq {cl_q} vs {}",
        old.cl_q
    );
}

#[test]
fn zero_rates_reproduce_the_plain_solve() {
    let w = WingPlanform::tapered(9.0, 1.6, 0.7, 0.6, 4f64.to_radians());
    let m = VlmModel::from_planform(&w, 10, 4, CamberModel::Naca4 { m: 0.03, p: 0.45 }).unwrap();
    let prepared = m.prepare().unwrap();
    let plain = m.solve(ALPHA, V, RHO, &[]).unwrap();
    let state = prepared
        .solve_state(
            &FlightState {
                alpha: ALPHA,
                ..FlightState::default()
            },
            V,
            RHO,
        )
        .unwrap();
    assert!((state.c_lift - plain.cl).abs() < 1e-9);
    assert!((state.c_drag - plain.cdi).abs() < 1e-9);
    assert!((state.c_pitch - plain.cm).abs() < 1e-9);
    assert!(state.c_side.abs() < 1e-9 && state.c_roll.abs() < 1e-9 && state.c_yaw.abs() < 1e-9);
}

#[test]
fn sideslip_response_has_mirror_symmetry() {
    let m = model(&dihedral_wing(6.0));
    let prepared = m.prepare().unwrap();
    let at = |beta: f64| {
        prepared
            .solve_state(
                &FlightState {
                    alpha: ALPHA,
                    beta,
                    ..FlightState::default()
                },
                V,
                RHO,
            )
            .unwrap()
    };
    let (pos, neg) = (at(0.08), at(-0.08));
    // Lift, drag and pitch are even in beta; side force, roll and yaw odd.
    assert!((pos.c_lift - neg.c_lift).abs() < 1e-9);
    assert!((pos.c_pitch - neg.c_pitch).abs() < 1e-9);
    assert!((pos.c_side + neg.c_side).abs() < 1e-9);
    assert!((pos.c_roll + neg.c_roll).abs() < 1e-9);
    assert!((pos.c_yaw + neg.c_yaw).abs() < 1e-9);
    assert!(pos.c_roll < 0.0, "wind from the right rolls the wing left");
}

#[test]
fn invalid_state_is_rejected() {
    let m = model(&WingPlanform::rectangular(8.0, 1.0));
    let prepared = m.prepare().unwrap();
    let bad = FlightState {
        beta: f64::NAN,
        ..FlightState::default()
    };
    assert!(matches!(
        prepared.solve_state(&bad, V, RHO),
        Err(VlmError::InvalidFlightCondition)
    ));
    assert!(matches!(
        prepared.solve_state(&FlightState::default(), 0.0, RHO),
        Err(VlmError::InvalidFlightCondition)
    ));
}
