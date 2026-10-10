//! Level-flight trim and lateral-directional modes of a complete aircraft.

use tpt_aero_core::{Vec3, WingPlanform, WingSection};
use tpt_aero_stability::lateral::{lateral_fd, LateralParameters, LateralPlant};
use tpt_aero_stability::trim::{trim_level_flight, TrimError, TrimRequest};
use tpt_aero_vlm::{ControlSymmetry, ReferenceGeometry, Surface, VlmModel};

const RHO: f64 = 1.225;
const MASS: f64 = 600.0;
const G: f64 = 9.80665;

/// Wing + tailplane (with elevator) + fin, CG at `cg` chords.
fn aircraft(dihedral_deg: f64, cg: f64) -> VlmModel {
    let section = |y: f64, chord: f64| WingSection {
        le: Vec3::new(0.0, y, 0.0),
        chord,
        twist: 0.0,
    };
    let fin = WingPlanform::from_sections(vec![section(0.0, 0.9), section(1.5, 0.5)]).unwrap();
    let wing = WingPlanform::tapered(10.0, 1.2, 1.2, 0.0, dihedral_deg.to_radians());
    VlmModel::from_surfaces(&[
        Surface::horizontal(wing.clone(), 20, 3),
        Surface::horizontal(WingPlanform::rectangular(4.0, 0.6), 10, 3)
            .at(4.0, 0.0, 0.0)
            .with_control("elevator", 0.6, 0.0, 2.0, ControlSymmetry::Symmetric),
        Surface::vertical(fin, 6, 3).at(3.8, 0.0, 0.0),
    ])
    .unwrap()
    .with_reference(ReferenceGeometry {
        area: wing.area(),
        chord: wing.mac(),
        span: wing.span(),
        moment_center: Vec3::new(cg * wing.mac(), 0.0, 0.0),
    })
}

fn request(v_inf: f64) -> TrimRequest<'static> {
    TrimRequest {
        weight_n: MASS * G,
        v_inf,
        density: RHO,
        elevator: "elevator",
        cd0: 0.03,
    }
}

#[test]
fn trim_balances_lift_and_pitching_moment() {
    let m = aircraft(0.0, 0.25);
    let t = trim_level_flight(&m, &request(50.0)).expect("trims");
    let q_s = 0.5 * RHO * 50.0 * 50.0 * m.reference().area;
    assert!((t.cl - MASS * G / q_s).abs() < 1e-9);

    // Re-solve at the trim point independently.
    let index = m.control_index("elevator").unwrap();
    let mut defl = vec![0.0; m.controls().len()];
    defl[index] = t.elevator_rad;
    let s = m
        .solve(t.alpha_rad, 50.0, RHO, &m.control_twists(&defl).unwrap())
        .unwrap();
    assert!((s.cl - t.cl).abs() < 1e-8, "CL {} vs {}", s.cl, t.cl);
    assert!(s.cm.abs() < 1e-8, "Cm = {}", s.cm);
    assert!((t.thrust_n - (s.cdi + 0.03) * q_s).abs() < 1e-6);
    assert!(t.iterations <= 6, "{} iterations", t.iterations);
    // A sensible cruise trim: a few degrees of alpha, a little up-elevator.
    assert!((1.0..8.0).contains(&t.alpha_rad.to_degrees()));
    assert!(t.elevator_rad < 0.0 && t.elevator_rad.to_degrees() > -10.0);
}

#[test]
fn slower_flight_needs_more_alpha_and_more_up_elevator() {
    let m = aircraft(0.0, 0.25);
    let slow = trim_level_flight(&m, &request(35.0)).unwrap();
    let mid = trim_level_flight(&m, &request(50.0)).unwrap();
    let fast = trim_level_flight(&m, &request(70.0)).unwrap();
    assert!(slow.alpha_rad > mid.alpha_rad && mid.alpha_rad > fast.alpha_rad);
    assert!(slow.elevator_rad < mid.elevator_rad && mid.elevator_rad < fast.elevator_rad);
    // Induced drag grows steeply at low speed.
    assert!(slow.cdi > 4.0 * mid.cdi);
}

#[test]
fn aft_cg_needs_less_up_elevator() {
    let fwd = trim_level_flight(&aircraft(0.0, 0.15), &request(50.0)).unwrap();
    let aft = trim_level_flight(&aircraft(0.0, 0.55), &request(50.0)).unwrap();
    assert!(
        fwd.elevator_rad < aft.elevator_rad,
        "forward CG {} vs aft CG {}",
        fwd.elevator_rad,
        aft.elevator_rad
    );
}

#[test]
fn trim_reports_bad_requests() {
    let m = aircraft(0.0, 0.25);
    let mut bad = request(50.0);
    bad.weight_n = -1.0;
    assert_eq!(trim_level_flight(&m, &bad), Err(TrimError::InvalidRequest));
    let mut bad = request(50.0);
    bad.elevator = "flaperon";
    assert_eq!(trim_level_flight(&m, &bad), Err(TrimError::NoSuchControl));
    let mut bad = request(50.0);
    bad.cd0 = f64::NAN;
    assert_eq!(trim_level_flight(&m, &bad), Err(TrimError::InvalidRequest));
}

fn lateral_modes(dihedral_deg: f64) -> tpt_aero_stability::lateral::LateralModes {
    let m = aircraft(dihedral_deg, 0.25);
    let d = lateral_fd(&m, 3.0_f64.to_radians(), 50.0, RHO).unwrap();
    let p = LateralParameters {
        v_inf: 50.0,
        span: m.reference().span,
        area: m.reference().area,
        density: RHO,
        mass: MASS,
        inertia_x: 1000.0,
        inertia_z: 1800.0,
        inertia_xz: 0.0,
    };
    LateralPlant::from_derivatives(&d, &p)
        .expect("valid parameters")
        .modes()
        .expect("standard spectrum")
}

#[test]
fn lateral_modes_have_the_classical_character() {
    let modes = lateral_modes(0.0);
    // Dutch roll: lightly damped oscillation of a few rad/s.
    assert!((2.0..6.0).contains(&modes.dutch_roll_wn), "{modes:?}");
    assert!((0.05..0.4).contains(&modes.dutch_roll_zeta), "{modes:?}");
    // Roll subsidence: fast and stable.
    let tau = modes.roll_time_constant().expect("stable roll mode");
    assert!(tau < 0.5, "roll time constant {tau}");
    // Spiral: far slower than the other modes.
    assert!(modes.spiral_root.abs() < 0.1, "{modes:?}");
}

#[test]
fn dihedral_stabilizes_the_spiral_mode() {
    let flat = lateral_modes(0.0);
    let dihedral = lateral_modes(6.0);
    // The classic spiral trade-off: effective dihedral (Cl_beta) helps the
    // spiral at the expense of Dutch-roll damping.
    assert!(flat.spiral_root > 0.0, "no dihedral: divergent spiral");
    assert!(dihedral.spiral_root < flat.spiral_root);
    assert!(dihedral.spiral_root < 0.002, "{dihedral:?}");
    assert!(dihedral.dutch_roll_zeta < flat.dutch_roll_zeta);
    assert!(flat.spiral_time_to_double().is_some());
}

#[test]
fn plant_structure_and_validation() {
    let m = aircraft(3.0, 0.25);
    let d = lateral_fd(&m, 3.0_f64.to_radians(), 50.0, RHO).unwrap();
    let params = LateralParameters {
        v_inf: 50.0,
        span: m.reference().span,
        area: m.reference().area,
        density: RHO,
        mass: MASS,
        inertia_x: 1000.0,
        inertia_z: 1800.0,
        inertia_xz: 120.0,
    };
    let plant = LateralPlant::from_derivatives(&d, &params).unwrap();
    // Kinematics rows: phi' = p, and gravity couples bank into beta'.
    assert_eq!(plant.a[3], [0.0, 1.0, 0.0, 0.0]);
    assert!((plant.a[0][3] - G / 50.0).abs() < 1e-12);
    // The trace equals the sum of eigenvalues.
    let trace: f64 = (0..4).map(|i| plant.a[i][i]).sum();
    let eig_sum: f64 = plant.eigenvalues().iter().map(|z| z.re).sum();
    assert!((trace - eig_sum).abs() < 1e-8, "{trace} vs {eig_sum}");
    // Invalid mass properties are rejected.
    for broken in [
        LateralParameters {
            mass: 0.0,
            ..params
        },
        LateralParameters {
            inertia_xz: 5000.0,
            ..params
        },
        LateralParameters {
            v_inf: f64::NAN,
            ..params
        },
    ] {
        assert!(LateralPlant::from_derivatives(&d, &broken).is_none());
    }
}
