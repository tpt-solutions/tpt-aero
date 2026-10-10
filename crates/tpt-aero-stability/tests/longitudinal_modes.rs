//! Full longitudinal plant (phugoid + short period) and control derivatives.

use tpt_aero_core::{Vec3, WingPlanform, WingSection};
use tpt_aero_stability::controls::{control_derivatives, ControlError};
use tpt_aero_stability::longitudinal::{longitudinal_fd, LongitudinalDerivatives};
use tpt_aero_stability::phugoid::{
    drag_alpha_fd, FullLongitudinalParameters, FullLongitudinalPlant,
};
use tpt_aero_stability::trim::{trim_level_flight, TrimRequest};
use tpt_aero_vlm::{ControlSymmetry, Surface, VlmError, VlmModel};

const RHO: f64 = 1.225;
const MASS: f64 = 600.0;
const G: f64 = 9.80665;
const CD0: f64 = 0.03;

fn aircraft() -> VlmModel {
    let section = |y: f64, chord: f64| WingSection {
        le: Vec3::new(0.0, y, 0.0),
        chord,
        twist: 0.0,
    };
    let fin = WingPlanform::from_sections(vec![section(0.0, 0.9), section(1.5, 0.5)]).unwrap();
    VlmModel::from_surfaces(&[
        Surface::horizontal(WingPlanform::rectangular(10.0, 1.2), 20, 3).with_control(
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

struct Case {
    d: LongitudinalDerivatives,
    p: FullLongitudinalParameters,
}

fn trimmed_case(m: &VlmModel, v: f64) -> Case {
    let t = trim_level_flight(
        m,
        &TrimRequest {
            weight_n: MASS * G,
            v_inf: v,
            density: RHO,
            elevator: "elevator",
            cd0: CD0,
        },
    )
    .expect("trims");
    let r = m.reference();
    Case {
        d: longitudinal_fd(m, t.alpha_rad, v, RHO, 1e-4),
        p: FullLongitudinalParameters {
            v_inf: v,
            chord: r.chord,
            area: r.area,
            density: RHO,
            mass: MASS,
            inertia_y: 900.0,
            cl: t.cl,
            cd: t.cdi + CD0,
            cd_alpha: drag_alpha_fd(m, t.alpha_rad, v, RHO).unwrap(),
        },
    }
}

#[test]
fn phugoid_follows_lanchester_scaling_and_is_much_slower_than_the_short_period() {
    let m = aircraft();
    let mut wn_times_v = Vec::new();
    for v in [35.0, 50.0, 70.0] {
        let c = trimmed_case(&m, v);
        let plant = FullLongitudinalPlant::from_derivatives(&c.d, &c.p).unwrap();
        let modes = plant.modes().expect("two oscillatory modes");

        // Lanchester: wn = sqrt(2) g / V, zeta = 1 / (sqrt(2) L/D). The full
        // model sits within ~15% in frequency and ~35% in damping.
        let lanchester_wn = 2.0_f64.sqrt() * G / v;
        let lanchester_zeta = c.p.cd / (2.0_f64.sqrt() * c.p.cl);
        assert!(
            (modes.phugoid_wn / lanchester_wn - 1.0).abs() < 0.2,
            "V={v}: wn {} vs Lanchester {lanchester_wn}",
            modes.phugoid_wn
        );
        assert!(
            (modes.phugoid_zeta / lanchester_zeta - 1.0).abs() < 0.35,
            "V={v}: zeta {} vs Lanchester {lanchester_zeta}",
            modes.phugoid_zeta
        );
        // Lightly damped but stable; short period far faster and well damped.
        assert!(modes.phugoid_zeta > 0.0 && modes.phugoid_zeta < 0.3);
        assert!(modes.short_period_wn > 15.0 * modes.phugoid_wn);
        assert!((0.3..1.0).contains(&modes.short_period_zeta));
        wn_times_v.push(modes.phugoid_wn * v);
    }
    // Phugoid frequency scales as 1/V.
    let (lo, hi) = wn_times_v
        .iter()
        .fold((f64::MAX, f64::MIN), |(a, b), x| (a.min(*x), b.max(*x)));
    assert!(hi / lo < 1.05, "wn*V should be constant: {wn_times_v:?}");
}

#[test]
fn phugoid_damping_drops_with_better_glide_ratio() {
    let m = aircraft();
    let zeta = |cd0_extra: f64| {
        let mut c = trimmed_case(&m, 50.0);
        c.p.cd += cd0_extra;
        FullLongitudinalPlant::from_derivatives(&c.d, &c.p)
            .unwrap()
            .modes()
            .unwrap()
            .phugoid_zeta
    };
    // More drag (worse L/D) damps the phugoid more.
    assert!(zeta(0.03) > zeta(0.0));
}

#[test]
fn statically_unstable_airframe_has_no_oscillatory_pair() {
    let m = aircraft();
    let mut c = trimmed_case(&m, 50.0);
    c.d.cm_alpha = 3.0; // clearly nose-up divergent
    let plant = FullLongitudinalPlant::from_derivatives(&c.d, &c.p).unwrap();
    assert!(plant.modes().is_none());
    assert!(plant
        .eigenvalues()
        .iter()
        .any(|z| z.im == 0.0 && z.re > 0.5));
}

#[test]
fn full_plant_rejects_bad_parameters() {
    let c = trimmed_case(&aircraft(), 50.0);
    for broken in [
        FullLongitudinalParameters { mass: 0.0, ..c.p },
        FullLongitudinalParameters {
            v_inf: f64::NAN,
            ..c.p
        },
        FullLongitudinalParameters {
            inertia_y: -1.0,
            ..c.p
        },
        FullLongitudinalParameters {
            cd_alpha: f64::INFINITY,
            ..c.p
        },
    ] {
        assert!(FullLongitudinalPlant::from_derivatives(&c.d, &broken).is_none());
    }
}

#[test]
fn control_derivatives_have_the_expected_character() {
    let m = aircraft();
    let alpha = 3.0_f64.to_radians();
    let e = control_derivatives(&m, "elevator", alpha, 50.0, RHO).unwrap();
    assert!(e.cl > 0.3 && e.cm < -1.0, "{e:?}");
    assert!(e.cy.abs() < 1e-6 && e.c_roll.abs() < 1e-6 && e.c_yaw.abs() < 1e-6);

    let a = control_derivatives(&m, "aileron", alpha, 50.0, RHO).unwrap();
    assert!(a.c_roll < -0.2, "{a:?}");
    assert!(a.cl.abs() < 1e-6 && a.cm.abs() < 1e-6);
    assert!(a.c_yaw > 0.0, "adverse yaw sense: {a:?}");

    let r = control_derivatives(&m, "rudder", alpha, 50.0, RHO).unwrap();
    assert!(r.cy < -0.1 && r.c_yaw > 0.05, "{r:?}");
    assert!(r.cl.abs() < 1e-6 && r.cm.abs() < 1e-6);

    // Control effectiveness is nearly independent of angle of attack here.
    let e0 = control_derivatives(&m, "elevator", 0.0, 50.0, RHO).unwrap();
    assert!((e0.cm / e.cm - 1.0).abs() < 0.1);
}

#[test]
fn control_derivatives_report_errors() {
    let m = aircraft();
    assert_eq!(
        control_derivatives(&m, "flaperon", 0.05, 50.0, RHO).unwrap_err(),
        ControlError::NoSuchControl
    );
    assert_eq!(
        control_derivatives(&m, "elevator", 0.05, 0.0, RHO).unwrap_err(),
        ControlError::Solve(VlmError::InvalidFlightCondition)
    );
}
