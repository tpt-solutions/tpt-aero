//! JSON round trips for the serializable result and parameter types.
#![cfg(feature = "serde")]

use tpt_aero::prelude::*;

fn round_trip<T>(value: &T) -> T
where
    T: serde::Serialize + serde::de::DeserializeOwned,
{
    let json = serde_json::to_string(value).expect("serializes");
    serde_json::from_str(&json).expect("deserializes")
}

fn wing_model() -> VlmModel {
    VlmModel::from_surfaces(&[
        Surface::horizontal(WingPlanform::rectangular(10.0, 1.2), 12, 3),
        Surface::horizontal(WingPlanform::rectangular(4.0, 0.6), 6, 3)
            .at(4.0, 0.0, 0.0)
            .with_control("elevator", 0.6, 0.0, 2.0, ControlSymmetry::Symmetric),
    ])
    .expect("valid configuration")
}

#[test]
fn derivative_sets_round_trip_exactly() {
    let m = wing_model();
    let alpha = 3.0_f64.to_radians();
    let long = longitudinal_fd(&m, alpha, 50.0, 1.225, 1e-4).expect("solve");
    assert_eq!(round_trip(&long), long);
    let lat = lateral_fd(&m, alpha, 50.0, 1.225).unwrap();
    assert_eq!(round_trip(&lat), lat);
    let ctl = control_derivatives(&m, "elevator", alpha, 50.0, 1.225).unwrap();
    assert_eq!(round_trip(&ctl), ctl);
}

#[test]
fn json_field_names_are_the_documented_ones() {
    let m = wing_model();
    let long = longitudinal_fd(&m, 0.05, 50.0, 1.225, 1e-4).expect("solve");
    let json = serde_json::to_value(long).unwrap();
    for key in [
        "cl_alpha", "cm_alpha", "cl_q", "cm_q", "cl_delta", "cm_delta",
    ] {
        assert!(
            json.get(key).is_some_and(|v| v.is_number()),
            "missing {key}"
        );
    }
}

#[test]
fn solutions_and_trim_round_trip() {
    let m = wing_model();
    let sol = m.solve(0.08, 50.0, 1.225, &[]).unwrap();
    let back = round_trip(&sol);
    assert_eq!(back.cl, sol.cl);
    assert_eq!(back.gamma, sol.gamma);
    assert_eq!(back.force, sol.force);
    assert_eq!(back.strip_cl.len(), sol.strip_cl.len());

    let state = m
        .prepare()
        .unwrap()
        .solve_state(
            &FlightState {
                alpha: 0.05,
                beta: 0.02,
                p: 0.1,
                ..FlightState::default()
            },
            50.0,
            1.225,
        )
        .unwrap();
    assert_eq!(round_trip(&state), state);

    let trim = trim_level_flight(
        &m,
        &TrimRequest {
            weight_n: 5886.0,
            v_inf: 50.0,
            density: 1.225,
            elevator: "elevator",
            cd0: 0.03,
        },
    )
    .unwrap();
    assert_eq!(round_trip(&trim), trim);
}

#[test]
fn geometry_and_atmosphere_serialize() {
    let v = Vec3::new(1.0, -2.5, 3.25);
    assert_eq!(round_trip(&v), v);
    let section = WingSection {
        le: v,
        chord: 1.1,
        twist: 0.02,
    };
    assert_eq!(round_trip(&section), section);
    // A planform serializes (for archiving a case) but is deliberately not
    // deserializable: loading must go through the validating constructors.
    let json = serde_json::to_value(WingPlanform::rectangular(8.0, 1.0)).unwrap();
    assert_eq!(json["sections"].as_array().map(Vec::len), Some(2));

    let isa = isa::at(5_000.0).unwrap();
    assert_eq!(round_trip(&isa), isa);
}

#[test]
fn rotor_results_serialize() {
    use tpt_aero::rotor::bemt::{BemtSolver, OperatingCondition};
    use tpt_aero::rotor::blade::{RotorGeometry, SectionModel};
    use tpt_aero::rotor::loss::{HubLoss, TipLoss};

    let geometry = RotorGeometry::rectangular(0.127, 0.0127, 2, 0.020, 12)
        .with_linear_twist(25.0_f64.to_radians(), 5.0_f64.to_radians());
    let solver = BemtSolver::new(
        geometry,
        SectionModel::Linear {
            cl_slope: 2.0 * core::f64::consts::PI,
            cl_max: 1.2,
            cd0: 0.012,
            cd_cl2: 0.02,
        },
        TipLoss::Prandtl,
        HubLoss::Prandtl { x_hub: 0.1 },
    )
    .expect("valid geometry");
    let op = OperatingCondition::hover(800.0, 1.225);
    assert_eq!(round_trip(&op), op);

    let sol = solver.solve(&op).expect("hover solved");
    let back = round_trip(&sol);
    assert_eq!(back.loads, sol.loads);
    assert_eq!(back.sections.len(), sol.sections.len());
    assert_eq!(back.iterations, sol.iterations);
}
