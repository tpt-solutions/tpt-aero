//! Property checks for the trim-request strategy: every generated request is
//! accepted by the trim solver, and every trim it returns satisfies the level
//! flight residuals it was asked for.

use proptest::prelude::*;
use tpt_aero_core::{Vec3, WingPlanform, WingSection};
use tpt_aero_stability::trim::{trim_level_flight, TrimError};
use tpt_aero_vlm::{ControlSymmetry, Surface, VlmModel};

/// The reference aircraft from the six-DOF tests: a wing with an aileron, and a
/// tailplane carrying a symmetric elevator.
fn aircraft() -> VlmModel {
    let wing_section = |y: f64, z: f64| WingSection {
        le: Vec3::new(0.0, y, z),
        chord: 1.2,
        twist: 0.0,
    };
    let dihedral_wing = WingPlanform::from_sections(vec![
        wing_section(-5.0, 0.5),
        wing_section(0.0, 0.0),
        wing_section(5.0, 0.5),
    ])
    .unwrap();
    VlmModel::from_surfaces(&[
        Surface::horizontal(dihedral_wing, 20, 3).with_control(
            "aileron",
            0.75,
            2.5,
            4.8,
            ControlSymmetry::Antisymmetric,
        ),
        Surface::horizontal(WingPlanform::rectangular(4.0, 0.6), 10, 3)
            .at(4.0, 0.0, 0.0)
            .with_control("elevator", 0.6, 0.0, 2.0, ControlSymmetry::Symmetric),
    ])
    .unwrap()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn generated_requests_are_never_rejected_as_invalid(p in tpt_aero_verify::trim_request()) {
        let model = aircraft();
        let req = p.request("elevator");
        // Any failure must be a solver outcome, never a rejected request.
        match trim_level_flight(&model, &req) {
            Err(TrimError::InvalidRequest) | Err(TrimError::NoSuchControl) => {
                prop_assert!(false, "generated request rejected: {p:?}");
            }
            _ => {}
        }
    }

    #[test]
    fn returned_trims_satisfy_level_flight(p in tpt_aero_verify::trim_request()) {
        let model = aircraft();
        let req = p.request("elevator");
        if let Ok(trim) = trim_level_flight(&model, &req) {
            let q_s = 0.5 * req.density * req.v_inf * req.v_inf * model.reference().area;
            let cl_required = req.weight_n / q_s;
            prop_assert!((trim.cl - cl_required).abs() < 1.0e-8 * cl_required.max(1.0),
                "cl {} vs required {cl_required} for {p:?}", trim.cl);
            prop_assert!(trim.alpha_rad.abs() <= tpt_aero_stability::trim::TRIM_ALPHA_LIMIT);
            prop_assert!(trim.thrust_n.is_finite() && trim.thrust_n >= 0.0);
        }
    }
}
