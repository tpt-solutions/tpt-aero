//! Steady coordinated turn and steady sideslip trim on the VLM.

use tpt_aero_core::{Vec3, WingPlanform, WingSection};
use tpt_aero_stability::trim::TrimError;
use tpt_aero_stability::turn::{coordinated_turn_rate, trim_steady_turn, TurnRequest, G};
use tpt_aero_vlm::{ControlSymmetry, Surface, VlmModel};

const RHO: f64 = 1.225;
const MASS: f64 = 600.0;
const V: f64 = 50.0;

/// Symmetric aircraft with aileron, elevator and rudder.
///
/// The wing has a dihedral so that its side force acts ahead of the fin's.
/// With a fin as the only lateral surface, sideslip and rudder produce side
/// force and yawing moment in the same ratio, and the lateral trim equations
/// are singular. The dihedral breaks that degeneracy.
fn aircraft() -> VlmModel {
    let section = |y: f64, chord: f64| WingSection {
        le: Vec3::new(0.0, y, 0.0),
        chord,
        twist: 0.0,
    };
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
    let fin = WingPlanform::from_sections(vec![section(0.0, 0.9), section(1.5, 0.5)]).unwrap();
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

fn request(bank_deg: f64, turn_rate: f64) -> TurnRequest<'static> {
    TurnRequest {
        weight_n: MASS * G,
        v_inf: V,
        density: RHO,
        bank_rad: bank_deg.to_radians(),
        turn_rate,
        elevator: "elevator",
        aileron: "aileron",
        rudder: "rudder",
    }
}

#[test]
fn coordinated_turn_is_nearly_sideslip_free_with_the_turn_load_factor() {
    let m = aircraft();
    let bank = 20.0_f64;
    let req = request(bank, coordinated_turn_rate(V, bank.to_radians()));
    let t = trim_steady_turn(&m, &req).expect("converges");
    // The rudder's side force is not exactly cancelled at zero sideslip, so
    // a small sideslip (a fraction of a degree) remains.
    assert!(t.beta_rad.abs() < 0.005, "β = {}", t.beta_rad);
    let expected_n = 1.0 / bank.to_radians().cos();
    // Agrees with 1/cos φ to the O(α²) effect of the pitch angle.
    assert!(
        (t.load_factor - expected_n).abs() < 1e-4,
        "n = {}",
        t.load_factor
    );
    // Aileron and rudder hold the bank and the yaw rate.
    assert!(
        t.aileron_rad.abs() > 1e-4 && t.rudder_rad.abs() > 1e-4,
        "{t:?}"
    );
    assert!(t.iterations <= 40, "iterations {}", t.iterations);
}

#[test]
fn turn_trim_is_mirror_symmetric_in_bank() {
    // A symmetric aircraft: reversing bank and yaw rate keeps α and the
    // elevator, and reverses the aileron and rudder.
    let m = aircraft();
    let right = trim_steady_turn(
        &m,
        &request(20.0, coordinated_turn_rate(V, 20f64.to_radians())),
    )
    .unwrap();
    let left = trim_steady_turn(
        &m,
        &request(-20.0, coordinated_turn_rate(V, -20f64.to_radians())),
    )
    .unwrap();
    assert!((right.alpha_rad - left.alpha_rad).abs() < 1e-7);
    assert!((right.elevator_rad - left.elevator_rad).abs() < 1e-7);
    assert!((right.aileron_rad + left.aileron_rad).abs() < 1e-7);
    assert!((right.rudder_rad + left.rudder_rad).abs() < 1e-7);
}

#[test]
fn steady_sideslip_with_bank_is_sideslipped_into_the_bank() {
    // Right bank (gravity pulls to the right) needs a leftward aerodynamic
    // side force, i.e. wind from the right: β > 0, since C_Yβ < 0.
    let m = aircraft();
    let t = trim_steady_turn(&m, &request(1.0, 0.0)).expect("converges");
    assert!(t.beta_rad > 0.0 && t.beta_rad < 1.0, "β = {}", t.beta_rad);
    // Without a turn the vertical force is W cos φ.
    assert!((t.load_factor - 1.0_f64.to_radians().cos()).abs() < 1e-4);
}

#[test]
fn turn_trim_reports_errors() {
    let m = aircraft();
    let mut bad = request(20.0, 0.05);
    bad.weight_n = -1.0;
    assert_eq!(
        trim_steady_turn(&m, &bad).unwrap_err(),
        TrimError::InvalidRequest
    );

    let mut bank_too_steep = request(20.0, 0.05);
    bank_too_steep.bank_rad = 1.6;
    assert_eq!(
        trim_steady_turn(&m, &bank_too_steep).unwrap_err(),
        TrimError::InvalidRequest
    );

    let mut missing = request(20.0, 0.05);
    missing.rudder = "flaperon";
    assert_eq!(
        trim_steady_turn(&m, &missing).unwrap_err(),
        TrimError::NoSuchControl
    );
}
