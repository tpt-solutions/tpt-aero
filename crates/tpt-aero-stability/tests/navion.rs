//! Navion-like wing and tail against flight-extracted longitudinal derivatives.
//!
//! Geometry: Table I and Figure 2 of Seckel & Morris / NASA TN D-6643 ("Aerodynamic
//! parameters of the Navion airplane extracted from flight data"): wing area
//! 17.112 m², AR 6.04, taper 0.54, 3.0° leading-edge sweep, 7.5° dihedral, MAC
//! 1.74 m, +2° root and −1° tip incidence; horizontal tail 4.0 m², AR 4.0, taper
//! 0.67, 6° sweep, set at −5° for the flight tests; moment centre at 0.25 MAC.
//! The tail arm is read from the Figure 2 three-view (about 4.7 m to the tail
//! quarter chord). It is an estimate from a drawing, not a tabulated value.
//!
//! Flight data (Table VII, clean configuration, 73.2 m/s): C_Zα = −4.33 (so
//! C_Lα ≈ 4.33) and C_mα = −0.77 about 25 % c (−0.63 at the test CG).
//!
//! This model has no fuselage, nacelle or propeller. The bands below are what
//! the model gives today, from the sources above, and the test fails if they
//! move. They are not agreement with the flight data.

use tpt_aero_core::{Vec3, WingPlanform, WingSection};
use tpt_aero_stability::longitudinal::longitudinal_fd;
use tpt_aero_vlm::{Surface, VlmModel};

const FLIGHT_CL_ALPHA: f64 = 4.33;
const FLIGHT_CM_ALPHA_25C: f64 = -0.77;

fn navion_like() -> VlmModel {
    let area = 17.112_f64;
    let aspect = 6.04_f64;
    let taper = 0.54_f64;
    let span = (aspect * area).sqrt();
    let root = 2.0 * area / (span * (1.0 + taper));
    let tip = taper * root;
    let half = span / 2.0;
    let dihedral = 7.5_f64.to_radians();
    let sweep = half * 2.996_f64.to_radians().tan();
    let station = |y: f64, twist_deg: f64, chord: f64, le_x: f64| WingSection {
        le: Vec3::new(le_x, y, y.abs() * dihedral.tan()),
        chord,
        twist: twist_deg.to_radians(),
    };
    let wing = WingPlanform::from_sections(vec![
        station(-half, -1.0, tip, sweep),
        station(0.0, 2.0, root, 0.0),
        station(half, -1.0, tip, sweep),
    ])
    .expect("valid wing");

    let tail_span = 4.0_f64;
    let tail_root = 2.0 * 4.0 / (tail_span * 1.67);
    let tail_sweep = (tail_span / 2.0) * 6.0_f64.to_radians().tan();
    let tail = WingPlanform::tapered(tail_span, tail_root, 0.67 * tail_root, tail_sweep, 0.0);
    VlmModel::from_surfaces(&[
        Surface::horizontal(wing, 20, 4),
        Surface::horizontal(tail, 8, 3)
            .at(5.0, 0.0, 0.0)
            .with_incidence((-5.0_f64).to_radians())
            .expect("valid tail incidence"),
    ])
    .expect("valid aircraft")
}

#[test]
fn wing_geometry_reproduces_table_i() {
    let m = navion_like();
    let r = m.reference();
    // Table I: S = 17.112 m², MAC = 1.74 m (the planform gives 1.733 m).
    assert!((r.area - 17.112).abs() < 1e-3, "S = {}", r.area);
    assert!((r.chord - 1.74).abs() < 0.02, "MAC = {}", r.chord);
}

#[test]
fn lift_slope_and_pitch_stiffness_against_flight_data() {
    let m = navion_like();
    let d = longitudinal_fd(&m, 2.0_f64.to_radians(), 73.2, 1.225, 1e-4).expect("solve");
    let cl_ratio = d.cl_alpha / FLIGHT_CL_ALPHA;
    let cm_ratio = d.cm_alpha / FLIGHT_CM_ALPHA_25C;
    // Observed: C_Lα ≈ 4.98 (+15 %), C_mα ≈ −1.37 (1.8× the flight value).
    assert!(
        (1.05..1.25).contains(&cl_ratio),
        "C_Lα = {} ({cl_ratio:.2}× flight)",
        d.cl_alpha
    );
    assert!(
        (1.4..2.2).contains(&cm_ratio),
        "C_mα = {} ({cm_ratio:.2}× flight)",
        d.cm_alpha
    );
    // Both have the right sign: statically stable and positive lift slope.
    assert!(d.cl_alpha > 0.0 && d.cm_alpha < 0.0);
}
