//! BEMT against the Caradonna–Tung model rotor in hover (transonic tip).
//!
//! Experimental thrust coefficient `C_T = 0.00473` at collective 8°, tip Mach
//! 0.877 (2500 rpm), from Caradonna & Tung, NASA TM-81232 (1981). The value is
//! as quoted in the CFD validation paper arXiv:2107.13895; the NASA table itself
//! was not checked. Rotor: R = 1.143 m, chord 0.1905 m, two untwisted blades.
//!
//! The BEMT here uses a linear lift-curve closure with no Mach correction and
//! no stall, so it is expected to over-predict at this tip Mach number. The
//! band below records the current discrepancy (about +23 %) and fails if it
//! grows, rather than asserting agreement that the model does not have.

use tpt_aero_rotor::{
    BemtSolver, HubLoss, OperatingCondition, RotorGeometry, SectionModel, TipLoss,
};

const CT_EXPERIMENT: f64 = 0.00473;

#[test]
fn hover_thrust_is_within_the_documented_band_of_the_measurement() {
    let geometry = RotorGeometry::rectangular(1.143, 0.01, 2, 0.1905, 40);
    let section = SectionModel::Linear {
        cl_slope: 2.0 * core::f64::consts::PI,
        cl_max: 1.4,
        cd0: 0.008,
        cd_cl2: 0.0,
    };
    let solver = BemtSolver::new(geometry, section, TipLoss::Prandtl, HubLoss::None).unwrap();
    let omega = 2500.0 * core::f64::consts::TAU / 60.0;
    let op = OperatingCondition {
        axial_speed: 0.0,
        omega,
        density: 1.2,
        collective: 8.0_f64.to_radians(),
    };
    let ct = solver.solve(&op).unwrap().loads.ct;
    let error = ct / CT_EXPERIMENT - 1.0;
    // Documented: about +23 % at the time of writing. Bound it at +30 %.
    assert!(
        error > 0.0,
        "BEMT is expected to over-predict here; C_T = {ct}"
    );
    assert!(
        error < 0.30,
        "C_T = {ct} is {:.0} % above the measurement",
        error * 100.0
    );
}
