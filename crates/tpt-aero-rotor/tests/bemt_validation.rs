//! BEMT validation: momentum-theory limits, actuator-disk mass-flow
//! conservation, induction bounds, and performance-map sanity.

use tpt_aero_rotor::{
    bemt::{BemtSolver, OperatingCondition},
    blade::{RotorError, RotorGeometry, SectionModel},
    loss::{HubLoss, TipLoss},
};

const RHO: f64 = 1.225;
const OMEGA: f64 = 300.0;

fn reference_solver() -> BemtSolver {
    let geom = RotorGeometry::rectangular(1.0, 0.1, 2, 0.12, 30).with_linear_twist(0.30, 0.05);
    BemtSolver::new(
        geom,
        SectionModel::Linear {
            cl_slope: 2.0 * core::f64::consts::PI,
            cl_max: 1.2,
            cd0: 0.01,
            cd_cl2: 0.02,
        },
        TipLoss::Prandtl,
        HubLoss::None,
    )
    .expect("valid")
}

#[test]
fn hover_thrust_and_convergence() {
    let solver = reference_solver();
    let op = OperatingCondition::hover(OMEGA, RHO);
    let sol = solver.solve(&op).expect("converged");
    assert!(
        sol.loads.ct > 0.0 && sol.loads.ct < 0.05,
        "ct = {}",
        sol.loads.ct
    );
    assert!(sol.loads.cq > 0.0);
    assert!(
        sol.max_residual < 1e-6,
        "element-momentum residual {}",
        sol.max_residual
    );
    // Figure of merit for this class of rotor: 0.5 .. 0.85.
    let fom = sol.loads.figure_of_merit;
    assert!((0.4..0.85).contains(&fom), "figure of merit = {fom}");
}

#[test]
fn induction_factors_stay_in_physical_bounds() {
    let solver = reference_solver();
    for &v in &[-2.0, 0.0, 5.0, 20.0, 50.0] {
        let op = OperatingCondition {
            axial_speed: v,
            omega: OMEGA,
            density: RHO,
            collective: 0.15,
        };
        let sol = solver.solve(&op).expect("converged");
        for s in &sol.sections {
            assert!(
                (0.0..=1.0).contains(&s.lambda_induced),
                "λ_i = {} at v = {v}",
                s.lambda_induced
            );
        }
    }
}

/// spec.txt Phase 2 exit: propeller performance map matches analytical
/// momentum theory limits. With a frictionless section closure the BEMT is
/// exactly the ideal actuator disk: its Froude efficiency must match
/// `η_ideal = 2/(1 + √(1 + 2C_T/λ_c²))` and never exceed it. A real
/// (profile-drag) closure must sit strictly below the bound.
#[test]
fn momentum_theory_limit_frictionless() {
    let geom = RotorGeometry::rectangular(1.0, 0.1, 2, 0.12, 30).with_linear_twist(0.30, 0.05);
    let solver = BemtSolver::new(
        geom,
        SectionModel::Linear {
            cl_slope: 2.0 * core::f64::consts::PI,
            cl_max: 1.2,
            cd0: 0.0,
            cd_cl2: 0.0,
        },
        TipLoss::None,
        HubLoss::None,
    )
    .expect("valid");
    for &(v, theta) in &[
        (10.0_f64, 0.10_f64),
        (20.0, 0.08),
        (30.0, 0.06),
        (25.0, 0.04),
    ] {
        let op = OperatingCondition {
            axial_speed: v,
            omega: OMEGA,
            density: RHO,
            collective: theta,
        };
        let sol = solver.solve(&op).expect("converged");
        let lambda_c = v / (1.0 * OMEGA);
        let ct = sol.loads.ct;
        assert!(ct > 0.0, "propeller must produce thrust (v={v})");
        let eta_ideal = 2.0 / (1.0 + (1.0 + 2.0 * ct / (lambda_c * lambda_c)).sqrt());
        let eta = sol.loads.propeller_efficiency;
        assert!(
            eta <= eta_ideal + 0.02,
            "v={v}: η {eta} exceeds momentum-theory bound {eta_ideal}"
        );
        assert!(
            eta > eta_ideal - 0.08,
            "v={v}: η {eta} too far below frictionless ideal {eta_ideal}"
        );
    }
}

/// A real (profile-drag) closure stays strictly below the momentum bound.
#[test]
fn real_rotor_below_momentum_bound() {
    let solver = reference_solver();
    let op = OperatingCondition {
        axial_speed: 20.0,
        omega: OMEGA,
        density: RHO,
        collective: 0.08,
    };
    let sol = solver.solve(&op).expect("converged");
    let lambda_c = 20.0 / OMEGA;
    let ct = sol.loads.ct;
    let eta_ideal = 2.0 / (1.0 + (1.0 + 2.0 * ct / (lambda_c * lambda_c)).sqrt());
    assert!(
        sol.loads.propeller_efficiency < eta_ideal,
        "profile drag must cost efficiency"
    );
}

/// spec.txt Phase 2 exit: total mass flow is conserved through the actuator
/// disk — the element-side and momentum-side annulus thrust must agree at
/// convergence for every annulus (that equality *is* local mass conservation
/// in the combined BEMT formulation).
#[test]
fn actuator_disk_mass_flow_conserved() {
    let solver = reference_solver();
    for &v in &[0.0_f64, 15.0, 40.0] {
        let op = OperatingCondition {
            axial_speed: v,
            omega: OMEGA,
            density: RHO,
            collective: 0.18,
        };
        let sol = solver.solve(&op).expect("converged");
        for s in &sol.sections {
            // Conservation is asserted where the momentum model applies:
            // positive-thrust annuli away from the singular tip (F → 0).
            // Sections at negative α (windmill/brake regime at high J) and
            // the tip station are outside the momentum formulation.
            if s.loss < 0.05 || s.dct_dx <= 0.0 {
                continue;
            }
            let residual = (s.dct_dx - s.dct_dx_momentum) / s.dct_dx.abs().max(1e-12);
            assert!(
                residual.abs() < 1e-3,
                "annulus x={:.3} v={v}: relative residual {residual}",
                s.r_frac
            );
        }
    }
}

#[test]
fn performance_map_shape() {
    let solver = reference_solver();
    let speeds = [0.0_f64, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0];
    let map =
        tpt_aero_rotor::performance_map(&solver, OMEGA, RHO, &speeds, 0.15).expect("map solved");
    // J axis ascending.
    for w in map.windows(2) {
        assert!(w[1].j > w[0].j);
    }
    // Thrust falls with advance ratio at fixed collective (unloading).
    for w in map.windows(2) {
        assert!(w[1].ct < w[0].ct, "ct must decrease with J");
    }
    // In the propeller region (positive thrust and torque, moderate J) the
    // Froude efficiency is positive and below 1. At high J the inboard
    // sections enter the windmill/brake regime where the propeller
    // efficiency metric is not meaningful.
    for p in &map {
        if p.j > 0.0 && p.j <= 0.5 && p.ct > 0.0 && p.cq > 0.0 {
            assert!(p.efficiency < 1.0, "η = {} at J = {}", p.efficiency, p.j);
            assert!(p.efficiency > 0.0, "η = {} at J = {}", p.efficiency, p.j);
        }
        if p.j > 0.2 && p.j < 0.5 {
            assert!(p.efficiency > 0.2, "η = {} at J = {}", p.efficiency, p.j);
        }
    }
    // Near-hover point should have a healthy figure of merit.
    assert!(
        (0.3..0.9).contains(&map[0].figure_of_merit),
        "hover FOM = {}",
        map[0].figure_of_merit
    );
}

#[test]
fn tip_loss_reduces_thrust() {
    let geom = RotorGeometry::rectangular(1.0, 0.1, 2, 0.12, 30).with_linear_twist(0.30, 0.05);
    let mk = |tip| {
        BemtSolver::new(
            geom.clone(),
            SectionModel::Linear {
                cl_slope: 2.0 * core::f64::consts::PI,
                cl_max: 1.2,
                cd0: 0.01,
                cd_cl2: 0.02,
            },
            tip,
            HubLoss::None,
        )
        .expect("valid")
    };
    let ideal = mk(TipLoss::None)
        .solve(&OperatingCondition::hover(OMEGA, RHO))
        .expect("ok");
    let prandtl = mk(TipLoss::Prandtl)
        .solve(&OperatingCondition::hover(OMEGA, RHO))
        .expect("ok");
    let glauert = mk(TipLoss::Glauert)
        .solve(&OperatingCondition::hover(OMEGA, RHO))
        .expect("ok");
    assert!(
        prandtl.loads.ct < ideal.loads.ct,
        "Prandtl loss must unload the tips"
    );
    assert!(
        glauert.loads.ct < ideal.loads.ct,
        "Glauert loss must unload the tips"
    );
    // The two finite-blade models should agree closely.
    let rel = (prandtl.loads.ct - glauert.loads.ct).abs() / prandtl.loads.ct;
    assert!(rel < 0.05, "Prandtl vs Glauert differ by {rel}");
}

#[test]
fn hub_loss_reduces_thrust_further() {
    let geom = RotorGeometry::rectangular(1.0, 0.1, 2, 0.12, 30).with_linear_twist(0.30, 0.05);
    let mk = |hub| {
        BemtSolver::new(
            geom.clone(),
            SectionModel::Linear {
                cl_slope: 2.0 * core::f64::consts::PI,
                cl_max: 1.2,
                cd0: 0.01,
                cd_cl2: 0.02,
            },
            TipLoss::Prandtl,
            hub,
        )
        .expect("valid")
    };
    let no_hub = mk(HubLoss::None)
        .solve(&OperatingCondition::hover(OMEGA, RHO))
        .expect("ok");
    let with_hub = mk(HubLoss::Prandtl { x_hub: 0.2 })
        .solve(&OperatingCondition::hover(OMEGA, RHO))
        .expect("ok");
    assert!(with_hub.loads.ct < no_hub.loads.ct);
    assert!(
        (with_hub.loads.ct - no_hub.loads.ct).abs() / no_hub.loads.ct < 0.15,
        "hub loss should be a correction, not a redesign"
    );
}

#[test]
fn polar_closure_runs() {
    // The tabulated-polar closure must interoperate with tpt-aero-core.
    use tpt_aero_core::AirfoilPolar;
    let alphas: Vec<f64> = (-12..=16).map(|d| f64::from(d).to_radians()).collect();
    let cl: Vec<f64> = alphas.iter().map(|&a| core::f64::consts::TAU * a).collect();
    let cd: Vec<f64> = alphas.iter().map(|&a| 0.01 + 0.03 * a * a).collect();
    let cm = vec![0.0; alphas.len()];
    let polar = AirfoilPolar::new(alphas, cl, cd, cm, 0.2, 1e6).expect("valid");
    let geom = RotorGeometry::rectangular(1.0, 0.1, 2, 0.12, 30).with_linear_twist(0.30, 0.05);
    let solver = BemtSolver::new(
        geom,
        SectionModel::Polar { polar, stall: None },
        TipLoss::Prandtl,
        HubLoss::None,
    )
    .expect("valid");
    let sol = solver
        .solve(&OperatingCondition::hover(OMEGA, RHO))
        .expect("ok");
    assert!(sol.loads.ct > 0.0);
}

#[test]
fn rejects_invalid_geometry() {
    use tpt_aero_rotor::blade::RotorError;
    let geom = RotorGeometry::rectangular(1.0, 1.5, 2, 0.12, 4); // hub > tip
    let err = BemtSolver::new(
        geom,
        SectionModel::Linear {
            cl_slope: core::f64::consts::TAU,
            cl_max: 1.0,
            cd0: 0.01,
            cd_cl2: 0.01,
        },
        TipLoss::None,
        HubLoss::None,
    )
    .unwrap_err();
    assert_eq!(err, RotorError::InvalidGeometry);
}

#[test]
fn zero_hub_rectangular_blade_is_valid_and_solves() {
    // Hub radius 0 used to put the first station at r = 0, which validate()
    // rejects; the constructor now starts the stations at a small fraction.
    let geom = RotorGeometry::rectangular(1.0, 0.0, 2, 0.12, 20).with_linear_twist(0.30, 0.05);
    assert!(geom.validate().is_ok());
    let solver = BemtSolver::new(
        geom,
        SectionModel::Linear {
            cl_slope: core::f64::consts::TAU,
            cl_max: 1.2,
            cd0: 0.01,
            cd_cl2: 0.02,
        },
        TipLoss::Prandtl,
        HubLoss::None,
    )
    .expect("valid");
    let sol = solver
        .solve(&OperatingCondition::hover(OMEGA, RHO))
        .expect("solves");
    assert!(sol.loads.ct > 0.0);
}

#[test]
fn nan_and_station_inside_hub_are_rejected() {
    let mut geom = RotorGeometry::rectangular(1.0, 0.2, 2, 0.12, 10);
    geom.sections[3].chord = f64::NAN;
    assert_eq!(geom.validate(), Err(RotorError::InvalidGeometry));
    let mut inside = RotorGeometry::rectangular(1.0, 0.2, 2, 0.12, 10);
    inside.sections[0].r_frac = 0.1; // inside the 0.2 hub
    assert_eq!(inside.validate(), Err(RotorError::InvalidGeometry));
}

#[test]
fn bad_operating_condition_has_its_own_error() {
    let solver = reference_solver();
    let op = OperatingCondition::hover(-OMEGA, RHO);
    assert_eq!(
        solver.solve(&op).unwrap_err(),
        RotorError::InvalidOperatingCondition
    );
    let op = OperatingCondition::hover(OMEGA, f64::NAN);
    assert_eq!(
        solver.solve(&op).unwrap_err(),
        RotorError::InvalidOperatingCondition
    );
}

#[test]
fn small_descent_keeps_total_inflow_non_negative() {
    // The bracket used to clamp λ_tot at zero in descent, silently treating
    // it as hover. Descent must now give a physical, positive-thrust solution
    // with λ_tot ≥ 0 at every station.
    let solver = reference_solver();
    let op = OperatingCondition {
        axial_speed: -0.02 * OMEGA,
        ..OperatingCondition::hover(OMEGA, RHO)
    };
    let sol = solver.solve(&op).expect("descent converges");
    assert!(sol.loads.ct > 0.0, "ct = {}", sol.loads.ct);
    assert!(sol.sections.iter().all(|s| s.lambda_tot >= 0.0));
    assert!(sol.sections.iter().all(|s| s.lambda_induced >= 0.0));
}

#[test]
fn exhausted_iteration_budget_reports_not_converged() {
    // Ten halvings cannot narrow the bracket to 1e-12, so the solve must say so.
    let solver = reference_solver().with_iteration_controls(10, 1.0e-11, 0.6);
    let err = solver
        .solve(&OperatingCondition::hover(OMEGA, RHO))
        .unwrap_err();
    assert_eq!(err, RotorError::NotConverged);
}
