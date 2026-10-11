//! Propeller-driven flight: a BEMT propeller supplies the thrust in a 6-DOF
//! simulation, and a fuselage and wing drag build-up sets the parasite drag.
//!
//! The propeller's shaft speed is found so that its thrust balances the level-
//! flight drag, the aircraft is trimmed there, and then the throttle is opened
//! by 1 % at t = 1 s. The response shows the speed build-up that the thrust
//! step produces.
//!
//! Further loads act on the airframe. The propeller's torque reaction rolls
//! the aircraft the other way, and the fuselage's first-order (Munk) lift and
//! pitching moment act about the centre of gravity. The propeller slipstream
//! raises the dynamic pressure at the elevator, which is behind the disk, and
//! that adds lift on the tail; see `tail_slipstream_loads`. The slipstream is
//! either the contracting momentum tube or the linearized vortex-cylinder field,
//! chosen by `SLIPSTREAM`; both are printed at the trim state.
//!
//! Near the trim point this propeller's thrust is steep in shaft speed (a 1 %
//! change in rpm gives a far larger change in thrust), so the step is small.
//!
//! Run with `cargo run --release -p tpt-aero --example propeller_flight`.

use tpt_aero::aero_core::viscous::{cd0_buildup, DragComponent};
use tpt_aero::aero_core::{Vec3, WingPlanform, WingSection};
use tpt_aero::panel::SlenderBody;
use tpt_aero::rotor::{
    BemtSolver, CylinderWake, HubLoss, OperatingCondition, RotorGeometry, SectionModel, Slipstream,
    SlipstreamModel, TipLoss,
};
use tpt_aero::stability::controls::control_derivatives;
use tpt_aero::stability::sixdof::{simulate_with_loads, Controls, RigidBody, SixDofConfig, State};
use tpt_aero::stability::trim::{trim_level_flight, TrimRequest};
use tpt_aero::vlm::{ControlSymmetry, Surface, VlmModel};

const MASS: f64 = 600.0;
const G: f64 = 9.806_65;
const V0: f64 = 50.0;
const RHO: f64 = 1.225;
const PROP_RADIUS: f64 = 1.0;
/// Distance from the fuselage nose to the centre of gravity (m).
const FUSELAGE_CG: f64 = 2.5;
/// Elevator: span (m), chord (m), and quarter-chord station from the nose (m).
const TAIL_SPAN: f64 = 4.0;
const TAIL_CHORD: f64 = 0.6;
const TAIL_X_QC: f64 = 4.15;

fn aircraft() -> VlmModel {
    let line = |x: f64, y: f64, chord: f64| WingSection {
        le: Vec3::new(x, y, 0.0),
        chord,
        twist: 0.0,
    };
    let wing = WingPlanform::from_sections(vec![
        line(0.0, -5.0, 1.2),
        line(0.0, 0.0, 1.2),
        line(0.0, 5.0, 1.2),
    ])
    .expect("valid wing");
    let fin = WingPlanform::from_sections(vec![line(0.0, 0.0, 0.9), line(0.0, 1.5, 0.5)])
        .expect("valid fin");
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
    .expect("valid aircraft")
}

/// Two-blade propeller with linear twist and the linear section closure.
fn propeller() -> BemtSolver {
    let geometry =
        RotorGeometry::rectangular(PROP_RADIUS, 0.1, 2, 0.12, 20).with_linear_twist(0.5, 0.15);
    let section = SectionModel::Linear {
        cl_slope: 2.0 * core::f64::consts::PI,
        cl_max: 1.4,
        cd0: 0.01,
        cd_cl2: 0.01,
    };
    BemtSolver::new(geometry, section, TipLoss::Prandtl, HubLoss::None).expect("valid propeller")
}

/// Fuselage: nose at x = 0, 4 m long, bluff base of radius 0.15 m.
fn fuselage_body() -> SlenderBody {
    let n = 41;
    let x: Vec<f64> = (0..n).map(|i| 4.0 * i as f64 / (n - 1) as f64).collect();
    let radius: Vec<f64> = x
        .iter()
        .map(|&xi| 0.15 * (xi / 4.0).sqrt().min(1.0))
        .collect();
    SlenderBody::new(x, radius).expect("valid fuselage profile")
}

/// Propeller thrust (N) and shaft torque (N·m) at shaft speed `rpm` and axial
/// inflow `speed` (m/s).
fn propeller_loads(solver: &BemtSolver, rpm: f64, speed: f64) -> (f64, f64) {
    let omega = rpm * core::f64::consts::TAU / 60.0;
    let sol = solver
        .solve(&OperatingCondition {
            axial_speed: speed,
            omega,
            density: RHO,
            collective: 0.0,
        })
        .expect("BEMT converges");
    let area = core::f64::consts::PI * PROP_RADIUS * PROP_RADIUS;
    let q_ref = RHO * area * (omega * PROP_RADIUS).powi(2);
    (sol.loads.ct * q_ref, sol.loads.cq * q_ref * PROP_RADIUS)
}

/// Which slipstream model the tail sees: the momentum-theory contracting tube,
/// or the linearized vortex-cylinder field (the full induced velocity).
#[derive(Debug, Clone, Copy, PartialEq)]
enum SlipstreamChoice {
    Tube,
    Cylinder,
}

/// The slipstream model used in the simulation.
const SLIPSTREAM: SlipstreamChoice = SlipstreamChoice::Cylinder;

/// Build the chosen slipstream for a disk producing `thrust` (N) at `speed`
/// (m/s). `None` if the inputs are invalid (for example, no thrust).
fn slipstream_model(choice: SlipstreamChoice, thrust: f64, speed: f64) -> Option<SlipstreamModel> {
    let tube = Slipstream::new(PROP_RADIUS, thrust.max(0.0), speed, RHO).ok()?;
    match choice {
        SlipstreamChoice::Tube => Some(SlipstreamModel::Tube(tube)),
        SlipstreamChoice::Cylinder => {
            CylinderWake::new(PROP_RADIUS, tube.disk_induced_velocity(), speed, RHO)
                .ok()
                .map(SlipstreamModel::Cylinder)
        }
    }
}

/// Extra tail lift (N, up positive) and nose-up pitching moment (N·m) from the
/// propeller slipstream at airspeed `speed` (m/s) and angle of attack `alpha`
/// (rad). The slipstream raises the elevator's dynamic pressure by the mean
/// `Δq` over its span, and the tail lift responds to that with the same lift
/// slope the elevator has at `q∞`. Without thrust there is no slipstream.
fn tail_slipstream_loads(
    choice: SlipstreamChoice,
    thrust: f64,
    speed: f64,
    alpha: f64,
) -> (f64, f64) {
    let Some(slip) = slipstream_model(choice, thrust, speed) else {
        return (0.0, 0.0);
    };
    let n = 9;
    let mean_dq = (0..n)
        .map(|i| {
            let y = -0.5 * TAIL_SPAN + TAIL_SPAN * i as f64 / (n - 1) as f64;
            slip.dynamic_pressure_increment(TAIL_X_QC, y)
        })
        .sum::<f64>()
        / n as f64;
    let area = TAIL_SPAN * TAIL_CHORD;
    let aspect = TAIL_SPAN * TAIL_SPAN / area;
    let cl_alpha = 2.0 * core::f64::consts::PI * aspect / (aspect + 2.0);
    let lift = mean_dq * area * cl_alpha * alpha;
    let arm = TAIL_X_QC - FUSELAGE_CG;
    (lift, -lift * arm)
}

/// Shaft speed (rpm) at which the propeller produces `target` thrust at `speed`,
/// by bisection (thrust is increasing in rpm over this range).
fn rpm_for_thrust(solver: &BemtSolver, target: f64, speed: f64) -> f64 {
    let (mut lo, mut hi) = (500.0_f64, 12_000.0_f64);
    for _ in 0..40 {
        let mid = 0.5 * (lo + hi);
        if propeller_loads(solver, mid, speed).0 < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

fn main() {
    let model = aircraft();
    // Parasite drag: wing skin friction plus a fuselage with its own form factor.
    // Fuselage wetted area ≈ 20 m² against the 12 m² reference area.
    let fuselage = DragComponent {
        reynolds: 1.0e7,
        wetted_ratio: 1.7,
        thickness_ratio: 0.0,
        form_override: Some(1.2),
        interference: 1.1,
    };
    let wing = DragComponent::lifting_surface(3.4e6, 2.0, 0.12);
    let cd0 = cd0_buildup(&[wing, fuselage]);
    println!("C_D0 (wing + fuselage) = {cd0:.5}");

    let trim = trim_level_flight(
        &model,
        &TrimRequest {
            weight_n: MASS * G,
            v_inf: V0,
            density: RHO,
            elevator: "elevator",
            cd0,
        },
    )
    .expect("level trim");
    println!(
        "level trim at {V0} m/s: α = {:.2}°, elevator = {:.2}°, drag = {:.1} N",
        trim.alpha_rad.to_degrees(),
        trim.elevator_rad.to_degrees(),
        trim.thrust_n
    );

    let prop = propeller();
    let rpm0 = rpm_for_thrust(&prop, trim.thrust_n, V0);
    let rpm1 = 1.01 * rpm0;
    // The level trim above knows nothing of the slipstream. Re-trim the elevator
    // so the slipstream's pitching moment at the trim state is balanced, using
    // the elevator's pitching-moment derivative from the VLM and the propeller's
    // thrust at the trim shaft speed.
    let thrust0 = propeller_loads(&prop, rpm0, V0 * trim.alpha_rad.cos()).0;
    let (_, slip_pitch) = tail_slipstream_loads(SLIPSTREAM, thrust0, V0, trim.alpha_rad);
    let cm_elevator = control_derivatives(&model, "elevator", trim.alpha_rad, V0, RHO)
        .expect("elevator derivative")
        .cm;
    let reference = model.reference();
    let q_s_c = 0.5 * RHO * V0 * V0 * reference.area * reference.chord;
    let elevator_trim = trim.elevator_rad - slip_pitch / (q_s_c * cm_elevator);
    println!(
        "slipstream re-trim: elevator {:.2}° -> {:.2}° (balances {slip_pitch:.1} N·m)",
        trim.elevator_rad.to_degrees(),
        elevator_trim.to_degrees()
    );
    println!("trim shaft speed {rpm0:.0} rpm; throttle step to {rpm1:.0} rpm (+1 %) at t = 1 s");

    let cfg = SixDofConfig {
        body: RigidBody {
            mass: MASS,
            ixx: 800.0,
            iyy: 900.0,
            izz: 1200.0,
            ixz: 0.0,
        },
        density: RHO,
        parasite_cd: cd0,
        elevator: "elevator",
        aileron: "aileron",
        rudder: "rudder",
        wind_ned: [0.0, 0.0, 0.0],
        thrust_body_fixed: false,
    };
    let (sa, ca) = trim.alpha_rad.sin_cos();
    let x0 = State {
        u: V0 * ca,
        v: 0.0,
        w: V0 * sa,
        p: 0.0,
        q: 0.0,
        r: 0.0,
        phi: 0.0,
        theta: trim.alpha_rad,
        psi: 0.0,
        north: 0.0,
        east: 0.0,
        altitude: 1000.0,
    };
    let controls = Controls {
        elevator: elevator_trim,
        ..Controls::default()
    };
    let throttle = |t: f64| if t < 1.0 { rpm0 } else { rpm1 };
    let thrust_of = |t: f64, s: &State| propeller_loads(&prop, throttle(t), s.u.max(1.0)).0;
    let torque_of = |t: f64, s: &State| propeller_loads(&prop, throttle(t), s.u.max(1.0)).1;
    let body = fuselage_body();
    // Additive loads in aircraft axes: the fuselage's Munk lift, side force and
    // moments about the CG, and the propeller's torque reaction on the airframe.
    // The trim is symmetric, so the sideslip terms are zero in this run.
    let extra = |t: f64, s: &State| {
        let alpha = s.w.atan2(s.u);
        let beta = (s.v / s.airspeed()).clamp(-1.0, 1.0).asin();
        let (lift, pitch) = body.munk_lift_and_moment(RHO, s.airspeed(), alpha, FUSELAGE_CG);
        let (side, yaw) = body.munk_side_force_and_yaw(RHO, s.airspeed(), beta, FUSELAGE_CG);
        let (tail_lift, tail_pitch) =
            tail_slipstream_loads(SLIPSTREAM, thrust_of(t, s), s.airspeed(), alpha);
        (
            [0.0, side, -(lift + tail_lift)],
            [-torque_of(t, s), pitch + tail_pitch, yaw],
        )
    };

    let dt = 0.02;
    let steps = 500;
    let history = simulate_with_loads(
        &model,
        &cfg,
        x0,
        dt,
        steps,
        &|_| controls,
        &thrust_of,
        &extra,
    )
    .expect("6-DOF run");

    let (_, s0) = history[0];
    let (fuselage_lift, fuselage_pitch) =
        body.munk_lift_and_moment(RHO, s0.airspeed(), s0.w.atan2(s0.u), FUSELAGE_CG);
    println!(
        "at trim: fuselage lift {fuselage_lift:.2} N, fuselage pitching moment {fuselage_pitch:.2} N·m, propeller torque reaction {:.2} N·m",
        -torque_of(0.0, &s0)
    );
    // The same trim state under both slipstream models, for comparison.
    for choice in [SlipstreamChoice::Tube, SlipstreamChoice::Cylinder] {
        let (tail_lift0, tail_pitch0) =
            tail_slipstream_loads(choice, thrust_of(0.0, &s0), s0.airspeed(), s0.w.atan2(s0.u));
        println!(
            "at trim, {choice:?} slipstream: tail lift {tail_lift0:.3} N, pitching moment {tail_pitch0:.3} N·m"
        );
    }

    println!("\n    t (s)   airspeed (m/s)   altitude (m)   roll rate (deg/s)   thrust (N)");
    for k in (0..=steps).step_by(50) {
        let (t, s) = history[k];
        println!(
            "{t:8.2}   {:14.3}   {:12.3}   {:17.3}   {:10.1}",
            s.airspeed(),
            s.altitude,
            s.p.to_degrees(),
            thrust_of(t, &s)
        );
    }
}
