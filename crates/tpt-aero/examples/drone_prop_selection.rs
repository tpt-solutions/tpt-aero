//! Drone propeller selection against a motor curve.
//!
//! A brushless motor is modelled by a linear torque-speed curve,
//! `Q_motor(n) = Q_stall (1 − n / n_0)`, where `n_0` is the no-load speed. For
//! each candidate propeller the operating speed is where the propeller's BEMT
//! torque equals the motor torque. Its thrust there is compared with the thrust
//! each rotor must produce, and the candidates are ranked by figure of merit.
//!
//! Scenario: a 1.5 kg quadrotor needs about 3.7 N of thrust per rotor in hover.
//!
//! Run with `cargo run --release -p tpt-aero --example drone_prop_selection`.

use tpt_aero::rotor::{
    BemtSolver, HubLoss, OperatingCondition, RotorGeometry, SectionModel, TipLoss,
};

const RHO: f64 = 1.225;
const THRUST_PER_ROTOR_N: f64 = 3.7;
/// Motor stall torque (N·m) and no-load speed (rpm).
const Q_STALL: f64 = 0.15;
const N_NO_LOAD_RPM: f64 = 11_000.0;

struct Candidate {
    name: &'static str,
    radius: f64,
    chord_root: f64,
    tip_ratio: f64,
    root_twist_deg: f64,
    tip_twist_deg: f64,
}

/// Thrust (N), power (W), torque (N·m) and figure of merit at shaft speed `rpm`.
fn perf(c: &Candidate, rpm: f64) -> (f64, f64, f64, f64) {
    let mut geometry = RotorGeometry::rectangular(c.radius, 0.08 * c.radius, 2, c.chord_root, 16)
        .with_linear_twist(c.root_twist_deg.to_radians(), c.tip_twist_deg.to_radians());
    let n = geometry.sections.len();
    for (i, sec) in geometry.sections.iter_mut().enumerate() {
        let t = i as f64 / (n - 1) as f64;
        sec.chord = c.chord_root * (1.0 + (c.tip_ratio - 1.0) * t);
    }
    let section = SectionModel::Linear {
        cl_slope: core::f64::consts::TAU,
        cl_max: 1.2,
        cd0: 0.012,
        cd_cl2: 0.012,
    };
    let solver =
        BemtSolver::new(geometry, section, TipLoss::Prandtl, HubLoss::None).expect("valid");
    let omega = rpm * core::f64::consts::TAU / 60.0;
    let sol = solver
        .solve(&OperatingCondition::hover(omega, RHO))
        .expect("converges");
    let area = core::f64::consts::PI * c.radius * c.radius;
    let q_ref = RHO * area * (omega * c.radius).powi(2);
    (
        sol.loads.ct * q_ref,
        sol.loads.cp * q_ref * omega * c.radius,
        sol.loads.cq * q_ref * c.radius,
        sol.loads.figure_of_merit,
    )
}

/// Operating speed (rpm) where the propeller torque equals the motor torque,
/// found by bisection on the difference.
fn operating_rpm(c: &Candidate) -> f64 {
    let motor_torque = |rpm: f64| Q_STALL * (1.0 - rpm / N_NO_LOAD_RPM).max(0.0);
    let excess = |rpm: f64| perf(c, rpm).2 - motor_torque(rpm);
    let (mut lo, mut hi) = (1000.0_f64, N_NO_LOAD_RPM);
    // Propeller torque grows with speed, the motor's falls; the root is bracketed.
    for _ in 0..50 {
        let mid = 0.5 * (lo + hi);
        if excess(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

fn main() {
    let candidates = [
        Candidate {
            name: "9x4.5 style",
            radius: 0.115,
            chord_root: 0.022,
            tip_ratio: 0.6,
            root_twist_deg: 30.0,
            tip_twist_deg: 12.0,
        },
        Candidate {
            name: "10x4.5 style",
            radius: 0.127,
            chord_root: 0.024,
            tip_ratio: 0.55,
            root_twist_deg: 30.0,
            tip_twist_deg: 11.0,
        },
        Candidate {
            name: "11x5 style",
            radius: 0.140,
            chord_root: 0.026,
            tip_ratio: 0.5,
            root_twist_deg: 32.0,
            tip_twist_deg: 10.0,
        },
    ];
    println!("=== Drone propeller selection ===");
    println!("motor: Q_stall {Q_STALL} N·m, no-load {N_NO_LOAD_RPM} rpm; need {THRUST_PER_ROTOR_N} N per rotor\n");
    println!(
        "{:<14} {:>8} {:>9} {:>9} {:>8} {:>6} note",
        "propeller", "rpm", "thrust[N]", "power[W]", "FM", "ok?"
    );
    let mut best: Option<(&str, f64)> = None;
    for c in &candidates {
        let rpm = operating_rpm(c);
        let (thrust, power, _torque, fm) = perf(c, rpm);
        let ok = thrust >= THRUST_PER_ROTOR_N;
        let note = if ok {
            "meets thrust"
        } else {
            "short of thrust"
        };
        println!(
            "{:<14} {rpm:>8.0} {thrust:>9.3} {power:>9.2} {fm:>8.3} {:>6} {note}",
            c.name,
            if ok { "yes" } else { "no" }
        );
        if ok && best.is_none_or(|(_, p)| power < p) {
            best = Some((c.name, power));
        }
    }
    match best {
        Some((name, power)) => println!(
            "\nLowest-power propeller that meets the thrust: {name} ({power:.2} W per rotor)."
        ),
        None => println!(
            "\nNo candidate meets the thrust: try a larger radius or a higher motor speed."
        ),
    }
}
