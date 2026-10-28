//! Nonlinear six-degree-of-freedom simulation: equilibrium, damping and a
//! steady coordinated turn, each checked against the trim solvers.

use tpt_aero_core::{Vec3, WingPlanform, WingSection};
use tpt_aero_stability::sixdof::{
    simulate, simulate_with_loads, Controls, RigidBody, SixDofConfig, SixDofError, State,
};
use tpt_aero_stability::trim::{trim_level_flight, TrimRequest};
use tpt_aero_stability::turn::{coordinated_turn_rate, trim_steady_turn, TurnRequest, G};
use tpt_aero_vlm::{ControlSymmetry, Surface, VlmModel};

const RHO: f64 = 1.225;
const MASS: f64 = 600.0;
const V: f64 = 50.0;
const CD0: f64 = 0.03;

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

fn config() -> SixDofConfig<'static> {
    SixDofConfig {
        body: RigidBody {
            mass: MASS,
            ixx: 800.0,
            iyy: 900.0,
            izz: 1200.0,
            ixz: 0.0,
        },
        density: RHO,
        parasite_cd: CD0,
        elevator: "elevator",
        aileron: "aileron",
        rudder: "rudder",
        wind_ned: [0.0, 0.0, 0.0],
        thrust_body_fixed: false,
        isa_density: false,
    }
}

#[test]
fn isa_density_matches_a_constant_run_at_the_same_altitude() {
    // One step from altitude 1000 m: the ISA density there must give the same
    // loads as passing that density in as the constant.
    let m = aircraft();
    let t = trim_level_flight(
        &m,
        &TrimRequest {
            weight_n: MASS * G,
            v_inf: V,
            density: RHO,
            elevator: "elevator",
            cd0: CD0,
        },
    )
    .unwrap();
    let x0 = steady_state(t.alpha_rad, 0.0, 0.0);
    let controls = Controls {
        elevator: t.elevator_rad,
        ..Controls::default()
    };
    let rho_1000 = tpt_aero_atmo::isa::at(1000.0).unwrap().density;
    let isa = SixDofConfig {
        isa_density: true,
        ..config()
    };
    let constant = SixDofConfig {
        density: rho_1000,
        ..config()
    };
    let run = |cfg: &SixDofConfig<'_>| {
        simulate(&m, cfg, x0, 1e-3, 1, &|_| controls, &|_, _| t.thrust_n).unwrap()
    };
    let (a, b) = (run(&isa), run(&constant));
    let (sa, sb) = (a[1].1, b[1].1);
    for (x, y) in [(sa.u, sb.u), (sa.w, sb.w), (sa.q, sb.q)] {
        assert!((x - y).abs() < 1e-12, "{x} vs {y}");
    }
}

#[test]
fn isa_density_refuses_altitude_outside_the_model() {
    let m = aircraft();
    let high = State {
        altitude: 90_000.0,
        ..steady_state(0.05, 0.0, 0.0)
    };
    let isa = SixDofConfig {
        isa_density: true,
        ..config()
    };
    let err = simulate(
        &m,
        &isa,
        high,
        1e-3,
        1,
        &|_| Controls::default(),
        &|_, _| 0.0,
    )
    .unwrap_err();
    assert_eq!(err, SixDofError::AltitudeOutOfRange);
}

/// Velocity and attitude of steady flight at `alpha`, `beta` and bank `phi`
/// with flight-path angle zero.
fn steady_state(alpha: f64, beta: f64, phi: f64) -> State {
    let (sa, ca) = alpha.sin_cos();
    let (sb, cb) = beta.sin_cos();
    let u = V * ca * cb;
    let v = V * sb;
    let w = V * sa * cb;
    // Level flight: u sin θ − v sin φ cos θ − w cos φ cos θ = 0.
    let theta = ((v * phi.sin() + w * phi.cos()) / u).atan();
    State {
        u,
        v,
        w,
        p: 0.0,
        q: 0.0,
        r: 0.0,
        phi,
        theta,
        psi: 0.0,
        north: 0.0,
        east: 0.0,
        altitude: 1000.0,
    }
}

#[test]
fn level_trim_is_an_equilibrium_of_the_six_dof_dynamics() {
    let m = aircraft();
    let t = trim_level_flight(
        &m,
        &TrimRequest {
            weight_n: MASS * G,
            v_inf: V,
            density: RHO,
            elevator: "elevator",
            cd0: CD0,
        },
    )
    .unwrap();
    let x0 = steady_state(t.alpha_rad, 0.0, 0.0);
    let controls = Controls {
        elevator: t.elevator_rad,
        ..Controls::default()
    };
    let hist = simulate(&m, &config(), x0, 0.01, 500, &|_| controls, &|_, _| {
        t.thrust_n
    })
    .unwrap();
    let end = hist.last().unwrap().1;
    assert!((end.u - x0.u).abs() < 0.05, "Δu = {}", end.u - x0.u);
    assert!((end.w - x0.w).abs() < 0.05, "Δw = {}", end.w - x0.w);
    assert!(
        (end.theta - x0.theta).abs() < 2e-3,
        "Δθ = {}",
        end.theta - x0.theta
    );
    assert!(
        (end.altitude - x0.altitude).abs() < 0.5,
        "Δh = {}",
        end.altitude - x0.altitude
    );
}

#[test]
fn a_pitch_rate_disturbance_damps_out() {
    // A statically stable aircraft damps a pitch-rate kick through the
    // short-period mode: the peak |q| in the second to third second is well
    // below the initial kick. This checks the sign of the pitch damping.
    let m = aircraft();
    let t = trim_level_flight(
        &m,
        &TrimRequest {
            weight_n: MASS * G,
            v_inf: V,
            density: RHO,
            elevator: "elevator",
            cd0: CD0,
        },
    )
    .unwrap();
    let mut x0 = steady_state(t.alpha_rad, 0.0, 0.0);
    x0.q = 0.05;
    let controls = Controls {
        elevator: t.elevator_rad,
        ..Controls::default()
    };
    let hist = simulate(&m, &config(), x0, 0.005, 600, &|_| controls, &|_, _| {
        t.thrust_n
    })
    .unwrap();
    let peak = hist
        .iter()
        .filter(|(time, _)| (2.0..3.0).contains(time))
        .map(|(_, s)| s.q.abs())
        .fold(0.0, f64::max);
    assert!(peak < 0.5 * x0.q, "peak |q| in [2,3) s = {peak}");
}

#[test]
fn coordinated_turn_trim_is_steady_under_six_dof_dynamics() {
    let m = aircraft();
    let bank = 20.0_f64.to_radians();
    let psi_dot = coordinated_turn_rate(V, bank);
    let turn = trim_steady_turn(
        &m,
        &TurnRequest {
            weight_n: MASS * G,
            v_inf: V,
            density: RHO,
            bank_rad: bank,
            turn_rate: psi_dot,
            elevator: "elevator",
            aileron: "aileron",
            rudder: "rudder",
        },
    )
    .unwrap();

    let mut x0 = steady_state(turn.alpha_rad, turn.beta_rad, bank);
    // Steady-turn body rates: p = −ψ̇ sin θ, q = ψ̇ sin φ cos θ, r = ψ̇ cos φ cos θ.
    let (sth, cth) = x0.theta.sin_cos();
    x0.p = -psi_dot * sth;
    x0.q = psi_dot * bank.sin() * cth;
    x0.r = psi_dot * bank.cos() * cth;
    let controls = Controls {
        elevator: turn.elevator_rad,
        aileron: turn.aileron_rad,
        rudder: turn.rudder_rad,
    };
    // Thrust acts along the airspeed vector. Its magnitude is the one that
    // cancels the acceleration along that vector at the start: T = −m (a₀ · V̂),
    // with a₀ the zero-thrust acceleration (one tiny RK4 step recovers it).
    let eps = 1.0e-5;
    let a0 = simulate(&m, &config(), x0, eps, 1, &|_| controls, &|_, _| 0.0).unwrap();
    let (s0, s1) = (a0[0].1, a0[1].1);
    let acc = [
        (s1.u - s0.u) / eps,
        (s1.v - s0.v) / eps,
        (s1.w - s0.w) / eps,
    ];
    let v_hat = [x0.u / V, x0.v / V, x0.w / V];
    let along = acc[0] * v_hat[0] + acc[1] * v_hat[1] + acc[2] * v_hat[2];
    let thrust = -MASS * along;

    let hist = simulate(&m, &config(), x0, 0.01, 300, &|_| controls, &|_, _| thrust).unwrap();
    let end = hist.last().unwrap().1;
    // Bank and yaw rate hold within a few percent over 3 s.
    assert!(
        (end.phi - bank).abs() < 0.03 * bank,
        "φ = {} deg",
        end.phi.to_degrees()
    );
    assert!((end.r - x0.r).abs() < 0.05 * x0.r, "r = {}", end.r);
    assert!((end.q - x0.q).abs() < 0.05 * x0.q.abs(), "q = {}", end.q);
    assert!(
        (end.altitude - x0.altitude).abs() < 1.0,
        "Δh = {}",
        end.altitude - x0.altitude
    );
}

#[test]
fn invalid_configurations_are_rejected() {
    use tpt_aero_stability::sixdof::SixDofError;
    let m = aircraft();
    let mut cfg = config();
    cfg.body.mass = -1.0;
    let x0 = steady_state(0.05, 0.0, 0.0);
    let err = simulate(&m, &cfg, x0, 0.01, 1, &|_| Controls::default(), &|_, _| 0.0).unwrap_err();
    assert_eq!(err, SixDofError::InvalidConfig);

    let mut bad_name = config();
    bad_name.rudder = "flaperon";
    let err = simulate(
        &m,
        &bad_name,
        x0,
        0.01,
        1,
        &|_| Controls::default(),
        &|_, _| 0.0,
    )
    .unwrap_err();
    assert_eq!(err, SixDofError::NoSuchControl);
}

#[test]
fn an_additive_roll_moment_accelerates_roll_as_m_over_ixx() {
    // A constant reaction torque (for example a propeller's) about the body x
    // axis gives roll acceleration M / I_xx from the trimmed state.
    let m = aircraft();
    let t = trim_level_flight(
        &m,
        &TrimRequest {
            weight_n: MASS * G,
            v_inf: V,
            density: RHO,
            elevator: "elevator",
            cd0: CD0,
        },
    )
    .unwrap();
    let x0 = steady_state(t.alpha_rad, 0.0, 0.0);
    let controls = Controls {
        elevator: t.elevator_rad,
        ..Controls::default()
    };
    let torque = 10.0;
    let dt = 0.001;
    let hist = simulate_with_loads(
        &m,
        &config(),
        x0,
        dt,
        1,
        &|_| controls,
        &|_, _| t.thrust_n,
        &|_, _| ([0.0; 3], [torque, 0.0, 0.0]),
    )
    .unwrap();
    let roll_accel = hist[1].1.p / dt;
    let expected = torque / config().body.ixx;
    assert!(
        (roll_accel / expected - 1.0).abs() < 0.02,
        "p' = {roll_accel}, M/Ixx = {expected}"
    );
}

#[test]
fn zero_gust_history_reproduces_the_steady_wind_run() {
    use tpt_aero_stability::sixdof::{simulate_partial, simulate_partial_with_wind};
    let m = aircraft();
    let x0 = steady_state(4.0_f64.to_radians(), 0.0, 0.0);
    let controls = Controls::default();
    let no_extra = |_: f64, _: &State| ([0.0; 3], [0.0; 3]);
    let (steady, e1) = simulate_partial(
        &m,
        &config(),
        x0,
        0.01,
        200,
        &|_| controls,
        &|_, _| 0.0,
        &no_extra,
    );
    let (still, e2) = simulate_partial_with_wind(
        &m,
        &config(),
        x0,
        0.01,
        200,
        &|_| controls,
        &|_, _| 0.0,
        &no_extra,
        &|_| [0.0; 3],
    );
    assert!(e1.is_none() && e2.is_none());
    assert_eq!(steady.last(), still.last());
}

#[test]
fn an_updraft_gust_changes_the_response() {
    use tpt_aero_stability::sixdof::{
        one_minus_cosine_gust, simulate_partial, simulate_partial_with_wind,
    };
    let m = aircraft();
    let x0 = steady_state(4.0_f64.to_radians(), 0.0, 0.0);
    let controls = Controls::default();
    let no_extra = |_: f64, _: &State| ([0.0; 3], [0.0; 3]);
    let (calm, _) = simulate_partial(
        &m,
        &config(),
        x0,
        0.01,
        300,
        &|_| controls,
        &|_, _| 0.0,
        &no_extra,
    );
    let gust = one_minus_cosine_gust(3.0, 10.0, 0.5, V).expect("valid gust");
    let (gusty, e) = simulate_partial_with_wind(
        &m,
        &config(),
        x0,
        0.01,
        300,
        &|_| controls,
        &|_, _| 0.0,
        &no_extra,
        &|t| gust(t),
    );
    assert!(e.is_none(), "{e:?}");
    // The gust is brief and the short period damps it, so compare the largest
    // deviation over the whole run, not the final sample.
    let peak_w_gap = calm
        .iter()
        .zip(&gusty)
        .map(|(a, b)| (a.1.w - b.1.w).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        peak_w_gap > 0.01,
        "gust had no effect: peak Δw = {peak_w_gap}"
    );
}

#[test]
fn one_minus_cosine_gust_peaks_at_u_ds_and_is_zero_outside_its_support() {
    use tpt_aero_stability::sixdof::one_minus_cosine_gust;
    let (u_ds, h, t0, v) = (3.0, 10.0, 0.5, 50.0);
    let gust = one_minus_cosine_gust(u_ds, h, t0, v).unwrap();
    // Penetration s = h, reached at t0 + h / v: NED down = −u_ds (updraft).
    assert!((gust(t0 + h / v)[2] + u_ds).abs() < 1e-12);
    // Zero before the gust starts and after penetration 2h.
    assert_eq!(gust(t0 - 0.1)[2], 0.0);
    assert!(gust(t0 + 2.0 * h / v + 0.01)[2].abs() < 1e-12);
    // Invalid inputs are refused.
    assert!(one_minus_cosine_gust(-1.0, h, t0, v).is_none());
    assert!(one_minus_cosine_gust(u_ds, 0.0, t0, v).is_none());
    assert!(one_minus_cosine_gust(u_ds, h, t0, 0.0).is_none());
}
