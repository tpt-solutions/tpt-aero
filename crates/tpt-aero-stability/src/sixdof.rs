//! Nonlinear six-degree-of-freedom rigid-body simulation with VLM aerodynamics.
//!
//! Body axes are the aircraft convention: `x` forward, `y` right, `z` down.
//! The VLM works in model axes (`x` aft, `y` right, `z` up); the two differ by
//! a 180° rotation about `y`, `R = diag(−1, 1, −1)`, which maps forces and
//! moments as well as vectors. The centre of gravity is taken to be the
//! model's moment reference point.
//!
//! Aerodynamic forces and moments come from
//! [`tpt_aero_vlm::PreparedVlm::solve_state`] at the instantaneous angle of
//! attack, sideslip and body rates, with the control deflections applied as
//! surface twists. The prepared (cached) model is rebuilt only when the
//! deflections change. Parasite drag is added as a constant `C_D0`, and
//! thrust is supplied by a closure, so a propeller model or any other
//! propulsion can plug in there. Thrust acts along the airspeed vector (the
//! flight path), which is the assumption behind the steady trim solvers in
//! this crate. A body-fixed thrust line would add a moment and a vertical
//! force that those solvers do not balance.
//!
//! Integration is fixed-step RK4 on the 12-state Euler-angle formulation
//! (Stevens & Lewis).

use alloc::vec;
use alloc::vec::Vec;

use tpt_aero_vlm::{FlightState, InfluenceKernel, PreparedVlm, VlmError, VlmModel};

use crate::turn::G;

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// Additional force (N) and moment (N·m) in aircraft axes at the centre of gravity.
pub type ExtraLoads = ([f64; 3], [f64; 3]);

/// Rigid-body mass properties (kg, kg·m²). `ixz` is the product of inertia
/// in the Stevens & Lewis sign convention, `I_xz = ∫ x z dm`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RigidBody {
    /// Mass (kg).
    pub mass: f64,
    /// Roll inertia `I_xx`.
    pub ixx: f64,
    /// Pitch inertia `I_yy`.
    pub iyy: f64,
    /// Yaw inertia `I_zz`.
    pub izz: f64,
    /// Product of inertia `I_xz`.
    pub ixz: f64,
}

/// Aircraft configuration for a simulation.
#[derive(Debug, Clone, Copy)]
pub struct SixDofConfig<'a> {
    /// Mass properties.
    pub body: RigidBody,
    /// Air density (kg/m³), constant over the run.
    pub density: f64,
    /// Parasite drag coefficient referenced to the model's reference area,
    /// added to the VLM's induced drag. Fuselage and other non-lifting drag
    /// enter here.
    pub parasite_cd: f64,
    /// Symmetric pitch control name.
    pub elevator: &'a str,
    /// Antisymmetric roll control name.
    pub aileron: &'a str,
    /// Yaw control name.
    pub rudder: &'a str,
    /// Steady wind of the air mass, in north–east–down axes (m/s). The zero
    /// vector is still air. Airspeed, α and β are taken relative to the air.
    pub wind_ned: [f64; 3],
    /// Thrust along the body x axis (`true`) rather than along the air-relative
    /// velocity (`false`, the default). Use `true` for a propeller on the nose.
    pub thrust_body_fixed: bool,
}

/// Control deflections (rad), VLM sign convention (see
/// [`tpt_aero_vlm::Control`]).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Controls {
    /// Elevator (positive = trailing edge down).
    pub elevator: f64,
    /// Aileron (positive = rolls left).
    pub aileron: f64,
    /// Rudder (positive = yaws right).
    pub rudder: f64,
}

/// The 12-state rigid-body vector in aircraft axes.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct State {
    /// Forward body velocity `u` (m/s).
    pub u: f64,
    /// Lateral body velocity `v` (m/s), positive right.
    pub v: f64,
    /// Vertical body velocity `w` (m/s), positive down.
    pub w: f64,
    /// Roll rate `p` (rad/s).
    pub p: f64,
    /// Pitch rate `q` (rad/s).
    pub q: f64,
    /// Yaw rate `r` (rad/s).
    pub r: f64,
    /// Bank angle `φ` (rad).
    pub phi: f64,
    /// Pitch angle `θ` (rad).
    pub theta: f64,
    /// Heading `ψ` (rad).
    pub psi: f64,
    /// North position (m).
    pub north: f64,
    /// East position (m).
    pub east: f64,
    /// Altitude (m), positive up.
    pub altitude: f64,
}

impl State {
    fn to_array(self) -> [f64; 12] {
        [
            self.u,
            self.v,
            self.w,
            self.p,
            self.q,
            self.r,
            self.phi,
            self.theta,
            self.psi,
            self.north,
            self.east,
            self.altitude,
        ]
    }

    fn from_array(a: [f64; 12]) -> Self {
        Self {
            u: a[0],
            v: a[1],
            w: a[2],
            p: a[3],
            q: a[4],
            r: a[5],
            phi: a[6],
            theta: a[7],
            psi: a[8],
            north: a[9],
            east: a[10],
            altitude: a[11],
        }
    }

    /// True airspeed (m/s).
    #[must_use]
    pub fn airspeed(&self) -> f64 {
        (self.u * self.u + self.v * self.v + self.w * self.w).sqrt()
    }
}

/// Why a simulation stopped or could not start.
#[derive(Debug, Clone, PartialEq)]
pub enum SixDofError {
    /// Non-positive mass or inertia, or a non-finite configuration value.
    InvalidConfig,
    /// The model has no control with one of the configured names.
    NoSuchControl,
    /// The airspeed is zero or the state is not finite.
    InvalidState,
    /// The pitch angle is within about 0.06° of ±90°, where the Euler-angle
    /// kinematics are singular. Stop the run or switch to quaternions.
    GimbalLock,
    /// The VLM solve failed.
    Vlm(VlmError),
}

impl core::fmt::Display for SixDofError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidConfig => {
                write!(f, "mass properties and density must be finite and positive")
            }
            Self::NoSuchControl => write!(f, "model has no control with a configured name"),
            Self::InvalidState => write!(f, "airspeed must be positive and the state finite"),
            Self::GimbalLock => write!(f, "pitch is near ±90°: Euler angles are singular"),
            Self::Vlm(e) => write!(f, "VLM solve failed: {e}"),
        }
    }
}

impl core::error::Error for SixDofError {}

impl From<VlmError> for SixDofError {
    fn from(e: VlmError) -> Self {
        Self::Vlm(e)
    }
}

/// A north–east–down vector expressed in body axes, for a ZYX (yaw, pitch,
/// roll) attitude `(phi, theta, psi)`: `v_b = R(φ,θ,ψ) · v_ned`.
fn ned_to_body(v: &[f64; 3], phi: f64, theta: f64, psi: f64) -> [f64; 3] {
    let (sph, cph) = phi.sin_cos();
    let (sth, cth) = theta.sin_cos();
    let (sps, cps) = psi.sin_cos();
    let r = [
        [cth * cps, cth * sps, -sth],
        [
            sph * sth * cps - cph * sps,
            sph * sth * sps + cph * cps,
            sph * cth,
        ],
        [
            cph * sth * cps + sph * sps,
            cph * sth * sps - sph * cps,
            cph * cth,
        ],
    ];
    core::array::from_fn(|i| r[i][0] * v[0] + r[i][1] * v[1] + r[i][2] * v[2])
}

/// Aircraft-axis dynamics with a prepared VLM, rebuilt when controls change.
struct Aero<'m> {
    model: &'m VlmModel,
    body: RigidBody,
    density: f64,
    parasite_cd: f64,
    ie: usize,
    ia: usize,
    ir: usize,
    key: Option<Controls>,
    prepared: Option<PreparedVlm<'m>>,
    /// Horseshoe geometry, built once. A control change then costs the normals
    /// and one LU, not a fresh N² kernel evaluation.
    kernel: InfluenceKernel,
    q_s: f64,
    /// Inverse of the inertia tensor, computed once from the body properties.
    inv_inertia: [[f64; 3]; 3],
    /// Thrust along body x rather than the relative wind.
    thrust_body_fixed: bool,
}

/// Pitch cosine below which the Euler-angle kinematics are refused (about
/// 0.06° from the gimbal-lock attitude).
const GIMBAL_COS_MIN: f64 = 1.0e-3;

/// Inverse of `I = [[Ixx, 0, −Ixz], [0, Iyy, 0], [−Ixz, 0, Izz]]`, in closed
/// form. Requires `Ixx·Izz > Ixz²` (checked by the caller).
fn inverse_inertia(b: &RigidBody) -> [[f64; 3]; 3] {
    let (a, bb, c, d) = (b.ixx, b.iyy, b.ixz, b.izz);
    let det_xz = a * d - c * c;
    [
        [d / det_xz, 0.0, c / det_xz],
        [0.0, 1.0 / bb, 0.0],
        [c / det_xz, 0.0, a / det_xz],
    ]
}

impl<'m> Aero<'m> {
    fn new(model: &'m VlmModel, cfg: &SixDofConfig<'_>) -> Result<Self, SixDofError> {
        let positive = |x: f64| x.is_finite() && x > 0.0;
        let b = cfg.body;
        if !(positive(b.mass)
            && positive(b.ixx)
            && positive(b.iyy)
            && positive(b.izz)
            && positive(cfg.density)
            && cfg.parasite_cd.is_finite()
            && cfg.parasite_cd >= 0.0
            && b.ixz.is_finite())
        {
            return Err(SixDofError::InvalidConfig);
        }
        // Positive definite inertia needs Ixx·Izz > Ixz² (the determinant of
        // the x–z block). Without it the inverse does not exist.
        if b.ixx * b.izz <= b.ixz * b.ixz {
            return Err(SixDofError::InvalidConfig);
        }
        Ok(Self {
            model,
            body: b,
            density: cfg.density,
            parasite_cd: cfg.parasite_cd,
            ie: model
                .control_index(cfg.elevator)
                .ok_or(SixDofError::NoSuchControl)?,
            ia: model
                .control_index(cfg.aileron)
                .ok_or(SixDofError::NoSuchControl)?,
            ir: model
                .control_index(cfg.rudder)
                .ok_or(SixDofError::NoSuchControl)?,
            key: None,
            prepared: None,
            kernel: model.influence_kernel()?,
            q_s: 0.0,
            inv_inertia: inverse_inertia(&b),
            thrust_body_fixed: cfg.thrust_body_fixed,
        })
    }

    /// Prepared VLM for `controls`, rebuilt when they differ from the last call.
    fn prepared(&mut self, controls: Controls) -> Result<&PreparedVlm<'m>, SixDofError> {
        if self.key != Some(controls) || self.prepared.is_none() {
            let mut deflections = vec![0.0; self.model.controls().len()];
            deflections[self.ie] = controls.elevator;
            deflections[self.ia] = controls.aileron;
            deflections[self.ir] = controls.rudder;
            let twists = self.model.control_twists(&deflections)?;
            self.prepared = Some(self.model.prepare_with_kernel(&self.kernel, &twists)?);
            self.key = Some(controls);
        }
        Ok(self.prepared.as_ref().expect("set above"))
    }

    /// Time derivative of the 12-state vector, with the air moving at `wind_ned`
    /// (NED, m/s) at this instant.
    fn derivative(
        &mut self,
        x: &[f64; 12],
        controls: Controls,
        thrust: f64,
        extra: ExtraLoads,
        wind_ned: [f64; 3],
    ) -> Result<[f64; 12], SixDofError> {
        let s = State::from_array(*x);
        // Air-relative body velocity: the aircraft velocity minus the wind,
        // which is given in NED and rotated into body axes (ZYX Euler).
        let wind_b = ned_to_body(&wind_ned, s.phi, s.theta, s.psi);
        let (u_r, v_r, w_r) = (s.u - wind_b[0], s.v - wind_b[1], s.w - wind_b[2]);
        let v_air = (u_r * u_r + v_r * v_r + w_r * w_r).sqrt();
        if !(v_air.is_finite() && v_air > 0.0 && x.iter().all(|a| a.is_finite())) {
            return Err(SixDofError::InvalidState);
        }
        let alpha = w_r.atan2(u_r);
        let beta = (v_r / v_air).clamp(-1.0, 1.0).asin();
        let (sa, ca) = alpha.sin_cos();
        let q_dyn = 0.5 * self.density * v_air * v_air;
        self.q_s = q_dyn * self.model.reference().area;

        // Model-frame angular velocity ω_m = R (p, q, r) = (−p, q, −r), projected
        // onto the stability axes that `solve_state` uses.
        let omega_m = [-s.p, s.q, -s.r];
        let forward = [-ca, 0.0, -sa];
        let down = [sa, 0.0, -ca];
        let dot = |a: &[f64; 3], b: &[f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let flight = FlightState {
            alpha,
            beta,
            p: dot(&forward, &omega_m),
            q: omega_m[1],
            r: dot(&down, &omega_m),
        };
        let rho = self.density;
        let prepared = self.prepared(controls)?;
        let loads = prepared.solve_state(&flight, v_air, rho)?;

        // Model → aircraft axes, R = diag(−1, 1, −1), for forces and moments.
        let [fx, fy, fz] = loads.force_n;
        let [mx, my, mz] = loads.moment_nm;
        let mass = self.body.mass;
        let drag_cd0 = self.parasite_cd * self.q_s;
        // Parasite drag acts along the air-relative velocity. Thrust acts along
        // that velocity too, unless the configuration fixes it to body x.
        let (ux, uy, uz) = (u_r / v_air, v_r / v_air, w_r / v_air);
        let (tx, ty, tz) = if self.thrust_body_fixed {
            (1.0, 0.0, 0.0)
        } else {
            (ux, uy, uz)
        };
        let force = [
            -fx - drag_cd0 * ux + thrust * tx,
            fy - drag_cd0 * uy + thrust * ty,
            -fz - drag_cd0 * uz + thrust * tz,
        ];
        let moment = [-mx, my, -mz];
        // Additive loads (fuselage, propeller reaction, ...), aircraft axes.
        let force: [f64; 3] = core::array::from_fn(|i| force[i] + extra.0[i]);
        let moment: [f64; 3] = core::array::from_fn(|i| moment[i] + extra.1[i]);

        // Gravity in body axes.
        let (sphi, cphi) = s.phi.sin_cos();
        let (sth, cth) = s.theta.sin_cos();
        if cth.abs() < GIMBAL_COS_MIN {
            return Err(SixDofError::GimbalLock);
        }
        let grav = [-G * sth, G * cth * sphi, G * cth * cphi];

        let du = s.r * s.v - s.q * s.w + force[0] / mass + grav[0];
        let dv = s.p * s.w - s.r * s.u + force[1] / mass + grav[1];
        let dw = s.q * s.u - s.p * s.v + force[2] / mass + grav[2];

        // I ω̇ = M − ω × (I ω), with I = [[Ixx,0,−Ixz],[0,Iyy,0],[−Ixz,0,Izz]].
        let b = self.body;
        let inertia = [
            [b.ixx, 0.0, -b.ixz],
            [0.0, b.iyy, 0.0],
            [-b.ixz, 0.0, b.izz],
        ];
        let omega = [s.p, s.q, s.r];
        let i_omega = [
            inertia[0][0] * omega[0] + inertia[0][2] * omega[2],
            inertia[1][1] * omega[1],
            inertia[2][0] * omega[0] + inertia[2][2] * omega[2],
        ];
        let gyro = [
            omega[1] * i_omega[2] - omega[2] * i_omega[1],
            omega[2] * i_omega[0] - omega[0] * i_omega[2],
            omega[0] * i_omega[1] - omega[1] * i_omega[0],
        ];
        let rhs: [f64; 3] = core::array::from_fn(|i| moment[i] - gyro[i]);
        let inv = &self.inv_inertia;
        let omega_dot: [f64; 3] =
            core::array::from_fn(|i| (0..3).map(|j| inv[i][j] * rhs[j]).sum());

        // Euler kinematics (ZYX).
        let tan_th = sth / cth;
        let dphi = s.p + tan_th * (s.q * sphi + s.r * cphi);
        let dtheta = s.q * cphi - s.r * sphi;
        let dpsi = (s.q * sphi + s.r * cphi) / cth;

        // Navigation: body velocity → NED, and altitude rate (up positive).
        let (spsi, cpsi) = s.psi.sin_cos();
        let dnorth = s.u * cth * cpsi
            + s.v * (sphi * sth * cpsi - cphi * spsi)
            + s.w * (cphi * sth * cpsi + sphi * spsi);
        let deast = s.u * cth * spsi
            + s.v * (sphi * sth * spsi + cphi * cpsi)
            + s.w * (cphi * sth * spsi - sphi * cpsi);
        let daltitude = s.u * sth - s.v * sphi * cth - s.w * cphi * cth;

        Ok([
            du,
            dv,
            dw,
            omega_dot[0],
            omega_dot[1],
            omega_dot[2],
            dphi,
            dtheta,
            dpsi,
            dnorth,
            deast,
            daltitude,
        ])
    }
}

/// Largest |α| and |β| (rad) for which the linear-in-state VLM loads are
/// trusted: the same 20° limit the trim solvers refuse beyond.
pub const LINEAR_ENVELOPE_LIMIT: f64 = 20.0 * core::f64::consts::PI / 180.0;

/// Air-relative angle of attack and sideslip at `state` in steady `wind_ned`
/// (rad), and whether both lie inside [`LINEAR_ENVELOPE_LIMIT`].
///
/// A six-DOF history sample outside the envelope is still integrated, but its
/// aerodynamic loads are extrapolated beyond the validated range. Check each
/// sample, or the whole run, before quoting the results.
#[must_use]
pub fn alpha_beta_envelope(state: &State, wind_ned: [f64; 3]) -> (f64, f64, bool) {
    let wind_b = ned_to_body(&wind_ned, state.phi, state.theta, state.psi);
    let (u_r, v_r, w_r) = (
        state.u - wind_b[0],
        state.v - wind_b[1],
        state.w - wind_b[2],
    );
    let v_air = (u_r * u_r + v_r * v_r + w_r * w_r).sqrt();
    let alpha = w_r.atan2(u_r);
    let beta = if v_air > 0.0 {
        (v_r / v_air).clamp(-1.0, 1.0).asin()
    } else {
        0.0
    };
    let inside = alpha.abs() <= LINEAR_ENVELOPE_LIMIT && beta.abs() <= LINEAR_ENVELOPE_LIMIT;
    (alpha, beta, inside)
}

/// A vertical one-minus-cosine gust for [`simulate_partial_with_wind`].
///
/// With penetration `s = v_air (t − t_start)`, the gust is zero before `s = 0`
/// and after `s = 2h`, and in between is an updraft of `u_ds/2 · (1 − cos(π s/h))`.
/// It peaks at `u_ds` when `s = h`. The NED down component is the negative of
/// that updraft.
///
/// Returns `None` unless `u_ds ≥ 0`, `h > 0`, `v_air > 0` and all inputs are finite.
#[must_use]
pub fn one_minus_cosine_gust(
    u_ds: f64,
    h: f64,
    t_start: f64,
    v_air: f64,
) -> Option<impl Fn(f64) -> [f64; 3]> {
    if !(u_ds.is_finite() && u_ds >= 0.0 && h.is_finite() && h > 0.0) {
        return None;
    }
    if !(t_start.is_finite() && v_air.is_finite() && v_air > 0.0) {
        return None;
    }
    Some(move |t: f64| {
        let s = v_air * (t - t_start);
        let down = if (0.0..=2.0 * h).contains(&s) {
            -0.5 * u_ds * (1.0 - (core::f64::consts::PI * s / h).cos())
        } else {
            0.0
        };
        [0.0, 0.0, down]
    })
}

/// Simulate `steps` RK4 steps of size `dt` from `x0` at `t = 0`, returning
/// `steps + 1` samples `(t, state)` including the initial state.
///
/// `controls(t)` gives the deflections and `thrust(t, state)` the thrust
/// magnitude (N, along the airspeed vector) at time `t`. Both are sampled at
/// the RK4 stages.
///
/// # Errors
/// [`SixDofError`] for an invalid configuration, unknown control names, a
/// zero or non-finite state during the run, or a VLM failure.
pub fn simulate(
    model: &VlmModel,
    cfg: &SixDofConfig<'_>,
    x0: State,
    dt: f64,
    steps: usize,
    controls: &dyn Fn(f64) -> Controls,
    thrust: &dyn Fn(f64, &State) -> f64,
) -> Result<Vec<(f64, State)>, SixDofError> {
    simulate_with_loads(model, cfg, x0, dt, steps, controls, thrust, &|_, _| {
        ([0.0; 3], [0.0; 3])
    })
}

/// [`simulate`] with additional loads: `extra(t, state)` returns a force (N)
/// and a moment (N·m), both in aircraft axes and acting at the centre of
/// gravity. Use it for a fuselage, a propeller's torque reaction, or any load
/// the VLM does not carry.
///
/// # Errors
/// As [`simulate`].
#[allow(clippy::too_many_arguments)]
pub fn simulate_with_loads(
    model: &VlmModel,
    cfg: &SixDofConfig<'_>,
    x0: State,
    dt: f64,
    steps: usize,
    controls: &dyn Fn(f64) -> Controls,
    thrust: &dyn Fn(f64, &State) -> f64,
    extra: &dyn Fn(f64, &State) -> ExtraLoads,
) -> Result<Vec<(f64, State)>, SixDofError> {
    let (history, error) = simulate_partial(model, cfg, x0, dt, steps, controls, thrust, extra);
    match error {
        None => Ok(history),
        Some(e) => Err(e),
    }
}

/// As [`simulate_with_loads`], but never discards the run: it returns every
/// sample computed before a failure, with the error alongside. The history
/// always starts with the initial state at `t = 0`. A run that stops early
/// (for example at [`SixDofError::GimbalLock`]) therefore keeps its useful part.
///
/// # Panics
/// Never. Invalid inputs are reported in the second element.
#[allow(clippy::too_many_arguments)]
pub fn simulate_partial(
    model: &VlmModel,
    cfg: &SixDofConfig<'_>,
    x0: State,
    dt: f64,
    steps: usize,
    controls: &dyn Fn(f64) -> Controls,
    thrust: &dyn Fn(f64, &State) -> f64,
    extra: &dyn Fn(f64, &State) -> ExtraLoads,
) -> (Vec<(f64, State)>, Option<SixDofError>) {
    simulate_partial_with_wind(model, cfg, x0, dt, steps, controls, thrust, extra, &|_| {
        cfg.wind_ned
    })
}

/// As [`simulate_partial`], with the NED air velocity given as a function of
/// time instead of the steady `cfg.wind_ned`. Use it for gusts: pass
/// [`one_minus_cosine_gust`] or any other wind history. `cfg.wind_ned` is
/// ignored.
///
/// # Panics
/// Never. Invalid inputs are reported in the second element.
#[allow(clippy::too_many_arguments)]
pub fn simulate_partial_with_wind(
    model: &VlmModel,
    cfg: &SixDofConfig<'_>,
    x0: State,
    dt: f64,
    steps: usize,
    controls: &dyn Fn(f64) -> Controls,
    thrust: &dyn Fn(f64, &State) -> f64,
    extra: &dyn Fn(f64, &State) -> ExtraLoads,
    wind: &dyn Fn(f64) -> [f64; 3],
) -> (Vec<(f64, State)>, Option<SixDofError>) {
    // Do not trust `steps` for the allocation; the vector grows as needed.
    let mut history = Vec::with_capacity((steps + 1).min(1 << 16));
    history.push((0.0, x0));
    if !(dt.is_finite() && dt > 0.0) {
        return (history, Some(SixDofError::InvalidState));
    }
    let mut aero = match Aero::new(model, cfg) {
        Ok(a) => a,
        Err(e) => return (history, Some(e)),
    };
    let mut x = x0.to_array();
    for step in 0..steps {
        let t = step as f64 * dt;
        let mut stage = |xs: &[f64; 12], ts: f64| -> Result<[f64; 12], SixDofError> {
            let c = controls(ts);
            let th = thrust(ts, &State::from_array(*xs));
            let ex = extra(ts, &State::from_array(*xs));
            let w = wind(ts);
            aero.derivative(xs, c, th, ex, w)
        };
        let stepped = (|| -> Result<[f64; 12], SixDofError> {
            let k1 = stage(&x, t)?;
            let x2: [f64; 12] = core::array::from_fn(|i| x[i] + 0.5 * dt * k1[i]);
            let k2 = stage(&x2, t + 0.5 * dt)?;
            let x3: [f64; 12] = core::array::from_fn(|i| x[i] + 0.5 * dt * k2[i]);
            let k3 = stage(&x3, t + 0.5 * dt)?;
            let x4: [f64; 12] = core::array::from_fn(|i| x[i] + dt * k3[i]);
            let k4 = stage(&x4, t + dt)?;
            Ok(core::array::from_fn(|i| {
                x[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
            }))
        })();
        match stepped {
            Ok(next) => {
                x = next;
                history.push(((step + 1) as f64 * dt, State::from_array(x)));
            }
            Err(e) => return (history, Some(e)),
        }
    }
    (history, None)
}

#[cfg(test)]
mod wind_tests {
    use super::ned_to_body;

    #[test]
    fn level_flight_maps_north_to_forward_and_heading_changes_it() {
        let north = [10.0, 0.0, 0.0];
        let level = ned_to_body(&north, 0.0, 0.0, 0.0);
        assert!((level[0] - 10.0).abs() < 1e-12 && level[1].abs() < 1e-12);
        // Heading east (ψ = 90°): north lies on the aircraft's left, i.e. −y.
        let east = ned_to_body(&north, 0.0, 0.0, core::f64::consts::FRAC_PI_2);
        assert!(
            east[0].abs() < 1e-12 && (east[1] + 10.0).abs() < 1e-12,
            "{east:?}"
        );
        // Down (+z in NED) is +z in body at level attitude.
        let down = ned_to_body(&[0.0, 0.0, 3.0], 0.0, 0.0, 0.0);
        assert!((down[2] - 3.0).abs() < 1e-12);
    }
}

#[cfg(test)]
mod envelope_tests {
    use super::*;

    fn state_with(u: f64, v: f64, w: f64) -> State {
        State {
            u,
            v,
            w,
            p: 0.0,
            q: 0.0,
            r: 0.0,
            phi: 0.0,
            theta: 0.0,
            psi: 0.0,
            north: 0.0,
            east: 0.0,
            altitude: 0.0,
        }
    }

    #[test]
    fn envelope_flags_alpha_and_beta_beyond_twenty_degrees() {
        let calm = [0.0; 3];
        let (a, b, ok) = alpha_beta_envelope(&state_with(30.0, 0.0, 0.0), calm);
        assert!(a.abs() < 1e-12 && b.abs() < 1e-12 && ok);
        // α = 25° with no sideslip is outside the envelope.
        let (a, _, ok) = alpha_beta_envelope(
            &state_with(30.0, 0.0, 30.0 * 25f64.to_radians().tan()),
            calm,
        );
        assert!((a - 25f64.to_radians()).abs() < 1e-9 && !ok);
        // β = 25° with α = 0 is outside the envelope.
        let (_, b, ok) = alpha_beta_envelope(
            &state_with(30.0, 30.0 * 25f64.to_radians().tan(), 0.0),
            calm,
        );
        assert!((b - 25f64.to_radians()).abs() < 1e-9 && !ok);
    }

    #[test]
    fn envelope_uses_air_relative_velocity_not_ground_velocity() {
        // Ground velocity is 20 m/s along body x, so a ground-speed α would be
        // zero. A 10 m/s updraft (NED down component −10) makes the airspeed
        // point 26.6° up, which is outside the envelope.
        let wind = [0.0, 0.0, -10.0];
        let (a, b, ok) = alpha_beta_envelope(&state_with(20.0, 0.0, 0.0), wind);
        assert!((a - 10f64.atan2(20.0)).abs() < 1e-9, "alpha = {a}");
        assert!(b.abs() < 1e-12 && !ok);
    }
}
