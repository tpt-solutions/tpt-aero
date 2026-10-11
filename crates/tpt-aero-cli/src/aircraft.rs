//! Aircraft case files: lifting surfaces with control surfaces, for steady
//! level and coordinated-turn trim.

use serde::{Deserialize, Serialize};
use tpt_aero::core::WingPlanform;
use tpt_aero::stability::trim::{trim_level_flight, TrimError, TrimRequest};
use tpt_aero::stability::turn::{coordinated_turn_rate, trim_steady_turn, TurnRequest, G};
use tpt_aero::vlm::{ControlSymmetry, Surface, VlmModel};

use crate::{finite, positive, CliError};

/// A fully commented starting point for `tpt-aero trim`.
pub const AIRCRAFT_TEMPLATE: &str = r#"# tpt-aero aircraft case. Run: tpt-aero trim aircraft.toml --speed 50 [--bank 20]

[mass]
mass = 600.0          # kg
cd0 = 0.03            # parasite drag coefficient (sets the thrust needed)

[[surface]]           # main wing
kind = "horizontal"   # or "vertical" (fins)
span = 10.0           # m
root_chord = 1.2      # m
tip_chord = 1.2       # m (default: root_chord)
sweep_le = 0.0        # m, aft offset of the tip leading edge
dihedral_deg = 5.0    # a little dihedral gives the wing its own lateral force, so
                      # sideslip and rudder can be told apart in a turn trim
n_span = 20
n_chord = 3
[[surface.control]]
name = "aileron"
hinge = 0.75          # control starts at this chord fraction and runs to the TE
inboard = 2.5         # m, spanwise extent from the centreline
outboard = 4.8
symmetry = "antisymmetric"

[[surface]]           # horizontal tail
kind = "horizontal"   # keep it clear of the fin: a tail at z = 0 meets the fin
                      # root, and the sideslip derivatives then go unreliable
span = 4.0
root_chord = 0.6
tip_chord = 0.6
x = 4.0               # m aft of the wing's reference point
n_span = 10
n_chord = 3
[[surface.control]]
name = "elevator"
hinge = 0.6
inboard = 0.0
outboard = 2.0
symmetry = "symmetric"

[[surface]]           # fin
kind = "vertical"
span = 1.5
root_chord = 0.9
tip_chord = 0.5
x = 3.8
n_span = 6
n_chord = 3
[[surface.control]]
name = "rudder"
hinge = 0.6
inboard = 0.0
outboard = 1.5
symmetry = "symmetric"
"#;

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct AircraftCase {
    mass: MassSpec,
    #[serde(rename = "surface")]
    surfaces: Vec<SurfaceSpec>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct MassSpec {
    mass: f64,
    #[serde(default = "default_cd0")]
    cd0: f64,
}

fn default_cd0() -> f64 {
    0.03
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct SurfaceSpec {
    #[serde(default = "default_kind")]
    kind: String,
    span: f64,
    root_chord: f64,
    tip_chord: Option<f64>,
    #[serde(default)]
    sweep_le: f64,
    #[serde(default)]
    dihedral_deg: f64,
    #[serde(default)]
    x: f64,
    #[serde(default)]
    y: f64,
    #[serde(default)]
    z: f64,
    #[serde(default = "default_n_span")]
    n_span: usize,
    #[serde(default = "default_n_chord")]
    n_chord: usize,
    #[serde(default, rename = "control")]
    controls: Vec<ControlSpec>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct ControlSpec {
    name: String,
    hinge: f64,
    inboard: f64,
    outboard: f64,
    symmetry: String,
}

fn default_kind() -> String {
    "horizontal".into()
}

fn default_n_span() -> usize {
    10
}

fn default_n_chord() -> usize {
    3
}

/// What to trim.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TrimKind {
    /// Steady level flight at zero bank.
    Level,
    /// Coordinated turn at the given bank angle (rad).
    Turn {
        /// Bank angle (rad), positive right wing down.
        bank_rad: f64,
    },
}

/// A steady trim point, in degrees and SI units.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct TrimReport {
    /// `"level"` or `"turn"`.
    pub kind: String,
    /// Airspeed (m/s).
    pub speed_ms: f64,
    /// Air density (kg/m³).
    pub density: f64,
    /// Bank angle (deg); zero for level flight.
    pub bank_deg: f64,
    /// Angle of attack (deg).
    pub alpha_deg: f64,
    /// Sideslip (deg); zero for level flight.
    pub beta_deg: f64,
    /// Elevator (deg).
    pub elevator_deg: f64,
    /// Aileron (deg); `None` for level flight.
    pub aileron_deg: Option<f64>,
    /// Rudder (deg); `None` for level flight.
    pub rudder_deg: Option<f64>,
    /// Thrust required (N); level flight only.
    pub thrust_n: Option<f64>,
    /// Load factor (lift over weight).
    pub load_factor: f64,
    /// Turn rate (deg/s); zero for level flight.
    pub turn_rate_deg_s: f64,
    /// Turn radius (m); `None` for level flight.
    pub turn_radius_m: Option<f64>,
}

fn build_model(case: &AircraftCase) -> Result<VlmModel, CliError> {
    let mut surfaces = Vec::with_capacity(case.surfaces.len());
    for (i, sp) in case.surfaces.iter().enumerate() {
        let label = format!("surface {}", i + 1);
        let span = positive(&format!("{label}.span"), sp.span)?;
        let root = positive(&format!("{label}.root_chord"), sp.root_chord)?;
        let tip = match sp.tip_chord {
            Some(t) => positive(&format!("{label}.tip_chord"), t)?,
            None => root,
        };
        let sweep = finite(&format!("{label}.sweep_le"), sp.sweep_le)?;
        let dihedral = finite(&format!("{label}.dihedral_deg"), sp.dihedral_deg)?.to_radians();
        if sp.n_span == 0 || sp.n_chord == 0 {
            return Err(CliError::Invalid(format!(
                "{label}: n_span and n_chord must be ≥ 1"
            )));
        }
        let planform = WingPlanform::tapered(span, root, tip, sweep, dihedral);
        let mut surface = match sp.kind.as_str() {
            "horizontal" => Surface::horizontal(planform, sp.n_span, sp.n_chord),
            "vertical" => Surface::vertical(planform, sp.n_span, sp.n_chord),
            other => {
                return Err(CliError::Invalid(format!(
                    "{label}.kind must be \"horizontal\" or \"vertical\", not \"{other}\""
                )))
            }
        }
        .at(
            finite(&format!("{label}.x"), sp.x)?,
            finite(&format!("{label}.y"), sp.y)?,
            finite(&format!("{label}.z"), sp.z)?,
        );
        for c in &sp.controls {
            let symmetry = match c.symmetry.as_str() {
                "symmetric" => ControlSymmetry::Symmetric,
                "antisymmetric" => ControlSymmetry::Antisymmetric,
                other => {
                    return Err(CliError::Invalid(format!(
                        "control {}: symmetry must be \"symmetric\" or \"antisymmetric\", not \"{other}\"",
                        c.name
                    )))
                }
            };
            let hinge = finite(&format!("control {}.hinge", c.name), c.hinge)?;
            if !(0.0..1.0).contains(&hinge) {
                return Err(CliError::Invalid(format!(
                    "control {}: hinge must be in [0, 1)",
                    c.name
                )));
            }
            surface = surface.with_control(
                &c.name,
                hinge,
                finite(&format!("control {}.inboard", c.name), c.inboard)?,
                finite(&format!("control {}.outboard", c.name), c.outboard)?,
                symmetry,
            );
        }
        surfaces.push(surface);
    }
    VlmModel::from_surfaces(&surfaces).map_err(|e| CliError::Invalid(e.to_string()))
}

/// Trim the aircraft described by `case_toml` at `speed` (m/s) and `altitude` (m).
///
/// # Errors
/// [`CliError::Parse`] for malformed TOML, [`CliError::Invalid`] for unusable
/// values or a missing control (level flight needs `elevator`; a turn needs
/// `elevator`, `aileron` and `rudder`), and [`CliError::Solve`] if the trim
/// does not converge.
pub fn run_trim(
    case_toml: &str,
    speed: f64,
    altitude: f64,
    kind: TrimKind,
) -> Result<TrimReport, CliError> {
    let case: AircraftCase =
        toml::from_str(case_toml).map_err(|e| CliError::Parse(e.to_string()))?;
    let speed = positive("speed", speed)?;
    let mass = positive("mass.mass", case.mass.mass)?;
    let cd0 = finite("mass.cd0", case.mass.cd0)?;
    if cd0 < 0.0 {
        return Err(CliError::Invalid("mass.cd0 must be ≥ 0".into()));
    }
    let density = tpt_aero::isa::at(altitude)
        .map_err(|e| CliError::Invalid(e.to_string()))?
        .density;
    let model = build_model(&case)?;
    let weight = mass * G;

    match kind {
        TrimKind::Level => {
            let t = trim_level_flight(
                &model,
                &TrimRequest {
                    weight_n: weight,
                    v_inf: speed,
                    density,
                    elevator: "elevator",
                    cd0,
                },
            )
            .map_err(trim_failed)?;
            Ok(TrimReport {
                kind: "level".into(),
                speed_ms: speed,
                density,
                bank_deg: 0.0,
                alpha_deg: t.alpha_rad.to_degrees(),
                beta_deg: 0.0,
                elevator_deg: t.elevator_rad.to_degrees(),
                aileron_deg: None,
                rudder_deg: None,
                thrust_n: Some(t.thrust_n),
                load_factor: 1.0,
                turn_rate_deg_s: 0.0,
                turn_radius_m: None,
            })
        }
        TrimKind::Turn { bank_rad } => {
            let bank = finite("bank", bank_rad)?;
            let psi_dot = coordinated_turn_rate(speed, bank);
            let t = trim_steady_turn(
                &model,
                &TurnRequest {
                    weight_n: weight,
                    v_inf: speed,
                    density,
                    bank_rad: bank,
                    turn_rate: psi_dot,
                    elevator: "elevator",
                    aileron: "aileron",
                    rudder: "rudder",
                },
            )
            .map_err(trim_failed)?;
            Ok(TrimReport {
                kind: "turn".into(),
                speed_ms: speed,
                density,
                bank_deg: bank.to_degrees(),
                alpha_deg: t.alpha_rad.to_degrees(),
                beta_deg: t.beta_rad.to_degrees(),
                elevator_deg: t.elevator_rad.to_degrees(),
                aileron_deg: Some(t.aileron_rad.to_degrees()),
                rudder_deg: Some(t.rudder_rad.to_degrees()),
                thrust_n: None,
                load_factor: t.load_factor,
                turn_rate_deg_s: psi_dot.to_degrees(),
                turn_radius_m: (psi_dot.abs() > 1.0e-12).then(|| speed / psi_dot.abs()),
            })
        }
    }
}

fn trim_failed(e: TrimError) -> CliError {
    match e {
        TrimError::NoSuchControl => CliError::Invalid(
            "the aircraft case lacks a control the trim needs (level flight: elevator;              turns: elevator, aileron and rudder)"
                .into(),
        ),
        other => CliError::Solve(other.to_string()),
    }
}

/// Human-readable summary of a trim report.
#[must_use]
pub fn trim_text(r: &TrimReport) -> String {
    let mut out = format!(
        "kind                {}\nspeed_ms            {:.3}\ndensity             {:.5}\nbank_deg            {:.4}\nalpha_deg           {:.4}\nbeta_deg            {:.4}\nelevator_deg        {:.4}\n",
        r.kind, r.speed_ms, r.density, r.bank_deg, r.alpha_deg, r.beta_deg, r.elevator_deg
    );
    if let Some(a) = r.aileron_deg {
        out.push_str(&format!("aileron_deg         {a:.4}\n"));
    }
    if let Some(d) = r.rudder_deg {
        out.push_str(&format!("rudder_deg          {d:.4}\n"));
    }
    if let Some(t) = r.thrust_n {
        out.push_str(&format!("thrust_n            {t:.3}\n"));
    }
    out.push_str(&format!("load_factor         {:.5}\n", r.load_factor));
    out.push_str(&format!("turn_rate_deg_s     {:.4}\n", r.turn_rate_deg_s));
    if let Some(rad) = r.turn_radius_m {
        out.push_str(&format!("turn_radius_m       {rad:.2}\n"));
    }
    out
}
