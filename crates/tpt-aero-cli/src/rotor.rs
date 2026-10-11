//! Rotor and propeller case files: one blade-element momentum operating point.

use serde::{Deserialize, Serialize};
use tpt_aero::rotor::{
    BemtSolver, HubLoss, OperatingCondition, RotorGeometry, SectionModel, TipLoss,
};

use crate::{finite, positive, CliError, MAX_MESH_COUNT};

/// A fully commented starting point for `tpt-aero rotor`.
pub const ROTOR_TEMPLATE: &str = r#"# tpt-aero rotor case. Run: tpt-aero rotor rotor.toml [--json]

[rotor]
radius = 0.5          # tip radius (m)
hub_radius = 0.05     # root cutout radius (m)
blades = 2
chord = 0.05          # constant chord (m)
stations = 20         # spanwise BEMT stations
root_twist_deg = 25.0 # geometric twist at the root cutout
tip_twist_deg = 8.0   # geometric twist at the tip
tip_loss = true       # Prandtl tip-loss correction

[section]             # linear 2-D closure: c_l = cl_slope·α (clamped), c_d = cd0 + cd_cl2·c_l²
cl_slope = 6.283
cl_max = 1.4
cd0 = 0.01
cd_cl2 = 0.01

[operating]
rpm = 6000.0          # shaft speed
speed = 0.0           # axial flow speed (m/s): 0 = hover, > 0 = propeller inflow
altitude = 0.0        # m (ISA density); or set `density` in kg/m^3 instead
# density = 1.225
collective_deg = 0.0  # added to the built-in twist
"#;

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RotorCase {
    rotor: RotorSpec,
    #[serde(default)]
    section: SectionSpec,
    operating: OperatingSpec,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RotorSpec {
    radius: f64,
    #[serde(default)]
    hub_radius: f64,
    blades: usize,
    chord: f64,
    #[serde(default = "default_stations")]
    stations: usize,
    #[serde(default)]
    root_twist_deg: f64,
    #[serde(default)]
    tip_twist_deg: f64,
    #[serde(default = "default_true")]
    tip_loss: bool,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct SectionSpec {
    cl_slope: f64,
    cl_max: f64,
    cd0: f64,
    cd_cl2: f64,
}

impl Default for SectionSpec {
    fn default() -> Self {
        Self {
            cl_slope: 2.0 * core::f64::consts::PI,
            cl_max: 1.4,
            cd0: 0.01,
            cd_cl2: 0.01,
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct OperatingSpec {
    rpm: f64,
    #[serde(default)]
    speed: f64,
    altitude: Option<f64>,
    density: Option<f64>,
    #[serde(default)]
    collective_deg: f64,
}

fn default_stations() -> usize {
    20
}

fn default_true() -> bool {
    true
}

/// Rotor performance at one operating point, in coefficient and dimensional form.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct RotorReport {
    /// Thrust (N).
    pub thrust_n: f64,
    /// Shaft torque (N·m).
    pub torque_nm: f64,
    /// Shaft power (W), including climb power.
    pub power_w: f64,
    /// Thrust coefficient `C_T`.
    pub ct: f64,
    /// Torque coefficient `C_Q`.
    pub cq: f64,
    /// Power coefficient `C_P`.
    pub cp: f64,
    /// Figure of merit (hover).
    pub figure_of_merit: f64,
    /// Propeller efficiency (axial flight).
    pub propeller_efficiency: f64,
}

/// Run one BEMT operating point described by a rotor case.
///
/// # Errors
/// [`CliError::Parse`] for malformed TOML, [`CliError::Invalid`] for unusable
/// values, and [`CliError::Solve`] if the BEMT iteration fails.
pub fn run_rotor(case_toml: &str) -> Result<RotorReport, CliError> {
    let case: RotorCase = toml::from_str(case_toml).map_err(|e| CliError::Parse(e.to_string()))?;
    let r = &case.rotor;
    let radius = positive("rotor.radius", r.radius)?;
    let chord = positive("rotor.chord", r.chord)?;
    let hub = finite("rotor.hub_radius", r.hub_radius)?;
    if !(0.0..radius).contains(&hub) {
        return Err(CliError::Invalid(
            "rotor.hub_radius must be in [0, radius)".into(),
        ));
    }
    if r.blades == 0 || r.stations < 2 {
        return Err(CliError::Invalid(
            "rotor.blades must be ≥ 1 and rotor.stations ≥ 2".into(),
        ));
    }
    if r.stations > MAX_MESH_COUNT {
        return Err(CliError::Invalid(format!(
            "rotor.stations must be <= {MAX_MESH_COUNT} (got {})",
            r.stations
        )));
    }
    let rpm = positive("operating.rpm", case.operating.rpm)?;
    let speed = finite("operating.speed", case.operating.speed)?;
    if speed < 0.0 {
        // Negative axial speed is a descent; BEMT accepts small negatives only.
        return Err(CliError::Invalid("operating.speed must be ≥ 0".into()));
    }
    let density = match (case.operating.density, case.operating.altitude) {
        (Some(_), Some(_)) => {
            return Err(CliError::Invalid(
                "set either operating.density or operating.altitude, not both".into(),
            ))
        }
        (Some(d), None) => positive("operating.density", d)?,
        (None, alt) => {
            let h = alt.unwrap_or(0.0);
            tpt_aero::isa::at(h)
                .map_err(|e| CliError::Invalid(e.to_string()))?
                .density
        }
    };
    let collective =
        finite("operating.collective_deg", case.operating.collective_deg)?.to_radians();

    let geometry = RotorGeometry::rectangular(radius, hub, r.blades, chord, r.stations)
        .with_linear_twist(
            finite("rotor.root_twist_deg", r.root_twist_deg)?.to_radians(),
            finite("rotor.tip_twist_deg", r.tip_twist_deg)?.to_radians(),
        );

    let s = &case.section;
    let section = SectionModel::Linear {
        cl_slope: positive("section.cl_slope", s.cl_slope)?,
        cl_max: positive("section.cl_max", s.cl_max)?,
        cd0: finite("section.cd0", s.cd0)?,
        cd_cl2: finite("section.cd_cl2", s.cd_cl2)?,
    };
    let tip_loss = if r.tip_loss {
        TipLoss::Prandtl
    } else {
        TipLoss::None
    };
    let solver = BemtSolver::new(geometry, section, tip_loss, HubLoss::None)
        .map_err(|e| CliError::Invalid(e.to_string()))?;

    let omega = rpm * core::f64::consts::TAU / 60.0;
    let sol = solver
        .solve(&OperatingCondition {
            axial_speed: speed,
            omega,
            density,
            collective,
        })
        .map_err(|e| CliError::Solve(e.to_string()))?;

    let area = core::f64::consts::PI * radius * radius;
    let tip_speed = omega * radius;
    let q_ref = density * area * tip_speed * tip_speed;
    Ok(RotorReport {
        thrust_n: sol.loads.ct * q_ref,
        torque_nm: sol.loads.cq * q_ref * radius,
        power_w: sol.loads.cp * q_ref * tip_speed,
        ct: sol.loads.ct,
        cq: sol.loads.cq,
        cp: sol.loads.cp,
        figure_of_merit: sol.loads.figure_of_merit,
        propeller_efficiency: sol.loads.propeller_efficiency,
    })
}

/// Human-readable summary of a rotor report.
#[must_use]
pub fn rotor_text(report: &RotorReport) -> String {
    format!(
        "thrust_n            {:.4}\ntorque_nm           {:.5}\npower_w             {:.3}\nct                  {:.6}\ncq                  {:.6}\ncp                  {:.6}\nfigure_of_merit     {:.4}\npropeller_efficiency {:.4}\n",
        report.thrust_n,
        report.torque_nm,
        report.power_w,
        report.ct,
        report.cq,
        report.cp,
        report.figure_of_merit,
        report.propeller_efficiency,
    )
}
