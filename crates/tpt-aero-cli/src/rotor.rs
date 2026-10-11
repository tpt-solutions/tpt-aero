//! Rotor and propeller case files: blade-element momentum operating points,
//! and sweeps over shaft speed.

use serde::{Deserialize, Serialize};
use tpt_aero::aero_core::AirfoilPolar;
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
chord = 0.05          # root chord (m)
chord_tip_ratio = 1.0 # tip chord / root chord (1.0 = constant chord)
stations = 20         # spanwise BEMT stations
root_twist_deg = 25.0 # geometric twist at the root cutout
tip_twist_deg = 8.0   # geometric twist at the tip
tip_loss = true       # Prandtl tip-loss correction
hub_loss = false      # Prandtl hub-loss correction

[section]             # linear 2-D closure: c_l = cl_slope·α (clamped), c_d = cd0 + cd_cl2·c_l²
cl_slope = 6.283
cl_max = 1.4
cd0 = 0.01
cd_cl2 = 0.01
# To use a tabulated polar instead, give all three arrays (degrees, -, -):
# alpha_deg = [-10.0, 0.0, 10.0]
# cl = [-0.8, 0.2, 1.2]
# cd = [0.03, 0.01, 0.02]

[operating]
rpm = 6000.0          # shaft speed
speed = 0.0           # axial flow speed (m/s): 0 = hover, > 0 = climb or propeller inflow, < 0 = descent
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
    #[serde(default = "default_one")]
    chord_tip_ratio: f64,
    #[serde(default = "default_stations")]
    stations: usize,
    #[serde(default)]
    root_twist_deg: f64,
    #[serde(default)]
    tip_twist_deg: f64,
    /// Optional chord per station (m), one per `stations`, hub side to tip.
    /// Overrides the linear taper from `chord` and `chord_tip_ratio`.
    chord_profile: Option<Vec<f64>>,
    /// Optional built-in twist per station (deg), one per `stations`, hub side
    /// to tip. Overrides the linear twist from the root and tip twists.
    twist_profile_deg: Option<Vec<f64>>,
    #[serde(default = "default_true")]
    tip_loss: bool,
    #[serde(default)]
    hub_loss: bool,
}

/// Section closure: linear by default, or a tabulated polar when all of
/// `alpha_deg`, `cl` and `cd` are given.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields, default)]
struct SectionSpec {
    cl_slope: f64,
    cl_max: f64,
    cd0: f64,
    cd_cl2: f64,
    alpha_deg: Option<Vec<f64>>,
    cl: Option<Vec<f64>>,
    cd: Option<Vec<f64>>,
}

impl Default for SectionSpec {
    fn default() -> Self {
        Self {
            cl_slope: 2.0 * core::f64::consts::PI,
            cl_max: 1.4,
            cd0: 0.01,
            cd_cl2: 0.01,
            alpha_deg: None,
            cl: None,
            cd: None,
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

fn default_one() -> f64 {
    1.0
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

/// Build the solver and the operating conditions from a parsed case. Returns
/// the solver, the air density, the radius and the collective (rad).
/// Reject a per-station array whose length is not the station count.
fn check_profile_len(name: &str, profile: &[f64], stations: usize) -> Result<(), CliError> {
    if profile.len() != stations {
        return Err(CliError::Invalid(format!(
            "{name} needs one value per station ({stations}), got {}",
            profile.len()
        )));
    }
    Ok(())
}

fn build(case: &RotorCase) -> Result<(BemtSolver, f64, f64, f64), CliError> {
    let r = &case.rotor;
    let radius = positive("rotor.radius", r.radius)?;
    let chord = positive("rotor.chord", r.chord)?;
    let taper = positive("rotor.chord_tip_ratio", r.chord_tip_ratio)?;
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

    let mut geometry = RotorGeometry::rectangular(radius, hub, r.blades, chord, r.stations)
        .with_linear_twist(
            finite("rotor.root_twist_deg", r.root_twist_deg)?.to_radians(),
            finite("rotor.tip_twist_deg", r.tip_twist_deg)?.to_radians(),
        );
    // Linear taper from the root chord to the tip chord.
    let n = geometry.sections.len();
    for (i, sec) in geometry.sections.iter_mut().enumerate() {
        let t = i as f64 / (n - 1).max(1) as f64;
        sec.chord = chord * (1.0 + (taper - 1.0) * t);
    }
    // Per-station arrays override the linear taper and twist, one value per station.
    if let Some(profile) = &r.chord_profile {
        check_profile_len("rotor.chord_profile", profile, n)?;
        for (sec, &c) in geometry.sections.iter_mut().zip(profile) {
            sec.chord = positive("rotor.chord_profile", c)?;
        }
    }
    if let Some(profile) = &r.twist_profile_deg {
        check_profile_len("rotor.twist_profile_deg", profile, n)?;
        for (sec, &t) in geometry.sections.iter_mut().zip(profile) {
            sec.twist = finite("rotor.twist_profile_deg", t)?.to_radians();
        }
    }

    let s = &case.section;
    let section = match (&s.alpha_deg, &s.cl, &s.cd) {
        (Some(alpha), Some(cl), Some(cd)) => {
            let alphas: Vec<f64> = alpha.iter().map(|a| a.to_radians()).collect();
            let cm = vec![0.0; alphas.len()];
            let polar = AirfoilPolar::new(alphas, cl.clone(), cd.clone(), cm, 0.0, 1.0e6)
                .map_err(|e| CliError::Invalid(format!("section polar: {e}")))?;
            SectionModel::Polar { polar, stall: None }
        }
        (None, None, None) => SectionModel::Linear {
            cl_slope: positive("section.cl_slope", s.cl_slope)?,
            cl_max: positive("section.cl_max", s.cl_max)?,
            cd0: finite("section.cd0", s.cd0)?,
            cd_cl2: finite("section.cd_cl2", s.cd_cl2)?,
        },
        _ => {
            return Err(CliError::Invalid(
                "section polar needs all of alpha_deg, cl and cd".into(),
            ))
        }
    };
    let tip_loss = if r.tip_loss {
        TipLoss::Prandtl
    } else {
        TipLoss::None
    };
    let hub_loss = if r.hub_loss {
        HubLoss::Prandtl {
            x_hub: hub / radius,
        }
    } else {
        HubLoss::None
    };
    let solver = BemtSolver::new(geometry, section, tip_loss, hub_loss)
        .map_err(|e| CliError::Invalid(e.to_string()))?;
    Ok((solver, density, radius, collective))
}

/// One operating point at shaft speed `rpm` and axial speed `speed` (m/s).
fn report_at(
    solver: &BemtSolver,
    speed: f64,
    density: f64,
    radius: f64,
    collective: f64,
    rpm: f64,
) -> Result<RotorReport, CliError> {
    let rpm = positive("operating.rpm", rpm)?;
    let speed = finite("operating.speed", speed)?;
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

fn parse(case_toml: &str) -> Result<RotorCase, CliError> {
    toml::from_str(case_toml).map_err(|e| CliError::Parse(e.to_string()))
}

/// Run one BEMT operating point described by a rotor case.
///
/// # Errors
/// [`CliError::Parse`] for malformed TOML, [`CliError::Invalid`] for unusable
/// values, and [`CliError::Solve`] if the BEMT iteration fails.
pub fn run_rotor(case_toml: &str) -> Result<RotorReport, CliError> {
    let case = parse(case_toml)?;
    let (solver, density, radius, collective) = build(&case)?;
    report_at(
        &solver,
        case.operating.speed,
        density,
        radius,
        collective,
        case.operating.rpm,
    )
}

/// The operating variable a sweep steps through; the others stay at the case
/// values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SweepVariable {
    /// Shaft speed in rpm.
    Rpm,
    /// Collective pitch added to the built-in twist, in degrees.
    CollectiveDeg,
    /// Axial (climb or inflow) speed in m/s. Also reports the advance ratio.
    AxialSpeed,
}

/// Sweep the shaft speed over `start..=stop` rpm in steps of `step`, and return
/// CSV: `rpm,thrust_n,torque_nm,power_w,ct,cq,cp`.
///
/// # Errors
/// As [`run_rotor_sweep_over`].
pub fn run_rotor_sweep(
    case_toml: &str,
    start: f64,
    stop: f64,
    step: f64,
) -> Result<String, CliError> {
    run_rotor_sweep_over(case_toml, SweepVariable::Rpm, start, stop, step)
}

/// Sweep `variable` over `start..=stop` in steps of `step` and return CSV.
///
/// The columns are the swept variable, then `thrust_n,torque_nm,power_w,ct,cq,cp`.
/// [`SweepVariable::AxialSpeed`] adds an `advance_ratio` column, `J = π λ_c`,
/// the propeller advance ratio with `D = 2R`.
///
/// # Errors
/// [`CliError::Invalid`] for a non-finite bound, a non-positive step, a reversed
/// range, more than 1000 points, or a non-positive `start` or `stop` in rpm;
/// otherwise as [`run_rotor`].
pub fn run_rotor_sweep_over(
    case_toml: &str,
    variable: SweepVariable,
    start: f64,
    stop: f64,
    step: f64,
) -> Result<String, CliError> {
    let start = finite("sweep start", start)?;
    let stop = finite("sweep stop", stop)?;
    let step = positive("sweep step", step)?;
    if stop < start {
        return Err(CliError::Invalid("sweep stop must be ≥ start".into()));
    }
    if variable == SweepVariable::Rpm {
        positive("sweep start", start)?;
        positive("sweep stop", stop)?;
    }
    let n = ((stop - start) / step + 1.0e-9).floor() as usize + 1;
    if n > 1000 {
        return Err(CliError::Invalid("sweep has more than 1000 points".into()));
    }
    let case = parse(case_toml)?;
    let (solver, density, radius, collective0) = build(&case)?;
    let mut out = match variable {
        SweepVariable::Rpm => String::from("rpm,thrust_n,torque_nm,power_w,ct,cq,cp\n"),
        SweepVariable::CollectiveDeg => {
            String::from("collective_deg,thrust_n,torque_nm,power_w,ct,cq,cp\n")
        }
        SweepVariable::AxialSpeed => {
            String::from("axial_speed_ms,advance_ratio,thrust_n,torque_nm,power_w,ct,cq,cp\n")
        }
    };
    for i in 0..n {
        let x = start + step * i as f64;
        let x = (x * 1.0e9).round() / 1.0e9;
        let (rpm, collective, speed) = match variable {
            SweepVariable::Rpm => (x, collective0, case.operating.speed),
            SweepVariable::CollectiveDeg => {
                (case.operating.rpm, x.to_radians(), case.operating.speed)
            }
            SweepVariable::AxialSpeed => (case.operating.rpm, collective0, x),
        };
        let r = report_at(&solver, speed, density, radius, collective, rpm)?;
        let label = match variable {
            SweepVariable::AxialSpeed => {
                // J = π λ_c with λ_c = V / (ΩR).
                let tip = rpm * core::f64::consts::TAU / 60.0 * radius;
                let j = core::f64::consts::PI * speed / tip;
                format!("{x},{j:.6}")
            }
            _ => format!("{x}"),
        };
        out.push_str(&format!(
            "{label},{:.6},{:.6},{:.4},{:.8},{:.8},{:.8}\n",
            r.thrust_n, r.torque_nm, r.power_w, r.ct, r.cq, r.cp
        ));
    }
    Ok(out)
}

/// Human-readable summary of a rotor report.
#[must_use]
pub fn rotor_text(report: &RotorReport) -> String {
    let rows: [(&str, String); 8] = [
        ("thrust_n", format!("{:.4}", report.thrust_n)),
        ("torque_nm", format!("{:.5}", report.torque_nm)),
        ("power_w", format!("{:.3}", report.power_w)),
        ("ct", format!("{:.6}", report.ct)),
        ("cq", format!("{:.6}", report.cq)),
        ("cp", format!("{:.6}", report.cp)),
        ("figure_of_merit", format!("{:.4}", report.figure_of_merit)),
        (
            "propeller_efficiency",
            format!("{:.4}", report.propeller_efficiency),
        ),
    ];
    rows.iter()
        .map(|(name, value)| format!("{name:<22}{value}\n"))
        .collect()
}
