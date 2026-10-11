//! Case-file driven wing analysis behind the `tpt-aero` command.
//!
//! A case is a TOML document ([`TEMPLATE`] shows every field). [`run_wing`]
//! sweeps the angle of attack through the VLM and returns the results as CSV,
//! so the binary stays a thin shell and everything here is unit-testable.

use std::fmt;

use serde::Deserialize;
use tpt_aero::prelude::*;
use tpt_aero::stability::lateral::lateral_fd;

mod aircraft;
mod rotor;

pub use aircraft::{run_trim, trim_text, TrimKind, TrimReport, AIRCRAFT_TEMPLATE};
pub use rotor::{rotor_text, run_rotor, RotorReport, ROTOR_TEMPLATE};

/// A fully commented starting point for `tpt-aero wing`.
pub const TEMPLATE: &str = r#"# tpt-aero wing case. Run: tpt-aero wing wing.toml --out polar.csv

[wing]
span = 10.0          # m
root_chord = 1.0     # m
tip_chord = 1.0      # m (default: root_chord)
sweep_le = 0.0       # m, aft offset of the tip leading edge
dihedral_deg = 0.0
# camber = { m = 0.02, p = 0.4 }   # NACA-style camber line (omit for flat plate)

[mesh]
n_span = 24
n_chord = 4

[conditions]
speed = 50.0         # m/s
altitude = 0.0       # m (ISA density); or set `density` in kg/m^3 instead
# density = 1.225
alpha = { start = -4.0, stop = 10.0, step = 2.0 }   # degrees
"#;

/// Everything that can go wrong running a case.
#[derive(Debug)]
pub enum CliError {
    /// The TOML could not be parsed or has missing/unknown fields.
    Parse(String),
    /// A value is physically or numerically unusable.
    Invalid(String),
    /// The solver failed.
    Solve(String),
}

impl fmt::Display for CliError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(m) => write!(f, "case file error: {m}"),
            Self::Invalid(m) => write!(f, "invalid case: {m}"),
            Self::Solve(m) => write!(f, "solver error: {m}"),
        }
    }
}

impl std::error::Error for CliError {}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Case {
    wing: WingSpec,
    #[serde(default)]
    mesh: MeshSpec,
    conditions: Conditions,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct WingSpec {
    span: f64,
    root_chord: f64,
    tip_chord: Option<f64>,
    #[serde(default)]
    sweep_le: f64,
    #[serde(default)]
    dihedral_deg: f64,
    camber: Option<CamberSpec>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct CamberSpec {
    m: f64,
    p: f64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct MeshSpec {
    #[serde(default = "default_n_span")]
    n_span: usize,
    #[serde(default = "default_n_chord")]
    n_chord: usize,
}

impl Default for MeshSpec {
    fn default() -> Self {
        Self {
            n_span: default_n_span(),
            n_chord: default_n_chord(),
        }
    }
}

fn default_n_span() -> usize {
    24
}
fn default_n_chord() -> usize {
    4
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Conditions {
    speed: f64,
    altitude: Option<f64>,
    density: Option<f64>,
    alpha: AlphaSweep,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct AlphaSweep {
    start: f64,
    stop: f64,
    step: f64,
}

fn positive(name: &str, v: f64) -> Result<f64, CliError> {
    if v.is_finite() && v > 0.0 {
        Ok(v)
    } else {
        Err(CliError::Invalid(format!(
            "{name} must be finite and > 0, got {v}"
        )))
    }
}

fn finite(name: &str, v: f64) -> Result<f64, CliError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(CliError::Invalid(format!("{name} must be finite, got {v}")))
    }
}

/// A validated case: the VLM model plus its flight condition and sweep.
struct Built {
    model: VlmModel,
    speed: f64,
    density: f64,
    alphas_deg: Vec<f64>,
}

fn build(case_toml: &str) -> Result<Built, CliError> {
    let case: Case = toml::from_str(case_toml).map_err(|e| CliError::Parse(e.to_string()))?;

    let span = positive("wing.span", case.wing.span)?;
    let root = positive("wing.root_chord", case.wing.root_chord)?;
    let tip = positive("wing.tip_chord", case.wing.tip_chord.unwrap_or(root))?;
    let sweep = finite("wing.sweep_le", case.wing.sweep_le)?;
    let dihedral = finite("wing.dihedral_deg", case.wing.dihedral_deg)?;
    if dihedral.abs() >= 90.0 {
        return Err(CliError::Invalid(
            "wing.dihedral_deg must be within (-90, 90)".into(),
        ));
    }
    let camber = match &case.wing.camber {
        None => CamberModel::None,
        Some(c) => {
            if !(c.m.is_finite() && c.m >= 0.0 && c.p.is_finite() && c.p > 0.0 && c.p < 1.0) {
                return Err(CliError::Invalid(
                    "wing.camber needs m >= 0 and 0 < p < 1".into(),
                ));
            }
            CamberModel::Naca4 { m: c.m, p: c.p }
        }
    };
    if case.mesh.n_span == 0 || case.mesh.n_chord == 0 {
        return Err(CliError::Invalid("mesh counts must be >= 1".into()));
    }

    let speed = positive("conditions.speed", case.conditions.speed)?;
    let density = match (case.conditions.density, case.conditions.altitude) {
        (Some(_), Some(_)) => {
            return Err(CliError::Invalid(
                "set either conditions.density or conditions.altitude, not both".into(),
            ))
        }
        (Some(rho), None) => positive("conditions.density", rho)?,
        (None, alt) => {
            isa::at(alt.unwrap_or(0.0))
                .map_err(|e| CliError::Invalid(format!("conditions.altitude: {e:?}")))?
                .density
        }
    };

    let sweep_spec = &case.conditions.alpha;
    let (a0, a1, da) = (
        finite("alpha.start", sweep_spec.start)?,
        finite("alpha.stop", sweep_spec.stop)?,
        finite("alpha.step", sweep_spec.step)?,
    );
    if da == 0.0 || (a1 - a0) * da < 0.0 {
        return Err(CliError::Invalid(
            "alpha.step must be non-zero and point from start towards stop".into(),
        ));
    }
    let n_points = ((a1 - a0) / da + 1e-9).floor() as usize + 1;
    if n_points > 10_000 {
        return Err(CliError::Invalid(
            "alpha sweep has more than 10000 points".into(),
        ));
    }

    let planform = WingPlanform::tapered(span, root, tip, sweep, dihedral.to_radians());
    let model = VlmModel::from_planform(&planform, case.mesh.n_span, case.mesh.n_chord, camber)
        .map_err(|e| CliError::Solve(e.to_string()))?
        .with_reference(ReferenceGeometry::from_planform(&planform));

    Ok(Built {
        model,
        speed,
        density,
        alphas_deg: (0..n_points).map(|i| a0 + da * i as f64).collect(),
    })
}

/// Run a wing case (TOML text) and return the polar as CSV.
///
/// Columns: `alpha_deg,cl,cdi,cm,cl_roll,cn,lift_n`.
///
/// # Errors
/// [`CliError`] for malformed or unphysical input and solver failures.
pub fn run_wing(case_toml: &str) -> Result<String, CliError> {
    let Built {
        model,
        speed,
        density,
        alphas_deg,
    } = build(case_toml)?;
    let mut out = String::from("alpha_deg,cl,cdi,cm,cl_roll,cn,lift_n\n");
    let alphas_rad: Vec<f64> = alphas_deg.iter().map(|a| a.to_radians()).collect();
    // One AIC factorization serves the whole sweep.
    let sols = model
        .solve_sweep(&alphas_rad, speed, density)
        .map_err(|e| CliError::Solve(e.to_string()))?;
    for ((&alpha_deg, &alpha), sol) in alphas_deg.iter().zip(&alphas_rad).zip(&sols) {
        // Lift is the force component normal to the freestream in the x-z plane.
        let lift = -sol.force.x * alpha.sin() + sol.force.z * alpha.cos();
        out.push_str(&format!(
            "{alpha_deg},{:.8},{:.8},{:.8},{:.8},{:.8},{lift:.3}\n",
            sol.cl, sol.cdi, sol.cm, sol.cl_roll, sol.cn
        ));
    }
    Ok(out)
}

/// Stability derivatives of the case's wing at one trim angle of attack, as
/// `derivative,value` CSV rows (per rad, or per non-dimensional rate
/// `q c/2V`, `p b/2V`, `r b/2V`; aircraft-standard signs).
///
/// Wing-alone results: no tail or fuselage contribution.
///
/// # Errors
/// [`CliError`] for a malformed case or a failed solve.
pub fn run_derivatives(case_toml: &str, alpha_deg: f64) -> Result<String, CliError> {
    let (_, _, long, lat) = derivative_sets(case_toml, alpha_deg)?;
    let rows = [
        ("cl_alpha", long.cl_alpha),
        ("cm_alpha", long.cm_alpha),
        ("cl_q", long.cl_q),
        ("cm_q", long.cm_q),
        ("cy_beta", lat.cy_beta),
        ("cl_beta", lat.cl_beta),
        ("cn_beta", lat.cn_beta),
        ("cy_p", lat.cy_p),
        ("cl_p", lat.cl_p),
        ("cn_p", lat.cn_p),
        ("cy_r", lat.cy_r),
        ("cl_r", lat.cl_r),
        ("cn_r", lat.cn_r),
    ];
    let mut out = String::from("derivative,value\n");
    for (name, value) in rows {
        out.push_str(&format!("{name},{value:.6}\n"));
    }
    Ok(out)
}

fn derivative_sets(
    case_toml: &str,
    alpha_deg: f64,
) -> Result<(f64, f64, LongitudinalDerivatives, LateralDerivatives), CliError> {
    let Built {
        model,
        speed,
        density,
        ..
    } = build(case_toml)?;
    let alpha = finite("alpha", alpha_deg)?.to_radians();
    let long = longitudinal_fd(&model, alpha, speed, density, 1.0e-4);
    let lat =
        lateral_fd(&model, alpha, speed, density).map_err(|e| CliError::Solve(e.to_string()))?;
    Ok((speed, density, long, lat))
}

/// As [`run_derivatives`], as pretty-printed JSON with the flight condition
/// and the `longitudinal` and `lateral` derivative sets (full precision).
///
/// # Errors
/// [`CliError`] for a malformed case or a failed solve.
pub fn run_derivatives_json(case_toml: &str, alpha_deg: f64) -> Result<String, CliError> {
    let (speed, density, long, lat) = derivative_sets(case_toml, alpha_deg)?;
    let value = serde_json::json!({
        "alpha_deg": alpha_deg,
        "speed": speed,
        "density": density,
        "longitudinal": long,
        "lateral": lat,
    });
    serde_json::to_string_pretty(&value).map_err(|e| CliError::Solve(e.to_string()))
}

/// Solve one angle of attack of the case's wing and render the lattice with
/// its circulation as legacy ASCII VTK (open in ParaView).
///
/// # Errors
/// [`CliError`] for a malformed case or a failed solve.
pub fn run_vtk(case_toml: &str, alpha_deg: f64) -> Result<String, CliError> {
    let Built {
        model,
        speed,
        density,
        ..
    } = build(case_toml)?;
    let alpha = finite("alpha", alpha_deg)?.to_radians();
    let sol = model
        .solve(alpha, speed, density, &[])
        .map_err(|e| CliError::Solve(e.to_string()))?;
    tpt_aero::vlm::export::to_vtk(
        &model,
        &format!("tpt-aero VLM alpha={alpha_deg} deg V={speed} m/s"),
        Some(&sol.gamma),
    )
    .map_err(|e| CliError::Solve(e.to_string()))
}

/// ISA table rows (`altitude_m,temperature_k,pressure_pa,density,speed_of_sound`).
///
/// # Errors
/// [`CliError::Invalid`] for an altitude outside the ISA range.
pub fn isa_csv(altitudes_m: &[f64]) -> Result<String, CliError> {
    let mut out = String::from("altitude_m,temperature_k,pressure_pa,density,speed_of_sound\n");
    for &h in altitudes_m {
        let s = isa::at(h).map_err(|e| CliError::Invalid(format!("altitude {h}: {e:?}")))?;
        out.push_str(&format!(
            "{h},{},{},{},{}\n",
            s.temperature, s.pressure, s.density, s.speed_of_sound
        ));
    }
    Ok(out)
}
