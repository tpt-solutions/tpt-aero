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

/// Largest accepted count of spanwise or chordwise panels (per surface, and
/// spanwise stations for the rotor).
pub(crate) const MAX_MESH_COUNT: usize = 500;

/// Largest accepted panel total per surface. The AIC matrix is dense, so its
/// memory and solve cost grow with the square and cube of this count.
pub(crate) const MAX_PANELS: usize = 2500;

/// Validate a spanwise/chordwise panel mesh: both counts ≥ 1, each at most
/// [`MAX_MESH_COUNT`], and their product at most [`MAX_PANELS`].
pub(crate) fn check_mesh(label: &str, n_span: usize, n_chord: usize) -> Result<(), CliError> {
    if n_span == 0 || n_chord == 0 {
        return Err(CliError::Invalid(format!(
            "{label}: n_span and n_chord must be >= 1"
        )));
    }
    if n_span > MAX_MESH_COUNT || n_chord > MAX_MESH_COUNT {
        return Err(CliError::Invalid(format!(
            "{label}: n_span and n_chord must be <= {MAX_MESH_COUNT}"
        )));
    }
    if n_span * n_chord > MAX_PANELS {
        return Err(CliError::Invalid(format!(
            "{label}: n_span × n_chord must be <= {MAX_PANELS} (got {})",
            n_span * n_chord
        )));
    }
    Ok(())
}

pub use aircraft::{
    aircraft_report, run_trim, run_trim_speed_sweep, trim_csv, trim_text, TrimKind, TrimReport,
    AIRCRAFT_TEMPLATE,
};
pub use rotor::{
    rotor_text, run_rotor, run_rotor_sweep, run_rotor_sweep_over, RotorReport, SweepVariable,
    ROTOR_TEMPLATE,
};

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

/// A small fixed-wing UAV (about 2 m span, 15 m/s cruise), ready to run.
pub const UAV_TEMPLATE: &str = r#"# tpt-aero UAV wing case: a 2 m, 15 m/s fixed-wing UAV. Run: tpt-aero wing uav.toml
[wing]
span = 2.0           # m
root_chord = 0.25    # m
tip_chord = 0.18     # m
sweep_le = 0.02      # m, aft offset of the tip leading edge
dihedral_deg = 3.0
camber = { m = 0.02, p = 0.4 }

[mesh]
n_span = 20
n_chord = 4

[conditions]
speed = 15.0         # m/s cruise
altitude = 100.0     # m (ISA)
alpha = { start = -2.0, stop = 10.0, step = 2.0 }   # degrees
cd0 = 0.02           # parasite drag, for the cd and l_over_d columns
"#;

/// A high-aspect-ratio glider (15 m span, 20 m/s), ready to run.
pub const GLIDER_TEMPLATE: &str = r#"# tpt-aero glider wing case: a 15 m high-aspect-ratio glider. Run: tpt-aero wing glider.toml
[wing]
span = 15.0          # m
root_chord = 1.0     # m
tip_chord = 0.5      # m
dihedral_deg = 3.0
camber = { m = 0.015, p = 0.4 }

[mesh]
n_span = 30
n_chord = 4

[conditions]
speed = 20.0         # m/s
altitude = 500.0     # m (ISA)
alpha = { start = -2.0, stop = 8.0, step = 1.0 }
cd0 = 0.008          # clean-wing parasite drag
"#;

/// A small quadrotor propeller (0.12 m radius, hover at 6000 rpm), ready to run.
pub const QUADROTOR_TEMPLATE: &str = r#"# tpt-aero quadrotor propeller case: 0.12 m rotor in hover. Run: tpt-aero rotor quadrotor.toml
[rotor]
radius = 0.12        # m
hub_radius = 0.01    # m
blades = 2
chord = 0.02         # m, root
chord_tip_ratio = 0.6
stations = 16
root_twist_deg = 28.0
tip_twist_deg = 10.0
tip_loss = true
hub_loss = false

[section]
cl_slope = 6.283
cl_max = 1.2
cd0 = 0.012
cd_cl2 = 0.012

[operating]
rpm = 6000.0
speed = 0.0          # hover
altitude = 0.0
collective_deg = 0.0
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
            Self::Invalid(m) => write!(f, "{m}"),
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
    /// Flight Mach number. Given: the Göthert compressible model is used.
    #[serde(default)]
    mach: Option<f64>,
    /// Zero-lift (parasite) drag coefficient added to the induced drag in the
    /// `cd` and `l_over_d` columns.
    #[serde(default)]
    cd0: f64,
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
    /// Parasite drag coefficient added to the induced drag.
    cd0: f64,
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
            // Beyond about 10 % the thin-camber line and the panel geometry
            // stop describing an airfoil; refuse rather than mislead.
            if c.m > 0.1 {
                return Err(CliError::Invalid(
                    "wing.camber m must be at most 0.1 (10 % camber)".into(),
                ));
            }
            CamberModel::Naca4 { m: c.m, p: c.p }
        }
    };
    check_mesh("mesh", case.mesh.n_span, case.mesh.n_chord)?;

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
                .map_err(|e| CliError::Invalid(format!("conditions.altitude: {e}")))?
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
    let mach = match case.conditions.mach {
        None => None,
        Some(m) => {
            if !(m.is_finite() && (0.0..1.0).contains(&m)) {
                return Err(CliError::Invalid(
                    "conditions.mach must satisfy 0 <= mach < 1 (subsonic)".into(),
                ));
            }
            Some(m)
        }
    };
    let model = match mach {
        None => VlmModel::from_planform(&planform, case.mesh.n_span, case.mesh.n_chord, camber),
        Some(m) => VlmModel::from_planform_compressible(
            &planform,
            case.mesh.n_span,
            case.mesh.n_chord,
            camber,
            m,
        ),
    }
    .map_err(|e| CliError::Solve(e.to_string()))?
    // The reference is physical; with_reference stretches it for compressible models.
    .with_reference(ReferenceGeometry::from_planform(&planform));
    let cd0 = finite("conditions.cd0", case.conditions.cd0)?;
    if cd0 < 0.0 {
        return Err(CliError::Invalid("conditions.cd0 must be >= 0".into()));
    }

    Ok(Built {
        cd0,
        model,
        speed,
        density,
        // Round to 1e-9° so that sums like 0.1 + 0.2 print as 0.3, not 0.30000000000000004.
        alphas_deg: (0..n_points)
            .map(|i| ((a0 + da * i as f64) * 1.0e9).round() / 1.0e9)
            .collect(),
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
        cd0,
        model,
        speed,
        density,
        alphas_deg,
    } = build(case_toml)?;
    let mut out = String::from("alpha_deg,cl,cdi,cm,cl_roll,cn,lift_n,cd,l_over_d\n");
    let alphas_rad: Vec<f64> = alphas_deg.iter().map(|a| a.to_radians()).collect();
    // One AIC factorization serves the whole sweep.
    let sols = model
        .solve_sweep(&alphas_rad, speed, density)
        .map_err(|e| CliError::Solve(e.to_string()))?;
    for ((&alpha_deg, &alpha), sol) in alphas_deg.iter().zip(&alphas_rad).zip(&sols) {
        // Lift is the force component normal to the freestream in the x-z plane.
        let lift = -sol.force.x * alpha.sin() + sol.force.z * alpha.cos();
        let cd = sol.cdi + cd0;
        // Near zero drag the ratio is not meaningful; leave it blank.
        let l_over_d = if cd > 1.0e-9 {
            format!("{:.4}", sol.cl / cd)
        } else {
            String::new()
        };
        out.push_str(&format!(
            "{alpha_deg},{:.8},{:.8},{:.8},{:.8},{:.8},{lift:.3},{cd:.8},{l_over_d}\n",
            sol.cl, sol.cdi, sol.cm, sol.cl_roll, sol.cn
        ));
    }
    Ok(out)
}

/// Spanwise loading of a wing case at one angle of attack, as CSV:
/// `strip,y_m,chord_m,cl_section`.
///
/// Each row is one spanwise strip, in the solver's strip order, which runs from
/// the `−y` tip to the `+y` tip across the full span. `y_m` is the strip's mean
/// collocation spanwise station, `chord_m` its mean chord, and `cl_section` its
/// section lift coefficient.
///
/// # Errors
/// As [`run_wing`], and [`CliError::Solve`] if the solve fails.
pub fn run_wing_spanwise(case_toml: &str, alpha_deg: f64) -> Result<String, CliError> {
    let alpha_deg = finite("alpha", alpha_deg)?;
    let Built {
        model,
        speed,
        density,
        ..
    } = build(case_toml)?;
    let sol = model
        .solve(alpha_deg.to_radians(), speed, density, &[])
        .map_err(|e| CliError::Solve(e.to_string()))?;
    let n = sol.strip_cl.len();
    let mut y_sum = vec![0.0; n];
    let mut chord_sum = vec![0.0; n];
    let mut count = vec![0usize; n];
    for e in model.lattice().elements() {
        if e.strip < n {
            y_sum[e.strip] += e.collocation.y;
            chord_sum[e.strip] += e.chord;
            count[e.strip] += 1;
        }
    }
    let mut out = String::from("strip,y_m,chord_m,cl_section\n");
    for s in 0..n {
        let c = count[s].max(1) as f64;
        out.push_str(&format!(
            "{s},{:.6},{:.6},{:.8}\n",
            y_sum[s] / c,
            chord_sum[s] / c,
            sol.strip_cl[s]
        ));
    }
    Ok(out)
}

/// Import an AVL deck and summarize it: the reference geometry, the surfaces,
/// and the lift coefficient at one angle of attack.
///
/// # Errors
/// [`CliError::Parse`] for a deck the importer cannot read, [`CliError::Invalid`]
/// for a flight condition or a model that is rejected.
pub fn avl_import_summary(deck_text: &str, alpha_deg: f64, speed: f64) -> Result<String, CliError> {
    let speed = positive("speed", speed)?;
    let alpha = finite("alpha", alpha_deg)?.to_radians();
    let deck = tpt_aero::vlm::parse_avl(deck_text).map_err(|e| CliError::Parse(e.to_string()))?;
    let model = deck
        .to_model()
        .map_err(|e| CliError::Invalid(e.to_string()))?;
    let sol = model
        .solve(alpha, speed, 1.225, &[])
        .map_err(|e| CliError::Solve(e.to_string()))?;
    let mut out = format!(
        "title       {}
mach        {}
S_ref       {}
c_ref       {}
b_ref       {}
",
        deck.title, deck.mach, deck.sref, deck.cref, deck.bref
    );
    for s in &deck.surfaces {
        out.push_str(&format!(
            "surface     {} ({} sections, {} chordwise x {} spanwise, {} controls)
",
            s.name,
            s.sections.len(),
            s.n_chord,
            s.n_span,
            s.controls.len()
        ));
    }
    out.push_str(&format!(
        "alpha_deg   {}
CL          {:.6}
CDi         {:.6}
Cm          {:.6}
",
        alpha_deg, sol.cl, sol.cdi, sol.cm
    ));
    Ok(out)
}

/// Export a wing case as an AVL deck: the same planform, one `SURFACE` with its
/// strip and chordwise counts, and the reference geometry the solver uses. The
/// camber line is not written; AVL decks use flat sections unless given `AFILE`.
///
/// # Errors
/// As [`run_wing`].
pub fn wing_to_avl(case_toml: &str) -> Result<String, CliError> {
    let case: Case = toml::from_str(case_toml).map_err(|e| CliError::Parse(e.to_string()))?;
    let w = &case.wing;
    let span = positive("wing.span", w.span)?;
    let root = positive("wing.root_chord", w.root_chord)?;
    let tip = positive("wing.tip_chord", w.tip_chord.unwrap_or(root))?;
    let sweep = finite("wing.sweep_le", w.sweep_le)?;
    let dihedral = finite("wing.dihedral_deg", w.dihedral_deg)?.to_radians();
    let planform = WingPlanform::tapered(span, root, tip, sweep, dihedral);
    let reference = ReferenceGeometry::from_planform(&planform);
    let surface = tpt_aero::vlm::Surface::horizontal(planform, case.mesh.n_span, case.mesh.n_chord);
    Ok(tpt_aero::vlm::avl::export(
        "tpt-aero wing export",
        &[surface],
        reference,
    ))
}

/// Warnings about where a wing case leaves the range the model handles well.
/// Empty for a case inside the envelope. Warnings do not stop a run.
///
/// # Errors
/// As [`run_wing`] for a malformed or invalid case.
pub fn wing_warnings(case_toml: &str) -> Result<Vec<String>, CliError> {
    let case: Case = toml::from_str(case_toml).map_err(|e| CliError::Parse(e.to_string()))?;
    let built = build(case_toml)?;
    let mut out = Vec::new();
    let w = &case.wing;
    let tip = w.tip_chord.unwrap_or(w.root_chord);
    // The stored area is stretched for compressible models; use the physical area.
    let area = built.model.reference().area / built.model.compressibility_beta();
    let aspect = w.span * w.span / area;
    if aspect > 30.0 {
        out.push(format!("aspect ratio {aspect:.1} is high: the induced-drag and lift-slope model loses accuracy above about 30"));
    }
    if aspect < 2.0 {
        out.push(format!(
            "aspect ratio {aspect:.2} is low: lifting-line trends are unreliable below about 2"
        ));
    }
    let taper = tip / w.root_chord;
    if !(0.2..=1.0).contains(&taper) {
        out.push(format!("taper ratio {taper:.2} is outside 0.2 to 1.0"));
    }
    let sweep_deg = (w.sweep_le / (0.5 * w.span)).atan().to_degrees();
    if sweep_deg.abs() > 45.0 {
        out.push(format!(
            "leading-edge sweep {sweep_deg:.1} deg is beyond 45 deg: check the planform"
        ));
    }
    let (a0, a1) = (case.conditions.alpha.start, case.conditions.alpha.stop);
    if a0.abs().max(a1.abs()) > 15.0 {
        out.push("the sweep reaches |alpha| above 15 deg, where attached flow and the linear VLM no longer hold".into());
    }
    if case.conditions.mach.is_none() && built.speed / 340.3 > 0.3 {
        out.push("the speed is above Mach 0.3 without conditions.mach: set mach for a compressible result".into());
    }
    Ok(out)
}

/// Mesh-convergence estimate for a wing case at its first sweep angle. Solves
/// the lift coefficient at `n_span`, `2·n_span`, `4·n_span` and `8·n_span`
/// strips and reports the first refinement whose change is below 1 %, with the
/// suggested `n_span`.
///
/// # Errors
/// As [`run_wing`].
pub fn wing_mesh_check(case_toml: &str) -> Result<String, CliError> {
    let case: Case = toml::from_str(case_toml).map_err(|e| CliError::Parse(e.to_string()))?;
    let base_n = case.mesh.n_span;
    if !case_toml.contains(&format!("n_span = {base_n}")) {
        return Err(CliError::Invalid(
            "mesh check needs an explicit mesh.n_span in the case file".into(),
        ));
    }
    let alpha = case.conditions.alpha.start.to_radians();
    let mut rows: Vec<(usize, f64)> = Vec::new();
    for k in 0..4 {
        let n = base_n << k;
        let edited = case_toml.replacen(&format!("n_span = {base_n}"), &format!("n_span = {n}"), 1);
        let built = build(&edited)?;
        let sol = built
            .model
            .solve(alpha, built.speed, built.density, &[])
            .map_err(|e| CliError::Solve(e.to_string()))?;
        rows.push((n, sol.cl));
    }
    let mut out = String::from("n_span,cl,rel_change\n");
    let mut suggested = rows[rows.len() - 1].0;
    let mut found = false;
    for (i, (n, cl)) in rows.iter().enumerate() {
        let change = if i == 0 {
            String::from("")
        } else {
            let prev = rows[i - 1].1;
            let rel = ((cl - prev) / cl.abs().max(1.0e-12)).abs();
            if !found && rel < 0.01 {
                suggested = *n;
                found = true;
            }
            format!("{rel:.5}")
        };
        out.push_str(&format!("{n},{cl:.8},{change}\n"));
    }
    if !found {
        out.push_str("# not converged to 1 % within 8x the requested n_span; refine further\n");
    } else {
        out.push_str(&format!(
            "# suggested n_span = {suggested} (first level changing cl by < 1 %)\n"
        ));
    }
    Ok(out)
}

/// The lift-drag polar of a wing case as a standalone SVG: `C_L` against `C_D`,
/// from the same sweep [`run_wing`] solves.
///
/// # Errors
/// As [`run_wing`].
pub fn wing_polar_svg(case_toml: &str) -> Result<String, CliError> {
    let csv = run_wing(case_toml)?;
    let mut pts: Vec<(f64, f64)> = Vec::new();
    for line in csv.lines().skip(1) {
        let cols: Vec<&str> = line.split(',').collect();
        if cols.len() < 8 {
            continue;
        }
        let cl: f64 = cols[1]
            .parse()
            .map_err(|_| CliError::Solve("bad cl".into()))?;
        let cd: f64 = cols[7]
            .parse()
            .map_err(|_| CliError::Solve("bad cd".into()))?;
        pts.push((cd, cl));
    }
    if pts.is_empty() {
        return Err(CliError::Invalid("the sweep has no points to plot".into()));
    }
    let (w, h, m) = (640.0_f64, 420.0_f64, 56.0_f64);
    let cd_max = pts.iter().map(|p| p.0).fold(0.0_f64, f64::max).max(1.0e-6);
    let cl_min = pts
        .iter()
        .map(|p| p.1)
        .fold(f64::INFINITY, f64::min)
        .min(0.0);
    let cl_max = pts
        .iter()
        .map(|p| p.1)
        .fold(f64::NEG_INFINITY, f64::max)
        .max(1.0e-6);
    let x = |cd: f64| m + (w - 2.0 * m) * cd / cd_max;
    let y = |cl: f64| h - m - (h - 2.0 * m) * (cl - cl_min) / (cl_max - cl_min);
    let path: String = pts
        .iter()
        .enumerate()
        .map(|(i, (cd, cl))| {
            format!(
                "{}{:.1},{:.1}",
                if i == 0 { "M" } else { " L" },
                x(*cd),
                y(*cl)
            )
        })
        .collect();
    Ok(format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {w} {h}\" width=\"{w}\" height=\"{h}\" font-family=\"sans-serif\" font-size=\"12\">\n\
<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n\
<line x1=\"{m}\" y1=\"{b}\" x2=\"{r}\" y2=\"{b}\" stroke=\"#444\"/>\n\
<line x1=\"{m}\" y1=\"{m}\" x2=\"{m}\" y2=\"{b}\" stroke=\"#444\"/>\n\
<text x=\"{mid}\" y=\"{tb}\" text-anchor=\"middle\">C_D</text>\n\
<text x=\"14\" y=\"{mid}\" transform=\"rotate(-90 14 {mid})\" text-anchor=\"middle\">C_L</text>\n\
<path d=\"{path}\" fill=\"none\" stroke=\"#1f6feb\" stroke-width=\"2\"/>\n\
</svg>\n",
        b = h - m,
        r = w - m,
        mid = w / 2.0,
        tb = h - 16.0,
    ))
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
    let long = longitudinal_fd(&model, alpha, speed, density, 1.0e-4)
        .map_err(|e| CliError::Solve(e.to_string()))?;
    let finite_long = [
        long.cl_alpha,
        long.cm_alpha,
        long.cl_q,
        long.cm_q,
        long.cl_delta,
        long.cm_delta,
    ]
    .iter()
    .all(|v| v.is_finite());
    if !finite_long {
        return Err(CliError::Solve(
            "longitudinal derivatives are not finite at this flight condition".into(),
        ));
    }
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
    isa_table(altitudes_m, 0.0, false)
}

/// ISA table with an optional temperature offset and imperial units.
///
/// `delta_isa_k` is the non-standard temperature in kelvin (see
/// [`tpt_aero::isa::at_isa_offset`]). With `imperial`, the altitudes are given in
/// feet and the rows are `altitude_ft,temperature_r,pressure_psf,density_slug_ft3,
/// speed_of_sound_fps`; otherwise they are the SI columns of [`isa_csv`].
///
/// # Errors
/// [`CliError::Invalid`] for an altitude outside the ISA range or an offset
/// beyond the supported ±100 K.
pub fn isa_table(altitudes: &[f64], delta_isa_k: f64, imperial: bool) -> Result<String, CliError> {
    const M_TO_FT: f64 = 3.280_84;
    const PA_TO_PSF: f64 = 0.020_885_43;
    const KG_M3_TO_SLUG_FT3: f64 = 0.001_940_32;
    const M_S_TO_FPS: f64 = 3.280_84;
    const K_TO_R: f64 = 1.8;
    let mut out = if imperial {
        String::from(
            "altitude_ft,temperature_r,pressure_psf,density_slug_ft3,speed_of_sound_fps
",
        )
    } else {
        String::from(
            "altitude_m,temperature_k,pressure_pa,density,speed_of_sound
",
        )
    };
    for &h in altitudes {
        let h_m = if imperial { h / M_TO_FT } else { h };
        let s = isa::at_isa_offset(h_m, delta_isa_k).map_err(|e| {
            CliError::Invalid(format!(
                "altitude {h} {}: {e}",
                if imperial { "ft" } else { "m" }
            ))
        })?;
        if imperial {
            out.push_str(&format!(
                "{h},{},{},{},{}
",
                s.temperature * K_TO_R,
                s.pressure * PA_TO_PSF,
                s.density * KG_M3_TO_SLUG_FT3,
                s.speed_of_sound * M_S_TO_FPS
            ));
        } else {
            out.push_str(&format!(
                "{h},{},{},{},{}
",
                s.temperature, s.pressure, s.density, s.speed_of_sound
            ));
        }
    }
    Ok(out)
}
