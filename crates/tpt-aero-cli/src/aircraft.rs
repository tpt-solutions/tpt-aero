//! Aircraft case files: lifting surfaces with control surfaces, for steady
//! level and coordinated-turn trim.

use serde::{Deserialize, Serialize};
use tpt_aero::aero_core::{Vec3, WingPlanform};
use tpt_aero::stability::longitudinal::longitudinal_fd;
use tpt_aero::stability::trim::{trim_level_flight, TrimError, TrimRequest};
use tpt_aero::stability::turn::{coordinated_turn_rate, trim_steady_turn, TurnRequest, G};
use tpt_aero::vlm::{ControlSymmetry, ReferenceGeometry, Surface, VlmModel};

use crate::{check_mesh, finite, positive, CliError};

/// A fully commented starting point for `tpt-aero trim`.
pub const AIRCRAFT_TEMPLATE: &str = r#"# tpt-aero aircraft case. Run: tpt-aero trim aircraft.toml --speed 50 [--bank 20]

[mass]
mass = 600.0          # kg
cd0 = 0.03            # parasite drag coefficient (sets the thrust needed)
cg_x = 0.25           # CG, streamwise from the wing root leading edge (m); sets the static margin
# cl_max = 1.5        # wing CL_max, if known; adds the CL_max margin to the trim report

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
outboard = 0.75      # |y| up to span/2: the rudder covers the whole fin
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
    /// Centre of gravity, streamwise position from the body origin (m). The
    /// moment reference for the pitch derivatives and the static margin.
    #[serde(default)]
    cg_x: f64,
    /// Maximum lift coefficient of the wing, if known. Enables the CL_max margin.
    #[serde(default)]
    cl_max: Option<f64>,
    /// Measured roll, pitch and yaw inertias (kg·m²). Without them the report
    /// uses labelled estimates for the mode frequencies.
    #[serde(default)]
    inertia_x: Option<f64>,
    #[serde(default)]
    inertia_y: Option<f64>,
    #[serde(default)]
    inertia_z: Option<f64>,
}

fn default_cd0() -> f64 {
    0.03
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct SurfaceSpec {
    #[serde(default = "default_kind")]
    kind: SurfaceKind,
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
    symmetry: ControlSymmetryKind,
}

/// How a lifting surface sits on the airframe. Unknown values are rejected
/// when the case is parsed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum SurfaceKind {
    /// A wing, tailplane or canard: spans the body y axis.
    Horizontal,
    /// A fin: rolled 90° to point up along body z.
    Vertical,
}

/// Whether a control deflects both halves the same way or oppositely.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum ControlSymmetryKind {
    /// Both halves deflect the same way (elevator, rudder).
    Symmetric,
    /// Opposite halves (aileron).
    Antisymmetric,
}

fn default_kind() -> SurfaceKind {
    SurfaceKind::Horizontal
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
    /// Lift coefficient `C_L = n·W/(q S)`.
    pub cl: f64,
    /// Drag coefficient `C_D = T/(q S)` (level flight, where thrust equals drag).
    pub cd: Option<f64>,
    /// Lift-to-drag ratio (level flight).
    pub lift_over_drag: Option<f64>,
    /// Static margin `−C_mα/C_Lα` about the CG, as a fraction of the reference
    /// chord (positive is statically stable). Level flight only.
    pub static_margin: Option<f64>,
    /// `C_L,max / C_L − 1` when the case gives `mass.cl_max`. Level flight only.
    pub cl_max_margin: Option<f64>,
    /// Load factor (lift over weight).
    pub load_factor: f64,
    /// Turn rate (deg/s); zero for level flight.
    pub turn_rate_deg_s: f64,
    /// Turn radius (m); `None` for level flight.
    pub turn_radius_m: Option<f64>,
}

fn build_model(case: &AircraftCase) -> Result<VlmModel, CliError> {
    let mut surfaces = Vec::with_capacity(case.surfaces.len());
    // Control names are how trim finds the elevator, aileron and rudder, so
    // they must be unique across the whole aircraft.
    let mut control_names: Vec<&str> = Vec::new();
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
        check_mesh(&label, sp.n_span, sp.n_chord)?;
        let planform = WingPlanform::tapered(span, root, tip, sweep, dihedral);
        let mut surface = match sp.kind {
            SurfaceKind::Horizontal => Surface::horizontal(planform, sp.n_span, sp.n_chord),
            SurfaceKind::Vertical => Surface::vertical(planform, sp.n_span, sp.n_chord),
        }
        .at(
            finite(&format!("{label}.x"), sp.x)?,
            finite(&format!("{label}.y"), sp.y)?,
            finite(&format!("{label}.z"), sp.z)?,
        );
        for c in &sp.controls {
            if c.name.trim().is_empty() {
                return Err(CliError::Invalid(format!(
                    "{label}: a control needs a name"
                )));
            }
            if control_names.contains(&c.name.as_str()) {
                return Err(CliError::Invalid(format!(
                    "control name \"{}\" is used more than once",
                    c.name
                )));
            }
            control_names.push(&c.name);
            let symmetry = match c.symmetry {
                ControlSymmetryKind::Symmetric => ControlSymmetry::Symmetric,
                ControlSymmetryKind::Antisymmetric => ControlSymmetry::Antisymmetric,
            };
            let hinge = finite(&format!("control {}.hinge", c.name), c.hinge)?;
            if !(0.0..1.0).contains(&hinge) {
                return Err(CliError::Invalid(format!(
                    "control {}: hinge must be in [0, 1)",
                    c.name
                )));
            }
            let inboard = finite(&format!("control {}.inboard", c.name), c.inboard)?;
            let outboard = finite(&format!("control {}.outboard", c.name), c.outboard)?;
            // inboard/outboard are |y| in metres; the planform runs ±span/2.
            let half_span = 0.5 * span;
            if !(0.0 <= inboard && inboard < outboard && outboard <= half_span) {
                return Err(CliError::Invalid(format!(
                    "control {}: need 0 <= inboard < outboard <= span/2 ({half_span} m), got {inboard} to {outboard}",
                    c.name
                )));
            }
            // A control must cover at least one spanwise strip, or it does nothing.
            let strip_width = span / sp.n_span as f64;
            let covered = (0..sp.n_span)
                .map(|k| ((k as f64 + 0.5) * strip_width - half_span).abs())
                .filter(|y| (inboard..=outboard).contains(y))
                .count();
            if covered == 0 {
                return Err(CliError::Invalid(format!(
                    "control {} covers no spanwise panel; widen it or raise n_span",
                    c.name
                )));
            }
            surface = surface.with_control(&c.name, hinge, inboard, outboard, symmetry);
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
    // The CG is the moment reference: it sets C_m and so the static margin.
    let r = *model.reference();
    let model = model.with_reference(ReferenceGeometry {
        area: r.area,
        chord: r.chord,
        span: r.span,
        moment_center: Vec3::new(case.mass.cg_x, 0.0, 0.0),
    });
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
            let q_s = 0.5 * density * speed * speed * model.reference().area;
            let cl = weight / q_s;
            let cd = t.thrust_n / q_s;
            let d = longitudinal_fd(&model, t.alpha_rad, speed, density, 1.0e-4)
                .map_err(|e| CliError::Solve(e.to_string()))?;
            let static_margin = -d.cm_alpha / d.cl_alpha;
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
                cl,
                cd: Some(cd),
                lift_over_drag: Some(cl / cd),
                static_margin: Some(static_margin),
                cl_max_margin: case.mass.cl_max.map(|m| m / cl - 1.0),
            })
        }
        TrimKind::Turn { bank_rad } => {
            let bank = finite("bank", bank_rad)?;
            // tan(bank) diverges at 90°, so a coordinated turn needs |bank| < 90°.
            if bank.abs() >= core::f64::consts::FRAC_PI_2 {
                return Err(CliError::Invalid(format!(
                    "--bank must satisfy |bank| < 90 deg (got {:.2} deg)",
                    bank.to_degrees()
                )));
            }
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
                cl: t.load_factor * weight
                    / (0.5 * density * speed * speed * model.reference().area),
                cd: None,
                lift_over_drag: None,
                static_margin: None,
                cl_max_margin: None,
            })
        }
    }
}

fn trim_failed(e: TrimError) -> CliError {
    match e {
        TrimError::NoSuchControl => CliError::Invalid(
            "the aircraft case lacks a control the trim needs (level flight: elevator; turns: elevator, aileron and rudder)"
                .into(),
        ),
        other => CliError::Solve(other.to_string()),
    }
}

/// The trim report as `(quantity, value, unit)` rows, in report order. Shared
/// by the text and CSV renderers so the two cannot drift apart.
fn trim_rows(r: &TrimReport) -> Vec<(&'static str, String, &'static str)> {
    let mut rows: Vec<(&'static str, String, &'static str)> = vec![
        ("kind", r.kind.clone(), "-"),
        ("speed_ms", format!("{:.3}", r.speed_ms), "m/s"),
        ("density", format!("{:.5}", r.density), "kg/m^3"),
        ("bank_deg", format!("{:.4}", r.bank_deg), "deg"),
        ("alpha_deg", format!("{:.4}", r.alpha_deg), "deg"),
        ("beta_deg", format!("{:.4}", r.beta_deg), "deg"),
        ("elevator_deg", format!("{:.4}", r.elevator_deg), "deg"),
    ];
    if let Some(a) = r.aileron_deg {
        rows.push(("aileron_deg", format!("{a:.4}"), "deg"));
    }
    if let Some(d) = r.rudder_deg {
        rows.push(("rudder_deg", format!("{d:.4}"), "deg"));
    }
    if let Some(t) = r.thrust_n {
        rows.push(("thrust_n", format!("{t:.3}"), "N"));
    }
    rows.push(("load_factor", format!("{:.5}", r.load_factor), "g"));
    rows.push(("cl", format!("{:.5}", r.cl), "-"));
    if let Some(cd) = r.cd {
        rows.push(("cd", format!("{cd:.5}"), "-"));
    }
    if let Some(ld) = r.lift_over_drag {
        rows.push(("lift_over_drag", format!("{ld:.3}"), "-"));
    }
    if let Some(sm) = r.static_margin {
        // Fraction of the reference chord; positive is stable.
        rows.push((
            "static_margin_mac",
            format!("{sm:.4}"),
            "fraction of MAC (+ = stable)",
        ));
    }
    if let Some(m) = r.cl_max_margin {
        rows.push(("cl_max_margin", format!("{m:.4}"), "C_Lmax/C_L - 1"));
    }
    rows.push((
        "turn_rate_deg_s",
        format!("{:.4}", r.turn_rate_deg_s),
        "deg/s",
    ));
    if let Some(rad) = r.turn_radius_m {
        rows.push(("turn_radius_m", format!("{rad:.2}"), "m"));
    }
    rows
}

/// Human-readable summary of a trim report: quantity, value and units columns.
#[must_use]
pub fn trim_text(r: &TrimReport) -> String {
    let rows = trim_rows(r);
    let name_w = rows.iter().map(|row| row.0.len()).max().unwrap_or(0);
    let value_w = rows.iter().map(|row| row.1.len()).max().unwrap_or(0);
    let mut out = String::new();
    for (name, value, unit) in rows {
        out.push_str(&format!("{name:<name_w$}  {value:>value_w$}  {unit}\n"));
    }
    out
}

/// The trim report as CSV with a `quantity,value,unit` header. Values are
/// unquoted; the one text value (`kind`) has no commas.
#[must_use]
pub fn trim_csv(r: &TrimReport) -> String {
    let mut out = String::from("quantity,value,unit\n");
    for (name, value, unit) in trim_rows(r) {
        out.push_str(&format!("{name},{value},{unit}\n"));
    }
    out
}

/// A Markdown stability report for an aircraft case: the level trim at `speed`
/// and `altitude`, the static margin, the longitudinal and lateral derivatives,
/// and the short-period and Dutch-roll modes. Inertias come from the case when
/// given, and otherwise from labelled estimates.
///
/// # Errors
/// As [`run_trim`], and [`CliError::Solve`] if a derivative solve fails.
pub fn aircraft_report(case_toml: &str, speed: f64, altitude: f64) -> Result<String, CliError> {
    use tpt_aero::stability::lateral::{lateral_fd, LateralParameters, LateralPlant};
    use tpt_aero::stability::plant::{LongitudinalPlant, PlantParameters};

    let trim = run_trim(case_toml, speed, altitude, TrimKind::Level)?;
    let case: AircraftCase =
        toml::from_str(case_toml).map_err(|e| CliError::Parse(e.to_string()))?;
    let model = build_model(&case)?;
    let r = *model.reference();
    let model = model.with_reference(ReferenceGeometry {
        area: r.area,
        chord: r.chord,
        span: r.span,
        moment_center: Vec3::new(case.mass.cg_x, 0.0, 0.0),
    });
    let mass = case.mass.mass;
    let alpha = trim.alpha_deg.to_radians();
    let rho = trim.density;
    let v = trim.speed_ms;
    let d = longitudinal_fd(&model, alpha, v, rho, 1.0e-4)
        .map_err(|e| CliError::Solve(e.to_string()))?;
    let lat = lateral_fd(&model, alpha, v, rho).map_err(|e| CliError::Solve(e.to_string()))?;
    let span = r.span;
    // (is_estimate, value): a measured inertia wins over the span estimate.
    let pick = |known: Option<f64>, fallback: f64| (known.is_none(), known.unwrap_or(fallback));
    let (est_y, iy) = pick(case.mass.inertia_y, mass * (span / 3.0).powi(2));
    let (est_x, ix) = pick(case.mass.inertia_x, mass * (span / 4.0).powi(2));
    let (est_z, iz) = pick(case.mass.inertia_z, mass * (span / 2.5).powi(2));
    let any_estimated = est_x || est_y || est_z;

    let mut out = String::new();
    out.push_str(
        "# Stability report

",
    );
    out.push_str(&format!(
        "Level flight at {v:.1} m/s, altitude {altitude:.0} m (ρ = {rho:.4} kg/m³), mass {mass:.0} kg, CG at x = {:.2} m.

",
        case.mass.cg_x
    ));
    out.push_str(
        "## Trim

| quantity | value |
|---|---|
",
    );
    out.push_str(&format!(
        "| angle of attack | {:.2} deg |
",
        trim.alpha_deg
    ));
    out.push_str(&format!(
        "| elevator | {:.2} deg |
",
        trim.elevator_deg
    ));
    out.push_str(&format!(
        "| C_L | {:.4} |
",
        trim.cl
    ));
    if let (Some(cd), Some(ld)) = (trim.cd, trim.lift_over_drag) {
        out.push_str(&format!(
            "| C_D | {cd:.5} |
| L/D | {ld:.2} |
"
        ));
    }
    if let Some(m) = trim.cl_max_margin {
        out.push_str(&format!(
            "| CL_max margin | {m:.3} |
"
        ));
    }
    out.push_str(
        "
## Static stability

",
    );
    let sm = -d.cm_alpha / d.cl_alpha;
    let verdict = if sm > 0.0 {
        "statically stable"
    } else {
        "statically unstable"
    };
    out.push_str(&format!(
        "Static margin {sm:.3} of the reference chord ({verdict}). Neutral point is {:.2} m aft of the CG.

",
        sm * r.chord
    ));
    out.push_str(
        "## Derivatives (per radian unless noted)

| derivative | value |
|---|---|
",
    );
    for (name, value) in [
        ("C_L,alpha", d.cl_alpha),
        ("C_m,alpha", d.cm_alpha),
        ("C_L,q", d.cl_q),
        ("C_m,q", d.cm_q),
        ("C_L,delta (uniform)", d.cl_delta),
        ("C_m,delta (uniform)", d.cm_delta),
        ("C_Y,beta", lat.cy_beta),
        ("C_l,beta", lat.cl_beta),
        ("C_n,beta", lat.cn_beta),
    ] {
        out.push_str(&format!(
            "| {name} | {value:.4} |
"
        ));
    }
    out.push_str(
        "
## Modes

",
    );
    let plant = LongitudinalPlant::from_derivatives(
        &d,
        &PlantParameters {
            v_inf: v,
            chord: r.chord,
            area: r.area,
            mass,
            inertia_y: iy,
            density: rho,
        },
        0.0,
    );
    match plant.short_period() {
        Some((wn, zeta)) => out.push_str(&format!(
            "- Short period: ω_n {wn:.3} rad/s, ζ {zeta:.3}
"
        )),
        None => out.push_str(
            "- Short period: not oscillatory (statically unstable)
",
        ),
    }
    let lateral = LateralPlant::from_derivatives(
        &lat,
        &LateralParameters {
            v_inf: v,
            span,
            area: r.area,
            density: rho,
            mass,
            inertia_x: ix,
            inertia_z: iz,
            inertia_xz: 0.0,
        },
    );
    match lateral.and_then(|p| p.modes()) {
        Some(m) => {
            out.push_str(&format!(
                "- Dutch roll: ω_n {:.3} rad/s, ζ {:.3}
- Roll subsidence: {:.3} 1/s
- Spiral: {:.4} 1/s ({})
",
                m.dutch_roll_wn,
                m.dutch_roll_zeta,
                m.roll_root,
                m.spiral_root,
                if m.spiral_root < 0.0 {
                    "stable"
                } else {
                    "divergent"
                }
            ));
        }
        None => out.push_str(
            "- Lateral modes: not available for these inputs
",
        ),
    }
    if any_estimated {
        out.push_str("
> Some inertias are estimates from the span (I_y ≈ m(b/3)², I_x ≈ m(b/4)², I_z ≈ m(b/2.5)²). Give `mass.inertia_x/y/z` in the case for measured values; the mode frequencies depend on them.
");
    }
    Ok(out)
}

/// Level trim swept over airspeed `start..=stop` m/s in steps of `step`, as CSV.
///
/// Columns: `speed_ms,alpha_deg,elevator_deg,cl,lift_over_drag,static_margin_mac,cl_max_margin`.
/// Each row is one [`run_trim`] at that speed, so a trim that fails at any speed
/// fails the sweep with that speed in the message.
///
/// # Errors
/// [`CliError::Invalid`] for a non-finite bound, a non-positive step, a reversed
/// range, or more than 1000 points; otherwise as [`run_trim`].
pub fn run_trim_speed_sweep(
    case_toml: &str,
    altitude: f64,
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
    let n = ((stop - start) / step + 1.0e-9).floor() as usize + 1;
    if n > 1000 {
        return Err(CliError::Invalid(
            "speed sweep has more than 1000 points".into(),
        ));
    }
    let mut out = String::from(
        "speed_ms,alpha_deg,elevator_deg,cl,lift_over_drag,static_margin_mac,cl_max_margin\n",
    );
    for i in 0..n {
        let v = ((start + step * i as f64) * 1.0e9).round() / 1.0e9;
        let r = run_trim(case_toml, v, altitude, TrimKind::Level)
            .map_err(|e| CliError::Solve(format!("trim failed at {v} m/s: {e}")))?;
        let opt = |x: Option<f64>| x.map_or(String::new(), |x| format!("{x:.6}"));
        out.push_str(&format!(
            "{v},{:.6},{:.6},{:.6},{},{},{}\n",
            r.alpha_deg,
            r.elevator_deg,
            r.cl,
            opt(r.lift_over_drag),
            opt(r.static_margin),
            opt(r.cl_max_margin),
        ));
    }
    Ok(out)
}
