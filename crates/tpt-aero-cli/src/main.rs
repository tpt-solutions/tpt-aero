//! `tpt-aero` command-line entry point.

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Parser, Subcommand, ValueEnum};
use tpt_aero_cli::{
    aircraft_report, avl_import_summary, isa_table, rotor_text, run_derivatives,
    run_derivatives_json, run_rotor, run_rotor_sweep_over, run_trim, run_trim_speed_sweep, run_vtk,
    run_wing, run_wing_spanwise, trim_csv, trim_text, wing_mesh_check, wing_mesh_sweep,
    wing_polar_svg, wing_to_avl, wing_warnings, SweepVariable, TrimKind, AIRCRAFT_TEMPLATE,
    GLIDER_TEMPLATE, QUADROTOR_TEMPLATE, ROTOR_TEMPLATE, TEMPLATE, UAV_TEMPLATE,
};

/// Aerodynamics from the command line: VLM wings, BEMT rotors, trim and
/// stability derivatives, from TOML case files.
#[derive(Debug, Parser)]
#[command(name = "tpt-aero", version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Run an alpha sweep of a wing; CSV to stdout unless `--out` is given.
    Wing {
        /// Wing case file (TOML).
        case: PathBuf,
        /// Write the CSV to this file instead of stdout.
        #[arg(short, long, value_name = "FILE")]
        out: Option<PathBuf>,
        /// Also write the lift-drag polar as an SVG file.
        #[arg(long, value_name = "FILE")]
        svg: Option<PathBuf>,
        /// Print a mesh-convergence estimate instead of the sweep.
        #[arg(long, conflicts_with_all = ["svg", "spanwise", "mesh_sweep"])]
        check_mesh: bool,
        /// Print C_L and C_Di at each of these spanwise strip counts (comma-separated)
        /// instead of the sweep, e.g. `--mesh-sweep 8,16,32`.
        #[arg(
            long,
            value_name = "N,N,...",
            value_delimiter = ',',
            conflicts_with_all = ["svg", "spanwise", "check_mesh"]
        )]
        mesh_sweep: Option<Vec<usize>>,
        /// Print the spanwise lift distribution at this angle (deg) instead of the sweep.
        #[arg(
            long,
            value_name = "DEG",
            allow_negative_numbers = true,
            conflicts_with_all = ["svg", "check_mesh"]
        )]
        spanwise: Option<f64>,
    },
    /// Stability derivatives at one angle of attack.
    Derivatives {
        /// Wing case file (TOML).
        case: PathBuf,
        /// Angle of attack in degrees.
        #[arg(long, default_value_t = 2.0, allow_negative_numbers = true)]
        alpha: f64,
        /// Print the derivatives as JSON (full precision).
        #[arg(long)]
        json: bool,
        /// Write the report to this file instead of stdout.
        #[arg(short, long, value_name = "FILE")]
        out: Option<PathBuf>,
    },
    /// Lattice and circulation as legacy ASCII VTK (open in ParaView).
    Vtk {
        /// Wing case file (TOML).
        case: PathBuf,
        /// Angle of attack in degrees.
        #[arg(long, default_value_t = 2.0, allow_negative_numbers = true)]
        alpha: f64,
        /// Write the VTK file here instead of stdout.
        #[arg(short, long, value_name = "FILE")]
        out: Option<PathBuf>,
    },
    /// One BEMT operating point: thrust, torque and power.
    Rotor {
        /// Rotor case file (TOML).
        case: PathBuf,
        /// Print the result as JSON.
        #[arg(long)]
        json: bool,
        /// Sweep the shaft speed instead: START:STOP:STEP in rpm, printed as CSV.
        #[arg(long, value_name = "START:STOP:STEP", conflicts_with_all = ["json", "collective_sweep", "speed_sweep"])]
        rpm_sweep: Option<String>,
        /// Sweep the collective pitch: START:STOP:STEP in degrees, printed as CSV.
        #[arg(long, value_name = "START:STOP:STEP", conflicts_with_all = ["json", "speed_sweep"])]
        collective_sweep: Option<String>,
        /// Sweep the axial speed (advance ratio J): START:STOP:STEP in m/s, CSV.
        #[arg(long, value_name = "START:STOP:STEP", conflicts_with = "json")]
        speed_sweep: Option<String>,
        /// Write the report to this file instead of stdout.
        #[arg(short, long, value_name = "FILE")]
        out: Option<PathBuf>,
    },
    /// Level trim, or a coordinated turn with `--bank`.
    Trim {
        /// Aircraft case file (TOML).
        aircraft: PathBuf,
        /// True airspeed in m/s (required unless `--speed-sweep` is given).
        #[arg(
            long,
            value_name = "M_S",
            required_unless_present = "speed_sweep",
            conflicts_with = "speed_sweep"
        )]
        speed: Option<f64>,
        /// Level trim over airspeed: START:STOP:STEP in m/s, printed as CSV.
        #[arg(long, value_name = "START:STOP:STEP", conflicts_with_all = ["bank", "json", "csv"])]
        speed_sweep: Option<String>,
        /// Geopotential altitude in metres (ISA).
        #[arg(long, default_value_t = 0.0, allow_negative_numbers = true)]
        altitude: f64,
        /// Bank angle in degrees for a coordinated turn (must be below 90).
        #[arg(long, value_name = "DEG", allow_negative_numbers = true)]
        bank: Option<f64>,
        /// Print the report as JSON.
        #[arg(long, conflicts_with = "csv")]
        json: bool,
        /// Print the report as CSV (quantity, value, unit).
        #[arg(long)]
        csv: bool,
        /// Write the report to this file instead of stdout.
        #[arg(short, long, value_name = "FILE")]
        out: Option<PathBuf>,
    },
    /// Markdown stability report for an aircraft case: trim, static margin,
    /// derivatives and modes.
    Report {
        /// Aircraft case file (TOML).
        aircraft: PathBuf,
        /// True airspeed in m/s.
        #[arg(long, value_name = "M_S")]
        speed: f64,
        /// Geopotential altitude in metres (ISA).
        #[arg(long, default_value_t = 0.0, allow_negative_numbers = true)]
        altitude: f64,
        /// Write the report to this file instead of stdout.
        #[arg(short, long, value_name = "FILE")]
        out: Option<PathBuf>,
    },
    /// Import or export AVL `.avl` decks.
    Avl {
        #[command(subcommand)]
        action: AvlCommand,
    },
    /// Print a commented case file to start from.
    Template {
        /// Which case file to print.
        #[arg(value_enum)]
        kind: TemplateKind,
        /// Write the template to this file instead of stdout.
        #[arg(short, long, value_name = "FILE")]
        out: Option<PathBuf>,
    },
    /// ISA properties at one or more altitudes, as CSV.
    Isa {
        /// Geopotential altitudes in metres (0 to 84852), or feet with `--imperial`.
        #[arg(
            required = true,
            allow_negative_numbers = true,
            value_name = "ALTITUDE"
        )]
        altitudes: Vec<f64>,
        /// Non-standard temperature offset in kelvin (within ±100).
        #[arg(
            long,
            default_value_t = 0.0,
            allow_negative_numbers = true,
            value_name = "K"
        )]
        delta_isa: f64,
        /// Imperial units: altitudes in feet; output in °R, psf, slug/ft³ and ft/s.
        #[arg(long)]
        imperial: bool,
        /// Write the CSV to this file instead of stdout.
        #[arg(short, long, value_name = "FILE")]
        out: Option<PathBuf>,
    },
}

#[derive(Debug, Subcommand)]
enum AvlCommand {
    /// Read an AVL deck, solve it at one angle of attack and summarize it.
    Import {
        /// AVL deck (`.avl`).
        deck: PathBuf,
        /// Angle of attack in degrees.
        #[arg(long, default_value_t = 2.0, allow_negative_numbers = true)]
        alpha: f64,
        /// Airspeed in m/s (sea-level density).
        #[arg(long, default_value_t = 50.0)]
        speed: f64,
        /// Write the summary to this file instead of stdout.
        #[arg(short, long, value_name = "FILE")]
        out: Option<PathBuf>,
    },
    /// Write a wing case as an AVL deck.
    Export {
        /// Wing case file (TOML).
        case: PathBuf,
        /// Write the deck to this file instead of stdout.
        #[arg(short, long, value_name = "FILE")]
        out: Option<PathBuf>,
    },
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum TemplateKind {
    Wing,
    Rotor,
    Aircraft,
    /// A small fixed-wing UAV wing.
    Uav,
    /// A high-aspect-ratio glider wing.
    Glider,
    /// A quadrotor propeller.
    Quadrotor,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(cli.command) {
        Ok(()) => ExitCode::SUCCESS,
        Err(msg) => {
            eprintln!("error: {msg}");
            ExitCode::FAILURE
        }
    }
}

/// Read a case file, naming the path in any error.
fn read_case(path: &Path) -> Result<String, String> {
    std::fs::read_to_string(path).map_err(|e| format!("cannot read {}: {e}", path.display()))
}

/// Write `text` to `out`, or to stdout when no file is given.
fn emit(text: &str, out: Option<&Path>) -> Result<(), String> {
    match out {
        Some(path) => {
            std::fs::write(path, text).map_err(|e| format!("cannot write {}: {e}", path.display()))
        }
        None => {
            print!("{text}");
            Ok(())
        }
    }
}

fn to_json<T: serde::Serialize>(value: &T) -> Result<String, String> {
    let mut json = serde_json::to_string_pretty(value).map_err(|e| e.to_string())?;
    json.push('\n');
    Ok(json)
}

fn run(command: Command) -> Result<(), String> {
    match command {
        Command::Wing {
            case,
            out,
            svg,
            check_mesh,
            mesh_sweep,
            spanwise,
        } => {
            let text = read_case(&case)?;
            for warning in wing_warnings(&text).map_err(|e| e.to_string())? {
                eprintln!("warning: {warning}");
            }
            if check_mesh {
                let report = wing_mesh_check(&text).map_err(|e| e.to_string())?;
                return emit(&report, out.as_deref());
            }
            if let Some(counts) = mesh_sweep {
                let csv = wing_mesh_sweep(&text, &counts).map_err(|e| e.to_string())?;
                return emit(&csv, out.as_deref());
            }
            if let Some(alpha) = spanwise {
                let csv = run_wing_spanwise(&text, alpha).map_err(|e| e.to_string())?;
                return emit(&csv, out.as_deref());
            }
            if let Some(path) = svg {
                let svg = wing_polar_svg(&text).map_err(|e| e.to_string())?;
                std::fs::write(&path, svg)
                    .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
            }
            let csv = run_wing(&text).map_err(|e| e.to_string())?;
            emit(&csv, out.as_deref())
        }
        Command::Derivatives {
            case,
            alpha,
            json,
            out,
        } => {
            let text = read_case(&case)?;
            let report = if json {
                run_derivatives_json(&text, alpha)
            } else {
                run_derivatives(&text, alpha)
            }
            .map_err(|e| e.to_string())?;
            emit(&format!("{}\n", report.trim_end()), out.as_deref())
        }
        Command::Vtk { case, alpha, out } => {
            let vtk = run_vtk(&read_case(&case)?, alpha).map_err(|e| e.to_string())?;
            emit(&vtk, out.as_deref())
        }
        Command::Rotor {
            case,
            json,
            rpm_sweep,
            collective_sweep,
            speed_sweep,
            out,
        } => {
            let sweep = match (rpm_sweep, collective_sweep, speed_sweep) {
                (Some(s), _, _) => Some((SweepVariable::Rpm, "--rpm-sweep", s)),
                (_, Some(s), _) => Some((SweepVariable::CollectiveDeg, "--collective-sweep", s)),
                (_, _, Some(s)) => Some((SweepVariable::AxialSpeed, "--speed-sweep", s)),
                _ => None,
            };
            if let Some((variable, flag, spec)) = sweep {
                let parts: Vec<f64> = spec
                    .split(':')
                    .map(|p| p.trim().parse::<f64>())
                    .collect::<Result<_, _>>()
                    .map_err(|_| format!("{flag} needs START:STOP:STEP numbers"))?;
                let [start, stop, step] = parts[..] else {
                    return Err(format!("{flag} needs exactly START:STOP:STEP"));
                };
                let csv = run_rotor_sweep_over(&read_case(&case)?, variable, start, stop, step)
                    .map_err(|e| e.to_string())?;
                return emit(&csv, out.as_deref());
            }
            let report = run_rotor(&read_case(&case)?).map_err(|e| e.to_string())?;
            let text = if json {
                to_json(&report)?
            } else {
                rotor_text(&report)
            };
            emit(&text, out.as_deref())
        }
        Command::Trim {
            aircraft,
            speed,
            speed_sweep,
            altitude,
            bank,
            json,
            csv,
            out,
        } => {
            if let Some(spec) = speed_sweep {
                let parts: Vec<f64> = spec
                    .split(':')
                    .map(|p| p.trim().parse::<f64>())
                    .collect::<Result<_, _>>()
                    .map_err(|_| "--speed-sweep needs START:STOP:STEP numbers".to_string())?;
                let [start, stop, step] = parts[..] else {
                    return Err("--speed-sweep needs exactly START:STOP:STEP".into());
                };
                let csv = run_trim_speed_sweep(&read_case(&aircraft)?, altitude, start, stop, step)
                    .map_err(|e| e.to_string())?;
                return emit(&csv, out.as_deref());
            }
            // Clap requires --speed unless a sweep is given.
            let speed = speed.ok_or("--speed is required")?;
            let kind = match bank {
                Some(deg) => TrimKind::Turn {
                    bank_rad: deg.to_radians(),
                },
                None => TrimKind::Level,
            };
            let report = run_trim(&read_case(&aircraft)?, speed, altitude, kind)
                .map_err(|e| e.to_string())?;
            let text = if json {
                to_json(&report)?
            } else if csv {
                trim_csv(&report)
            } else {
                trim_text(&report)
            };
            emit(&text, out.as_deref())
        }
        Command::Report {
            aircraft,
            speed,
            altitude,
            out,
        } => {
            let md = aircraft_report(&read_case(&aircraft)?, speed, altitude)
                .map_err(|e| e.to_string())?;
            emit(&md, out.as_deref())
        }
        Command::Avl { action } => match action {
            AvlCommand::Import {
                deck,
                alpha,
                speed,
                out,
            } => {
                let summary = avl_import_summary(&read_case(&deck)?, alpha, speed)
                    .map_err(|e| e.to_string())?;
                emit(&summary, out.as_deref())
            }
            AvlCommand::Export { case, out } => {
                let deck = wing_to_avl(&read_case(&case)?).map_err(|e| e.to_string())?;
                emit(&deck, out.as_deref())
            }
        },
        Command::Template { kind, out } => {
            let text = match kind {
                TemplateKind::Wing => TEMPLATE,
                TemplateKind::Rotor => ROTOR_TEMPLATE,
                TemplateKind::Aircraft => AIRCRAFT_TEMPLATE,
                TemplateKind::Uav => UAV_TEMPLATE,
                TemplateKind::Glider => GLIDER_TEMPLATE,
                TemplateKind::Quadrotor => QUADROTOR_TEMPLATE,
            };
            emit(text, out.as_deref())
        }
        Command::Isa {
            altitudes,
            delta_isa,
            imperial,
            out,
        } => {
            let csv = isa_table(&altitudes, delta_isa, imperial).map_err(|e| e.to_string())?;
            emit(&csv, out.as_deref())
        }
    }
}
