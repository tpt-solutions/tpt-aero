//! `tpt-aero` command-line entry point.

use std::process::ExitCode;

use tpt_aero_cli::{isa_csv, run_derivatives, run_wing, TEMPLATE};

const USAGE: &str = "tpt-aero - aerodynamics from the command line

USAGE:
    tpt-aero wing <case.toml> [--out <file.csv>]   run an alpha sweep (CSV to stdout by default)
    tpt-aero derivatives <case.toml> [--alpha <deg>]   stability derivatives (default alpha 2)
    tpt-aero template wing                         print a commented case file
    tpt-aero isa <altitude_m>...                   ISA properties as CSV
    tpt-aero --help
";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(msg) => {
            eprintln!("error: {msg}");
            ExitCode::FAILURE
        }
    }
}

fn run(args: &[String]) -> Result<(), String> {
    match args.first().map(String::as_str) {
        None | Some("--help" | "-h" | "help") => {
            print!("{USAGE}");
            Ok(())
        }
        Some("template") => match args.get(1).map(String::as_str) {
            Some("wing") => {
                print!("{TEMPLATE}");
                Ok(())
            }
            _ => Err("usage: tpt-aero template wing".into()),
        },
        Some("derivatives") => {
            let path = args
                .get(1)
                .ok_or("usage: tpt-aero derivatives <case.toml> [--alpha <deg>]")?;
            let alpha = match args.get(2).map(String::as_str) {
                None => 2.0,
                Some("--alpha") => args
                    .get(3)
                    .ok_or("--alpha needs a value in degrees")?
                    .parse::<f64>()
                    .map_err(|_| "--alpha must be a number".to_string())?,
                Some(other) => return Err(format!("unexpected argument: {other}")),
            };
            let text =
                std::fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))?;
            print!(
                "{}",
                run_derivatives(&text, alpha).map_err(|e| e.to_string())?
            );
            Ok(())
        }
        Some("isa") => {
            let alts = args[1..]
                .iter()
                .map(|a| a.parse::<f64>().map_err(|_| format!("not a number: {a}")))
                .collect::<Result<Vec<_>, _>>()?;
            if alts.is_empty() {
                return Err("usage: tpt-aero isa <altitude_m>...".into());
            }
            print!("{}", isa_csv(&alts).map_err(|e| e.to_string())?);
            Ok(())
        }
        Some("wing") => {
            let path = args
                .get(1)
                .ok_or("usage: tpt-aero wing <case.toml> [--out <file.csv>]")?;
            let out_path = match args.get(2).map(String::as_str) {
                None => None,
                Some("--out") => Some(args.get(3).ok_or("--out needs a file name")?),
                Some(other) => return Err(format!("unexpected argument: {other}")),
            };
            let text =
                std::fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))?;
            let csv = run_wing(&text).map_err(|e| e.to_string())?;
            match out_path {
                Some(p) => std::fs::write(p, csv).map_err(|e| format!("cannot write {p}: {e}")),
                None => {
                    print!("{csv}");
                    Ok(())
                }
            }
        }
        Some(other) => Err(format!("unknown command: {other}\n\n{USAGE}")),
    }
}
