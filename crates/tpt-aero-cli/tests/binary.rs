//! Binary-level checks: argument parsing, exit codes and file output of the
//! `tpt-aero` executable itself.

use std::path::PathBuf;
use std::process::{Command, Output};

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_tpt-aero"))
}

fn run(args: &[&str]) -> Output {
    bin().args(args).output().expect("binary runs")
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("tpt-aero-cli-test-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("temp dir");
    dir.join(name)
}

#[test]
fn version_flag_reports_the_crate_version() {
    let out = run(&["--version"]);
    assert!(out.status.success());
    let text = String::from_utf8_lossy(&out.stdout);
    assert!(text.contains(env!("CARGO_PKG_VERSION")), "{text}");
}

#[test]
fn every_subcommand_has_help() {
    for cmd in [
        "wing",
        "derivatives",
        "vtk",
        "rotor",
        "trim",
        "template",
        "isa",
    ] {
        let out = run(&[cmd, "--help"]);
        assert!(out.status.success(), "{cmd} --help failed");
        assert!(!out.stdout.is_empty(), "{cmd} --help printed nothing");
    }
}

#[test]
fn unknown_subcommand_and_bad_flags_exit_with_usage_error() {
    assert_eq!(run(&["nonsense"]).status.code(), Some(2));
    assert_eq!(run(&["isa", "1000", "--bogus"]).status.code(), Some(2));
    // Flags work before the positional arguments as well as after.
    let path = scratch("flag-first.csv");
    let p = path.to_str().expect("utf-8 path");
    assert!(run(&["isa", "--out", p, "1000"]).status.success());
    assert!(path.exists());
}

#[test]
fn runtime_errors_exit_one_with_a_message() {
    let out = run(&["isa", "90000"]);
    assert_eq!(out.status.code(), Some(1));
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.starts_with("error:"), "{err}");
    assert!(err.contains("84852"), "range not reported: {err}");
    assert!(
        !err.contains("AltitudeOutOfRange"),
        "debug output leaked: {err}"
    );
}

#[test]
fn missing_case_file_names_the_path() {
    let out = run(&["wing", "no-such-file.toml"]);
    assert_eq!(out.status.code(), Some(1));
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("no-such-file.toml"), "{err}");
}

#[test]
fn isa_writes_to_the_out_file() {
    let path = scratch("isa.csv");
    let out = run(&["isa", "0", "5000", "-o", path.to_str().expect("utf-8 path")]);
    assert!(out.status.success());
    assert!(out.stdout.is_empty(), "nothing should go to stdout");
    let csv = std::fs::read_to_string(&path).expect("file written");
    assert!(csv.lines().count() >= 3, "{csv}");
}

#[test]
fn bank_past_ninety_degrees_is_refused() {
    let aircraft = scratch("aircraft.toml");
    std::fs::write(&aircraft, tpt_aero_cli::AIRCRAFT_TEMPLATE).expect("write");
    let out = run(&[
        "trim",
        aircraft.to_str().expect("utf-8 path"),
        "--speed",
        "50",
        "--bank",
        "95",
    ]);
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("90 deg"));
}

#[test]
fn negative_numeric_values_are_parsed_as_values_not_flags() {
    // -100 m is below the ISA range: it must reach the solver and fail there
    // (exit 1), not be rejected by the argument parser (exit 2).
    let out = run(&["isa", "-100"]);
    assert_eq!(
        out.status.code(),
        Some(1),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
}
