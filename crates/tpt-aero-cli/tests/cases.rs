//! Runs the ready case files in `cases/` and compares every number with the
//! expected output stored beside each one (relative tolerance 1e-6).

use std::path::PathBuf;

fn cases_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../cases")
}

fn read(name: &str) -> String {
    std::fs::read_to_string(cases_dir().join(name))
        .unwrap_or_else(|e| panic!("cannot read {name}: {e}"))
}

/// Compare two texts token by token: numbers to a relative 1e-6, other tokens exactly.
fn assert_close(actual: &str, expected: &str, what: &str) {
    let split = |s: &str| -> Vec<String> {
        s.split(|c: char| c == ',' || c.is_whitespace())
            .filter(|t| !t.is_empty())
            .map(str::to_owned)
            .collect()
    };
    let (a, e) = (split(actual), split(expected));
    assert_eq!(a.len(), e.len(), "{what}: token count differs");
    for (x, y) in a.iter().zip(&e) {
        match (x.parse::<f64>(), y.parse::<f64>()) {
            (Ok(xv), Ok(yv)) => {
                let scale = yv.abs().max(1.0e-9);
                assert!(
                    (xv - yv).abs() <= 1.0e-6 * scale,
                    "{what}: {x} vs expected {y}"
                );
            }
            _ => assert_eq!(x, y, "{what}"),
        }
    }
}

#[test]
fn uav_case_matches_expected() {
    let out = tpt_aero_cli::run_wing(&read("uav.toml")).expect("runs");
    assert_close(&out, &read("uav.expected.csv"), "uav wing sweep");
}

#[test]
fn glider_case_matches_expected() {
    let out = tpt_aero_cli::run_wing(&read("glider.toml")).expect("runs");
    assert_close(&out, &read("glider.expected.csv"), "glider wing sweep");
}

#[test]
fn quadrotor_case_matches_expected() {
    let report = tpt_aero_cli::run_rotor(&read("quadrotor.toml")).expect("runs");
    let got = serde_json::to_string_pretty(&report).expect("json");
    let expected: serde_json::Value =
        serde_json::from_str(&read("quadrotor.expected.json")).expect("json");
    let got: serde_json::Value = serde_json::from_str(&got).expect("json");
    for (key, value) in expected.as_object().expect("object") {
        let g = got[key].as_f64().expect("number");
        let e = value.as_f64().expect("number");
        assert!(
            (g - e).abs() <= 1.0e-6 * e.abs().max(1.0e-9),
            "{key}: {g} vs {e}"
        );
    }
}

#[test]
fn aircraft_trim_matches_expected() {
    let report = tpt_aero_cli::run_trim(
        &read("aircraft.toml"),
        50.0,
        0.0,
        tpt_aero_cli::TrimKind::Level,
    )
    .expect("trims");
    assert_close(
        &tpt_aero_cli::trim_text(&report),
        &read("aircraft.trim-50ms.expected.txt"),
        "aircraft trim",
    );
}
