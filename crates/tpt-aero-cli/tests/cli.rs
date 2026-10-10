use tpt_aero_cli::{isa_csv, run_wing, CliError, TEMPLATE};

fn rows(csv: &str) -> Vec<Vec<f64>> {
    csv.lines()
        .skip(1)
        .map(|l| l.split(',').map(|v| v.parse().unwrap()).collect())
        .collect()
}

#[test]
fn template_runs_and_matches_vlm() {
    let csv = run_wing(TEMPLATE).expect("template is a valid case");
    let r = rows(&csv);
    assert_eq!(r.len(), 8, "-4..10 step 2 is 8 points");
    // alpha = 0 row: flat plate has zero lift; positive alpha lifts.
    let zero = r.iter().find(|row| row[0] == 0.0).unwrap();
    assert!(zero[1].abs() < 1e-6, "CL(0) = {}", zero[1]);
    let a6 = r.iter().find(|row| row[0] == 6.0).unwrap();
    assert!(a6[1] > 0.5 && a6[1] < 0.6, "CL(6) = {}", a6[1]);
    // Dimensional lift is consistent with q S CL.
    let q_s = 0.5 * 1.225 * 50.0 * 50.0 * 10.0;
    assert!((a6[6] / (q_s * a6[1]) - 1.0).abs() < 0.02);
}

#[test]
fn rejects_bad_input() {
    let bad_step = TEMPLATE.replace("step = 2.0", "step = -2.0");
    assert!(matches!(run_wing(&bad_step), Err(CliError::Invalid(_))));
    let both = TEMPLATE.replace("# density = 1.225", "density = 1.225");
    assert!(matches!(run_wing(&both), Err(CliError::Invalid(_))));
    let typo = TEMPLATE.replace("span = 10.0", "spn = 10.0");
    assert!(matches!(run_wing(&typo), Err(CliError::Parse(_))));
    let zero_speed = TEMPLATE.replace("speed = 50.0", "speed = 0.0");
    assert!(matches!(run_wing(&zero_speed), Err(CliError::Invalid(_))));
    assert!(matches!(run_wing("not toml ["), Err(CliError::Parse(_))));
}

#[test]
fn isa_rows() {
    let csv = isa_csv(&[0.0, 11_000.0]).unwrap();
    let r = rows(&csv);
    assert!((r[0][3] - 1.225).abs() < 1e-3);
    assert!((r[1][1] - 216.65).abs() < 1e-6);
    assert!(isa_csv(&[-5.0]).is_err());
}
