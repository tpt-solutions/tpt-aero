use tpt_aero_cli::{
    isa_csv, run_derivatives, run_derivatives_json, run_rotor, run_trim, run_wing, CliError,
    TrimKind, AIRCRAFT_TEMPLATE, ROTOR_TEMPLATE, TEMPLATE,
};

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

fn asymmetry_case() -> String {
    // Tapered, swept, dihedral and cambered: every parameter that used to
    // produce a skewed (non-mirror-symmetric) wing.
    TEMPLATE
        .replace("tip_chord = 1.0 ", "tip_chord = 0.5 ")
        .replace("sweep_le = 0.0 ", "sweep_le = 1.2 ")
        .replace("dihedral_deg = 0.0", "dihedral_deg = 5.0")
        .replace(
            "alpha = { start = -4.0, stop = 10.0, step = 2.0 }",
            "alpha = { start = 0.0, stop = 6.0, step = 3.0 }",
        )
}

#[test]
fn tapered_swept_wing_is_mirror_symmetric() {
    let csv = run_wing(&asymmetry_case()).expect("valid case");
    for row in rows(&csv) {
        // cl_roll and cn must vanish for a mirror-symmetric wing in
        // symmetric flight (columns: alpha, cl, cdi, cm, cl_roll, cn, lift).
        assert!(row[4].abs() < 1e-7, "cl_roll = {} at {}", row[4], row[0]);
        assert!(row[5].abs() < 1e-7, "cn = {} at {}", row[5], row[0]);
    }
}

#[test]
fn derivatives_report_expected_signs() {
    let csv = run_derivatives(TEMPLATE, 4.0).expect("valid case");
    let get = |name: &str| -> f64 {
        csv.lines()
            .find_map(|l| l.strip_prefix(&format!("{name},")))
            .unwrap_or_else(|| panic!("missing {name}"))
            .parse()
            .unwrap()
    };
    assert!((4.5..5.5).contains(&get("cl_alpha")));
    assert!(get("cm_q") < 0.0, "pitch damping");
    assert!(get("cl_p") < 0.0, "roll damping");
    assert!(run_derivatives(TEMPLATE, f64::NAN).is_err());
}

#[test]
fn vtk_export_contains_the_lattice_and_circulation() {
    let vtk = tpt_aero_cli::run_vtk(TEMPLATE, 4.0).expect("valid case");
    assert!(vtk.starts_with("# vtk DataFile Version 3.0"));
    // Template mesh: 24 x 4 = 96 panels.
    assert!(vtk.contains("POINTS 384 float"));
    assert!(vtk.contains("SCALARS gamma float 1"));
    assert!(tpt_aero_cli::run_vtk(TEMPLATE, f64::NAN).is_err());
}

#[test]
fn derivatives_json_matches_the_csv_values() {
    let json: serde_json::Value =
        serde_json::from_str(&run_derivatives_json(TEMPLATE, 4.0).expect("valid case"))
            .expect("valid JSON");
    assert_eq!(json["alpha_deg"], 4.0);
    assert_eq!(json["speed"], 50.0);
    let csv = run_derivatives(TEMPLATE, 4.0).unwrap();
    let csv_value = |name: &str| -> f64 {
        csv.lines()
            .find_map(|l| l.strip_prefix(&format!("{name},")))
            .unwrap()
            .parse()
            .unwrap()
    };
    let cl_alpha = json["longitudinal"]["cl_alpha"].as_f64().unwrap();
    assert!((cl_alpha - csv_value("cl_alpha")).abs() < 1e-5);
    let cl_p = json["lateral"]["cl_p"].as_f64().unwrap();
    assert!((cl_p - csv_value("cl_p")).abs() < 1e-5);
    assert!(run_derivatives_json(TEMPLATE, f64::NAN).is_err());
}

#[test]
fn rotor_template_gives_consistent_hover_loads() {
    let r = run_rotor(ROTOR_TEMPLATE).expect("template is a valid rotor case");
    assert!(r.thrust_n > 0.0 && r.torque_nm > 0.0 && r.power_w > 0.0);
    assert!(
        r.figure_of_merit > 0.0 && r.figure_of_merit <= 1.0,
        "FM {}",
        r.figure_of_merit
    );
    // Dimensional forms follow the coefficients: T = C_T ρ A (ΩR)².
    let omega = 6000.0 * core::f64::consts::TAU / 60.0;
    let area = core::f64::consts::PI * 0.5 * 0.5;
    let t = r.ct * 1.225 * area * (omega * 0.5).powi(2);
    assert!((r.thrust_n / t - 1.0).abs() < 1e-6);
}

#[test]
fn rotor_rejects_bad_input() {
    let no_rpm = ROTOR_TEMPLATE.replace("rpm = 6000.0", "rpm = 0.0");
    assert!(matches!(run_rotor(&no_rpm), Err(CliError::Invalid(_))));
    let big_hub = ROTOR_TEMPLATE.replace("hub_radius = 0.05", "hub_radius = 0.6");
    assert!(matches!(run_rotor(&big_hub), Err(CliError::Invalid(_))));
    let both_density = ROTOR_TEMPLATE.replace("# density = 1.225", "density = 1.225");
    assert!(matches!(
        run_rotor(&both_density),
        Err(CliError::Invalid(_))
    ));
    let typo = ROTOR_TEMPLATE.replace("radius = 0.5", "radus = 0.5");
    assert!(matches!(run_rotor(&typo), Err(CliError::Parse(_))));
}

#[test]
fn level_trim_balances_weight_and_drag() {
    let t = run_trim(AIRCRAFT_TEMPLATE, 50.0, 0.0, TrimKind::Level).expect("template trims");
    assert_eq!(t.kind, "level");
    assert!(
        t.alpha_deg > 0.0 && t.alpha_deg < 10.0,
        "α = {}",
        t.alpha_deg
    );
    // Thrust equals drag: T = (C_Di + C_D0) q S, with C_L = W/(qS), so the
    // thrust is at least the parasite drag q S C_D0 = 0.03 q S.
    let thrust = t.thrust_n.expect("level trim reports thrust");
    let q_s = 0.5 * 1.225 * 50.0 * 50.0 * 12.0;
    assert!(thrust > 0.03 * q_s && thrust < 0.06 * q_s, "T = {thrust}");
    assert_eq!(t.load_factor, 1.0);
}

#[test]
fn coordinated_turn_trim_has_the_turn_load_factor_and_mirror_symmetry() {
    let bank = 20.0_f64.to_radians();
    let right = run_trim(
        AIRCRAFT_TEMPLATE,
        50.0,
        0.0,
        TrimKind::Turn { bank_rad: bank },
    )
    .unwrap();
    let left = run_trim(
        AIRCRAFT_TEMPLATE,
        50.0,
        0.0,
        TrimKind::Turn { bank_rad: -bank },
    )
    .unwrap();
    assert!(
        (right.load_factor - 1.0 / bank.cos()).abs() < 1e-4,
        "n = {}",
        right.load_factor
    );
    assert!((right.alpha_deg - left.alpha_deg).abs() < 1e-6);
    assert!((right.elevator_deg - left.elevator_deg).abs() < 1e-6);
    assert!((right.beta_deg + left.beta_deg).abs() < 1e-6);
    assert!((right.aileron_deg.unwrap() + left.aileron_deg.unwrap()).abs() < 1e-6);
    assert!((right.rudder_deg.unwrap() + left.rudder_deg.unwrap()).abs() < 1e-6);
    assert!(right.turn_radius_m.unwrap() > 100.0);
}

#[test]
fn trim_rejects_bad_input() {
    assert!(matches!(
        run_trim(AIRCRAFT_TEMPLATE, 0.0, 0.0, TrimKind::Level),
        Err(CliError::Invalid(_))
    ));
    let bad_kind = AIRCRAFT_TEMPLATE.replace("kind = \"vertical\"", "kind = \"canard\"");
    assert!(matches!(
        run_trim(&bad_kind, 50.0, 0.0, TrimKind::Level),
        Err(CliError::Invalid(_))
    ));
    // A turn needs an aileron and a rudder.
    let no_aileron = AIRCRAFT_TEMPLATE.replace("name = \"aileron\"", "name = \"flap\"");
    assert!(matches!(
        run_trim(&no_aileron, 50.0, 0.0, TrimKind::Turn { bank_rad: 0.3 }),
        Err(CliError::Invalid(_))
    ));
    assert!(matches!(
        run_trim(AIRCRAFT_TEMPLATE, 50.0, 200_000.0, TrimKind::Level),
        Err(CliError::Invalid(_))
    ));
}
