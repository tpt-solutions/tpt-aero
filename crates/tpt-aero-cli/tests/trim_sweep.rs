//! Level trim over airspeed: one CSV row per speed, and a reversed range refused.

use tpt_aero_cli::run_trim_speed_sweep;

const CASE: &str = include_str!("../../../cases/aircraft.toml");

#[test]
fn speed_sweep_gives_one_row_per_speed() {
    let csv = run_trim_speed_sweep(CASE, 0.0, 30.0, 50.0, 10.0).expect("sweep");
    let rows: Vec<&str> = csv.lines().collect();
    assert_eq!(rows.len(), 1 + 3, "{csv}");
    assert!(rows[0].starts_with("speed_ms,alpha_deg,elevator_deg"));
    assert!(rows[1].starts_with("30,"), "{}", rows[1]);
    // Lift must fall as the speed rises for a level trim.
    let cl = |row: &str| row.split(',').nth(3).unwrap().parse::<f64>().unwrap();
    assert!(cl(rows[1]) > cl(rows[3]));
}

#[test]
fn speed_sweep_refuses_a_reversed_range_and_zero_step() {
    assert!(run_trim_speed_sweep(CASE, 0.0, 50.0, 30.0, 10.0).is_err());
    assert!(run_trim_speed_sweep(CASE, 0.0, 30.0, 50.0, 0.0).is_err());
}
