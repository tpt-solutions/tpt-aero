//! The `wing --mesh-sweep` table: one CSV row per spanwise strip count.

use tpt_aero_cli::{wing_mesh_sweep, CliError};

const GLIDER: &str = include_str!("../../../cases/glider.toml");

#[test]
fn mesh_sweep_has_one_row_per_strip_count() {
    let csv = wing_mesh_sweep(GLIDER, &[8, 16, 32]).expect("sweep solves");
    let rows: Vec<&str> = csv.lines().collect();
    assert_eq!(rows[0], "n_span,cl,cdi");
    assert_eq!(rows.len(), 4);
    for (row, n) in rows[1..].iter().zip([8, 16, 32]) {
        let fields: Vec<&str> = row.split(',').collect();
        assert_eq!(fields[0], n.to_string());
        let cl: f64 = fields[1].parse().expect("cl parses");
        let cdi: f64 = fields[2].parse().expect("cdi parses");
        assert!(cl.is_finite() && cdi.is_finite() && cdi >= 0.0, "{row}");
    }
}

#[test]
fn mesh_sweep_rejects_empty_and_zero_counts() {
    assert!(matches!(
        wing_mesh_sweep(GLIDER, &[]),
        Err(CliError::Invalid(_))
    ));
    assert!(matches!(
        wing_mesh_sweep(GLIDER, &[8, 0]),
        Err(CliError::Invalid(_))
    ));
}
