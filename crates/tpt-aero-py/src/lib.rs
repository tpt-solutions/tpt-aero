//! Python bindings for tpt-aero.
//!
//! Each function takes the TOML case text the `tpt-aero` command reads and
//! returns the text it prints, so Python and the command line agree. Errors
//! raise `ValueError` with the same message the command line shows.
//!
//! ```python
//! import tpt_aero
//! csv = tpt_aero.wing_csv(open("uav.toml").read())
//! ```

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

fn err<E: core::fmt::Display>(e: E) -> PyErr {
    PyValueError::new_err(e.to_string())
}

/// The wing alpha sweep as CSV (`tpt-aero wing`).
#[pyfunction]
fn wing_csv(case: &str) -> PyResult<String> {
    tpt_aero_cli::run_wing(case).map_err(err)
}

/// The lift-drag polar of a wing case as an SVG document.
#[pyfunction]
fn wing_svg(case: &str) -> PyResult<String> {
    tpt_aero_cli::wing_polar_svg(case).map_err(err)
}

/// Validity warnings for a wing case (a list; empty when in the envelope).
#[pyfunction]
fn wing_warnings(case: &str) -> PyResult<Vec<String>> {
    tpt_aero_cli::wing_warnings(case).map_err(err)
}

/// The mesh-convergence estimate of a wing case (`tpt-aero wing --check-mesh`).
#[pyfunction]
fn wing_mesh_check(case: &str) -> PyResult<String> {
    tpt_aero_cli::wing_mesh_check(case).map_err(err)
}

/// The level trim of an aircraft case, as the `tpt-aero trim` text report.
#[pyfunction]
fn trim_text(case: &str, speed_ms: f64, altitude_m: f64) -> PyResult<String> {
    let report = tpt_aero_cli::run_trim(case, speed_ms, altitude_m, tpt_aero_cli::TrimKind::Level)
        .map_err(err)?;
    Ok(tpt_aero_cli::trim_text(&report))
}

/// The Markdown stability report for an aircraft case (`tpt-aero report`).
#[pyfunction]
fn report_markdown(case: &str, speed_ms: f64, altitude_m: f64) -> PyResult<String> {
    tpt_aero_cli::aircraft_report(case, speed_ms, altitude_m).map_err(err)
}

/// One BEMT operating point of a rotor case, as a dict of the report fields.
/// A dict loads straight into a pandas DataFrame.
#[pyfunction]
fn rotor<'py>(py: Python<'py>, case: &str) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    let r = tpt_aero_cli::run_rotor(case).map_err(err)?;
    let d = pyo3::types::PyDict::new(py);
    d.set_item("thrust_n", r.thrust_n)?;
    d.set_item("torque_nm", r.torque_nm)?;
    d.set_item("power_w", r.power_w)?;
    d.set_item("ct", r.ct)?;
    d.set_item("cq", r.cq)?;
    d.set_item("cp", r.cp)?;
    d.set_item("figure_of_merit", r.figure_of_merit)?;
    d.set_item("propeller_efficiency", r.propeller_efficiency)?;
    Ok(d)
}

/// The rotor RPM sweep as CSV (`tpt-aero rotor --rpm-sweep`).
#[pyfunction]
fn rotor_sweep_csv(case: &str, start_rpm: f64, stop_rpm: f64, step_rpm: f64) -> PyResult<String> {
    tpt_aero_cli::run_rotor_sweep(case, start_rpm, stop_rpm, step_rpm).map_err(err)
}

/// The standard atmosphere at each altitude (m), as CSV.
#[pyfunction]
fn isa_csv(altitudes_m: Vec<f64>) -> PyResult<String> {
    tpt_aero_cli::isa_table(&altitudes_m, 0.0, false).map_err(err)
}

/// The module: `import tpt_aero`.
#[pymodule]
fn tpt_aero(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(wing_csv, m)?)?;
    m.add_function(wrap_pyfunction!(wing_svg, m)?)?;
    m.add_function(wrap_pyfunction!(wing_warnings, m)?)?;
    m.add_function(wrap_pyfunction!(wing_mesh_check, m)?)?;
    m.add_function(wrap_pyfunction!(trim_text, m)?)?;
    m.add_function(wrap_pyfunction!(report_markdown, m)?)?;
    m.add_function(wrap_pyfunction!(rotor, m)?)?;
    m.add_function(wrap_pyfunction!(rotor_sweep_csv, m)?)?;
    m.add_function(wrap_pyfunction!(isa_csv, m)?)?;
    Ok(())
}
