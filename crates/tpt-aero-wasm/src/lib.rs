//! WebAssembly bindings for tpt-aero.
//!
//! Each function takes the same TOML case text the command line reads, and
//! returns the same text the command line prints, so a browser page and
//! `tpt-aero` give identical numbers. Errors come back as a JavaScript string.
//!
//! ```js
//! import init, { wing_svg, trim_text, isa_table } from "./tpt_aero_wasm.js";
//! await init();
//! document.body.innerHTML = wing_svg(caseText);
//! ```

use wasm_bindgen::prelude::*;

fn to_js<E: core::fmt::Display>(e: E) -> JsValue {
    JsValue::from_str(&e.to_string())
}

/// The wing alpha sweep as CSV (the `tpt-aero wing` output).
///
/// # Errors
/// A JavaScript string with the message for a malformed or invalid case.
#[wasm_bindgen]
pub fn wing_csv(case: &str) -> Result<String, JsValue> {
    tpt_aero_cli::run_wing(case).map_err(to_js)
}

/// The lift-drag polar of a wing case as an SVG document.
///
/// # Errors
/// As [`wing_csv`].
#[wasm_bindgen]
pub fn wing_svg(case: &str) -> Result<String, JsValue> {
    tpt_aero_cli::wing_polar_svg(case).map_err(to_js)
}

/// Validity warnings for a wing case, one per line (empty when in envelope).
///
/// # Errors
/// As [`wing_csv`].
#[wasm_bindgen]
pub fn wing_warnings(case: &str) -> Result<String, JsValue> {
    tpt_aero_cli::wing_warnings(case)
        .map(|w| w.join("\n"))
        .map_err(to_js)
}

/// The level trim of an aircraft case, as the `tpt-aero trim` text report.
///
/// # Errors
/// As [`wing_csv`].
#[wasm_bindgen]
pub fn trim_text(case: &str, speed_ms: f64, altitude_m: f64) -> Result<String, JsValue> {
    let report = tpt_aero_cli::run_trim(case, speed_ms, altitude_m, tpt_aero_cli::TrimKind::Level)
        .map_err(to_js)?;
    Ok(tpt_aero_cli::trim_text(&report))
}

/// The Markdown stability report for an aircraft case.
///
/// # Errors
/// As [`wing_csv`].
#[wasm_bindgen]
pub fn report_markdown(case: &str, speed_ms: f64, altitude_m: f64) -> Result<String, JsValue> {
    tpt_aero_cli::aircraft_report(case, speed_ms, altitude_m).map_err(to_js)
}

/// One BEMT operating point of a rotor case, as the `tpt-aero rotor` text.
///
/// # Errors
/// As [`wing_csv`].
#[wasm_bindgen]
pub fn rotor_text(case: &str) -> Result<String, JsValue> {
    let report = tpt_aero_cli::run_rotor(case).map_err(to_js)?;
    Ok(tpt_aero_cli::rotor_text(&report))
}

/// The standard atmosphere at each altitude (metres), as CSV.
///
/// # Errors
/// As [`wing_csv`], for an altitude outside 0 to 84 852 m.
#[wasm_bindgen]
pub fn isa_table(altitudes_m: &[f64]) -> Result<String, JsValue> {
    tpt_aero_cli::isa_table(altitudes_m, 0.0, false).map_err(to_js)
}

/// The starting wing case, for a page to show in its editor.
#[wasm_bindgen]
pub fn wing_template() -> String {
    tpt_aero_cli::UAV_TEMPLATE.to_owned()
}
