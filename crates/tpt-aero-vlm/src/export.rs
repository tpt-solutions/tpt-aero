//! Export of the lattice and its solution for visualization.
//!
//! [`to_vtk`] writes legacy ASCII VTK polydata (one quad per panel, with the
//! circulation as cell data), readable by ParaView, VisIt, PyVista and most
//! other viewers. Panels are drawn from the bound-vortex line (quarter chord)
//! with the panel chord, so tapered panels are shown with a constant chord
//! per panel: a picture, not exact geometry.

use alloc::{format, string::String, vec::Vec};
use core::fmt::Write as _;

use tpt_aero_core::Vec3;

use crate::solver::{VlmError, VlmModel};

/// One panel quad: leading-edge and trailing-edge corners at both bound-vortex ends.
fn panel_corners(model: &VlmModel) -> Vec<[Vec3<f64>; 4]> {
    model
        .lattice()
        .elements()
        .iter()
        .map(|e| {
            // The bound vortex sits at the panel quarter chord; chords run
            // along body x for every surface.
            let le = Vec3::new(-0.25 * e.chord, 0.0, 0.0);
            let te = Vec3::new(0.75 * e.chord, 0.0, 0.0);
            [
                e.bound_a + le,
                e.bound_b + le,
                e.bound_b + te,
                e.bound_a + te,
            ]
        })
        .collect()
}

/// Render `model` as legacy ASCII VTK polydata.
///
/// Every panel carries its `surface` and `strip` index as cell data. With
/// `gamma` (one circulation per element, as in `VlmSolution::gamma`) the file
/// also carries `gamma` (m²/s), ready to colour by in the viewer.
///
/// # Errors
/// [`VlmError::DimensionMismatch`] if `gamma` does not have one entry per element.
pub fn to_vtk(model: &VlmModel, title: &str, gamma: Option<&[f64]>) -> Result<String, VlmError> {
    let elements = model.lattice().elements();
    if gamma.is_some_and(|g| g.len() != elements.len()) {
        return Err(VlmError::DimensionMismatch);
    }
    let quads = panel_corners(model);
    let n = quads.len();
    let one_line_title: String = title
        .chars()
        .map(|c| if c.is_control() { ' ' } else { c })
        .collect();

    let mut out = String::new();
    // Writing to a String cannot fail.
    let _ = writeln!(out, "# vtk DataFile Version 3.0");
    let _ = writeln!(out, "{one_line_title}");
    let _ = writeln!(out, "ASCII");
    let _ = writeln!(out, "DATASET POLYDATA");
    let _ = writeln!(out, "POINTS {} float", 4 * n);
    for quad in &quads {
        for p in quad {
            let _ = writeln!(out, "{:.6} {:.6} {:.6}", p.x, p.y, p.z);
        }
    }
    let _ = writeln!(out, "POLYGONS {} {}", n, 5 * n);
    for i in 0..n {
        let b = 4 * i;
        let _ = writeln!(out, "4 {} {} {} {}", b, b + 1, b + 2, b + 3);
    }
    let _ = writeln!(out, "CELL_DATA {n}");
    let mut scalar = |name: &str, kind: &str, values: &mut dyn Iterator<Item = String>| {
        let _ = writeln!(out, "SCALARS {name} {kind} 1");
        let _ = writeln!(out, "LOOKUP_TABLE default");
        for v in values {
            let _ = writeln!(out, "{v}");
        }
    };
    scalar(
        "surface",
        "int",
        &mut elements.iter().map(|e| format!("{}", e.surface)),
    );
    scalar(
        "strip",
        "int",
        &mut elements.iter().map(|e| format!("{}", e.strip)),
    );
    if let Some(g) = gamma {
        scalar("gamma", "float", &mut g.iter().map(|v| format!("{v:.8e}")));
    }
    Ok(out)
}
