use tpt_aero_core::{Vec3, WingPlanform, WingSection};
use tpt_aero_vlm::export::to_vtk;
use tpt_aero_vlm::{CamberModel, Surface, VlmError, VlmModel};

fn parse_points(vtk: &str) -> Vec<[f64; 3]> {
    let mut lines = vtk.lines();
    let header = lines.find(|l| l.starts_with("POINTS")).expect("POINTS");
    let n: usize = header.split_whitespace().nth(1).unwrap().parse().unwrap();
    lines
        .take(n)
        .map(|l| {
            let v: Vec<f64> = l.split_whitespace().map(|x| x.parse().unwrap()).collect();
            [v[0], v[1], v[2]]
        })
        .collect()
}

/// The values of the named cell-data scalar.
fn scalar_values(vtk: &str, name: &str) -> Vec<String> {
    let mut lines = vtk.lines();
    lines
        .find(|l| l.starts_with(&format!("SCALARS {name} ")))
        .unwrap_or_else(|| panic!("no scalar {name}"));
    lines.next(); // LOOKUP_TABLE
    lines
        .take_while(|l| !l.starts_with("SCALARS") && !l.is_empty())
        .map(str::to_string)
        .collect()
}

#[test]
fn vtk_is_a_well_formed_polydata_file_with_circulation() {
    let wing = WingPlanform::rectangular(8.0, 1.0);
    let model = VlmModel::from_planform(&wing, 8, 3, CamberModel::None).unwrap();
    let sol = model.solve(0.08, 40.0, 1.225, &[]).unwrap();
    let vtk = to_vtk(&model, "test\nwing", Some(&sol.gamma)).unwrap();

    let mut lines = vtk.lines();
    assert_eq!(lines.next(), Some("# vtk DataFile Version 3.0"));
    assert_eq!(
        lines.next(),
        Some("test wing"),
        "title is one sanitized line"
    );
    assert_eq!(lines.next(), Some("ASCII"));
    assert_eq!(lines.next(), Some("DATASET POLYDATA"));

    let n = model.lattice().len();
    let pts = parse_points(&vtk);
    assert_eq!(pts.len(), 4 * n);
    assert!(pts.iter().flatten().all(|v| v.is_finite()));
    // Wing planform: |y| <= 4, chord from the leading to the trailing edge.
    assert!(pts.iter().all(|p| p[1].abs() <= 4.0 + 1e-9));
    assert!(pts.iter().all(|p| p[0] >= -1e-9 && p[0] <= 1.0 + 1e-9));

    assert!(vtk.contains(&format!("POLYGONS {} {}", n, 5 * n)));
    assert!(vtk.contains(&format!("CELL_DATA {n}")));
    let gamma: Vec<f64> = scalar_values(&vtk, "gamma")
        .iter()
        .map(|v| v.parse().unwrap())
        .collect();
    assert_eq!(gamma.len(), n);
    for (a, b) in gamma.iter().zip(&sol.gamma) {
        assert!((a - b).abs() <= 1e-7 * b.abs().max(1.0));
    }
    assert_eq!(scalar_values(&vtk, "strip").len(), n);
}

#[test]
fn surfaces_are_tagged_and_fin_panels_stand_up() {
    let fin = WingPlanform::from_sections(vec![
        WingSection {
            le: Vec3::new(0.0, 0.0, 0.0),
            chord: 0.8,
            twist: 0.0,
        },
        WingSection {
            le: Vec3::new(0.0, 1.5, 0.0),
            chord: 0.5,
            twist: 0.0,
        },
    ])
    .unwrap();
    let model = VlmModel::from_surfaces(&[
        Surface::horizontal(WingPlanform::rectangular(6.0, 1.0), 6, 2),
        Surface::vertical(fin, 3, 2).at(3.5, 0.0, 0.0),
    ])
    .unwrap();
    let vtk = to_vtk(&model, "aircraft", None).unwrap();
    assert!(!vtk.contains("SCALARS gamma"));
    let surface: Vec<i32> = scalar_values(&vtk, "surface")
        .iter()
        .map(|v| v.parse().unwrap())
        .collect();
    assert_eq!(surface.iter().filter(|s| **s == 0).count(), 6 * 2);
    assert_eq!(surface.iter().filter(|s| **s == 1).count(), 3 * 2);
    // The fin panels (last in the file) reach up to z = 1.5 and stay on y = 0.
    let pts = parse_points(&vtk);
    let fin_pts = &pts[4 * 12..];
    assert!(fin_pts.iter().all(|p| p[1].abs() < 1e-9));
    assert!(fin_pts.iter().map(|p| p[2]).fold(0.0, f64::max) > 1.0);
}

#[test]
fn mismatched_gamma_is_rejected() {
    let model = VlmModel::from_planform(
        &WingPlanform::rectangular(4.0, 1.0),
        4,
        2,
        CamberModel::None,
    )
    .unwrap();
    assert!(matches!(
        to_vtk(&model, "x", Some(&[0.0, 1.0])),
        Err(VlmError::DimensionMismatch)
    ));
}
