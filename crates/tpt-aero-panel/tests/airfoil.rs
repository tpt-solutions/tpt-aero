//! Panel solvers on a closed-trailing-edge NACA section.

use tpt_aero_core::airfoil::naca4_closed;
use tpt_aero_panel::{solve_doublet_panels, solve_source_panels, PanelError};

fn contour(code: &str, n: usize) -> Vec<(f64, f64)> {
    let a = naca4_closed(code, n).expect("valid code");
    a.x.iter().zip(&a.y).map(|(&x, &y)| (x, y)).collect()
}

#[test]
fn naca0012_source_panels_match_reference_suction() {
    let pts = contour("0012", 80);
    let sol = solve_source_panels(&pts, 0.0, 1.0).expect("solvable");
    // No gap panel: every Cp is physical, even at the trailing edge.
    assert!(sol.cp.iter().all(|&c| c.is_finite() && c <= 1.0 + 1e-9));
    let (imin, cp_min) =
        sol.cp.iter().enumerate().fold(
            (0, f64::INFINITY),
            |b, (i, &c)| if c < b.1 { (i, c) } else { b },
        );
    // Thin-airfoil/XFOIL-inviscid reference: Cp_min about -0.43 near x/c = 0.1.
    assert!((cp_min + 0.43).abs() < 0.04, "Cp_min = {cp_min}");
    let x_min = pts[imin].0;
    assert!(
        (0.05..=0.2).contains(&x_min),
        "suction peak at x/c = {x_min}"
    );
    // Trailing-edge Cp stays small and positive (no spurious -0.59 spike).
    assert!(
        sol.cp[0] > -0.05 && sol.cp[0] < 0.5,
        "TE Cp = {}",
        sol.cp[0]
    );
}

#[test]
fn naca0012_lift_slope_close_to_thick_airfoil_theory() {
    let pts = contour("0012", 80);
    let alpha = 5.0_f64.to_radians();
    let sol = solve_doublet_panels(&pts, alpha, 1.0, 1.0).expect("solvable");
    // Thick-airfoil theory: dCL/dalpha = 2 pi (1 + 0.77 t/c) = 6.86 /rad.
    let expected = 2.0 * std::f64::consts::PI * (1.0 + 0.77 * 0.12) * alpha;
    assert!(
        (sol.cl - expected).abs() < 0.04,
        "CL = {} vs {expected}",
        sol.cl
    );
}

#[test]
fn clockwise_contours_are_rejected() {
    let mut pts = contour("0012", 40);
    pts.reverse();
    assert_eq!(
        solve_source_panels(&pts, 0.0, 1.0).unwrap_err(),
        PanelError::ClockwiseContour
    );
}
