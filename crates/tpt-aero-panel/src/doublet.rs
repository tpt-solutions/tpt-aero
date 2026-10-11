//! 2-D lifting flow: constant-strength source panels plus a uniform vortex
//! sheet on the same contour, closed by a trailing-edge Kutta condition.
//!
//! Unknowns: one source strength per panel (thickness and incidence response)
//! and one uniform vortex-sheet strength `g` (circulation per unit length)
//! over the whole contour. Equations: flow tangency at every panel midpoint,
//! plus the Kutta condition that the tangential velocities on the two panels
//! meeting at the trailing edge are equal and opposite.
//!
//! Panel influences are the exact integrals for a constant-strength segment
//! (Hess & Smith 1967; Katz & Plotkin, *Low-Speed Aerodynamics*, ch. 11), not
//! a Gauss-quadrature approximation. The bound circulation sits on the body
//! surface rather than at an interior point, so the surface speed near the
//! trailing edge and mid-chord is not dominated by a nearby vortex core.
//!
//! The module keeps its historical name; it is not a doublet-panel method.

use alloc::vec;
use alloc::vec::Vec;

use tpt_aero_core::linalg::lu_solve;

use crate::source::{signed_area2, PanelError, Point2};

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// Result of a lifting panel analysis.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DoubletSolution {
    /// Source strength per unit length on each panel, `σ_j` (velocity units).
    pub strengths: Vec<f64>,
    /// Total circulation `Γ = −g·P` (m²/s), with `P` the contour length. Sign
    /// is chosen so that positive lift gives positive `Γ`.
    pub circulation: f64,
    /// Tangential velocity (normalized by V∞) at each panel midpoint.
    pub tangential_velocity: Vec<f64>,
    /// Pressure coefficient at each panel midpoint.
    pub cp: Vec<f64>,
    /// Section lift coefficient `C_L = 2Γ/(V∞·c)`.
    pub cl: f64,
}

/// Solve the lifting source/vortex-sheet panel system on a closed contour.
///
/// * `points`: contour vertices, counter-clockwise: start at the trailing
///   edge, run along the upper surface to the leading edge, then back along
///   the lower surface. The contour closes implicitly from the last vertex to
///   the first, so do **not** repeat the first vertex at the end; the Kutta
///   condition uses the first and last panels, which meet at the TE.
/// * `alpha_rad`: angle of attack.
/// * `v_inf`: freestream speed; `chord`: reference chord.
///
/// # Errors
/// [`PanelError::InvalidGeometry`] for degenerate or non-finite contours, and
/// for a non-finite `alpha`, or a non-positive or non-finite `v_inf` or `chord`;
/// [`PanelError::ClockwiseContour`] if the contour is ordered clockwise,
/// [`PanelError::Singular`] if the influence system is singular.
pub fn solve_doublet_panels(
    points: &[Point2],
    alpha_rad: f64,
    v_inf: f64,
    chord: f64,
) -> Result<DoubletSolution, PanelError> {
    let n = points.len();
    // Reject NaN and zero flight conditions before any division.
    let finite_input = points
        .iter()
        .all(|&(px, py)| px.is_finite() && py.is_finite());
    if !(finite_input
        && alpha_rad.is_finite()
        && v_inf.is_finite()
        && v_inf > 0.0
        && chord.is_finite()
        && chord > 0.0)
    {
        return Err(PanelError::InvalidGeometry);
    }
    if n < 4 {
        return Err(PanelError::InvalidGeometry);
    }
    if signed_area2(points) <= 0.0 {
        return Err(PanelError::ClockwiseContour);
    }

    // Panel geometry: start point, unit tangent, outward unit normal, length,
    // midpoint.
    let mut len = vec![0.0; n];
    let mut tx = vec![0.0; n];
    let mut ty = vec![0.0; n];
    let mut xc = vec![0.0; n];
    let mut yc = vec![0.0; n];
    for i in 0..n {
        let (x1, y1) = points[i];
        let (x2, y2) = points[(i + 1) % n];
        let (dx, dy) = (x2 - x1, y2 - y1);
        len[i] = (dx * dx + dy * dy).sqrt();
        if len[i].is_nan() || len[i] <= 0.0 {
            return Err(PanelError::InvalidGeometry);
        }
        tx[i] = dx / len[i];
        ty[i] = dy / len[i];
        xc[i] = 0.5 * (x1 + x2);
        yc[i] = 0.5 * (y1 + y2);
    }
    // Outward normal of a CCW contour is (ty, -tx); the interior lies to the
    // left of each panel, along (-ty, tx).
    let nx = |i: usize| ty[i];
    let ny = |i: usize| -tx[i];

    // Influence of panel j at point (px, py), per unit source strength and per
    // unit vortex-sheet strength, as global (x, y) velocity vectors.
    //
    // In the panel frame (origin at the start point, x along the panel, y to
    // the left) a constant source of strength 1 induces
    //   u = ln(r1^2/r2^2) / (4 pi),   v = (theta2 - theta1) / (2 pi),
    // and a constant vortex sheet of strength 1 (counter-clockwise) is the
    // same kernel rotated by 90 degrees: (u, v) -> (-v, u).
    let panel_influence = |j: usize, px: f64, py: f64| -> ((f64, f64), (f64, f64)) {
        let (x0, y0) = points[j];
        let (dx, dy) = (px - x0, py - y0);
        let xl = dx * tx[j] + dy * ty[j];
        let yl = -dx * ty[j] + dy * tx[j];
        let l = len[j];
        let r1_sq = (xl * xl + yl * yl).max(1.0e-30);
        let r2_sq = ((xl - l) * (xl - l) + yl * yl).max(1.0e-30);
        let theta1 = yl.atan2(xl);
        let theta2 = yl.atan2(xl - l);
        let src_u = (r1_sq.ln() - r2_sq.ln()) / (4.0 * core::f64::consts::PI);
        let src_v = (theta2 - theta1) / (2.0 * core::f64::consts::PI);
        // Local to global: a local vector (u, v) maps to u*t + v*left.
        let to_global = |u: f64, v: f64| (u * tx[j] - v * ty[j], u * ty[j] + v * tx[j]);
        (to_global(src_u, src_v), to_global(-src_v, src_u))
    };

    // Influence coefficients at midpoint i from panel j:
    //   kn_s, kt_s: source, normal and tangential to panel i
    //   kn_v, kt_v: vortex sheet, normal and tangential to panel i
    // Self-terms: a source contributes sigma/2 along the outward normal and a
    // vortex sheet contributes g/2 along the tangent on the outer side; the
    // principal-value parts vanish, so the self terms are set explicitly.
    let mut kn_s = vec![vec![0.0; n]; n];
    let mut kt_s = vec![vec![0.0; n]; n];
    let mut kn_v = vec![vec![0.0; n]; n];
    let mut kt_v = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                kn_s[i][j] = 0.5;
                continue;
            }
            let (s_vec, v_vec) = panel_influence(j, xc[i], yc[i]);
            kn_s[i][j] = s_vec.0 * nx(i) + s_vec.1 * ny(i);
            kt_s[i][j] = s_vec.0 * tx[i] + s_vec.1 * ty[i];
            kn_v[i][j] = v_vec.0 * nx(i) + v_vec.1 * ny(i);
            kt_v[i][j] = v_vec.0 * tx[i] + v_vec.1 * ty[i];
        }
    }
    // Own-panel vortex contribution to the outer tangential velocity.
    let self_tangential_vortex = 0.5;

    // Assemble: unknowns [sigma_0 .. sigma_{n-1}, g].
    let size = n + 1;
    let mut a = vec![vec![0.0; size]; size];
    let mut rhs = vec![0.0; size];
    let (cos_a, sin_a) = (alpha_rad.cos(), alpha_rad.sin());
    for i in 0..n {
        for j in 0..n {
            a[i][j] = kn_s[i][j];
            if j != i {
                a[i][n] += kn_v[i][j];
            }
        }
        rhs[i] = -v_inf * (cos_a * nx(i) + sin_a * ny(i));
    }

    // Kutta: V.t0 + V.t_{n-1} = 0 (equal speeds leaving the trailing edge).
    let k = n - 1;
    for j in 0..n {
        a[n][j] = kt_s[0][j] + kt_s[k][j];
    }
    let vortex_t0: f64 = (1..n).map(|j| kt_v[0][j]).sum::<f64>() + self_tangential_vortex;
    let vortex_tk: f64 =
        (0..n).filter(|&j| j != k).map(|j| kt_v[k][j]).sum::<f64>() + self_tangential_vortex;
    a[n][n] = vortex_t0 + vortex_tk;
    rhs[n] = -v_inf * (cos_a * (tx[0] + tx[k]) + sin_a * (ty[0] + ty[k]));

    let sol = lu_solve(&a, &rhs).map_err(|_| PanelError::Singular)?;
    let sigma: Vec<f64> = sol[..n].to_vec();
    let g = sol[n];

    // Circulation (counter-clockwise g·P) and lift. A positive lift needs the
    // clockwise-positive circulation Γ = -g·P, with C_L = 2Γ/(V∞ c).
    let perimeter: f64 = len.iter().sum();
    let circulation = -g * perimeter;
    let cl = 2.0 * circulation / (v_inf * chord);

    // Surface speed and Cp at the midpoints.
    let mut tangential = Vec::with_capacity(n);
    let mut cp = Vec::with_capacity(n);
    for i in 0..n {
        let mut vt = v_inf * (cos_a * tx[i] + sin_a * ty[i]);
        for j in 0..n {
            vt += sigma[j] * kt_s[i][j];
        }
        let mut vortex_t = self_tangential_vortex;
        for j in 0..n {
            if j != i {
                vortex_t += kt_v[i][j];
            }
        }
        vt += g * vortex_t;
        tangential.push(vt / v_inf);
        cp.push(1.0 - (vt / v_inf).powi(2));
    }

    Ok(DoubletSolution {
        strengths: sigma,
        circulation,
        tangential_velocity: tangential,
        cp,
        cl,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn symmetric_airfoil(n: usize) -> Vec<Point2> {
        // Thin symmetric NACA-like contour: TE → upper → LE → lower, closing
        // back to the exact TE point. The closed-te coefficient keeps the
        // surface points strictly off the TE so no panel degenerates.
        let thickness = |x: f64| {
            5.0 * 0.12
                * (0.2969 * x.sqrt() - 0.126_0 * x - 0.351_6 * x * x + 0.284_3 * x * x * x
                    - 0.103_6 * x * x * x * x)
        };
        let stations: Vec<f64> = (0..=n)
            .map(|i| 0.5 * (1.0 - (core::f64::consts::PI * i as f64 / n as f64).cos()))
            .collect();
        let mut pts: Vec<Point2> = vec![(1.0, 0.0)];
        for i in (1..n).rev() {
            let x = stations[i];
            pts.push((x, thickness(x)));
        }
        pts.push((0.0, 0.0));
        for i in 1..n {
            let x = stations[i];
            pts.push((x, -thickness(x)));
        }
        pts
    }

    #[test]
    fn flat_plate_lift_slope() {
        // Symmetric thin airfoil: C_L ≈ 2πα.
        let pts = symmetric_airfoil(36);
        let alpha = 5.0_f64.to_radians();
        let sol = solve_doublet_panels(&pts, alpha, 1.0, 1.0).expect("solvable");
        let expected = 2.0 * core::f64::consts::PI * alpha;
        assert!(
            (sol.cl - expected).abs() < 0.15,
            "cl = {} vs thin-airfoil {expected}",
            sol.cl
        );
    }

    #[test]
    fn zero_alpha_zero_lift_symmetric() {
        let pts = symmetric_airfoil(36);
        let sol = solve_doublet_panels(&pts, 0.0, 1.0, 1.0).expect("solvable");
        assert!(sol.cl.abs() < 1e-6, "cl = {}", sol.cl);
    }

    #[test]
    fn stagnation_cp_positive() {
        let pts = symmetric_airfoil(48);
        let sol = solve_doublet_panels(&pts, 3.0_f64.to_radians(), 1.0, 1.0).expect("ok");
        let max_cp = sol.cp.iter().fold(f64::NEG_INFINITY, |m, &c| m.max(c));
        assert!(max_cp > 0.5, "max cp {max_cp}");
        assert!(sol.cp.iter().all(|&c| c.is_finite() && c < 1.0 + 1e-9));
    }

    #[test]
    fn rejects_clockwise_contour() {
        let mut pts = symmetric_airfoil(24);
        pts.reverse();
        assert_eq!(
            solve_doublet_panels(&pts, 0.0, 1.0, 1.0).unwrap_err(),
            PanelError::ClockwiseContour
        );
    }

    #[test]
    fn rejects_degenerate() {
        assert_eq!(
            solve_doublet_panels(&[(0.0, 0.0), (0.5, 0.1), (1.0, 0.0)], 0.0, 1.0, 1.0).unwrap_err(),
            PanelError::InvalidGeometry
        );
    }

    #[test]
    fn zero_incidence_pressure_matches_source_solver() {
        // No circulation at zero incidence, so the lifting solver's surface
        // pressure must equal the non-lifting source solver's. The source
        // solver uses Gauss quadrature for its influences, so the two agree to
        // quadrature accuracy (about 3e-3 at this resolution), not exactly.
        let pts = symmetric_airfoil(60);
        let lifting = solve_doublet_panels(&pts, 0.0, 1.0, 1.0).expect("solvable");
        let source = crate::source::solve_source_panels(&pts, 0.0, 1.0).expect("solvable");
        for (a, b) in lifting.cp.iter().zip(&source.cp) {
            assert!((a - b).abs() < 5e-3, "lifting cp {a} vs source cp {b}");
        }
    }

    #[test]
    fn aft_chord_pressure_jump_matches_thin_airfoil_theory() {
        // Upper-minus-lower Cp of a thin section tends to the flat-plate value
        // 4 alpha sqrt((1 - x)/x) away from the leading edge. Pair each upper
        // midpoint with the lower midpoint nearest in x.
        let foil = tpt_aero_core::airfoil::naca4_closed("0012", 200).expect("valid NACA code");
        let pts: Vec<Point2> = foil.x.iter().zip(&foil.y).map(|(&x, &y)| (x, y)).collect();
        let alpha = 4.0_f64.to_radians();
        let sol = solve_doublet_panels(&pts, alpha, 50.0, 1.0).expect("solvable");
        let n = pts.len();
        let mid = |i: usize| 0.5 * (pts[i].0 + pts[(i + 1) % n].0);
        let mid_y = |i: usize| 0.5 * (pts[i].1 + pts[(i + 1) % n].1);
        let upper: Vec<usize> = (0..n).filter(|&i| mid_y(i) > 0.0).collect();
        let lower: Vec<usize> = (0..n).filter(|&i| mid_y(i) < 0.0).collect();
        let mut checked = 0;
        for &iu in &upper {
            let x = mid(iu);
            if !(0.6..=0.9).contains(&x) {
                continue;
            }
            let il = *lower
                .iter()
                .min_by(|&&a, &&b| (mid(a) - x).abs().total_cmp(&(mid(b) - x).abs()))
                .expect("lower surface");
            let jump = sol.cp[il] - sol.cp[iu];
            let theory = 4.0 * alpha * ((1.0 - x) / x).sqrt();
            assert!(
                (jump - theory).abs() < 0.15 * theory,
                "x/c {x:.3}: jump {jump:.3} vs thin-airfoil {theory:.3}"
            );
            checked += 1;
        }
        assert!(checked > 10, "only {checked} stations checked");
    }

    #[test]
    fn lifting_pressure_is_smooth_near_mid_chord() {
        // A bound vortex close to the surface used to leave a spike in Cp
        // near mid-chord. The surface Cp must now vary smoothly.
        let pts = symmetric_airfoil(100);
        let sol = solve_doublet_panels(&pts, 4.0_f64.to_radians(), 50.0, 1.0).expect("solvable");
        let worst = (1..sol.cp.len() - 1)
            .map(|i| (sol.cp[i] - 0.5 * (sol.cp[i - 1] + sol.cp[i + 1])).abs())
            .fold(0.0_f64, f64::max);
        assert!(worst < 0.2, "largest local Cp kink {worst}");
    }

    #[test]
    fn kutta_condition_holds_at_trailing_edge() {
        // Equal speeds leave the trailing edge on both sides.
        let pts = symmetric_airfoil(80);
        let sol = solve_doublet_panels(&pts, 5.0_f64.to_radians(), 1.0, 1.0).expect("solvable");
        let n = sol.tangential_velocity.len();
        let (vt_first, vt_last) = (sol.tangential_velocity[0], sol.tangential_velocity[n - 1]);
        assert!(
            (vt_first + vt_last).abs() < 1e-6,
            "Kutta residual {}",
            vt_first + vt_last
        );
    }
}

#[cfg(test)]
mod input_tests {
    use super::*;
    use alloc::vec;

    fn square() -> Vec<Point2> {
        vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    }

    #[test]
    fn non_finite_or_zero_flight_conditions_are_rejected() {
        assert!(solve_doublet_panels(&square(), 0.0, 10.0, 1.0).is_ok());
        assert!(solve_doublet_panels(&square(), f64::NAN, 10.0, 1.0).is_err());
        assert!(solve_doublet_panels(&square(), 0.0, 0.0, 1.0).is_err());
        assert!(solve_doublet_panels(&square(), 0.0, 10.0, 0.0).is_err());
        let mut nan_point = square();
        nan_point[2].0 = f64::NAN;
        assert!(solve_doublet_panels(&nan_point, 0.0, 10.0, 1.0).is_err());
    }
}
