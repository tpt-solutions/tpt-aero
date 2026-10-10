//! Two-dimensional source-panel method (Hess–Smith) for non-lifting bodies.
//!
//! Constant-strength source panels over a closed polygonal contour, Neumann
//! flow-tangency at the panel midpoints, solved with the native dense LU.
//! Produces the tangential surface velocity and pressure coefficient.

use alloc::{vec, vec::Vec};

use tpt_aero_vlm::lu_solve;

/// Errors from panel solves.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PanelError {
    /// The influence matrix was numerically singular.
    Singular,
    /// Geometry needs at least 3 distinct points forming a closed contour.
    InvalidGeometry,
}

impl core::fmt::Display for PanelError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            PanelError::Singular => write!(f, "panel influence matrix is singular"),
            PanelError::InvalidGeometry => {
                write!(f, "need ≥3 distinct points forming a closed contour")
            }
        }
    }
}

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// A 2-D point (x, y).
pub type Point2 = (f64, f64);

/// Result of a source-panel analysis.
#[derive(Debug, Clone)]
pub struct SourceSolution {
    /// Source strengths per panel.
    pub strengths: Vec<f64>,
    /// Tangential velocity at each panel midpoint (normalized by V∞).
    pub tangential_velocity: Vec<f64>,
    /// Pressure coefficient at each panel midpoint.
    pub cp: Vec<f64>,
}

/// Constant-strength source panels on a closed contour.
///
/// * `points`: contour vertices, ordered counter-clockwise for an exterior
///   flow (the panel normal points out of the body).
/// * `alpha_rad`: freestream angle of attack relative to the x axis.
///
/// The boundary condition is zero normal velocity at each panel midpoint;
/// the tangential velocity follows from the total induced field.
///
/// # Errors
/// [`PanelError::InvalidGeometry`] for degenerate contours,
/// [`PanelError::Singular`] if the influence system is singular.
pub fn solve_source_panels(
    points: &[Point2],
    alpha_rad: f64,
    v_inf: f64,
) -> Result<SourceSolution, PanelError> {
    let n = points.len();
    if n < 3 {
        return Err(PanelError::InvalidGeometry);
    }
    // Panel geometry: endpoints, midpoint, outward normal, length.
    let mut mid = Vec::with_capacity(n);
    let mut nrm = Vec::with_capacity(n);
    let mut len = Vec::with_capacity(n);
    let mut tangent = Vec::with_capacity(n);
    for i in 0..n {
        let (x1, y1) = points[i];
        let (x2, y2) = points[(i + 1) % n];
        let dx = x2 - x1;
        let dy = y2 - y1;
        let l = (dx * dx + dy * dy).sqrt();
        if l <= 0.0 {
            return Err(PanelError::InvalidGeometry);
        }
        len.push(l);
        mid.push(((x1 + x2) * 0.5, (y1 + y2) * 0.5));
        // CCW contour: outward normal is (+dy, −dx)/l.
        nrm.push((dy / l, -dx / l));
        tangent.push((dx / l, dy / l));
    }

    // Influence coefficients by 8-point Gauss–Legendre quadrature of the
    // point-source kernel v = (σ l w / 2π)·R/|R|² along each panel — immune
    // to the branch-cut pitfalls of the closed-form angle expressions. The
    // self-panel normal influence is the analytic ½ (Hess–Smith).
    const GAUSS: [(f64, f64); 4] = [
        (-0.861_136_311_594_052_6, 0.347_854_845_137_453_85),
        (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.861_136_311_594_052_6, 0.347_854_845_137_453_85),
    ];
    let panel_velocity =
        |i_target: usize, (xi, yi): Point2, j: usize, sigma_j: f64| -> (f64, f64) {
            if j == i_target {
                return (0.0, 0.0);
            }
            let (x1, y1) = points[j];
            let (x2, y2) = points[(j + 1) % n];
            let half = len[j] * 0.5;
            let mut ug = 0.0;
            let mut vg = 0.0;
            for &(g, w) in GAUSS.iter() {
                let s = half * (1.0 + g);
                let sx = x1 + (x2 - x1) * (s / len[j]);
                let sy = y1 + (y2 - y1) * (s / len[j]);
                let rx = xi - sx;
                let ry = yi - sy;
                let r2 = (rx * rx + ry * ry).max(1.0e-20);
                let f = sigma_j * half * w / (2.0 * core::f64::consts::PI);
                ug += f * rx / r2;
                vg += f * ry / r2;
            }
            (ug, vg)
        };

    let mut a = vec![vec![0.0; n]; n];
    let mut rhs = vec![0.0; n];
    for i in 0..n {
        let (xi, yi) = mid[i];
        let (nxi, nyi) = nrm[i];
        for j in 0..n {
            a[i][j] = if i == j {
                0.5
            } else {
                let (ug, vg) = panel_velocity(i, (xi, yi), j, 1.0);
                ug * nxi + vg * nyi
            };
        }
        // Freestream through the outward normal.
        rhs[i] = -v_inf * (alpha_rad.cos() * nxi + alpha_rad.sin() * nyi);
    }

    let sigma = lu_solve(&a, &rhs).map_err(|_| PanelError::Singular)?;

    // Surface speed: tangential component of (freestream + induction); the
    // self panel has zero tangential contribution at its own midpoint, and
    // the normal components cancel by the boundary condition, so the speed
    // is |V_t|.
    let mut tangential = Vec::with_capacity(n);
    let mut cp = Vec::with_capacity(n);
    for i in 0..n {
        let (xi, yi) = mid[i];
        let (txi, tyi) = tangent[i];
        let mut ug = alpha_rad.cos() * v_inf;
        let mut vg = alpha_rad.sin() * v_inf;
        for (j, sigma_j) in sigma.iter().enumerate() {
            let (du, dv) = panel_velocity(i, (xi, yi), j, *sigma_j);
            ug += du;
            vg += dv;
        }
        let vt = ug * txi + vg * tyi;
        tangential.push(vt / v_inf);
        cp.push(1.0 - (vt / v_inf).powi(2));
    }
    Ok(SourceSolution {
        strengths: sigma,
        tangential_velocity: tangential,
        cp,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn circle_points(n: usize, radius: f64) -> Vec<Point2> {
        (0..n)
            .map(|i| {
                let t = 2.0 * core::f64::consts::PI * i as f64 / n as f64;
                (radius * t.cos(), radius * t.sin())
            })
            .collect()
    }

    #[test]
    fn cylinder_cp_matches_analytic() {
        // Circular cylinder, α = 0: Cp(θ) = 1 − 4 sin²θ on the surface.
        let pts = circle_points(64, 1.0);
        let sol = solve_source_panels(&pts, 0.0, 1.0).expect("solvable");
        for (i, &cp) in sol.cp.iter().enumerate() {
            let t = 2.0 * core::f64::consts::PI * (i as f64 + 0.5) / 64.0;
            let expected = 1.0 - 4.0 * t.sin().powi(2);
            assert!(
                (cp - expected).abs() < 0.08,
                "panel {i} (θ={t:.2}): {cp} vs {expected}"
            );
        }
        // Stagnation Cp = 1 at nose, min ≈ −3 at shoulders.
        assert!(sol.cp[0] > 0.9);
        let min_cp = sol.cp.iter().fold(f64::INFINITY, |m, &c| m.min(c));
        assert!((min_cp + 3.0).abs() < 0.3, "min cp {min_cp}");
    }

    #[test]
    fn dalembert_zero_drag() {
        // Closed body in potential flow: pressure-integrated drag ≈ 0.
        let pts = circle_points(64, 1.0);
        let sol = solve_source_panels(&pts, 0.0, 1.0).expect("solvable");
        let mut drag = 0.0;
        for (i, &cp) in sol.cp.iter().enumerate() {
            let (x1, y1) = pts[i];
            let (x2, y2) = pts[(i + 1) % 64];
            let dx = x2 - x1;
            let dy = y2 - y1;
            let l = (dx * dx + dy * dy).sqrt();
            // Outward normal × Cp × panel length, x-component.
            let (nx, _) = (dy / l, -dx / l);
            drag += cp * nx * l;
        }
        assert!(drag.abs() < 1e-3, "d'Alembert drag {drag}");
    }

    #[test]
    fn rejects_degenerate() {
        assert_eq!(
            solve_source_panels(&[(0.0, 0.0), (1.0, 0.0)], 0.0, 1.0).unwrap_err(),
            PanelError::InvalidGeometry
        );
        let dup = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
        assert_eq!(
            solve_source_panels(&dup, 0.0, 1.0).unwrap_err(),
            PanelError::InvalidGeometry
        );
    }
}
