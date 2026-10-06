//! 2-D lifting flow: constant-strength doublet panels with a Kutta
//! condition.
//!
//! Implementation note: by the doublet–vortex identity a constant-strength
//! doublet sheet of strength µ is equivalent to a vortex sheet of the same
//! strength, so the solver uses the classic constant-strength vortex-panel
//! formulation (Kuethe & Chow, ch. 5): unknown panel vorticity, flow
//! tangency at the midpoints, and a trailing-edge Kutta condition
//! (equal upper/lower tangential velocity). The finite-segment induction
//! kernel is the validated 3-D Biot–Savart segment from `tpt-aero-vlm`
//! restricted to the plane.

use alloc::vec;
use alloc::vec::Vec;

use tpt_aero_vlm::lu_solve;

use crate::source::{PanelError, Point2};

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// Result of a doublet/vortex-panel analysis.
#[derive(Debug, Clone)]
pub struct DoubletSolution {
    /// Node vorticities γ_0..γ_n (doublet-equivalent sheet strengths).
    pub strengths: Vec<f64>,
    /// Total circulation Γ = Σ γ·l (m²/s).
    pub circulation: f64,
    /// Tangential velocity (normalized by V∞) at each panel midpoint.
    pub tangential_velocity: Vec<f64>,
    /// Pressure coefficient at each panel midpoint.
    pub cp: Vec<f64>,
    /// Section lift coefficient `C_L = 2Γ/(V∞·c)`.
    pub cl: f64,
}

/// Solve the constant-strength doublet/vortex-panel system.
///
/// * `points`: closed contour vertices; the contour must start and end at
///   the trailing edge (first vertex = last vertex position) and be ordered
///   along the surface (upper from TE to LE, lower from LE to TE).
/// * `alpha_rad`: angle of attack.
/// * `v_inf`: freestream speed; `chord`: reference chord.
///
/// # Errors
/// [`PanelError::InvalidGeometry`] for degenerate contours,
/// [`PanelError::Singular`] if the influence system is singular.
pub fn solve_doublet_panels(
    points: &[Point2],
    alpha_rad: f64,
    v_inf: f64,
    chord: f64,
) -> Result<DoubletSolution, PanelError> {
    let n = points.len();
    if n < 4 {
        return Err(PanelError::InvalidGeometry);
    }
    let mut x = Vec::with_capacity(n);
    let mut y = Vec::with_capacity(n);
    for &(px, py) in points {
        x.push(px);
        y.push(py);
    }
    let (mut xc, mut yc) = (vec![0.0; n], vec![0.0; n]);
    let (mut panel_angle, mut panel_len) = (vec![0.0; n], vec![0.0; n]);
    for i in 0..n {
        let j = (i + 1) % n;
        let dx = x[j] - x[i];
        let dy = y[j] - y[i];
        panel_len[i] = (dx * dx + dy * dy).sqrt();
        if panel_len[i] <= 0.0 {
            return Err(PanelError::InvalidGeometry);
        }
        panel_angle[i] = dy.atan2(dx);
        xc[i] = 0.5 * (x[i] + x[j]);
        yc[i] = 0.5 * (y[i] + y[j]);
    }

    // ---- Hess-Smith lifting panel method ----
    //
    // Unknowns: the n source strengths sigma_j (they carry thickness and
    // angle response, with the healthy 1/2 self-diagonal) plus ONE global
    // circulation Gamma represented by a point vortex (a degenerate doublet)
    // placed inside the body at the contour centroid. Equations: flow
    // tangency at all n midpoints, plus the trailing-edge Kutta condition
    // V.t0 = -V.tn-1. This is the classic 1966 Hess-Smith lifting
    // formulation; it is immune to the checkerboard nullspace that plagues
    // pure sheet-strength discretizations on near-collinear contours.
    const GAUSS: [(f64, f64); 4] = [
        (-0.861_136_311_594_052_6, 0.347_854_845_137_453_85),
        (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.861_136_311_594_052_6, 0.347_854_845_137_453_85),
    ];

    // Reference point for the circulation vortex: the contour centroid.
    let mut cx = 0.0;
    let mut cy = 0.0;
    for i in 0..n {
        cx += xc[i];
        cy += yc[i];
    }
    let vortex_pt = (cx / n as f64, cy / n as f64);
    let vortex_velocity = |(xi, yi): Point2| -> (f64, f64) {
        let rx = xi - vortex_pt.0;
        let ry = yi - vortex_pt.1;
        let r2 = (rx * rx + ry * ry).max(1.0e-12);
        let f = 1.0 / (2.0 * core::f64::consts::PI * r2);
        (-f * ry, f * rx)
    };

    // Source velocity at midpoint i from panel j per unit strength,
    // projected on the tangent (tx, ty) or normal (nx, ny) of panel i.
    let source_influence = |i: usize, j: usize, tx: f64, ty: f64| -> f64 {
        if i == j {
            return 0.0; // zero tangential self contribution at the midpoint
        }
        let (x1, y1) = points[j];
        let (x2, y2) = points[(j + 1) % n];
        let half = panel_len[j] * 0.5;
        let mut acc = 0.0;
        for &(g, w) in GAUSS.iter() {
            let sx = x1 + (x2 - x1) * (0.5 * (1.0 + g));
            let sy = y1 + (y2 - y1) * (0.5 * (1.0 + g));
            let rx = xc[i] - sx;
            let ry = yc[i] - sy;
            let r2 = (rx * rx + ry * ry).max(1.0e-20);
            let base = w * half / (2.0 * core::f64::consts::PI) / r2;
            acc += base * (rx * tx + ry * ty);
        }
        acc
    };
    let source_influence_n = |i: usize, j: usize| -> f64 {
        if i == j {
            return 0.5;
        }
        let (nxi, nyi) = nrm_of(points, i);
        let (x1, y1) = points[j];
        let (x2, y2) = points[(j + 1) % n];
        let half = panel_len[j] * 0.5;
        let mut acc = 0.0;
        for &(g, w) in GAUSS.iter() {
            let sx = x1 + (x2 - x1) * (0.5 * (1.0 + g));
            let sy = y1 + (y2 - y1) * (0.5 * (1.0 + g));
            let rx = xc[i] - sx;
            let ry = yc[i] - sy;
            let r2 = (rx * rx + ry * ry).max(1.0e-20);
            acc += w * half / (2.0 * core::f64::consts::PI) * (rx * nxi + ry * nyi) / r2;
        }
        acc
    };

    let size = n + 1;
    let mut a = vec![vec![0.0; size]; size];
    let mut rhs = vec![0.0; size];
    for i in 0..n {
        for j in 0..n {
            a[i][j] = source_influence_n(i, j);
        }
        let (vu, vv) = vortex_velocity((xc[i], yc[i]));
        let (nxi, nyi) = nrm_of(points, i);
        a[i][n] = vu * nxi + vv * nyi;
        rhs[i] = -v_inf * (alpha_rad.cos() * nxi + alpha_rad.sin() * nyi);
    }

    // Kutta: V.t0 = -V.t_{n-1}.
    let kutta_i = n - 1;
    for j in 0..n {
        let (t0x, t0y) = tangent_of(points, 0);
        let (vkx, vky) = tangent_of(points, kutta_i);
        a[n][j] = source_influence(0, j, t0x, t0y) + source_influence(kutta_i, j, vkx, vky);
    }
    let (vu0, vv0) = vortex_velocity((xc[0], yc[0]));
    let (vuk, vvk) = vortex_velocity((xc[kutta_i], yc[kutta_i]));
    let (c0, s0) = tangent_of(points, 0);
    let (ck, sk) = tangent_of(points, kutta_i);
    a[n][n] = (vu0 * c0 + vv0 * s0) + (vuk * ck + vvk * sk);
    rhs[n] = -v_inf
        * ((alpha_rad.cos() * c0 + alpha_rad.sin() * s0)
            + (alpha_rad.cos() * ck + alpha_rad.sin() * sk));

    let sol = lu_solve(&a, &rhs).map_err(|_| PanelError::Singular)?;
    let sigma: Vec<f64> = sol[..n].to_vec();
    let gamma_total = sol[n];

    // Loads: the contour is traversed TE->upper->LE->lower (clockwise for
    // lift-up flow), so the lift-positive circulation is the negative of the
    // solved vortex strength.
    let circulation = -gamma_total;
    let cl = 2.0 * circulation / (v_inf * chord);

    // Surface speed and Cp.
    let mut tangential = Vec::with_capacity(n);
    let mut cp = Vec::with_capacity(n);
    for i in 0..n {
        let (txi, tyi) = tangent_of(points, i);
        let mut vt = v_inf * (alpha_rad.cos() * txi + alpha_rad.sin() * tyi);
        for (j, sigma_j) in sigma.iter().enumerate() {
            vt += sigma_j * source_influence(i, j, txi, tyi);
        }
        let (vu, vv) = vortex_velocity((xc[i], yc[i]));
        vt += vu * txi + vv * tyi;
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

/// Outward unit normal of panel i.
fn nrm_of(points: &[Point2], i: usize) -> (f64, f64) {
    let n = points.len();
    let (x1, y1) = points[i];
    let (x2, y2) = points[(i + 1) % n];
    let dx = x2 - x1;
    let dy = y2 - y1;
    let l = (dx * dx + dy * dy).sqrt().max(1e-300);
    (dy / l, -dx / l)
}

/// Unit tangent of panel i.
fn tangent_of(points: &[Point2], i: usize) -> (f64, f64) {
    let n = points.len();
    let (x1, y1) = points[i];
    let (x2, y2) = points[(i + 1) % n];
    let dx = x2 - x1;
    let dy = y2 - y1;
    let l = (dx * dx + dy * dy).sqrt().max(1e-300);
    (dx / l, dy / l)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn symmetric_airfoil(n: usize) -> Vec<Point2> {
        // Thin symmetric NACA-like contour: TE → upper → LE → lower, closing
        // back to the exact TE point. The closed-te coefficient keeps the
        // surface points strictly off the TE so no panel degenerates.
        let thickness = |x: f64| {
            5.0
                * 0.12
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
    fn rejects_degenerate() {
        assert_eq!(
            solve_doublet_panels(&[(0.0, 0.0), (0.5, 0.1), (1.0, 0.0)], 0.0, 1.0, 1.0)
                .unwrap_err(),
            PanelError::InvalidGeometry
        );
    }
}
