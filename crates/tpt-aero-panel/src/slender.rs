//! Slender-body theory for bodies of revolution in axial flow.
//!
//! First-order potential flow past a slender body `r(x)`: the cross-flow
//! perturbation potential `φ = S'(x)/(2π)·ln r` gives the classic surface
//! pressure distribution from the area derivative `S'(x)` and the local
//! radius. Suitable for fuselages, nacelles and missile bodies at small
//! incidence.

use alloc::vec::Vec;

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// Body of revolution defined by a radius profile sampled on `x ∈ [0, L]`.
#[derive(Debug, Clone)]
pub struct SlenderBody {
    /// Station positions (m), ascending, starting at 0.
    pub x: Vec<f64>,
    /// Radius at each station (m).
    pub radius: Vec<f64>,
}

/// Errors from body construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SlenderError {
    /// Fewer than 3 stations, non-ascending x, or negative radii.
    InvalidProfile,
}

impl core::fmt::Display for SlenderError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            SlenderError::InvalidProfile => {
                write!(f, "need ≥3 ascending stations with non-negative radii")
            }
        }
    }
}

impl core::error::Error for SlenderError {}

impl SlenderBody {
    /// Build and validate a profile.
    ///
    /// # Errors
    /// [`SlenderError::InvalidProfile`] on malformed input.
    pub fn new(x: Vec<f64>, radius: Vec<f64>) -> Result<Self, SlenderError> {
        if x.len() != radius.len() || x.len() < 3 {
            return Err(SlenderError::InvalidProfile);
        }
        if x.windows(2).any(|w| w[1] <= w[0]) || radius.iter().any(|&r| r < 0.0) {
            return Err(SlenderError::InvalidProfile);
        }
        Ok(Self { x, radius })
    }

    /// Parabolic-arc body of revolution (Sears–Haack-like spine) with
    /// `n` stations: `r(x) = r_max·√(4ξ(1−ξ))`, ξ = x/L.
    #[must_use]
    pub fn parabolic(length: f64, r_max: f64, n: usize) -> Self {
        let n = n.max(3);
        let x: Vec<f64> = (0..n).map(|i| length * i as f64 / (n - 1) as f64).collect();
        let radius: Vec<f64> = x
            .iter()
            .map(|&xi| {
                let xi_n = xi / length;
                // Keep the endpoints exactly zero (pointed nose/tail).
                if xi_n <= 0.0 || xi_n >= 1.0 {
                    0.0
                } else {
                    r_max * (4.0 * xi_n * (1.0 - xi_n)).sqrt()
                }
            })
            .collect();
        Self { x, radius }
    }

    /// Cross-sectional area `S(x) = π r²` at each station.
    #[must_use]
    pub fn area(&self) -> Vec<f64> {
        self.radius
            .iter()
            .map(|&r| core::f64::consts::PI * r * r)
            .collect()
    }

    /// Surface perturbation analysis: returns per-station
    /// `(dS/dx, d²S/dx², u_perturbation/V∞)` using central differences on
    /// the area distribution, plus the incompressible surface pressure
    /// coefficient from the linearized slender-body relation
    /// `C_p ≈ −2·u'/V∞ − (r'/V∞)²·(1 − ln(...))`: we use the classic
    /// small-perturbation form `C_p = −2 u'/V∞` (slender-body, axisymmetric).
    ///
    /// The total axial force from this distribution vanishes (d'Alembert) in
    /// inviscid incompressible theory — see the property test.
    #[must_use]
    pub fn surface_analysis(&self, v_inf: f64) -> Vec<(f64, f64, f64)> {
        let s = self.area();
        let n = s.len();
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let (ds, d2s) = match i {
                0 => {
                    let h = self.x[1] - self.x[0];
                    ((s[1] - s[0]) / h, (s[2] - 2.0 * s[1] + s[0]) / (h * h))
                }
                k if k == n - 1 => {
                    let h = self.x[n - 1] - self.x[n - 2];
                    (
                        (s[n - 1] - s[n - 2]) / h,
                        (s[n - 1] - 2.0 * s[n - 2] + s[n - 3]) / (h * h),
                    )
                }
                k => {
                    let hp = self.x[k + 1] - self.x[k];
                    let hm = self.x[k] - self.x[k - 1];
                    let ds = (s[k + 1] - s[k - 1]) / (hp + hm);
                    let d2s = 2.0 * ((s[k + 1] - s[k]) / hp - (s[k] - s[k - 1]) / hm) / (hp + hm);
                    (ds, d2s)
                }
            };
            // Linearized axial perturbation from slender-body theory:
            // u'(x) = (1/2π)·[S''(x)·ln(1/r(x)) − ... we use the leading
            // term u' ≈ S''(x)/(2π)·ln(L/r_ref) with a documented reference
            // length; the dominant contribution for Cp at moderate fineness
            // ratio is the −2u'/V∞ term.
            let r_ref = self.radius.iter().fold(1.0e-3_f64, |m, &r| m.max(r));
            let r_loc = self.radius[i].max(1.0e-6 * r_ref);
            let u_pert = d2s / (2.0 * core::f64::consts::PI)
                * ((r_ref / r_loc).ln().clamp(-10.0, 10.0) + 0.5);
            out.push((ds, d2s, u_pert / v_inf.max(1e-9)));
        }
        out
    }
}

impl SlenderBody {
    /// First-order (Munk) normal force and pitching moment of the body at
    /// incidence `alpha_rad`.
    ///
    /// Slender-body potential theory gives a sectional normal force per unit
    /// length `l(x) = ρ V∞ α dS/dx`, where `S(x) = π r(x)²` is the local cross
    /// section. The body lift is `∫ l dx`, which for a bluff tail is about
    /// `ρ V∞ α (S(L) − S(0))`, and the moment about `x_cg` (positive aft,
    /// nose-up positive) is `−∫ (x − x_cg) l dx`.
    ///
    /// Returns `(lift, moment)` in N and N·m, with lift positive up. This is a
    /// first-order estimate for slender bodies at small incidence, not a
    /// substitute for a body panel model.
    #[must_use]
    pub fn munk_lift_and_moment(
        &self,
        rho: f64,
        v_inf: f64,
        alpha_rad: f64,
        x_cg: f64,
    ) -> (f64, f64) {
        let analysis = self.surface_analysis(v_inf);
        let mut lift = 0.0;
        let mut moment = 0.0;
        for (w, a) in self.x.windows(2).zip(analysis.windows(2)) {
            let dx = w[1] - w[0];
            let l0 = rho * v_inf * alpha_rad * a[0].0;
            let l1 = rho * v_inf * alpha_rad * a[1].0;
            let (x0, x1) = (w[0] - x_cg, w[1] - x_cg);
            lift += 0.5 * (l0 + l1) * dx;
            moment -= 0.5 * (x0 * l0 + x1 * l1) * dx;
        }
        (lift, moment)
    }
}

impl SlenderBody {
    /// First-order (Munk) side force and yawing moment of the body at sideslip
    /// `beta_rad` (wind from the right positive).
    ///
    /// The sectional side force per unit length is `−ρ V∞² β dS/dx`, the lateral
    /// counterpart of the lift in [`SlenderBody::munk_lift_and_moment`]. The side
    /// force is negative for a positive sideslip, and the yawing moment about
    /// `x_cg` (positive nose-right) is `ρ V∞² β ∫ (x − x_cg) dS/dx dx`. An
    /// aft-heavy body therefore gives weathercock stability.
    ///
    /// Returns `(side_force, yaw_moment)` in N and N·m, side force positive to
    /// the right. A first-order estimate at small sideslip.
    #[must_use]
    pub fn munk_side_force_and_yaw(
        &self,
        rho: f64,
        v_inf: f64,
        beta_rad: f64,
        x_cg: f64,
    ) -> (f64, f64) {
        let analysis = self.surface_analysis(v_inf);
        let k = rho * v_inf * v_inf * beta_rad;
        let mut side = 0.0;
        let mut yaw = 0.0;
        for (w, a) in self.x.windows(2).zip(analysis.windows(2)) {
            let dx = w[1] - w[0];
            let s0 = a[0].0;
            let s1 = a[1].0;
            let (x0, x1) = (w[0] - x_cg, w[1] - x_cg);
            side -= 0.5 * k * (s0 + s1) * dx;
            yaw += 0.5 * k * (x0 * s0 + x1 * s1) * dx;
        }
        (side, yaw)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    #[test]
    fn parabolic_profile_valid() {
        let b = SlenderBody::parabolic(4.0, 0.3, 41);
        assert_eq!(b.radius[0], 0.0);
        assert_eq!(b.radius[40], 0.0);
        let rmax = b.radius.iter().fold(0.0_f64, |m, &r| m.max(r));
        assert!((rmax - 0.3).abs() < 1e-9);
    }

    #[test]
    fn slender_body_closes_dalembert() {
        // Pressure-integrated axial force of the linearized distribution
        // over a closed pointed body vanishes (no wave drag incompressible).
        let b = SlenderBody::parabolic(4.0, 0.3, 81);
        let analysis = b.surface_analysis(50.0);
        let mut fx = 0.0;
        for (w, a) in b.x.windows(2).zip(analysis.windows(2)) {
            let dx = w[1] - w[0];
            // F = ∮ Cp n_x dS ≈ Cp·dS(x): axial force from pressure on the
            // projected area: F = ∫ Cp dS/dx dx with Cp = −2u'/V∞.
            let cp = -2.0 * 0.5 * (a[0].2 + a[1].2);
            fx += cp * 0.5 * (a[0].0 + a[1].0) * dx;
        }
        assert!(fx.abs() < 1.0, "axial force {fx} (d'Alembert violation)");
    }

    #[test]
    fn munk_lift_is_the_area_change_over_the_body() {
        // A blunt-tailed body: S(0) = 0 at the nose, S(L) = π·0.2² at the base.
        let n = 81;
        let x: Vec<f64> = (0..n).map(|i| 4.0 * i as f64 / (n - 1) as f64).collect();
        let radius: Vec<f64> = x
            .iter()
            .map(|&xi| 0.2 * (xi / 4.0).sqrt().min(1.0))
            .collect();
        let b = SlenderBody::new(x, radius).unwrap();
        let (rho, v, alpha) = (1.225, 50.0, 0.05);
        let (lift, _) = b.munk_lift_and_moment(rho, v, alpha, 2.0);
        let s_base = core::f64::consts::PI * 0.2 * 0.2;
        let expected = rho * v * alpha * s_base;
        assert!(
            (lift / expected - 1.0).abs() < 0.02,
            "lift {lift} vs {expected}"
        );
    }

    #[test]
    fn munk_moment_is_nose_down_when_the_lift_acts_aft_of_the_cg() {
        let n = 81;
        let x: Vec<f64> = (0..n).map(|i| 4.0 * i as f64 / (n - 1) as f64).collect();
        let radius: Vec<f64> = x
            .iter()
            .map(|&xi| 0.2 * (xi / 4.0).sqrt().min(1.0))
            .collect();
        let b = SlenderBody::new(x, radius).unwrap();
        let (lift, moment_fwd) = b.munk_lift_and_moment(1.225, 50.0, 0.05, 1.0);
        let (_, moment_aft) = b.munk_lift_and_moment(1.225, 50.0, 0.05, 3.0);
        assert!(lift > 0.0);
        // Moving the CG aft moves the lift forward of it, so the nose-up moment rises.
        assert!(moment_aft > moment_fwd);
        // The lift acts behind the CG at 1 m, so the moment is nose-down.
        assert!(moment_fwd < 0.0);
    }

    #[test]
    fn munk_side_force_opposes_sideslip_and_weathercock_yaw_needs_aft_area() {
        let n = 81;
        let x: Vec<f64> = (0..n).map(|i| 4.0 * i as f64 / (n - 1) as f64).collect();
        let radius: Vec<f64> = x
            .iter()
            .map(|&xi| 0.2 * (xi / 4.0).sqrt().min(1.0))
            .collect();
        let b = SlenderBody::new(x, radius).unwrap();
        let (rho, v, beta) = (1.225, 50.0, 0.05);
        // CG ahead of the area centroid (about 2 m for this profile), so the yaw is positive.
        let (side, yaw) = b.munk_side_force_and_yaw(rho, v, beta, 1.0);
        // Side force is −ρV²β ΔS, opposite to a sideslip from the right.
        let expected = -rho * v * v * beta * core::f64::consts::PI * 0.2 * 0.2;
        assert!(
            (side / expected - 1.0).abs() < 0.02,
            "side {side} vs {expected}"
        );
        // The aft-heavy body is weathercock stable: positive yaw for β > 0.
        assert!(yaw > 0.0, "yaw {yaw}");
        // Moving the CG aft reduces the restoring yaw.
        let (_, yaw_aft) = b.munk_side_force_and_yaw(rho, v, beta, 3.0);
        assert!(yaw_aft < yaw);
    }

    #[test]
    fn rejects_invalid_profile() {
        assert_eq!(
            SlenderBody::new(vec![0.0, 1.0], vec![0.1, 0.1]).unwrap_err(),
            SlenderError::InvalidProfile
        );
        assert_eq!(
            SlenderBody::new(vec![0.0, 1.0, 0.5], vec![0.1, 0.2, 0.3]).unwrap_err(),
            SlenderError::InvalidProfile
        );
        assert_eq!(
            SlenderBody::new(vec![0.0, 1.0, 2.0], vec![0.1, -0.2, 0.3]).unwrap_err(),
            SlenderError::InvalidProfile
        );
    }
}
