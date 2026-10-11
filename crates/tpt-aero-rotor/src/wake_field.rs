//! Full induced-velocity field of a uniformly loaded actuator disk: the
//! linearized vortex-cylinder wake.
//!
//! For uniform loading the bound circulation is constant along the disk, so
//! the only trailing vorticity is shed at the rim. The wake is therefore a
//! single cylinder of radius `R` carrying azimuthal vorticity `γ = 2 v_d` per
//! unit length, from the disk downstream (Heyson & Katzoff; Jewel & Heyson,
//! NASA Memo 4-15-59L, which describes the uniformly loaded rotor's wake as "a
//! single vortex cylinder"). A vortex sheet of strength `γ` jumps the axial
//! velocity by `γ`, from `0` outside the wake to `2 v_d` inside far downstream.
//!
//! Each element of the cylinder is a ring of radius `R` and circulation `γ dz'`.
//! The induced velocity is the integral of the exact ring kernel over the
//! cylinder:
//!
//! ```text
//! u(z, r) = 2 v_d ∫₀^∞ ring(R, z'; z, r) dz'
//! ```
//!
//! On the axis this integral is exact, `u_z = v_d (1 + z/√(z² + R²))`, so the
//! model reproduces the classical on-axis result and extends it off the axis.
//!
//! The model is linear: it does not include the wake contraction of momentum
//! theory, which is what [`crate::Slipstream`] adds on the axis. The ring
//! kernel is evaluated in closed form with complete elliptic integrals; the
//! `z'` integral is by composite Simpson after the substitution `z' = R tan t`.

use alloc::vec;

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

use crate::blade::RotorError;

/// Largest Simpson panel count (in `t`) the wake integral may use. Must be even.
const WAKE_PANELS: usize = 4096;
/// Panel count at which the adaptive wake integral starts.
const WAKE_PANELS_START: usize = 64;
/// Relative agreement between successive panel counts that ends the integral.
const WAKE_RTOL: f64 = 1.0e-9;

/// Complete elliptic integrals `(K(k), E(k))` for modulus `k` (`0 ≤ k < 1`),
/// by the arithmetic-geometric mean. `kp = √(1 − k²)` is passed separately so
/// that it stays accurate as `k → 1`.
pub(crate) fn elliptic_ke(k: f64, kp: f64) -> (f64, f64) {
    let (mut a, mut b, mut c) = (1.0_f64, kp, k);
    // E = K (1 − Σ 2^{n−1} c_n²), with c_0 = k.
    let mut sum = 0.5 * c * c;
    let mut power = 1.0;
    for _ in 0..64 {
        let an = 0.5 * (a + b);
        let bn = (a * b).sqrt();
        c = 0.5 * (a - b);
        a = an;
        b = bn;
        sum += power * c * c;
        power *= 2.0;
        if c.abs() < 1.0e-17 {
            break;
        }
    }
    let kk = core::f64::consts::FRAC_PI_2 / a;
    (kk, kk * (1.0 - sum))
}

/// Induced velocity `(u_z, u_r)` of a unit-circulation vortex ring of radius
/// `a` centred at axial station `z0`, at the point `(z, r)`. Positive
/// circulation gives `u_z > 0` on the axis. Singular on the ring itself.
pub(crate) fn ring_velocity(a: f64, z0: f64, z: f64, r: f64) -> (f64, f64) {
    let zeta = z - z0;
    let d = (a + r) * (a + r) + zeta * zeta;
    let q = (a - r) * (a - r) + zeta * zeta;
    // The point lies on this ring's filament (only reachable when a quadrature
    // node hits the field point exactly, e.g. r = R at the disk plane). The
    // filament's self-induced velocity is singular there, so this node
    // contributes nothing; the neighbouring nodes carry the integral.
    if q <= f64::MIN_POSITIVE {
        return (0.0, 0.0);
    }
    let k = (4.0 * a * r / d).sqrt();
    let kp = (q / d).sqrt();
    let (kk, ee) = elliptic_ke(k, kp);
    let sqrt_d = d.sqrt();
    let two_pi = 2.0 * core::f64::consts::PI;
    let u_z = (kk + (a * a - r * r - zeta * zeta) / q * ee) / (two_pi * sqrt_d);
    let u_r = if r.abs() < 1.0e-12 {
        0.0
    } else {
        zeta * (-kk + (a * a + r * r + zeta * zeta) / q * ee) / (two_pi * r * sqrt_d)
    };
    (u_z, u_r)
}

/// Linearized vortex-cylinder wake of a uniformly loaded actuator disk.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CylinderWake {
    radius: f64,
    v_inf: f64,
    v_d: f64,
    density: f64,
}

impl CylinderWake {
    /// Field of a disk of radius `radius` (m) with induced velocity `v_d`
    /// (m/s) in freestream `v_inf` (m/s) at air density `density` (kg/m³).
    ///
    /// # Errors
    /// [`RotorError::InvalidGeometry`] if `radius` or `density` is not
    /// positive, `v_inf` is negative, or any input is non-finite.
    pub fn new(radius: f64, v_d: f64, v_inf: f64, density: f64) -> Result<Self, RotorError> {
        let finite = [radius, v_d, v_inf, density].iter().all(|v| v.is_finite());
        if !finite || radius <= 0.0 || density <= 0.0 || v_inf < 0.0 {
            return Err(RotorError::InvalidGeometry);
        }
        Ok(Self {
            radius,
            v_inf,
            v_d,
            density,
        })
    }

    /// Disk radius (m).
    #[must_use]
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Induced velocity at the disk (m/s).
    #[must_use]
    pub fn disk_induced_velocity(&self) -> f64 {
        self.v_d
    }

    /// Induced velocity `(u_z, u_r)` at axial station `z` (downstream of the
    /// disk) and radius `r` (m/s).
    ///
    /// The integrand is smooth except at the wake itself (`r = R`), where the
    /// sheet's velocity jump makes the quadrature less accurate near the surface.
    #[must_use]
    pub fn induced(&self, z: f64, r: f64) -> (f64, f64) {
        // Axisymmetric: evaluate on the positive side and flip u_r for r < 0.
        let sign = if r < 0.0 { -1.0 } else { 1.0 };
        let r = r.abs();
        let gamma = 2.0 * self.v_d;
        // Double the Simpson panel count until successive estimates agree.
        // Smooth field points stop after a few hundred ring evaluations rather
        // than 4001; points near the wake edge refine up to the cap.
        let mut panels = WAKE_PANELS_START;
        let (mut acc_z, mut acc_r) = self.simpson_sums(z, r, panels);
        while panels < WAKE_PANELS {
            let next = (panels * 2).min(WAKE_PANELS);
            let (nz, nr) = self.simpson_sums(z, r, next);
            let scale = nz.abs().max(nr.abs()).max(1.0e-12);
            let converged = (nz - acc_z).abs().max((nr - acc_r).abs()) <= WAKE_RTOL * scale;
            (acc_z, acc_r) = (nz, nr);
            panels = next;
            if converged {
                break;
            }
        }
        let h = core::f64::consts::FRAC_PI_2 / panels as f64;
        let scale = gamma * h / 3.0;
        (scale * acc_z, sign * scale * acc_r)
    }

    /// Unscaled Simpson sums of the ring kernel over `panels` panels in `t`,
    /// with `z' = R tan t` and `dz' = R sec² t dt`. The induced velocity is
    /// `2 v_d · h / 3` times these sums, where `h = (π/2)/panels`.
    fn simpson_sums(&self, z: f64, r: f64, panels: usize) -> (f64, f64) {
        let rad = self.radius;
        let half = core::f64::consts::FRAC_PI_2;
        let h = half / panels as f64;
        let mut acc_z = 0.0;
        let mut acc_r = 0.0;
        for i in 0..=panels {
            let t = h * i as f64;
            // Stop just short of the endpoint, where tan t diverges.
            let t = t.min(half * (1.0 - 1.0e-12));
            let zp = rad * t.tan();
            let jac = rad / (t.cos() * t.cos());
            let w = if i == 0 || i == panels {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            };
            let (uz, ur) = ring_velocity(rad, zp, z, r);
            acc_z += w * uz * jac;
            acc_r += w * ur * jac;
        }
        (acc_z, acc_r)
    }

    /// Total flow at `(z, r)`: axial speed `V + u_z` and radial speed `u_r`
    /// (m/s).
    #[must_use]
    pub fn local_flow(&self, z: f64, r: f64) -> (f64, f64) {
        let (uz, ur) = self.induced(z, r);
        (self.v_inf + uz, ur)
    }

    /// Dynamic-pressure increment `½ρ[(V + u_z)² + u_r² − V²]` at `(z, r)` (Pa).
    #[must_use]
    pub fn dynamic_pressure_increment(&self, z: f64, r: f64) -> f64 {
        let (u, w) = self.local_flow(z, r);
        0.5 * self.density * (u * u + w * w - self.v_inf * self.v_inf)
    }
}

/// Induced axial velocities of the cylinder wake sampled along the disk plane,
/// for checks and plotting. Returns `(r, u_z)` at `n` evenly spaced radii on
/// `[0, R]`.
#[must_use]
pub fn disk_plane_profile(wake: &CylinderWake, n: usize) -> alloc::vec::Vec<(f64, f64)> {
    let mut out = vec![];
    for i in 0..n {
        let r = wake.radius * i as f64 / (n.max(2) - 1) as f64;
        // Sit a hair above the disk so the sheet's own jump is not sampled.
        let (uz, _) = wake.induced(1.0e-9 * wake.radius, r);
        out.push((r, uz));
    }
    out
}

/// Choice of slipstream model: the momentum-theory contracting tube, or the
/// linearized vortex-cylinder field. Both give the same on-axis velocity; they
/// differ in the wake radius (the tube contracts, the cylinder does not) and
/// the cylinder adds the crossflow.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum SlipstreamModel {
    /// Momentum-theory contracting stream tube with a blended edge.
    Tube(crate::slipstream::Slipstream),
    /// Linearized vortex-cylinder wake with the full induced field.
    Cylinder(CylinderWake),
}

impl SlipstreamModel {
    /// Induced velocity at the disk (m/s).
    #[must_use]
    pub fn disk_induced_velocity(&self) -> f64 {
        match self {
            Self::Tube(t) => t.disk_induced_velocity(),
            Self::Cylinder(c) => c.disk_induced_velocity(),
        }
    }

    /// Total flow `(axial, radial)` at station `z` and radius `r` (m/s). The
    /// tube has no crossflow, so its radial component is zero.
    #[must_use]
    pub fn local_flow(&self, z: f64, r: f64) -> (f64, f64) {
        match self {
            Self::Tube(t) => (t.freestream() + t.velocity_increment(z, r), 0.0),
            Self::Cylinder(c) => c.local_flow(z, r),
        }
    }

    /// Dynamic-pressure increment at `(z, r)` (Pa).
    #[must_use]
    pub fn dynamic_pressure_increment(&self, z: f64, r: f64) -> f64 {
        match self {
            Self::Tube(t) => t.dynamic_pressure_increment(z, r),
            Self::Cylinder(c) => c.dynamic_pressure_increment(z, r),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn exact_on_axis(v_d: f64, radius: f64, z: f64) -> f64 {
        v_d * (1.0 + z / (z * z + radius * radius).sqrt())
    }

    #[test]
    fn on_axis_matches_closed_form_to_quadrature_accuracy() {
        let w = CylinderWake::new(1.0, 0.75, 25.0, 1.225).expect("valid");
        for &z in &[0.1, 0.5, 1.0, 2.0, 4.0, 10.0] {
            let (uz, ur) = w.induced(z, 0.0);
            let want = exact_on_axis(0.75, 1.0, z);
            assert!((uz - want).abs() < 1e-9 * want, "z = {z}: {uz} vs {want}");
            assert!(ur.abs() < 1e-12);
        }
    }

    #[test]
    fn disk_plane_induced_velocity_is_uniform_inside_the_disk() {
        let w = CylinderWake::new(1.0, 0.75, 25.0, 1.225).expect("valid");
        for (r, uz) in disk_plane_profile(&w, 10) {
            if r < 0.95 {
                assert!((uz - 0.75).abs() < 1e-6 * 0.75, "r = {r}: {uz}");
            }
        }
    }

    #[test]
    fn flow_is_divergence_free_off_the_wake() {
        // ∂u_z/∂z + (1/r) ∂(r u_r)/∂r = 0 for an incompressible potential field.
        let w = CylinderWake::new(1.0, 0.75, 25.0, 1.225).expect("valid");
        let (z, r, h) = (0.6, 0.4, 1.0e-4);
        let duz_dz = (w.induced(z + h, r).0 - w.induced(z - h, r).0) / (2.0 * h);
        let rur = |rr: f64| rr * w.induced(z, rr).1;
        let div_r = (rur(r + h) - rur(r - h)) / (2.0 * h) / r;
        assert!(
            (duz_dz + div_r).abs() < 1e-4,
            "divergence = {}",
            duz_dz + div_r
        );
    }

    #[test]
    fn off_axis_field_converges_with_wake_resolution() {
        let w = CylinderWake::new(1.0, 0.75, 25.0, 1.225).expect("valid");
        let (zp_, rp_) = (0.5, 0.5);
        let fine = {
            let rad = w.radius;
            let gamma = 2.0 * w.v_d;
            let half = core::f64::consts::FRAC_PI_2;
            let n = 2 * WAKE_PANELS;
            let h = half / n as f64;
            let (mut az, mut ar) = (0.0, 0.0);
            for i in 0..=n {
                let t = (h * i as f64).min(half * (1.0 - 1.0e-12));
                let zp = rad * t.tan();
                let jac = rad / (t.cos() * t.cos());
                let wgt = if i == 0 || i == n {
                    1.0
                } else if i % 2 == 1 {
                    4.0
                } else {
                    2.0
                };
                let (uz, ur) = ring_velocity(rad, zp, zp_, rp_);
                az += wgt * uz * jac;
                ar += wgt * ur * jac;
            }
            (gamma * h / 3.0 * az, gamma * h / 3.0 * ar)
        };
        let (uz, ur) = w.induced(zp_, rp_);
        assert!((uz - fine.0).abs() < 1e-6, "u_z {uz} vs fine {}", fine.0);
        assert!((ur - fine.1).abs() < 1e-6, "u_r {ur} vs fine {}", fine.1);
    }

    #[test]
    fn far_wake_inside_the_cylinder_approaches_twice_the_disk_value() {
        let w = CylinderWake::new(1.0, 0.75, 25.0, 1.225).expect("valid");
        let (uz, _) = w.induced(10.0, 0.0);
        assert!((uz - 2.0 * 0.75).abs() / 1.5 < 0.01, "u_z(10R) = {uz}");
    }

    #[test]
    fn wake_edge_region_converges_with_resolution() {
        // Points just inside and outside the wake, where the sheet's jump makes
        // the integrand sharpest. Compare with a 2× finer quadrature.
        let w = CylinderWake::new(1.0, 0.75, 25.0, 1.225).expect("valid");
        for &(z, r) in &[(4.15, 0.9), (4.15, 1.1)] {
            let fine = {
                let rad = w.radius;
                let half = core::f64::consts::FRAC_PI_2;
                let n = 2 * WAKE_PANELS;
                let h = half / n as f64;
                let (mut az, mut ar) = (0.0, 0.0);
                for i in 0..=n {
                    let t = (h * i as f64).min(half * (1.0 - 1.0e-12));
                    let zp = rad * t.tan();
                    let jac = rad / (t.cos() * t.cos());
                    let wgt = if i == 0 || i == n {
                        1.0
                    } else if i % 2 == 1 {
                        4.0
                    } else {
                        2.0
                    };
                    let (uz, ur) = ring_velocity(rad, zp, z, r);
                    az += wgt * uz * jac;
                    ar += wgt * ur * jac;
                }
                (2.0 * w.v_d * h / 3.0 * az, 2.0 * w.v_d * h / 3.0 * ar)
            };
            let (uz, ur) = w.induced(z, r);
            assert!(
                (uz - fine.0).abs() < 1e-4,
                "u_z at r = {r}: {uz} vs {}",
                fine.0
            );
            assert!(
                (ur - fine.1).abs() < 1e-4,
                "u_r at r = {r}: {ur} vs {}",
                fine.1
            );
        }
    }

    #[test]
    fn field_is_axisymmetric_about_the_axis() {
        // u_z is even in r, u_r is odd: the field on the −y side of the tail.
        let w = CylinderWake::new(1.0, 0.75, 25.0, 1.225).expect("valid");
        for &(z, r) in &[(4.15, 0.5), (4.15, 2.0), (1.0, 1.4)] {
            let (uz_p, ur_p) = w.induced(z, r);
            let (uz_m, ur_m) = w.induced(z, -r);
            assert!((uz_p - uz_m).abs() < 1e-12, "u_z at ({z}, ±{r})");
            assert!((ur_p + ur_m).abs() < 1e-12, "u_r at ({z}, ±{r})");
        }
    }

    #[test]
    fn rejects_invalid_inputs() {
        assert!(CylinderWake::new(0.0, 1.0, 0.0, 1.2).is_err());
        assert!(CylinderWake::new(1.0, 1.0, -1.0, 1.2).is_err());
        assert!(CylinderWake::new(1.0, f64::NAN, 0.0, 1.2).is_err());
    }
}

#[cfg(test)]
mod wall_tests {
    use super::*;

    #[test]
    fn field_on_the_wake_edge_is_finite() {
        // r = R at the disk plane used to hit a quadrature node exactly on the
        // ring and return NaN.
        let w = CylinderWake::new(1.0, 5.0, 10.0, 1.225).expect("valid");
        for z in [0.0, 0.5, 2.0] {
            let (uz, ur) = w.induced(z, 1.0);
            assert!(uz.is_finite() && ur.is_finite(), "z = {z}: ({uz}, {ur})");
            assert!(w.dynamic_pressure_increment(z, 1.0).is_finite());
        }
    }

    #[test]
    fn adaptive_integral_matches_the_fixed_fine_grid() {
        // The adaptive count must agree with the 4096-panel reference it caps at,
        // off axis and close to the wake edge where refinement is needed.
        let w = CylinderWake::new(1.0, 0.8, 0.0, 1.225).expect("valid");
        let h = core::f64::consts::FRAC_PI_2 / WAKE_PANELS as f64;
        let scale = 2.0 * 0.8 * h / 3.0;
        for &(z, r) in &[(0.2, 0.3), (0.5, 0.9), (1.0, 0.99), (2.0, 1.4), (0.05, 0.7)] {
            let (rz, rr) = w.simpson_sums(z, r, WAKE_PANELS);
            let (uz, ur) = w.induced(z, r);
            let (ref_z, ref_r) = (scale * rz, scale * rr);
            let tol = 1e-7 * ref_z.abs().max(ref_r.abs()).max(1e-9);
            assert!((uz - ref_z).abs() < tol, "({z}, {r}) u_z {uz} vs {ref_z}");
            assert!((ur - ref_r).abs() < tol, "({z}, {r}) u_r {ur} vs {ref_r}");
        }
    }
}
