//! Propeller slipstream: a momentum-theory stream tube behind an actuator disk.
//!
//! The disk carries thrust `T` at radius `R` in a freestream `V`. Momentum
//! theory gives the induced velocity at the disk, `v_d`, from
//!
//! ```text
//! T = 2 ρ A v_d (V + v_d),   A = π R²,
//! ```
//!
//! and the far wake moves at `V + 2 v_d`. On the axis, a uniformly loaded disk
//! has the classical closed form
//!
//! ```text
//! Δu(x) = v_d (1 + x / √(x² + R²)),
//! ```
//!
//! with `x` measured downstream from the disk. This is `v_d` at the disk, `2 v_d`
//! far downstream, and tends to zero far upstream. Mass conservation contracts
//! the tube: `R_w(x) = R √((V + v_d) / (V + Δu(x)))`.
//!
//! Off the axis the model is a **tube**, not a field. Inside the tube the axial
//! increment is taken equal to its on-axis value, outside it is zero, and the
//! edge is blended linearly across ±10 % of the tube radius. This is the
//! standard preliminary-design treatment of a tail in a slipstream. It has no
//! swirl, no induced crossflow, and no off-axis velocity profile, so it is a
//! first-order estimate of the tail's dynamic-pressure change, not a
//! replacement for a full induced-velocity field.

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

use crate::blade::RotorError;

/// Half-width of the tube edge blend, as a fraction of the tube radius.
const EDGE_BLEND: f64 = 0.1;

/// A momentum-theory slipstream behind a propeller or rotor disk.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Slipstream {
    /// Disk radius `R` (m).
    radius: f64,
    /// Freestream speed `V` (m/s).
    v_inf: f64,
    /// Air density `ρ` (kg/m³).
    density: f64,
    /// Induced velocity at the disk `v_d` (m/s).
    induced: f64,
}

impl Slipstream {
    /// Build the slipstream for a disk of radius `radius` (m) producing thrust
    /// `thrust` (N) in freestream `v_inf` (m/s) at air density `density`
    /// (kg/m³).
    ///
    /// # Errors
    /// [`RotorError::InvalidGeometry`] if `radius` or `density` is not positive,
    /// if `v_inf` is negative, or if `thrust` is negative (windmilling is not
    /// modelled). Non-finite inputs are rejected too.
    pub fn new(radius: f64, thrust: f64, v_inf: f64, density: f64) -> Result<Self, RotorError> {
        let finite = [radius, thrust, v_inf, density]
            .iter()
            .all(|v| v.is_finite());
        if !finite || radius <= 0.0 || density <= 0.0 || v_inf < 0.0 || thrust < 0.0 {
            return Err(RotorError::InvalidGeometry);
        }
        let area = core::f64::consts::PI * radius * radius;
        // Positive root of v² + V v − T/(2ρA) = 0.
        let disk_loading = thrust / (2.0 * density * area);
        let half_v = 0.5 * v_inf;
        let induced = (half_v * half_v + disk_loading).sqrt() - half_v;
        Ok(Self {
            radius,
            v_inf,
            density,
            induced,
        })
    }

    /// Disk radius `R` (m).
    #[must_use]
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Freestream speed `V` (m/s).
    #[must_use]
    pub fn freestream(&self) -> f64 {
        self.v_inf
    }

    /// Induced velocity at the disk, `v_d` (m/s).
    #[must_use]
    pub fn disk_induced_velocity(&self) -> f64 {
        self.induced
    }

    /// Axial velocity increment on the axis at distance `x` downstream of the
    /// disk (m/s). Negative `x` is upstream, where the increment is smaller.
    #[must_use]
    pub fn axial_increment(&self, x: f64) -> f64 {
        let r = self.radius;
        self.induced * (1.0 + x / (x * x + r * r).sqrt())
    }

    /// Tube radius at axial station `x` (m), from mass conservation.
    ///
    /// With no thrust and no freestream there is no contraction to speak of,
    /// so the tube keeps the disk radius rather than dividing by zero.
    #[must_use]
    pub fn tube_radius(&self, x: f64) -> f64 {
        let upstream = self.v_inf + self.induced;
        let total = self.v_inf + self.axial_increment(x);
        if upstream <= 0.0 || total <= 0.0 {
            return self.radius;
        }
        self.radius * (upstream / total).sqrt()
    }

    /// Axial velocity increment at station `x` and radial distance `r` from the
    /// axis (m/s). Inside the tube it equals the on-axis value, outside it is
    /// zero, and it blends linearly across the tube edge.
    #[must_use]
    pub fn velocity_increment(&self, x: f64, r: f64) -> f64 {
        let s = r.abs() / self.tube_radius(x);
        // 1 well inside, 0 well outside, linear across 1 ± EDGE_BLEND.
        let weight = ((1.0 + EDGE_BLEND - s) / (2.0 * EDGE_BLEND)).clamp(0.0, 1.0);
        weight * self.axial_increment(x)
    }

    /// Increase in dynamic pressure `Δq = ½ρ[(V + Δu)² − V²]` (Pa) at station
    /// `x` and radial distance `r`. Defined at zero freestream speed as well.
    #[must_use]
    pub fn dynamic_pressure_increment(&self, x: f64, r: f64) -> f64 {
        let du = self.velocity_increment(x, r);
        0.5 * self.density * (2.0 * self.v_inf * du + du * du)
    }

    /// Change in lift on a surface of area `area` (m²) with lift-curve slope
    /// `cl_alpha` (per rad) at incidence `alpha` (rad), at station `x` and radial
    /// distance `r` (N). Only the dynamic-pressure change is applied: the
    /// slipstream is axial, so it raises the local speed without changing the
    /// local incidence in this model.
    #[must_use]
    pub fn lift_increment(&self, x: f64, r: f64, area: f64, cl_alpha: f64, alpha: f64) -> f64 {
        self.dynamic_pressure_increment(x, r) * area * cl_alpha * alpha
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * b.abs().max(1.0)
    }

    #[test]
    fn disk_velocity_satisfies_momentum_theory() {
        for &(t, v) in &[(50.0, 0.0), (120.0, 20.0), (400.0, 60.0)] {
            let s = Slipstream::new(0.9, t, v, 1.225).expect("valid");
            let a = core::f64::consts::PI * 0.81;
            let vd = s.disk_induced_velocity();
            let thrust = 2.0 * 1.225 * a * vd * (v + vd);
            assert!(close(thrust, t, 1e-9), "T = {thrust}, want {t}");
        }
    }

    #[test]
    fn static_disk_velocity_is_classical_hover_value() {
        let (t, r, rho) = (200.0, 1.0, 1.225);
        let s = Slipstream::new(r, t, 0.0, rho).expect("valid");
        let vh = (t / (2.0 * rho * core::f64::consts::PI * r * r)).sqrt();
        assert!(close(s.disk_induced_velocity(), vh, 1e-12));
    }

    #[test]
    fn axis_increment_runs_from_disk_value_to_twice_it() {
        let s = Slipstream::new(0.9, 150.0, 25.0, 1.225).expect("valid");
        let vd = s.disk_induced_velocity();
        assert!(close(s.axial_increment(0.0), vd, 1e-12));
        assert!(close(s.axial_increment(1.0e4), 2.0 * vd, 1e-3));
        assert!(s.axial_increment(-1.0e4).abs() < 1e-3 * vd);
    }

    #[test]
    fn tube_conserves_mass_flow() {
        let s = Slipstream::new(0.9, 150.0, 25.0, 1.225).expect("valid");
        let a0 = core::f64::consts::PI * 0.81;
        let q0 = a0 * (25.0 + s.disk_induced_velocity());
        for &x in &[0.0, 0.5, 2.0, 10.0, 100.0] {
            let rw = s.tube_radius(x);
            let q = core::f64::consts::PI * rw * rw * (25.0 + s.axial_increment(x));
            assert!(close(q, q0, 1e-9), "x = {x}: Q = {q}, Q0 = {q0}");
        }
        // Far-wake contraction ratio is (V + v_d)/(V + 2 v_d) in area.
        let vd = s.disk_induced_velocity();
        let far = s.tube_radius(1.0e6);
        let want = 0.9 * ((25.0 + vd) / (25.0 + 2.0 * vd)).sqrt();
        assert!(close(far, want, 1e-6));
    }

    #[test]
    fn velocity_increment_is_one_sided_at_the_tube_edge() {
        let s = Slipstream::new(0.9, 150.0, 25.0, 1.225).expect("valid");
        let x = 1.5;
        let rw = s.tube_radius(x);
        let du = s.axial_increment(x);
        assert!(close(s.velocity_increment(x, 0.0), du, 1e-12));
        assert!(close(s.velocity_increment(x, 0.5 * rw), du, 1e-12));
        assert!(s.velocity_increment(x, 2.0 * rw).abs() < 1e-12);
        // Monotone decreasing across the edge.
        let mut prev = f64::INFINITY;
        for i in 0..=40 {
            let r = rw * (0.8 + 0.01 * i as f64);
            let v = s.velocity_increment(x, r);
            assert!(v <= prev + 1e-12);
            prev = v;
        }
    }

    #[test]
    fn dynamic_pressure_increment_is_exact_and_nonnegative_for_positive_thrust() {
        let s = Slipstream::new(0.9, 150.0, 25.0, 1.225).expect("valid");
        let du = s.velocity_increment(0.0, 0.0);
        let want = 0.5 * 1.225 * ((25.0 + du).powi(2) - 25.0_f64.powi(2));
        assert!(close(s.dynamic_pressure_increment(0.0, 0.0), want, 1e-12));
        assert!(s.dynamic_pressure_increment(0.5, 0.2) > 0.0);
        // Static case: Δq = ½ρ Δu² (no freestream term).
        let st = Slipstream::new(0.9, 150.0, 0.0, 1.225).expect("valid");
        let du0 = st.velocity_increment(0.0, 0.0);
        assert!(close(
            st.dynamic_pressure_increment(0.0, 0.0),
            0.5 * 1.225 * du0 * du0,
            1e-12
        ));
    }

    #[test]
    fn rejects_invalid_inputs() {
        assert!(Slipstream::new(0.0, 10.0, 0.0, 1.225).is_err());
        assert!(Slipstream::new(1.0, -5.0, 0.0, 1.225).is_err());
        assert!(Slipstream::new(1.0, 5.0, -1.0, 1.225).is_err());
        assert!(Slipstream::new(1.0, 5.0, 0.0, 0.0).is_err());
        assert!(Slipstream::new(f64::NAN, 5.0, 0.0, 1.225).is_err());
    }
}

#[cfg(test)]
mod zero_thrust_tests {
    use super::*;

    #[test]
    fn zero_thrust_and_zero_speed_are_finite() {
        let s = Slipstream::new(1.0, 0.0, 0.0, 1.225).expect("valid");
        for &(x, r) in &[(0.0, 0.0), (0.5, 0.0), (0.5, 1.0), (-1.0, 0.3)] {
            assert!(s.velocity_increment(x, r).is_finite());
            assert!(s.dynamic_pressure_increment(x, r).is_finite());
            assert!(s.tube_radius(x).is_finite());
        }
        assert_eq!(s.dynamic_pressure_increment(0.5, 0.0), 0.0);
    }
}
