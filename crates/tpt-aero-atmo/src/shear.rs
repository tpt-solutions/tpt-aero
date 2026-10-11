//! Wind shear profiles and discrete gust models.

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// Mean wind speed at height `z` (m) by the power law
/// `u(z) = u_ref (z/z_ref)^p`.
///
/// `p` is the terrain exponent: 0.10 (open sea/flat), 0.14 (open terrain),
/// 0.20 (suburban), 0.33 (urban).
#[must_use]
pub fn power_law(u_ref: f64, z_ref: f64, z: f64, p: f64) -> f64 {
    if z <= 0.0 || z_ref <= 0.0 {
        return 0.0;
    }
    u_ref * (z / z_ref).powf(p)
}

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// Mean wind speed at height `z` by the logarithmic law
/// `u(z) = u*/κ · ln(z/z0)` with von Kármán constant κ = 0.4.
#[must_use]
pub fn log_law(friction_velocity: f64, z0: f64, z: f64) -> f64 {
    const KAPPA: f64 = 0.4;
    if z <= z0 || z0 <= 0.0 {
        return 0.0;
    }
    friction_velocity / KAPPA * (z / z0).ln()
}

/// Discrete gust shape.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GustShape {
    /// `1 − cos` ramp-up/ramp-down over the gradient distance `2H`
    /// (FAR/JAR 25-ish discrete gust).
    OneMinusCosine,
    /// Sharp-edged step gust.
    SharpEdge,
}

/// A discrete gust: peak speed `u_ds` (m/s) with gradient distance `h`
/// (m) reached over the penetration `0 .. 2h`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiscreteGust {
    /// Design gust speed at the peak (m/s).
    pub u_ds: f64,
    /// Gust gradient half-distance (m).
    pub h: f64,
    /// Shape.
    pub shape: GustShape,
}

impl DiscreteGust {
    /// Build a gust, validating `u_ds ≥ 0`, `h > 0`.
    ///
    /// # Errors
    /// Panics are avoided; returns `None` on invalid parameters.
    #[must_use]
    pub fn new(u_ds: f64, h: f64, shape: GustShape) -> Option<Self> {
        if !(u_ds >= 0.0 && h > 0.0 && u_ds.is_finite() && h.is_finite()) {
            return None;
        }
        Some(Self { u_ds, h, shape })
    }

    /// Gust speed at penetration distance `x` (m) into the gust.
    ///
    /// * `1 − cos`: `u = u_ds/2 · (1 − cos(π x/h))` for `0 ≤ x ≤ 2h`, then 0.
    /// * sharp edge: `u = u_ds` for `x ≥ 0`.
    #[must_use]
    pub fn speed_at(&self, x: f64) -> f64 {
        match self.shape {
            GustShape::SharpEdge => {
                if x >= 0.0 {
                    self.u_ds
                } else {
                    0.0
                }
            }
            GustShape::OneMinusCosine => {
                if x <= 0.0 || x >= 2.0 * self.h {
                    0.0
                } else {
                    0.5 * self.u_ds * (1.0 - (core::f64::consts::PI * x / self.h).cos())
                }
            }
        }
    }

    /// The spatial gradient du/dx at penetration `x` (1/s when multiplied by
    /// the aircraft speed). Peak gradient of the 1−cos shape is `u_ds π/(2h)`
    /// at `x = h/2`, where `sin(π x/h)` reaches 1.
    #[must_use]
    pub fn gradient_at(&self, x: f64) -> f64 {
        match self.shape {
            GustShape::SharpEdge => 0.0, // discontinuous by construction
            GustShape::OneMinusCosine => {
                if x <= 0.0 || x >= 2.0 * self.h {
                    0.0
                } else {
                    0.5 * self.u_ds * core::f64::consts::PI / self.h
                        * (core::f64::consts::PI * x / self.h).sin()
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn power_law_scales() {
        let u = power_law(10.0, 10.0, 100.0, 0.14);
        let expected = 10.0 * 10.0_f64.powf(0.14);
        assert!((u - expected).abs() < 1e-12);
        assert_eq!(power_law(10.0, 10.0, 0.0, 0.14), 0.0);
    }

    #[test]
    fn log_law_shape() {
        let u = log_law(0.5, 0.03, 10.0);
        let expected = 0.5 / 0.4 * (10.0_f64 / 0.03).ln();
        assert!((u - expected).abs() < 1e-12);
        assert_eq!(log_law(0.5, 0.03, 0.03), 0.0);
    }

    #[test]
    fn one_minus_cos_shape() {
        let g = DiscreteGust::new(10.0, 30.0, GustShape::OneMinusCosine).expect("valid");
        assert_eq!(g.speed_at(0.0), 0.0);
        // Peak at the gradient distance x = h.
        assert!((g.speed_at(30.0) - 10.0).abs() < 1e-9);
        // Zero again at 2h.
        assert_eq!(g.speed_at(60.0), 0.0);
        // Half-way points at h/2 and 3h/2.
        assert!((g.speed_at(15.0) - 5.0).abs() < 1e-9);
        // Peak gradient at x = h/2 (du/dx ∝ sin(πx/h)).
        let peak = 0.5 * 10.0 * core::f64::consts::PI / 30.0;
        assert!((g.gradient_at(15.0) - peak).abs() < 1e-9);
        assert!(g.gradient_at(30.0).abs() < 1e-9);
    }

    #[test]
    fn sharp_edge_shape() {
        let g = DiscreteGust::new(5.0, 30.0, GustShape::SharpEdge).expect("valid");
        assert_eq!(g.speed_at(-1.0), 0.0);
        assert_eq!(g.speed_at(0.0), 5.0);
        assert_eq!(g.speed_at(100.0), 5.0);
    }

    #[test]
    fn rejects_invalid() {
        assert!(DiscreteGust::new(-1.0, 30.0, GustShape::SharpEdge).is_none());
        assert!(DiscreteGust::new(5.0, 0.0, GustShape::SharpEdge).is_none());
    }
}
