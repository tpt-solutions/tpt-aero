//! Unit-safe aerodynamic quantities.
//!
//! Phantom-typed roles distinguish *where* a velocity lives: the unperturbed
//! [`Freestream`] far from the body versus the [`Local`] flow at a surface
//! point. Mixing them up is a compile error, not a silent wrong answer.
//!
//! Non-dimensional coefficients live in [`crate::coefficients`]; dimensional
//! forces never convert to coefficients implicitly.

use core::marker::PhantomData;

use tpt_math_units::si::f64::{Angle, MassDensity, Pressure, Velocity};
use tpt_math_units::si::angle::radian;
use tpt_math_units::si::mass_density::kilogram_per_cubic_meter;
use tpt_math_units::si::pressure::pascal;
use tpt_math_units::si::velocity::meter_per_second;

/// Marker role: the unperturbed flow far upstream of the body.
#[derive(Debug, Clone, Copy)]
pub struct Freestream;

/// Marker role: the local flow velocity at a point on (or induced by) the body.
#[derive(Debug, Clone, Copy)]
pub struct Local;

/// A velocity tagged with its aerodynamic role.
///
/// `FlowVelocity<Freestream>` and `FlowVelocity<Local>` are distinct types:
/// adding them, assigning one to the other, or passing one where the other is
/// expected is a compile error. Re-tagging requires the explicit
/// [`FlowVelocity::retag`].
#[derive(Debug, Clone, Copy)]
pub struct FlowVelocity<Role> {
    v: Velocity,
    _role: PhantomData<fn(Role)>,
}

impl<Role> FlowVelocity<Role> {
    /// Construct from a raw SI velocity magnitude (m/s).
    #[must_use]
    pub fn from_si(meters_per_second: f64) -> Self {
        Self {
            v: Velocity::new::<meter_per_second>(meters_per_second),
            _role: PhantomData,
        }
    }

    /// Construct from a `uom` velocity.
    #[must_use]
    pub const fn from_quantity(v: Velocity) -> Self {
        Self {
            v,
            _role: PhantomData,
        }
    }

    /// Magnitude in m/s.
    #[must_use]
    pub fn si(&self) -> f64 {
        self.v.get::<meter_per_second>()
    }

    /// The underlying `uom` quantity.
    #[must_use]
    pub const fn quantity(&self) -> &Velocity {
        &self.v
    }

    /// Change the role tag. Only for the rare, deliberate cross-role reuse —
    /// the compiler now cannot catch a mistake you make here.
    #[must_use]
    pub fn retag<Role2>(self) -> FlowVelocity<Role2> {
        FlowVelocity {
            v: self.v,
            _role: PhantomData,
        }
    }
}

/// Alias: freestream (flight) velocity.
pub type FreestreamVelocity = FlowVelocity<Freestream>;
/// Alias: local flow velocity at a surface/induction point.
pub type LocalVelocity = FlowVelocity<Local>;

/// Dynamic pressure *q* = ½ ρ V² (Pa).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DynamicPressure {
    q: Pressure,
}

impl DynamicPressure {
    /// From the freestream speed and density.
    #[must_use]
    pub fn from_flow(speed: &FreestreamVelocity, density: &MassDensity) -> Self {
        let v = speed.si();
        Self {
            q: Pressure::new::<pascal>(0.5 * density.get::<kilogram_per_cubic_meter>() * v * v),
        }
    }

    /// From a raw SI value in Pa.
    #[must_use]
    pub fn from_si(pascals: f64) -> Self {
        Self {
            q: Pressure::new::<pascal>(pascals),
        }
    }

    /// Value in Pa.
    #[must_use]
    pub fn si(&self) -> f64 {
        self.q.get::<pascal>()
    }

    /// The underlying `uom` pressure.
    #[must_use]
    pub const fn quantity(&self) -> &Pressure {
        &self.q
    }
}

/// Mach number *M* = V / a. Dimensionless, but a distinct type so it cannot be
/// confused with a Reynolds number or a coefficient.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct MachNumber(f64);

impl MachNumber {
    /// Construct from a raw value. Panics if negative or non-finite — there is
    /// no physical subsonic model for M < 0 in this library.
    #[must_use]
    pub fn new(mach: f64) -> Self {
        assert!(mach >= 0.0 && mach.is_finite(), "Mach must be >= 0: {mach}");
        Self(mach)
    }

    /// Raw value.
    #[must_use]
    pub const fn value(&self) -> f64 {
        self.0
    }

    /// True in the subsonic regime this library models (M < 1).
    #[must_use]
    pub const fn is_subsonic(&self) -> bool {
        self.0 < 1.0
    }
}

/// Reynolds number *Re* = ρ V L / μ. Distinct type; strictly positive.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct ReynoldsNumber(f64);

impl ReynoldsNumber {
    /// Construct from a raw value. Panics if ≤ 0 or non-finite.
    #[must_use]
    pub fn new(re: f64) -> Self {
        assert!(re > 0.0 && re.is_finite(), "Reynolds must be > 0: {re}");
        Self(re)
    }

    /// Raw value.
    #[must_use]
    pub const fn value(&self) -> f64 {
        self.0
    }
}

/// Angle of attack α (rad). Positive nose-up relative to the freestream.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct AngleOfAttack {
    alpha: Angle,
}

impl AngleOfAttack {
    /// From radians.
    #[must_use]
    pub fn from_radians(rad: f64) -> Self {
        Self {
            alpha: Angle::new::<radian>(rad),
        }
    }

    /// From degrees.
    #[must_use]
    pub fn from_degrees(deg: f64) -> Self {
        Self::from_radians(deg.to_radians())
    }

    /// Value in radians.
    #[must_use]
    pub fn radians(&self) -> f64 {
        self.alpha.get::<radian>()
    }

    /// Value in degrees.
    #[must_use]
    pub fn degrees(&self) -> f64 {
        self.radians().to_degrees()
    }

    /// The underlying `uom` angle.
    #[must_use]
    pub const fn quantity(&self) -> &Angle {
        &self.alpha
    }
}

/// Angle of sideslip β (rad). Positive with the relative wind from the right.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct AngleOfSideslip {
    beta: Angle,
}

impl AngleOfSideslip {
    /// From radians.
    #[must_use]
    pub fn from_radians(rad: f64) -> Self {
        Self {
            beta: Angle::new::<radian>(rad),
        }
    }

    /// From degrees.
    #[must_use]
    pub fn from_degrees(deg: f64) -> Self {
        Self::from_radians(deg.to_radians())
    }

    /// Value in radians.
    #[must_use]
    pub fn radians(&self) -> f64 {
        self.beta.get::<radian>()
    }

    /// Value in degrees.
    #[must_use]
    pub fn degrees(&self) -> f64 {
        self.radians().to_degrees()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_math_units::si::f64::Velocity as UomVelocity;
    use tpt_math_units::si::velocity::meter_per_second as mps;

    #[test]
    fn dynamic_pressure_from_flow() {
        let v = FreestreamVelocity::from_quantity(UomVelocity::new::<mps>(50.0));
        let rho = MassDensity::new::<kilogram_per_cubic_meter>(1.225);
        let q = DynamicPressure::from_flow(&v, &rho);
        let expected = 0.5 * 1.225 * 50.0_f64.powi(2);
        assert!((q.si() - expected).abs() < 1e-9);
    }

    #[test]
    fn roles_are_distinct_types() {
        let freestream = FreestreamVelocity::from_si(10.0);
        // Compile-time separation: only via explicit retag.
        let _local: LocalVelocity = freestream.retag();
        assert_eq!(freestream.si(), 10.0);
    }

    #[test]
    fn mach_and_re_validation() {
        assert_eq!(MachNumber::new(0.3).value(), 0.3);
        assert!(MachNumber::new(0.3).is_subsonic());
        assert!(!MachNumber::new(1.2).is_subsonic());
        assert_eq!(ReynoldsNumber::new(1e6).value(), 1e6);
    }

    #[test]
    #[should_panic]
    fn mach_negative_panics() {
        let _ = MachNumber::new(-0.1);
    }

    #[test]
    fn angle_roundtrip() {
        let a = AngleOfAttack::from_degrees(5.0);
        assert!((a.degrees() - 5.0).abs() < 1e-12);
        assert!((a.radians() - 5.0_f64.to_radians()).abs() < 1e-12);
    }
}
