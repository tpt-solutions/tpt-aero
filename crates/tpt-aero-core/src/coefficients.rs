//! Non-dimensional aerodynamic coefficients and dimensional forces.
//!
//! Coefficients (`C_L`, `C_D`, `C_m`, …) and dimensional quantities
//! (`newton`, `newton·metre`) are distinct types. Converting between them
//! requires the reference dynamic pressure and area explicitly — there is no
//! `From`/`Into` shortcut that could hide a missing *q·S* factor.

use tpt_math_units::si::f64::{Area, Force, Torque};
use tpt_math_units::si::area::square_meter;
use tpt_math_units::si::force::newton;
use tpt_math_units::si::length::meter;
use tpt_math_units::si::torque::newton_meter;

use crate::num::Real as _;

use crate::quantities::DynamicPressure;

macro_rules! coefficient {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
        pub struct $name(f64);

        impl $name {
            /// Wrap a raw non-dimensional value.
            #[must_use]
            pub const fn new(value: f64) -> Self {
                Self(value)
            }

            /// Raw non-dimensional value.
            #[must_use]
            pub const fn value(&self) -> f64 {
                self.0
            }

            /// Multiply by *q·S* to obtain the dimensional result.
            ///
            /// For force coefficients the result is a force in N; for moment
            /// coefficients use [`CoefficientWithArea::to_moment`] with the
            /// reference chord instead.
            #[must_use]
            pub fn times_q_s(&self, q: &DynamicPressure, s: &Area) -> Force {
                Force::new::<newton>(self.0 * q.si() * s.get::<square_meter>())
            }
        }
    };
}

coefficient!(
    /// Lift coefficient `C_L`.
    LiftCoefficient
);
coefficient!(
    /// Drag coefficient `C_D`.
    DragCoefficient
);
coefficient!(
    /// Side-force coefficient `C_Y`.
    SideForceCoefficient
);

macro_rules! moment_coefficient {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
        pub struct $name(f64);

        impl $name {
            /// Wrap a raw non-dimensional value.
            #[must_use]
            pub const fn new(value: f64) -> Self {
                Self(value)
            }

            /// Raw non-dimensional value.
            #[must_use]
            pub const fn value(&self) -> f64 {
                self.0
            }

            /// Multiply by *q·S·c* to obtain the dimensional moment in N·m.
            #[must_use]
            pub fn times_q_s_c(
                &self,
                q: &DynamicPressure,
                s: &Area,
                c: &tpt_math_units::si::f64::Length,
            ) -> Torque {
                Torque::new::<newton_meter>(
                    self.0 * q.si() * s.get::<square_meter>() * c.get::<meter>(),
                )
            }
        }
    };
}

moment_coefficient!(
    /// Pitching-moment coefficient `C_m` (about the y/body-pitch axis).
    PitchingMomentCoefficient
);
moment_coefficient!(
    /// Rolling-moment coefficient `C_l` (about the x/body-roll axis).
    RollingMomentCoefficient
);
moment_coefficient!(
    /// Yawing-moment coefficient `C_n` (about the z/body-yaw axis).
    YawingMomentCoefficient
);

/// Dimensional force vector in body axes (N).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ForceVector {
    /// Force along the body x axis.
    pub x: Force,
    /// Force along the body y axis.
    pub y: Force,
    /// Force along the body z axis.
    pub z: Force,
}

impl ForceVector {
    /// Zero force vector.
    #[must_use]
    pub fn zero() -> Self {
        Self {
            x: Force::new::<newton>(0.0),
            y: Force::new::<newton>(0.0),
            z: Force::new::<newton>(0.0),
        }
    }

    /// Magnitude in N.
    #[must_use]
    pub fn magnitude(&self) -> f64 {
        let (x, y, z) = (
            self.x.get::<newton>(),
            self.y.get::<newton>(),
            self.z.get::<newton>(),
        );
        (x * x + y * y + z * z).sqrt()
    }
}

/// Dimensional moment vector in body axes (N·m).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MomentVector {
    /// Rolling moment about x.
    pub l: Torque,
    /// Pitching moment about y.
    pub m: Torque,
    /// Yawing moment about z.
    pub n: Torque,
}

impl MomentVector {
    /// Zero moment vector.
    #[must_use]
    pub fn zero() -> Self {
        Self {
            l: Torque::new::<newton_meter>(0.0),
            m: Torque::new::<newton_meter>(0.0),
            n: Torque::new::<newton_meter>(0.0),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_math_units::si::area::square_meter as sqm;
    use tpt_math_units::si::f64::Length;
    use tpt_math_units::si::length::meter;

    #[test]
    fn coefficient_to_force() {
        let cl = LiftCoefficient::new(0.5);
        let q = DynamicPressure::from_si(1000.0);
        let s = Area::new::<sqm>(2.0);
        let lift = cl.times_q_s(&q, &s);
        assert!((lift.get::<newton>() - 1000.0).abs() < 1e-9);
    }

    #[test]
    fn moment_coefficient_to_torque() {
        let cm = PitchingMomentCoefficient::new(-0.05);
        let q = DynamicPressure::from_si(2000.0);
        let s = Area::new::<sqm>(1.0);
        let c = Length::new::<meter>(0.5);
        let m = cm.times_q_s_c(&q, &s, &c);
        assert!((m.get::<newton_meter>() + 50.0).abs() < 1e-9);
    }
}
