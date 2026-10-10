//! # tpt-aero-core
//!
//! The foundational layer of tpt-aero: unit-safe aerodynamic quantities,
//! non-dimensional coefficients, airfoil polars with stall modeling, and
//! lifting-surface geometry primitives.
//!
//! Everything in this crate is `no_std` + `alloc` compatible and forbidden
//! from `unsafe` (workspace-wide `unsafe_code = "forbid"`).
//!
//! ## Organisation
//!
//! * [`quantities`] — phantom-typed freestream vs local velocities, dynamic
//!   pressure, Mach/Reynolds numbers, α/β angles.
//! * [`coefficients`] — `C_L`/`C_D`/`C_m`/… strictly separated from Newtons.
//! * [`polar`] — tabulated airfoil polars, (M, Re) blending, stall model.
//! * [`airfoil`] — NACA 4-digit generation and coordinate transforms.
//! * [`geometry`] — `Vec3`, lifting-line panels, quad surface panels, wing
//!   planforms with refinement.
//! * [`num`] — the [`Real`] abstraction that lets solver kernels run on `f64`
//!   or on forward-mode dual numbers (exact first derivatives).
//!
//! # Example
//!
//! ```
//! use tpt_aero_core::{LiftCoefficient, MachNumber, WingPlanform};
//!
//! let wing = WingPlanform::rectangular(10.0, 1.0);
//! let mach = MachNumber::new(0.3);
//! let cl = LiftCoefficient::new(0.5);
//! assert!(cl.value() > 0.0 && mach.value() < 0.5);
//! ```

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub mod airfoil;
pub mod coefficients;
pub mod geometry;
pub mod num;
pub mod polar;
pub mod quantities;

pub use airfoil::{AirfoilCoordinates, AirfoilError};
pub use coefficients::{
    DragCoefficient, ForceVector, LiftCoefficient, MomentVector, PitchingMomentCoefficient,
    RollingMomentCoefficient, SideForceCoefficient, YawingMomentCoefficient,
};
pub use geometry::{
    GeometryError, LiftingLinePanel, SurfacePanel, Vec3, WingPlanform, WingSection,
};
pub use num::Real;
pub use polar::{AirfoilPolar, PolarError, PolarSet, StallModel};
pub use quantities::{
    AngleOfAttack, AngleOfSideslip, DynamicPressure, Freestream, FreestreamVelocity, LocalVelocity,
    MachNumber, ReynoldsNumber,
};
