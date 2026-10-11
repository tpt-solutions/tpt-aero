//! # tpt-aero-atmo
//!
//! Atmospheric models for flight simulation:
//!
//! * [`isa`] — International Standard Atmosphere (0..86 km geopotential),
//!   with pressure-altitude inversion, ΔISA temperature offsets and
//!   Sutherland viscosity.
//! * [`shear`] — power-law / logarithmic wind shear and discrete gusts
//!   (`1−cos`, sharp-edge).
//! * [`turbulence`] — Dryden and von Kármán continuous turbulence with
//!   discretized second-order shaping filters, MIL-HDBK-1797 scale lengths,
//!   and a `no_std`-friendly deterministic noise source.
//!
//! # Example
//!
//! ```
//! use tpt_aero_atmo::isa;
//!
//! let sea_level = isa::at(0.0).expect("in range");
//! assert!((sea_level.density - 1.225).abs() < 2e-4);
//! let cruise = isa::at(11_000.0).expect("in range");
//! assert!((cruise.temperature - 216.65).abs() < 1e-9);
//! ```

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub mod isa;
pub mod shear;
pub mod turbulence;

pub use isa::{
    altitude_from_pressure, at, at_isa_offset, dynamic_viscosity, kinematic_viscosity, IsaError,
    IsaState, G0, R_AIR,
};
pub use shear::{log_law, power_law, DiscreteGust, GustShape};
pub use turbulence::{simulate_axis, TurbulenceAxis, TurbulenceField, TurbulenceModel, WhiteNoise};
