//! # tpt-aero
//!
//! Umbrella crate: feature-gated re-exports of the tpt-aero aerodynamics and
//! flight-mechanics stack. Enable the sub-crates you need:
//!
//! * `core` — unit-safe quantities, coefficients, airfoil polars, geometry.
//! * `vlm` — horseshoe-vortex Vortex Lattice Method (differentiable).
//! * `panel` — source/doublet panel methods, slender bodies, compressibility.
//! * `rotor` — Blade Element Momentum Theory, inflow models, propeller maps.
//! * `atmo` — ISA atmosphere, wind shear, gusts, Dryden/von Karman turbulence.
//! * `stability` — stability & control derivatives, longitudinal plant.
//! * `verify` — proptest strategies and Kani harnesses (std-only, opt-in).
//!
//! All solver crates are `no_std` + `alloc`; `std` is on by default for
//! convenience and can be dropped for embedded targets.
//!
//! # Example
//!
//! ```
//! use tpt_aero::{WingPlanform, VlmModel, CamberModel};
//!
//! let wing = WingPlanform::rectangular(10.0, 1.0);
//! let model = VlmModel::from_planform(&wing, 24, 4, CamberModel::None).expect("valid");
//! let sol = model.solve(5.0_f64.to_radians(), 50.0, 1.225, &[]).expect("solvable");
//! assert!(sol.cl > 0.4 && sol.cl < 0.5, "C_L = {}", sol.cl);
//! ```

#![no_std]
#![forbid(unsafe_code)]

#[cfg(feature = "std")]
extern crate std;

/// Unit-safe quantities, coefficients, airfoil polars, geometry.
#[cfg(feature = "core")]
pub use tpt_aero_core as core;

/// Vortex Lattice Method.
#[cfg(feature = "vlm")]
pub use tpt_aero_vlm as vlm;

/// Source/doublet panel methods and compressibility.
#[cfg(feature = "panel")]
pub use tpt_aero_panel as panel;

/// Blade Element Momentum Theory and rotor inflow models.
#[cfg(feature = "rotor")]
pub use tpt_aero_rotor as rotor;

/// ISA atmosphere, gusts, and turbulence.
#[cfg(feature = "atmo")]
pub use tpt_aero_atmo as atmo;

/// Stability & control derivatives and plant matrices.
#[cfg(feature = "stability")]
pub use tpt_aero_stability as stability;

/// Verification strategies and harnesses (std-only, opt-in).
#[cfg(all(feature = "verify", feature = "std"))]
pub use tpt_aero_verify as verify;

// Convenience re-exports of the most-used types (default features).
#[cfg(feature = "atmo")]
pub use tpt_aero_atmo::isa;
#[cfg(feature = "core")]
pub use tpt_aero_core::{
    AngleOfAttack, DynamicPressure, MachNumber, ReynoldsNumber, WingPlanform, WingSection,
};
#[cfg(feature = "rotor")]
pub use tpt_aero_rotor::{BemtSolver, OperatingCondition, RotorGeometry};
#[cfg(feature = "vlm")]
pub use tpt_aero_vlm::{CamberModel, VlmModel, VlmSolution};

/// One-line import of the most common types: `use tpt_aero::prelude::*;`.
pub mod prelude {
    #[cfg(feature = "atmo")]
    pub use tpt_aero_atmo::{isa, IsaState};
    #[cfg(feature = "core")]
    pub use tpt_aero_core::{
        airfoil::{naca4, naca4_closed},
        AngleOfAttack, DynamicPressure, MachNumber, ReynoldsNumber, Vec3, WingPlanform,
        WingSection,
    };
    #[cfg(feature = "panel")]
    pub use tpt_aero_panel::{compress, solve_source_panels};
    #[cfg(feature = "rotor")]
    pub use tpt_aero_rotor::{BemtSolver, OperatingCondition, RotorGeometry};
    #[cfg(feature = "stability")]
    pub use tpt_aero_stability::{
        longitudinal::{longitudinal_fd, LongitudinalDerivatives},
        plant::{LongitudinalPlant, PlantParameters},
    };
    #[cfg(feature = "vlm")]
    pub use tpt_aero_vlm::{CamberModel, ReferenceGeometry, VlmModel, VlmSolution};
}
