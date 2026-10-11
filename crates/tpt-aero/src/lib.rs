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
//!
//! # A whole aircraft in a few lines
//!
//! ```
//! use tpt_aero::prelude::*;
//!
//! // Wing + tailplane with an elevator; trim for level flight at 50 m/s and
//! // read off the static stability.
//! let aircraft = VlmModel::from_surfaces(&[
//!     Surface::horizontal(WingPlanform::rectangular(10.0, 1.2), 20, 3),
//!     Surface::horizontal(WingPlanform::rectangular(4.0, 0.6), 10, 3)
//!         .at(4.0, 0.0, 0.0)
//!         .with_control("elevator", 0.6, 0.0, 2.0, ControlSymmetry::Symmetric),
//! ])
//! .expect("valid configuration");
//!
//! let trim = trim_level_flight(
//!     &aircraft,
//!     &TrimRequest {
//!         weight_n: 600.0 * 9.80665,
//!         v_inf: 50.0,
//!         density: 1.225,
//!         elevator: "elevator",
//!         cd0: 0.03,
//!     },
//! )
//! .expect("trims");
//! assert!((1.0..8.0).contains(&trim.alpha_rad.to_degrees()));
//!
//! let d = longitudinal_fd(&aircraft, trim.alpha_rad, 50.0, 1.225, 1e-4).expect("solve");
//! assert!(d.cm_alpha < 0.0, "statically stable");
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
        controls::{control_derivatives, ControlDerivatives},
        lateral::{lateral_fd, LateralDerivatives, LateralParameters, LateralPlant},
        longitudinal::{longitudinal_fd, LongitudinalDerivatives},
        phugoid::{drag_alpha_fd, FullLongitudinalParameters, FullLongitudinalPlant},
        plant::{LongitudinalPlant, PlantParameters},
        trim::{trim_level_flight, TrimRequest},
    };
    #[cfg(feature = "vlm")]
    pub use tpt_aero_vlm::{
        CamberModel, ControlSymmetry, FlightState, ReferenceGeometry, Surface, VlmModel,
        VlmSolution,
    };
}
