//! # tpt-aero-rotor
//!
//! Blade Element Momentum Theory (BEMT) for propellers, rotors, and turbines:
//! iterative induction solve with physical-bound clamping, Prandtl/Glauert
//! tip losses, hub-loss corrections, skewed-wake and first-order dynamic
//! inflow, and advance-ratio performance maps.
//!
//! # Example
//!
//! ```
//! use tpt_aero_rotor::{
//!     blade::{BladeSection, RotorGeometry, SectionModel},
//!     bemt::{BemtSolver, OperatingCondition},
//!     loss::{HubLoss, TipLoss},
//! };
//!
//! let geom = RotorGeometry::rectangular(1.0, 0.1, 2, 0.12, 20)
//!     .with_linear_twist(0.30, 0.05);
//! let solver = BemtSolver::new(
//!     geom,
//!     SectionModel::Linear { cl_slope: 6.2832, cl_max: 1.2, cd0: 0.01, cd_cl2: 0.01 },
//!     TipLoss::Prandtl,
//!     HubLoss::None,
//! )
//! .expect("valid");
//! let op = OperatingCondition::hover(300.0, 1.225);
//! let sol = solver.solve(&op).expect("converged");
//! assert!(sol.loads.ct > 0.0);
//! assert!(sol.sections.iter().all(|s| s.lambda_induced >= 0.0 && s.lambda_induced <= 1.0));
//! ```

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub mod bemt;
pub mod wake_field;
pub mod blade;
pub mod inflow;
pub mod loss;
pub mod maps;
pub mod slipstream;

pub use bemt::{
    solve_induction_bounded, BemtSolution, BemtSolver, OperatingCondition, RotorLoads,
    SectionResult,
};
pub use blade::{BladeSection, RotorError, RotorGeometry, SectionModel};
pub use inflow::{DynamicInflow, SkewedWake};
pub use loss::{HubLoss, TipLoss};
pub use maps::{performance_map, MapPoint};
pub use slipstream::Slipstream;
pub use wake_field::{CylinderWake, SlipstreamModel};
