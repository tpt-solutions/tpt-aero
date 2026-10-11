//! # tpt-aero-vlm
//!
//! From-scratch Vortex Lattice Method: horseshoe-vortex lattices on
//! [`tpt_aero_core::WingPlanform`] geometry, aerodynamic influence coefficient
//! (AIC) assembly, and native linear solves (tpt-math dense partial-pivot LU
//! for `f64`; a from-scratch generic LU so the identical code path runs on
//! forward-mode dual numbers — differentiable AIC assembly for gradient-based
//! wing optimization).
//!
//! # Example
//!
//! ```
//! use tpt_aero_core::WingPlanform;
//! use tpt_aero_vlm::{CamberModel, VlmModel};
//!
//! // Rectangular AR = 10 wing at 5° angle of attack.
//! let wing = WingPlanform::rectangular(10.0, 1.0);
//! let model = VlmModel::from_planform(&wing, 24, 4, CamberModel::None).expect("valid");
//! let sol = model.solve(5.0_f64.to_radians(), 50.0, 1.225, &[]).expect("solvable");
//!
//! // Lifting-line estimate: C_Lα = 2π·AR/(AR+2) ⇒ C_L ≈ 0.457.
//! assert!((sol.cl - 0.457).abs() < 0.05, "cl = {}", sol.cl);
//! assert!(sol.cl > 0.0);
//! ```

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub mod aircraft;
pub mod avl;
pub mod export;
pub mod lattice;
pub mod optimize;
pub mod solver;
pub mod vortex;

pub use aircraft::{Control, ControlRegion, ControlSymmetry, Surface};
pub use avl::{parse as parse_avl, AvlDeck, AvlError};
pub use lattice::{CamberModel, Lattice, LatticeElement, LatticeError, SpanSpacing};
pub use optimize::{optimize_linear_twist, optimize_twist_shape, TwistOptimum, TwistShapeOptimum};
pub use solver::{
    validity_flags, FlightState, InfluenceKernel, PreparedVlm, ReferenceGeometry, StateLoads,
    ValidityFlags, VlmError, VlmModel, VlmSolution, VlmSolutionFields, WakeModel,
    VALIDITY_ALPHA_LIMIT, VALIDITY_MACH_LIMIT, VALIDITY_REYNOLDS_MIN,
};
/// Native linear solves, re-exported from `tpt-aero-core` so existing paths keep working.
pub use tpt_aero_core::linalg;
pub use tpt_aero_core::linalg::{cholesky_spd, lu_solve, lu_solve_generic, LuFactors, SolveError};
pub use vortex::{horseshoe_velocity, segment_velocity, semi_infinite_velocity};
