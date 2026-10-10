//! # tpt-aero-panel
//!
//! Panel methods for higher-fidelity 3-D/2-D aerodynamics:
//!
//! * [`source`] — 2-D constant-strength source panels (Hess–Smith) for
//!   non-lifting bodies; validated against the circular-cylinder analytic
//!   pressure distribution.
//! * [`doublet`] — constant-strength doublet (vortex-panel equivalent)
//!   method with a Kutta condition for 2-D lifting flow.
//! * [`slender`] — slender-body theory for bodies of revolution.
//! * [`compress`] — Prandtl–Glauert and von Kármán–Tsien compressibility
//!   corrections.
//!
//! All linear solves go through the in-house `tpt-math` dense LU.
//!
//! # Example
//!
//! ```
//! use tpt_aero_panel::{compress, source};
//!
//! // Cylinder Cp amplification under Prandtl–Glauert at M = 0.6.
//! let cp_i = -3.0;               // shoulder suction, incompressible
//! let cp_c = compress::prandtl_glauert(cp_i, 0.6).unwrap();
//! assert!(cp_c < cp_i);
//!
//! // Source panels on a circle reproduce Cp = 1 − 4 sin²θ.
//! let pts: Vec<(f64, f64)> = (0..48).map(|i| {
//!     let t = 2.0 * std::f64::consts::PI * i as f64 / 48.0;
//!     (t.cos(), t.sin())
//! }).collect();
//! let sol = source::solve_source_panels(&pts, 0.0, 1.0).expect("solvable");
//! assert!(sol.cp[0] > 0.9);
//! ```

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub mod compress;
pub mod doublet;
pub mod slender;
pub mod source;

pub use compress::{beta, is_applicable, karman_tsien, prandtl_glauert, CompressibilityError};
pub use doublet::{solve_doublet_panels, DoubletSolution};
pub use slender::{SlenderBody, SlenderError};
pub use source::{solve_source_panels, PanelError, SourceSolution};
