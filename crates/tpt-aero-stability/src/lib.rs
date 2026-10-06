//! # tpt-aero-stability
//!
//! Flight-mechanics derivative extraction and plant assembly:
//!
//! * [`longitudinal`] — finite-difference and dual-number extraction of
//!   `C_{Lα}`, `C_{mα}`, `C_{Lq}`, `C_{mq}`, `C_{Lδ}`, `C_{mδ}` from a VLM
//!   plant, plus thin-airfoil analytical references.
//! * [`plant`] — linearized longitudinal state-space plant (`A`, `B`) with
//!   short-period characteristics: the handoff hook for a future
//!   `tpt-control` crate.
//!
//! # Example
//!
//! ```
//! use tpt_aero_core::WingPlanform;
//! use tpt_aero_stability::{longitudinal::longitudinal_fd, plant::{LongitudinalPlant, PlantParameters}};
//! use tpt_aero_vlm::{CamberModel, VlmModel};
//!
//! let wing = WingPlanform::rectangular(10.0, 1.0);
//! // Reference the moment at the leading edge for static stability.
//! let model = VlmModel::from_planform(&wing, 16, 4, CamberModel::None)
//!     .expect("valid")
//!     .with_reference(tpt_aero_vlm::ReferenceGeometry {
//!         area: wing.area(),
//!         chord: wing.mac(),
//!         span: wing.span(),
//!         moment_center: tpt_aero_core::Vec3::new(0.0, 0.0, 0.0),
//!     });
//! let d = longitudinal_fd(&model, 2.0_f64.to_radians(), 50.0, 1.225, 1e-4);
//! assert!((4.7..=5.2).contains(&d.cl_alpha));
//! assert!(d.cm_alpha < 0.0, "static stability ahead of the AC");
//!
//! let p = PlantParameters { v_inf: 50.0, chord: 1.0, area: 10.0, mass: 500.0, inertia_y: 400.0, density: 1.225 };
//! let plant = LongitudinalPlant::from_derivatives(&d, &p, 0.05);
//! let (wn, zeta) = plant.short_period();
//! assert!(wn > 0.0 && zeta > 0.0);
//! ```

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub mod longitudinal;
pub mod plant;
