# Changelog

All notable changes to `tpt-aero-panel` are documented here. Workspace-wide
notes are in the root [CHANGELOG.md](../../CHANGELOG.md).

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Breaking
- `compress::{beta, prandtl_glauert, karman_tsien}` return `Result<f64, CompressibilityError>` instead of panicking.
- The crate no longer depends on `tpt-aero-vlm`. `lu_solve` moved to `tpt-aero-core::linalg`.

### Added
- Example `airfoil_pressure`: doublet lift, zero-incidence source pressure and Prandtl-Glauert correction.
- `doublet::solve_doublet_panels`: lifting panel solver with a Kutta condition (Hess-Smith style), returning `DoubletSolution`.
- `source::solve_source_panels`: constant-strength source panels for non-lifting bodies.
- `airfoil::naca4_closed` input path, for sharp trailing edges with no gap panel.
- `PanelError::ClockwiseContour`: clockwise contours are rejected.
- `slender::SlenderBody`: `munk_lift_and_moment` (first-order lift and moment in incidence, `rho V alpha dS/dx`) and `munk_side_force_and_yaw` (side force and yawing moment at sideslip).
- Criterion benchmark (`benches/panel.rs`).
- Property test: a symmetric section at zero incidence has zero lift.

### Fixed
- **`solve_doublet_panels` surface pressure at incidence was not physical.** The bound circulation was a point vortex at the contour centroid, which sits close to the surface near mid-chord. The surface point nearest it saw a huge induced speed, leaving a suction spike (NACA 0012 at 4 degrees: Cp -2.1 at x/c 0.5, growing with panel count). The circulation is now a uniform vortex sheet on the contour, with one extra unknown and the same Kutta condition. Lift is unchanged to within 0.1 % and now converges with panel count.
- Panel influences in the lifting solver are exact constant-panel integrals, not 4-point Gauss quadrature, which was inaccurate for adjacent panels.
- The lifting solver's zero-incidence pressure equals the source solver's to quadrature accuracy.
- `doublet.rs` documentation corrected to match the implementation.
- NACA 0012 validated against published inviscid pressure distributions and thick-airfoil theory (`Cp_min` within 0.04, lift slope within 0.04).
- Contour winding is checked, so counter-clockwise input is required and documented.
