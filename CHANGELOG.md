# Changelog

All notable changes to this workspace are documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Breaking
- `tpt_aero_stability::plant::LongitudinalPlant::short_period` returns `Option<(f64, f64)>` (`None` when statically unstable).
- `tpt_aero_rotor::DynamicInflow::step` takes an extra `tip_speed` argument (m/s).
- `tpt_aero_panel::{beta, prandtl_glauert, karman_tsien}` return `Result<f64, CompressibilityError>` instead of panicking.
- `naca4` returns `AirfoilError::InvalidResolution` for fewer than 2 stations; accepts thickness up to 40 %.

### Added
- Control surfaces: `Surface::with_control`, `VlmModel::controls`/`control_index`/`control_twists` (elevator, aileron, rudder as incidence changes aft of a hinge), `VlmSolution::stability_roll_yaw`. Aileron roll authority and adverse yaw validated against lifting line.
- Multi-surface aircraft: `Surface` (`horizontal`, `vertical`, `.at`, `.with_incidence`, `.with_camber`) and `VlmModel::from_surfaces` for wing + tailplane + fin with shared downwash/sidewash; `Lattice::combine` and `Lattice::transformed`.
- Lifting-line validation for tapered wings (CL within 4%, loading within 7%, Oswald e within 5%).
- `aircraft_stability` example (neutral point, derivatives, CG sweep with short-period mode).
- Sideslip and body rates in the VLM: `PreparedVlm::solve_state(FlightState)` returns wind-axis forces and body moments in aircraft-standard signs.
- `tpt_aero_stability::lateral::lateral_fd`: `C_Yb, C_lb, C_nb, C_Yp, C_lp, C_np, C_Yr, C_lr, C_nr`.
- `tpt-aero derivatives case.toml --alpha 4` CLI subcommand.
- `VlmModel::prepare` / `PreparedVlm` / `VlmModel::solve_sweep`: factor the AIC and precompute induced-velocity influences once, then solve each angle in O(N^2) (about 20x faster for a 17-angle sweep). `LuFactors` for reusable LU solves.
- criterion benches for the VLM (`cargo bench -p tpt-aero-vlm`).
- `naca4_closed`: closed-trailing-edge NACA 4-digit generator for panel solvers.
- `PanelError::ClockwiseContour`: panel solvers reject clockwise contours.
- `tpt-aero-cli` crate: `tpt-aero wing case.toml --out polar.csv`, `tpt-aero template wing`, `tpt-aero isa <alt>...`.
- `VlmError::InvalidFlightCondition`; `VlmModel::solve` validates angle, speed and density.
- `MachNumber::try_new`, `ReynoldsNumber::try_new`.
- `tpt_aero::prelude`; examples `isa_table`, `airfoil_cp`, `stability_derivatives`.

### Fixed
- **`solve_state` moments and rates are now about stability axes** (body axes rotated by alpha), the convention of stability derivatives. Roll and yaw previously mixed by alpha; found by comparing aileron adverse yaw against lifting-line theory (the sign was wrong in body axes). Lateral derivatives change slightly at non-zero alpha.
- **`WingSection::twist` was ignored by the VLM solvers** (only the per-solve `twists` slice was applied). Effective twist is now section twist plus the slice.
- **`WingPlanform::tapered` built a skewed, non-symmetric wing** (one segment from the left tip to the right tip), so any taper, sweep or dihedral broke mirror symmetry. It now builds three sections (tip, root, tip) with V-dihedral. Results for tapered, swept or dihedral wings change.
- Longitudinal plant: control-matrix area double count, lift control row, missing `C_Lq` and gravity terms.
- BEMT: reported iteration count; thrust no longer integrated inside the hub cutout.
- Dynamic inflow time constant was off by a factor of the tip speed.
- Frame-convention doc in `tpt-aero-core` now matches the solver.
- CI checks out `tpt-math`.
