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
- `StateSpace::with_feedback` / `eigenvalues` for closed-loop analysis; `control_design` example (yaw damper cuts the ringing yaw rate ~9x and stabilizes the spiral; pitch-attitude hold gives ~1 s rise with no overshoot).
- `sim`: RK4 `StateSpace` simulation; `FullLongitudinalPlant::state_space` (elevator + vertical gust inputs) and `LateralPlant::state_space` (aileron, rudder, sidewash). Phugoid period/decay, Dutch roll, steady updraft alignment and roll-rate balance verified in the time domain; 1-cos gust alleviation matches FAR 25 Kg within 5% at 12.5 chords.
- `gust_response` example: discrete-gust gradient sweep and Dryden turbulence load factors, using the atmosphere crate.
- `tpt_aero_vlm::export::to_vtk` and `tpt-aero vtk case.toml --alpha 5 --out lattice.vtk`: legacy ASCII VTK of the lattice with circulation, surface and strip cell data (ParaView, PyVista).
- `phugoid::FullLongitudinalPlant`: four-state longitudinal plant with speed/drag couplings; short-period and phugoid modes (phugoid checked against Lanchester: omega*V constant, damping tracks 1/(sqrt2 L/D)). `drag_alpha_fd` helper.
- `controls::control_derivatives`: lift, drag, side force and the three moments per radian of any named control.
- `tpt_aero_stability::trim::trim_level_flight`: Newton trim of angle of attack and elevator for level flight (plus thrust required).
- `LateralPlant`/`LateralParameters`/`LateralModes`: lateral-directional state-space plant `[beta, p, r, phi]` with Dutch-roll, roll-subsidence and spiral modes; `modes::eigenvalues4` (no-dependency 4x4 eigenvalues).
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
