# Changelog

All notable changes to this workspace are documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Breaking
- `tpt_aero_stability::plant::LongitudinalPlant::short_period` returns `Option<(f64, f64)>` (`None` when statically unstable).
- `tpt_aero_rotor::DynamicInflow::step` takes an extra `tip_speed` argument (m/s).
- `tpt_aero_panel::{beta, prandtl_glauert, karman_tsien}` return `Result<f64, CompressibilityError>` instead of panicking.
- `naca4` returns `AirfoilError::InvalidResolution` for fewer than 2 stations; accepts thickness up to 40 %.
- `tpt_aero_vlm::VlmError` has two new variants (`FreestreamWakeNotPrepared`, `InvalidGroundHeight`); exhaustive matches need the new arms.
- `tpt_aero_vlm::StateLoads` has two new fields (`force_n`, `moment_nm`); struct literals need them.
- `tpt_aero_panel` no longer depends on `tpt_aero_vlm`; `lu_solve` moved to `tpt_aero_core::linalg` (still re-exported from `tpt_aero_vlm`).
- Error types (`VlmError`, `LatticeError`, `PanelError`, `RotorError`, `IsaError`, `SolveError`, `AirfoilError`, `GeometryError`, `PolarError`, `SlenderError`) implement `core::error::Error`, so `?` converts them into `Box<dyn Error>`.

### Added
- `tpt_aero_rotor::CylinderWake`: full induced-velocity field of a uniformly loaded actuator disk, linearized vortex-cylinder wake (the wake is a cylinder of ring strength `2 v_d`, per NASA Memo 4-15-59L). Closed-form ring kernel with complete elliptic integrals; the wake integral is by composite Simpson. Exact on the axis, divergence-free, uniform over the disk.
- `tpt_aero_rotor::SlipstreamModel`: choice between the contracting tube (`Slipstream`) and the cylinder field, with common `local_flow` and `dynamic_pressure_increment`. `propeller_flight` selects it with `SLIPSTREAM`.
- `tpt_aero_rotor::Slipstream`: momentum-theory propeller slipstream (induced velocity at the disk, on-axis axial increment, contracting tube with an edge blend, dynamic-pressure increment at a point). Example `slipstream_tail` in `tpt-aero-rotor`.
- `propeller_flight` example: the slipstream raises the elevator's dynamic pressure, adding tail lift and a nose-down moment; the elevator is re-trimmed to balance it at the start state.
- AVL 3.52 comparison: decks and run script in `docs/validation/avl/`; `tpt-aero-stability/tests/avl_reference.rs` checks C_Lα, C_Di, C_Yβ, C_lβ, C_nβ, C_lp, C_lr, C_nr against recorded AVL values.
- Optional `serde` feature (`tpt-aero/serde`, per crate) on results, parameters and derivative sets; JSON round trips are exact with `serde_json`'s `float_roundtrip`. `WingPlanform` is serialize-only (deserialize through the validating constructors).
- `tpt-aero derivatives case.toml --json`.
- Documentation: rewritten README, per-crate READMEs, `docs/CONVENTIONS.md` (frames, signs, control sense) and `docs/VALIDATION.md` (what is checked against what, and known weaknesses); whole-aircraft doctest in the umbrella crate.
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
- VLM options (all opt-in, defaults unchanged): cosine spanwise spacing (`SpanSpacing::Cosine`, `from_planform_spaced`, `WingPlanform::refine_cosine`); freestream-aligned trailing legs (`WakeModel::Freestream`, per-angle solves only); ground effect by image horseshoes (`with_ground_height`); subsonic compressibility by Göthert's rule (`from_planform_compressible`, lift slope within 2 % of Helmbold-Diederich).
- `prepare_with_twists`: prepared solves with control deflections, for control-deflected stability loads.
- `tpt_aero_core::viscous`: flat-plate skin friction (laminar, Schlichting mixed), Hoerner form factor and a `C_D0` component build-up.
- `tpt_aero_atmo`: `at_isa_offset` (non-standard temperature, hydrostatically consistent) and Sutherland `dynamic_viscosity` / `kinematic_viscosity`.
- `tpt_aero_stability::turn`: nonlinear steady turn and steady sideslip trim (`trim_steady_turn`, `coordinated_turn_rate`), solved by continuation from level flight. Checked for mirror symmetry and load factor 1/cos φ, and steady under 6-DOF integration.
- `tpt_aero_stability::sixdof`: nonlinear 6-DOF simulation (aircraft axes, Euler angles, VLM aerodynamics, thrust along the airspeed vector, parasite `C_D0`). Level trim is an equilibrium; a pitch-rate kick damps.
- `tpt_aero_stability::sim::StateSpace::with_integrator`: integral action (autotrim) for attitude holds; the pitch hold tracks a step with no steady offset.
- `control_derivatives`: `C_Lδa`/`C_nδr` checked; `C_nδr` matches the fin side force times its moment arm to 1 %.
- CLI: `tpt-aero rotor <rotor.toml> [--json]` (BEMT operating point) and `tpt-aero trim <aircraft.toml> --speed V [--altitude h] [--bank deg] [--json]` (level trim or coordinated turn); `tpt-aero template rotor|aircraft`.
- `propeller_flight` example: BEMT propeller thrust and torque reaction in the 6-DOF run, with a fuselage and wing parasite drag build-up and the fuselage's Munk lift and pitching moment.
- `sixdof::simulate_with_loads`: additive force and moment (aircraft axes) at the centre of gravity, for loads the VLM does not carry. Checked: a constant roll moment gives `M / I_xx` acceleration.
- `SlenderBody::munk_lift_and_moment`: first-order slender-body lift and moment in incidence (`ρ V α dS/dx`).
- `SlenderBody::munk_side_force_and_yaw`: first-order side force and yawing moment at sideslip.
- Published-data comparisons (results recorded in `docs/VALIDATION.md`): Caradonna-Tung hover thrust (BEMT over-predicts by about 23 %) and Navion wing-tail lift slope and pitch stiffness from NASA TN D-6643 (C_Lα +15 %; C_mα 1.8× the flight value, attributed to the missing fuselage).
- Benches for the rotor (`bemt`) and panel (`panel`) crates.
- mdBook sources in `docs/book` and a cargo-generate template in `templates/tpt-aero-analysis`.
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
- `tpt-aero-rotor` did not build without `std` (float methods were not in scope); it now builds for thumbv6m.
- Feature wiring: `libm`, `serde` and the umbrella sub-crate features now imply `alloc`; umbrella `verify` enables the verify crate's `std`. The full feature powerset (52 combinations) now compiles.
- **`solve_state` moments and rates are now about stability axes** (body axes rotated by alpha), the convention of stability derivatives. Roll and yaw previously mixed by alpha; found by comparing aileron adverse yaw against lifting-line theory (the sign was wrong in body axes). Lateral derivatives change slightly at non-zero alpha.
- **`WingSection::twist` was ignored by the VLM solvers** (only the per-solve `twists` slice was applied). Effective twist is now section twist plus the slice.
- **`WingPlanform::tapered` built a skewed, non-symmetric wing** (one segment from the left tip to the right tip), so any taper, sweep or dihedral broke mirror symmetry. It now builds three sections (tip, root, tip) with V-dihedral. Results for tapered, swept or dihedral wings change.
- Longitudinal plant: control-matrix area double count, lift control row, missing `C_Lq` and gravity terms.
- BEMT: reported iteration count; thrust no longer integrated inside the hub cutout.
- Dynamic inflow time constant was off by a factor of the tip speed.
- Frame-convention doc in `tpt-aero-core` now matches the solver.
- CI checks out `tpt-math`.
