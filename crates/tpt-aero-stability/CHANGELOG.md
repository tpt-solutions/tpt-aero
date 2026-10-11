# Changelog

All notable changes to `tpt-aero-stability` are documented here. Workspace-wide
notes are in the root [CHANGELOG.md](../../CHANGELOG.md).

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `alpha_beta_envelope` and `LINEAR_ENVELOPE_LIMIT` (20°): flag air-relative |α| or |β| outside the linear VLM envelope.
- `simulate_partial_with_wind`: the six-DOF run takes the NED air velocity as a function of time, for gusts. `simulate_partial` delegates to it with the steady `cfg.wind_ned`.
- `one_minus_cosine_gust`: a vertical one-minus-cosine updraft for that function, returning `None` for invalid inputs.

### Changed
- Six-DOF builds the VLM influence kernel once and reuses it on every control change.
- Turn trim tests its force residuals relative to W/(qS) rather than absolutely.
- `trim` Newton step uses a squared norm for its damping comparison (same decisions; `no_std`).

### Breaking
- `LongitudinalPlant::short_period` returns `Option<(f64, f64)>`: `None` when the aircraft is statically unstable.

### Added
- Example `short_period`: wing to longitudinal derivatives to short-period frequency and damping.
- `longitudinal_fd`: `C_Lα`, `C_mα`, `C_Lq`, `C_mq`, `C_Lδ`, `C_mδ` by finite difference, plus `drag_alpha_fd` and a thin-airfoil reference helper.
- `cl_alpha_dual`: lift slope from dual numbers.
- `lateral_fd`: `C_Yβ`, `C_lβ`, `C_nβ`, `C_Yp`, `C_lp`, `C_np`, `C_Yr`, `C_lr`, `C_nr`.
- `LateralPlant`, `LateralParameters`, `LateralModes`: lateral-directional state-space plant `[β, p, r, φ]` with Dutch-roll, roll-subsidence and spiral modes. `modes::eigenvalues4` is a dependency-free 4x4 eigenvalue routine.
- `controls::control_derivatives`: lift, drag, side force and the three moments per radian of any named control.
- `trim::trim_level_flight` with `TrimRequest`: Newton trim of angle of attack and elevator, reporting thrust required.
- `turn`: nonlinear steady turn and steady sideslip trim (`trim_steady_turn`, `coordinated_turn_rate`), solved by continuation from level flight.
- `sixdof`: nonlinear six-DOF simulation in aircraft axes with Euler angles, VLM aerodynamics, thrust along the airspeed vector and parasite `C_D0`. `simulate_with_loads` adds external force and moment at the centre of gravity.
- `phugoid::FullLongitudinalPlant`: four-state longitudinal plant with speed and drag couplings, giving short-period and phugoid modes.
- `sim`: RK4 `StateSpace` simulation, with `FullLongitudinalPlant::state_space` (elevator and vertical gust inputs) and `LateralPlant::state_space` (aileron, rudder and sidewash inputs).
- `StateSpace::with_feedback` and `eigenvalues` for closed-loop analysis.
- `StateSpace::with_integrator`: integral action (autotrim) for attitude holds.
- AVL 3.52 comparison test (`tests/avl_reference.rs`) for `C_Lα`, `C_Dᵢ`, `C_Yβ`, `C_lβ`, `C_nβ`, `C_lp`, `C_lr` and `C_nr`.
- Navion comparison against NASA TN D-6643 (`tests/navion.rs`).
- Optional `serde` feature on derivative sets, trim results and related types.

### Fixed
- Longitudinal plant: control-matrix area double count, a wrong lift control row, and missing `C_Lq` and gravity terms. The plant matrix is rebuilt.
- Stability derivatives propagate solve errors from the VLM instead of returning meaningless values.
- `C_nδr` matches the fin side force times its moment arm to 1 %.
