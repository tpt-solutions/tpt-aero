# Changelog

All notable changes to `tpt-aero-rotor` are documented here. Workspace-wide
notes are in the root [CHANGELOG.md](../../CHANGELOG.md).

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `BemtSolver::with_swirl(bool)` and `swirl()`: solve the tangential induction `a'` per annulus. Off by default so the validated hover numbers do not move.

### Changed
- `CylinderWake` uses an adaptive Simpson rule (to 1e-9 relative, capped at 4096 panels) rather than a fixed 4001-point rule.

### Breaking
- `SectionResult` has a new `swirl` field (zero unless swirl is enabled). Struct-literal construction needs it.

### Breaking
- `DynamicInflow::step` takes an extra `tip_speed` argument (m/s). The time constant was previously off by a factor of the tip speed.

### Added
- Example `hover_bemt`: hover point, thrust in newtons, and advance-ratio map.
- `CylinderWake`: full induced-velocity field of a uniformly loaded actuator disk, from a linearised vortex-cylinder wake (ring strength `2 v_d`). Closed-form ring kernel with complete elliptic integrals; exact on the axis, divergence-free and uniform over the disk.
- `Slipstream`: momentum-theory propeller slipstream (induced velocity at the disk, on-axis axial increment, contracting tube with an edge blend, dynamic-pressure increment at a point).
- `SlipstreamModel`: choice between the contracting tube and the cylinder field, with a common `local_flow` and `dynamic_pressure_increment`.
- Example `slipstream_tail`: propeller slipstream over a horizontal tail.
- `performance_map`: advance-ratio sweep of thrust, torque and efficiency.
- `BemtSolver::with_iteration_controls` for tuning the induction iteration.
- Criterion benchmark (`benches/bemt.rs`).
- Published-data comparison against Caradonna-Tung hover (results in `docs/VALIDATION.md`).
- Optional `serde` feature on parameter and result types.

### Fixed
- The crate did not build without `std` (float methods were not in scope). It now builds for `thumbv6m`.
- BEMT reported the wrong iteration count. It now reports the real count.
- Thrust is no longer integrated inside the hub cutout.
- Efficiency and figure of merit are zero once thrust reverses, rather than meaningless.
- `DynamicInflow` time constant corrected (unit factor of the tip speed).
