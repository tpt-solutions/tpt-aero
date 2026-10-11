# Changelog

All notable changes to `tpt-aero-verify` are documented here. Workspace-wide
notes are in the root [CHANGELOG.md](../../CHANGELOG.md).

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `trim_request` strategy (type `TrimRequestParams`, with `density()` and `request()`) for level-flight trim, with property tests: generated requests are never rejected as invalid, and any returned trim satisfies the lift residual and the ±20° envelope.
- `section_polar`, `flight_condition` and `control_deflection` strategies, with property tests: polar lookups stay within the table's range and match the linear slope, ISA density falls with altitude, and sampled flight conditions are inside the ISA range.

### Added
- Example `property_sweep`: using the proptest strategies outside a test suite.
- proptest strategies: `symmetric_wing`, `angle_of_attack`, `freestream_speed`, `lattice_resolution`.
- Physical-invariant property tests: symmetry gives zero roll and yaw, freestream reversal reverses lift, and dynamic-pressure scaling.
- Panel symmetry property (zero lift for a symmetric section at zero incidence).
- Kani proof harnesses (`cfg(kani)`) for BEMT induction bounds and linear-solve soundness. Run with `cargo kani`.
- A persisted proptest regression file (`tests/properties.proptest-regressions`) so previously found failures are replayed.

### Fixed
- Initial release candidate; no fixes recorded against this crate yet.
