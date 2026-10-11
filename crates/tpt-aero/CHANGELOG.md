# Changelog

All notable changes to the `tpt-aero` umbrella crate are documented here.
Changes to the sub-crates are recorded in their own changelogs; workspace-wide
notes are in the root [CHANGELOG.md](../../CHANGELOG.md).

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `prelude` module: one-line import of the common types (`use tpt_aero::prelude::*;`).
- Optional `serde` feature forwarded to the enabled sub-crates, with round-trip tests.
- `verify` feature (opt-in): re-exports `tpt-aero-verify` as `tpt_aero::verify`.
- Runnable examples: `aircraft_stability`, `gust_response`, `control_design`, `wing_vlm`, `airfoil_cp`, `propeller_sweep`, `propeller_flight`, `stability_derivatives`, `isa_table`.
- Whole-aircraft doctest and single-wing doctest.
- Sub-crates re-exported as `aero_core`, `vlm`, `panel`, `rotor`, `atmo`, `stability`.

### Changed
- Per-crate examples now cover every solver crate; see each crate's README for its `examples/` programs.
- `std`, `alloc` and `libm` are forwarded to every enabled sub-crate.

### Fixed
- Feature wiring: `libm`, `serde` and the umbrella sub-crate features imply `alloc`. `verify` enables the verify crate's `std`. The full feature powerset compiles.
