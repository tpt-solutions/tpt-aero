# Changelog

All notable changes to `tpt-aero-core` are documented here. Workspace-wide
notes are in the root [CHANGELOG.md](../../CHANGELOG.md).

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- `naca4` rejects nonzero camber with camber position `p = 0` (`InvalidNacaCode`). Such codes were silently flattened to symmetric sections.
- `PolarSet` blends Reynolds numbers in ln(Re) rather than linearly in Re, so tables are weighted by decades. Non-positive Re returns NaN.
- `linalg` tests pin that the cached and generic LU paths share one singularity threshold.

### Breaking
- `naca4` returns `AirfoilError::InvalidResolution` for fewer than 2 stations and accepts thickness up to 40 %.
- Error types (`AirfoilError`, `PolarError`, `GeometryError`, `SolveError`) implement `core::error::Error`, so `?` converts them into `Box<dyn Error>`.
- `lu_solve` moved here from `tpt-aero-panel`. It is still re-exported from `tpt-aero-vlm`.

### Added
- Example `quickstart`: geometry, airfoil, checked quantities, drag build-up and linear solve.
- `linalg` module: `LuFactors` (reusable LU factors), `lu_solve`, `lu_solve_generic`, `cholesky_spd`, `SolveError`.
- `airfoil::naca4_closed`: closed-trailing-edge NACA 4-digit generator for panel solvers.
- `quantities`: `MachNumber::try_new` and `ReynoldsNumber::try_new` return `Option` for invalid input.
- `viscous`: flat-plate skin friction (laminar, Schlichting mixed), Hoerner form factor and a `C_D0` component build-up.
- `polar`: `AirfoilPolar`, `PolarSet`, `StallModel` for tabulated section data.
- `num::Real` trait, so geometry and coefficient code runs over dual numbers.
- Optional `serde` feature on result and parameter types, valid with `no_std` + `alloc`.
- `libm` feature for hosts without hardware floats.

### Changed
- Frame-convention documentation now matches the solver: see `docs/CONVENTIONS.md`.

### Fixed
- `WingPlanform::tapered` built a skewed, non-symmetric wing. It now builds three sections (tip, root, tip) with V-dihedral, so any taper, sweep or dihedral keeps mirror symmetry.
- Feature wiring: `libm` and `serde` imply `alloc`. The full feature powerset compiles.
