# Changelog

All notable changes to `tpt-aero-atmo` are documented here. Workspace-wide
notes are in the root [CHANGELOG.md](../../CHANGELOG.md).

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Example `atmosphere_gusts`: standard and non-standard atmosphere, power-law shear and 1-cos gust profile.
- `isa::at`: International Standard Atmosphere from 0 to 84 852 m geopotential, returning `IsaState`.
- `isa::at_isa_offset`: non-standard temperature day, hydrostatically consistent with the standard atmosphere.
- Sutherland `dynamic_viscosity` and `kinematic_viscosity`.
- `shear::power_law` and `shear::log_law` wind profiles.
- `shear::DiscreteGust` with `GustShape` (1-cos and sharp-edged), `speed_at` and `gradient_at`.
- `turbulence`: Dryden and von Kármán `TurbulenceModel`, `TurbulenceAxis`, `TurbulenceField` with low- and medium-altitude presets, and `simulate_axis`.
- `turbulence::WhiteNoise`: a reproducible `no_std` noise source (uniform and Gaussian), so simulations repeat exactly for a given seed.
- Optional `serde` feature on state and parameter types.

### Fixed
- Initial release candidate; no fixes recorded against this crate yet.
