# Changelog

All notable changes to `tpt-aero-vlm` are documented here. Workspace-wide
notes are in the root [CHANGELOG.md](../../CHANGELOG.md).

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Breaking
- `VlmError` has two new variants (`FreestreamWakeNotPrepared`, `InvalidGroundHeight`). Exhaustive matches need the new arms.
- `StateLoads` has two new fields (`force_n`, `moment_nm`). Struct literals need them.
- `lu_solve` now lives in `tpt-aero-core::linalg`. It is still re-exported here, so existing paths keep working.

### Added
- Example `wing_polar`: alpha-sweep polar, compared with elliptic-wing lift slope.
- `Surface` (`horizontal`, `vertical`, `.at`, `.with_incidence`, `.with_camber`) and `VlmModel::from_surfaces` for whole-aircraft models with shared downwash and sidewash.
- `Lattice::combine` and `Lattice::transformed`.
- Control surfaces: `Surface::with_control`, `VlmModel::controls`, `control_index`, `control_twists`, and `VlmSolution::stability_roll_yaw`.
- `PreparedVlm` (from `VlmModel::prepare`) and `VlmModel::solve_sweep`: the AIC and LU factors are computed once and each angle costs O(N²). About 20x faster for a 17-angle sweep. `LuFactors` for reusable solves.
- `FlightState` and `PreparedVlm::solve_state`: sideslip and body rates, wind-axis forces, body moments in aircraft-standard signs.
- `prepare_with_twists`: prepared solves with control deflections applied.
- Cosine spanwise spacing (`SpanSpacing::Cosine`, `from_planform_spaced`, `WingPlanform::refine_cosine`).
- `WakeModel::Freestream`: trailing legs aligned with the freestream (per-angle solves only).
- Ground effect by image horseshoes (`with_ground_height`).
- Subsonic compressibility by Göthert's rule (`from_planform_compressible`). Lift slope within 2 % of Helmbold-Diederich.
- `VlmError::InvalidFlightCondition`. `VlmModel::solve` validates angle, speed and density.
- `solve_generic` over dual numbers for gradients.
- `export::to_vtk`: legacy ASCII VTK of the lattice with circulation, surface and strip cell data (ParaView, PyVista).
- Criterion benches (`cargo bench -p tpt-aero-vlm`).
- Optional `serde` feature on solution types.

### Changed
- Local velocities are computed once per solve rather than per control point.

### Fixed
- **`solve_state` moments and rates are about stability axes** (body axes rotated by alpha), the convention of stability derivatives. Roll and yaw previously mixed by alpha. Lateral derivatives change slightly at non-zero alpha.
- **`WingSection::twist` was ignored by the solvers**; only the per-solve `twists` slice was applied. The effective twist is now the section twist plus the slice.
- Non-finite and non-positive flight conditions are rejected rather than producing NaN results.
