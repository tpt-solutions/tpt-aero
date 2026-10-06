# tpt-aero — Phased Build Tracking

Tracks all tasks for the tpt-aero pure-Rust aerodynamics/flight-mechanics library, per [spec.txt](spec.txt). Dual-licensed MIT OR Apache-2.0 (TPT Solutions); no Apache-2.0-only runtime dependencies (Kani is an approved dev-only exception).

## Phase 0 — Project Setup & Licensing

- [x] Initialize git repository
- [x] Add `LICENSE-MIT` and `LICENSE-APACHE` (copyright TPT Solutions 2026)
- [x] Workspace `Cargo.toml` (resolver "2"; workspace package fields: `edition = "2021"`, `rust-version = "1.84"`, `license = "MIT OR Apache-2.0"`, `authors = ["TPT Solutions"]`)
- [x] `rust-toolchain.toml` (stable channel, `rustfmt` + `clippy` components, no_std targets e.g. `thumbv6m-none-eabi`, `wasm32-unknown-unknown`)
- [x] `deny.toml` (cargo-deny license policy: allow MIT / MIT-OR-Apache-2.0 / dual-licensed only; reject Apache-2.0-only crates)
- [x] Workspace lints: `unsafe_code = "forbid"`, `clippy.all = "warn"`
- [x] CI workflow `.github/workflows/ci.yml`: fmt check, clippy `-D warnings`, tests via `cargo-nextest` + doctests, no_std build check, `cargo-deny check`, feature-powerset via `cargo-hack`, bench-smoke (`cargo bench --no-run`)
- [x] README skeleton (overview, crate map, license badges)
- [x] Scaffold all 8 crate directories under `crates/` with empty `lib.rs` + per-crate `Cargo.toml` (`.workspace = true` inheritance):
  - [x] `tpt-aero-core`
  - [x] `tpt-aero-vlm`
  - [x] `tpt-aero-panel`
  - [x] `tpt-aero-rotor`
  - [x] `tpt-aero-atmo`
  - [x] `tpt-aero-stability`
  - [x] `tpt-aero-verify`
  - [x] `tpt-aero` (umbrella)
- [x] Wire up `tpt-math` as local path dependency for `tpt-aero-core`

## Phase 1 — Foundation & VLM (Weeks 1-3)

### tpt-aero-core
- [x] Unit-safe wrapper types: Dynamic Pressure (Pa), Mach number, Reynolds number, Angle of Attack (rad)
- [x] Phantom-type distinction between freestream conditions, local flow velocities, aerodynamic forces
- [x] Non-dimensional coefficient types (CL, CD, Cm) strictly separated from dimensional forces (N)
- [x] Airfoil polar interpolation (Cl/Cd/Cm vs Alpha/Mach/Re)
- [x] Stall modeling
- [x] Geometry primitives: lifting-line panels
- [x] Geometry primitives: quadrilateral surface panels
- [x] Airfoil coordinate transformations
- [x] `no_std` + `alloc` compatibility verified (CI no_std build target)

### tpt-aero-vlm
- [x] From-scratch horseshoe vortex lattice implementation
- [x] AIC (Aerodynamic Influence Coefficient) matrix assembly
- [x] Native LU solve (via tpt-math)
- [x] Native Cholesky solve (via tpt-math)
- [x] Camber modeling
- [x] Twist (washout) modeling
- [x] Dihedral modeling
- [x] Differentiable AIC assembly (gradient-based wing shape optimization, e.g. induced-drag minimization via twist distribution)

### tpt-aero-verify (initial)
- [x] proptest strategies for generating valid 3D lifting geometries
- [x] Property test: symmetric geometries at zero AoA produce zero rolling/yawing moments (VLM)
- [x] Property test: reversing freestream velocity vector reverses the lift vector

### Phase 1 exit criteria
- [x] VLM symmetry and zero-lift conditions verified via proptest

## Phase 2 — Rotors & Atmosphere (Weeks 4-6)

### tpt-aero-rotor
- [x] BEMT (Blade Element Momentum Theory) core model
- [x] Iterative induction factor solver
- [x] Prandtl tip loss model
- [x] Glauert tip loss model
- [x] Hub loss corrections
- [x] Skewed wake model
- [x] Dynamic inflow model (helicopter/VTOL UAV rotors)
- [x] Propeller performance mapping: thrust vs advance ratio
- [x] Propeller performance mapping: power vs advance ratio
- [x] Propeller performance mapping: efficiency vs advance ratio

### tpt-aero-atmo
- [x] ISA (International Standard Atmosphere) model: density, temperature, pressure, speed of sound vs altitude
- [x] Wind shear model
- [x] Gust models: 1-cos
- [x] Gust models: sharp-edge
- [x] Dryden continuous turbulence model
- [x] von Kármán continuous turbulence model

### Phase 2 exit criteria
- [x] Integration test: propeller performance map matches analytical momentum theory limits
- [x] Property test: total mass flow conserved through actuator disk in BEMT

## Phase 3 — Panel Methods & Stability (Weeks 7-9)

### tpt-aero-panel
- [x] Source panel method
- [x] Doublet panel method
- [x] Slender body theory for rapid fuselage aerodynamics
- [x] Prandtl-Glauert compressibility correction
- [x] Karman-Tsien compressibility correction

### tpt-aero-stability
- [x] Finite-difference stability derivative extraction
- [x] Analytical stability derivative extraction
- [x] Derivatives: CLα, Cmq, Cnβ (and other standard longitudinal/lateral-directional derivatives)
- [x] Control surface effectiveness derivatives: CLδe, Cmδa (and other control derivatives)
- [x] Plant-matrix coupling hook for future `tpt-control` integration

### tpt-aero-verify (complete)
- [x] Kani harness: BEMT induction factors remain within physical bounds (0 ≤ a ≤ 1)
- [x] Kani harness: BEMT induction factors do not diverge / satisfy Glauert's condition
- [x] Kani proof: VLM AIC matrices do not contain singularities for valid geometries
- [x] Invariant check: symmetric geometries at zero AoA produce zero rolling/yawing moments (panel methods)

## Phase 4 — Umbrella Crate, Docs & Release

- [ ] `tpt-aero` umbrella crate: feature-gated re-exports of all sub-crates
- [ ] Full rustdoc pass across all crates
- [ ] Doctested usage examples per crate
- [ ] `examples/` directory: sample wing VLM analysis
- [ ] `examples/` directory: sample propeller BEMT sweep
- [ ] Final `cargo-deny` / license audit pass (confirm zero Apache-2.0-only runtime deps)
- [ ] Version 0.1.0 tagging
- [ ] Publish readiness check (crates.io or internal registry, per house convention)

## Future Integration (tracked, not scheduled)

- [ ] `tpt-flight-control` / `tpt-multibody-dynamics`: consume `tpt-aero-stability` linearized plant matrices + `tpt-aero-atmo` disturbances for 6-DOF equations of motion
- [ ] `tpt-physics` (CFD): define handoff boundary for high-AoA/separated-flow regimes where VLM/Panel methods fail
- [ ] `tpt-systems-optimisation`: gradient-based MDO using differentiable VLM/BEMT
- [ ] `tpt-fem-modal`: aeroelastic coupling (flutter analysis) via unsteady aerodynamics + structural mode shapes
