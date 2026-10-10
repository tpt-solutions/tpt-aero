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
- [x] Derivatives: CLα, Cmα, CLq, Cmq (longitudinal)
- [x] Lateral-directional derivatives: Cyβ, Clβ, Cnβ, Cyp, Clp, Cnp, Cyr, Clr, Cnr (wing-alone, via sideslip/rates in `solve_state`; Cmq cross-checked against the incidence method to 0.1%)
- [x] Control surface effectiveness derivatives: CLδ, Cmδ (longitudinal flap/elevator)
- [ ] Aileron/rudder derivatives: Clδa, Cnδr
- [x] Plant-matrix coupling hook for future `tpt-control` integration

### tpt-aero-verify (complete)
- [x] Kani harness: BEMT induction factors remain within physical bounds (0 ≤ a ≤ 1)
- [x] Kani harness: BEMT induction factors do not diverge / satisfy Glauert's condition
- [x] Kani proof: VLM AIC matrices do not contain singularities for valid geometries
- [x] Invariant check: symmetric geometries at zero AoA produce zero rolling/yawing moments (panel methods)

## Phase 4 — Umbrella Crate, Docs & Release

- [x] `tpt-aero` umbrella crate: feature-gated re-exports of all sub-crates
- [x] Full rustdoc pass across all crates (`RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps` clean)
- [x] Doctested usage examples per crate (all 8 crates have at least one passing doctest)
- [x] `examples/` directory: sample wing VLM analysis
- [x] `examples/` directory: sample propeller BEMT sweep
- [x] Final `cargo-deny` / license audit pass (confirm zero Apache-2.0-only runtime deps) — `cargo deny check` passes all categories
- [ ] Version 0.1.0 tagging (version is already `0.1.0` in `Cargo.toml`; tag not yet created)
- [ ] Publish readiness check (crates.io or internal registry, per house convention) — `tpt-aero-core` dry-run passes; umbrella `tpt-aero` requires the sub-crates published first, in dependency order, with the umbrella last

## Follow-up — VLM Oswald factor near zero lift

Context: `e = CL²/(π·AR·CDi)` is unstable where `CL → 0` for the cambered wing (NACA m=0.02, p=0.4; AR 10). Flat-wing `CDi/CL²` is constant (~0.032) across the sweep.

- [x] `examples/wing_vlm.rs`: report Oswald factor only for `|CL| ≥ 0.1`, otherwise `n/a`
- [x] Tried Trefftz-style induced drag (trailing legs only) — no change in output (bound-segment term cancels by reciprocity); reverted
- [x] Tried freestream-aligned wake for the drag term only — inconsistent with the body-aligned solve; reverted
- [x] Fourier-fit the spanwise `strip_cl` at −2° / 0° / 4°: odd modes a3, a5 carry the induced drag at near-zero CL (mode-based e ≈ 0.53 at −2°, ≈ 0.91 at 0°)
- [x] Thin-airfoil zero-lift angle for NACA m=0.02, p=0.4: α_L0 = −2.077° (corrected integral)
- [x] VLM zero-lift angle: −2.171° (32×4), −2.124° (64×8) — the 4-chord-panel value is the outlier; see the convergence items below
- [x] Measure VLM zero-lift angle vs chordwise panel count (n_span = 32): n_chord 4 → −2.171°, 8 → −2.123°, 16 → −2.124°. Chordwise error converges by 8 panels
- [x] Measure VLM zero-lift angle vs spanwise panel count (n_chord = 8): n_span 16 → −2.119°, 32 → −2.123°, 64 → −2.124°. Spanwise resolution is not the driver
- [x] Test freestream-aligned wake for the zero-lift angle: −2.1232° (unchanged). Wake alignment ruled out; solver reverted
- [x] Residual −0.047° offset from thin-airfoil α_L0 (−2.077°): decays ~1/AR (AR 5 → −0.101°, 10 → −0.046°, 20 → −0.021°, 40 → −0.010°). A finite-aspect-ratio property of the horseshoe VLM, not a camber-treatment error. Lifting-line theory predicts no shift for an untwisted planar wing, so this is characterised by trend, not validated against an external reference
- [x] Library remedy: documented (`VlmSolution` near-zero-lift caveat: Oswald ratio ill-conditioned as CL → 0; zero-lift offset ~1/AR). Example guard kept. No formulation change — the Trefftz-style and freestream-wake experiments did not alter results
- [x] Regression test: flat-wing `CDi/CL²` constant across an alpha sweep (`flat_wing_induced_drag_scales_with_cl_squared`, 1% bound; measured spread ≈ 0.2%)

## Review follow-ups (2026-10)

Done:
- [x] CI clones sibling `tpt-math`; examples run in CI
- [x] BEMT: real iteration count, no thrust inside hub cutout
- [x] `DynamicInflow::step` takes tip speed (unit fix)
- [x] Longitudinal plant matrix rebuilt; `short_period` returns `Option`
- [x] VLM rejects non-finite / non-positive flight conditions
- [x] `compress::*` return `Result`; `MachNumber`/`ReynoldsNumber::try_new`; `naca4` no longer panics
- [x] `tpt-aero-cli`: TOML wing cases to CSV, `template`, `isa` subcommands
- [x] Panel: `naca4_closed` (sharp TE, no gap panel), CCW winding check, `doublet.rs` docs corrected, NACA 0012 validated (Cp_min, lift slope)
- [x] VLM: `prepare`/`solve_sweep` (AIC + LU cached across angles, ~20x), local velocities computed once per solve, criterion benches
- [x] VTK export of lattice + circulation (library and CLI)
- [x] Phugoid-capable longitudinal plant (Lanchester cross-check) and control derivatives for any named control
- [x] Level-flight trim solver; lateral plant with Dutch roll / roll subsidence / spiral modes (dihedral stabilizes spiral, trims vs speed/CG behave classically)
- [x] Control surfaces (elevator/aileron/rudder via hinge-line incidence); aileron roll and adverse yaw validated against lifting line (roll within 15%); stability-axis moments
- [x] Multi-surface aircraft (`Surface`, `from_surfaces`): wing + tail + fin; downwash reduces tail lift slope to ~60-70%, NP/Cn_beta/Cn_r physically sensible
- [x] Planform twist now applied by the solvers (was silently ignored)
- [x] Lifting-line validation of tapered wings (found near-field Oswald e up to ~3% above LL; CDi is optimistic, e can exceed 1 at low resolution)
- [x] Fixed `WingPlanform::tapered` (was a skewed asymmetric wing); symmetry regression tests in core and CLI
- [x] Prelude + `isa_table`, `airfoil_cp`, `stability_derivatives` examples

Open:
- [ ] Panel imports `lu_solve` from vlm (layering); move to core / tpt-math
- [ ] VLM: opt-in cosine spanwise spacing; wake along freestream
- [ ] Turn/steady-sideslip trim; full 6-DOF nonlinear time simulation; ailerons/rudder in the lateral plant (B matrices); fuselage and propeller effects
- [ ] Viscous drag, VLM compressibility, ground effect, ISA offset/viscosity
- [ ] `serde` feature on result types, CLI: rotor/stability subcommands, mdBook, cargo-generate templates
- [ ] Validation: AVL / Caradonna-Tung / Navion reference cases; benches for rotor and panel
- [ ] Tag v0.1.0 and publish in dependency order

