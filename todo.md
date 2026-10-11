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
- [x] Aileron/rudder derivatives: Clδa, Cnδr (Cnδr matches the fin side force times its moment arm to 1 %, test `rudder_yaw_moment_matches_fin_side_force_times_arm`; Clδa validated against lifting line, see Review follow-ups)
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
- [ ] Version 0.1.0 tagging (version is already `0.1.0` in `Cargo.toml`; tag not yet created) *Not in this pass:* publishing, excluded by request.
- [x] Publish readiness check: `cargo publish --workspace --dry-run --allow-dirty` passes for all 9 crates (checked after the 2026-10 changes; nothing uploaded). The actual publish stays open, in dependency order with the umbrella last. If the dry run reports a missing `linalg` or similar, clear the stale verification build with `cargo clean -p <crate>` and remove the `-<hash>` directory under `~/.cargo/registry/src`

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
- [x] Closed-loop state feedback (`with_feedback`): yaw damper and pitch-attitude hold examples
- [x] Linear time-domain simulation (RK4) with elevator/aileron/rudder and gust inputs; FAR 25 gust alleviation check
- [x] `serde` feature on result and parameter types (no_std + alloc); `derivatives --json`; fixed rotor no_std build and the feature powerset
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
- [x] Panel imports `lu_solve` from vlm (layering); moved to `tpt-aero-core::linalg`, re-exported from vlm
- [x] VLM: opt-in cosine spanwise spacing (`SpanSpacing::Cosine`); opt-in freestream wake (`WakeModel::Freestream`, per-angle solves only)
- [x] Turn/steady-sideslip trim (`trim_steady_turn`: coordinated turn and steady sideslip, continuation from level flight)
- [x] Nonlinear 6-DOF simulation (`sixdof::simulate`, aircraft axes, VLM aerodynamics, Euler-angle RK4); checked against the level trim and the turn trim
- [x] Integrator/autotrim in the attitude hold (`with_integrator`, integral error test)
- [x] Fuselage and propeller effects, part 1. Done: fuselage parasite drag, Munk lift, side force and moments (`SlenderBody::munk_lift_and_moment`, `munk_side_force_and_yaw`), propeller thrust and torque reaction from BEMT, all through `sixdof::simulate_with_loads` (`propeller_flight` example).
- [x] Propeller slipstream and its effect on the tail: two models, selectable. `SlipstreamModel::Tube` (momentum contracting tube, closed-form on-axis) and `SlipstreamModel::Cylinder` (`CylinderWake`, linearized vortex-cylinder field with full crossflow, exact on the axis and divergence-free). The `propeller_flight` example uses the cylinder by default and prints both at trim. Limits: linear (no wake contraction), tail angle of attack taken as the wing's (no downwash), and the VLM does not see the slipstream.
- [x] Viscous drag (`tpt-aero-core::viscous`: flat-plate skin friction, Hoerner form factor, C_D0 build-up)
- [x] VLM compressibility (Göthert rule, `from_planform_compressible`; matches Helmbold-Diederich lift slope ratio within 2 %)
- [x] Ground effect (image horseshoes, `with_ground_height`)
- [x] ISA offset and Sutherland viscosity (`at_isa_offset`, `dynamic_viscosity`, `kinematic_viscosity`)
- [x] CLI: `rotor` and `trim` subcommands (with `template rotor|aircraft`)
- [x] mdBook sources in `docs/book` (builds with mdbook 0.5)
- [x] cargo-generate template in `templates/tpt-aero-analysis` (verified by substituting the placeholder and building; cargo-generate itself was not run, and the template depends on a published `tpt-aero`)
- [x] Benches for rotor (`bemt`) and panel (`panel`)
- [x] Validation against published data: Caradonna-Tung hover `C_T` (BEMT +23 %, recorded and bounded) and the Navion wing-tail `C_Lα` and `C_mα` from NASA TN D-6643 (recorded with bands; the moment gap is the missing fuselage). The tail arm is read from the report's figure. Sources and results are in `docs/VALIDATION.md`
- [x] AVL 3.52 comparison on matched geometry (built from source with gfortran; decks and run script in `docs/validation/avl/`; `tpt-aero-stability/tests/avl_reference.rs`): wing alone, dihedral wing, refined T-tail. Agreement 1–4 % on the converged grid. Found and recorded: the VLM converges first order in spanwise strips, and the CLI template's tail-on-fin-root layout gives sideslip derivatives ~30 % off AVL
- [ ] Validation still open: XFLR5 comparison (not run), swept / cambered / control-surface configurations against AVL, and measured data for the full six-DOF and compressible results
- [x] `cargo publish --workspace --dry-run` passes for all 9 crates (dependency order handled by cargo)
- [x] README rewritten; real per-crate READMEs; `docs/CONVENTIONS.md` and `docs/VALIDATION.md`
- Note: repeated `cargo publish --workspace --dry-run` runs reuse an unpacked copy of the packaged crates in Cargo's registry cache (a `-<hash>` directory under `registry/src`, named after the temporary registry). After changing a crate without bumping its version, delete that directory first or verification builds against stale sources (this produced a false CLI failure once).
- [ ] Tag v0.1.0 and publish (needs a decision: bump to 0.1.0 final changelog entry, push, confirm GitHub CI is green first) *Not in this pass:* publishing, excluded by request.

## Platform review (2026-10-11)

Source: three read-only audits of all crates. Items marked (confirmed) were spot-checked in the code; the rest are auditor-reported and need a failing test before the fix.

### 1. Bugs (each fix lands with a regression test)
- [x] Panel `slender.rs`: Munk lift was `rho*V*alpha*dS/dx`, now `rho*V^2*alpha*dS/dx` (confirmed); test and docs corrected
- [x] Rotor `bemt.rs`: thrust element used `cl*cos(phi) + cd*sin(phi)`, now `- cd*sin(phi)` at both sites (confirmed). Still to add: a high-advance-ratio test against the momentum limit
- [x] CLI `aircraft.rs`: garbled run of spaces in the missing-control error (confirmed)
- [x] Stability `longitudinal.rs`: seven `unwrap_or(f64::NAN)` sites turn failed solves into NaN; now return `Result<_, VlmError>` (`longitudinal_fd`, `cl_alpha_dual`), and CLI `derivatives` rejects non-finite output (test `failed_solve_is_an_error_not_nan`)
- [x] CLI `trim --bank`: check `|bank| < 90 deg` before `coordinated_turn_rate`; clear message (test `turn_rejects_bank_at_or_past_90_degrees`)
- [x] CLI control geometry: validate inboard < outboard, within span, no duplicate names, and that the control covers at least one element (now a misleading `Singular`); warn on hinge quantisation to chord-panel boundaries (Inboard < outboard ≤ span/2, at least one strip covered, unique names; template rudder corrected to the fin half-span.)
- [x] CLI mesh-size cap on `n_span`, `n_chord`, `stations` (each ≤ 500, panel product ≤ 2500 per surface; test `mesh_size_is_capped`)
- [x] 6-DOF: guard gimbal lock near theta = +-90 deg; return partial history on error like the linear sim; check `ixx*izz > ixz^2` in `Aero::new`; invert the inertia tensor once; cap `steps` allocation **Partly done:** `GimbalLock` error near θ = ±90°, closed-form inertia inverse, `Ixx·Izz > Ixz²` check, `simulate_partial` keeps history on error, capped allocation. Still open: AIC refactored on every control change (performance), thrust along airspeed (documented only), and the alpha/beta envelope flag. (GimbalLock error, closed-form inertia inverse, Ixx·Izz > Ixz² check, simulate_partial keeps history, allocation capped.)
- [ ] 6-DOF: AIC is refactored whenever `Controls` differs bitwise, up to 4 LU factorisations per RK4 step; quantise or linearise the control effect **Partly done:** the horseshoe kernel is now built once per model (`VlmModel::influence_kernel` and `prepare_with_kernel`), so a control change costs the normals and one LU instead of N² horseshoe evaluations. Still open: quantising or linearising the control effect, which would also remove the LU.
- [ ] 6-DOF: thrust is forced along the airspeed vector, not the body x-axis; document or add a body-fixed thrust option; flag alpha/beta outside the linear-VLM envelope **Partly done:** `thrust_body_fixed` and a steady `wind_ned` are in. **Update:** `alpha_beta_envelope` flags |α| or |β| beyond 20° for each sample. Still open: the flag as a per-sample history field, and atmosphere and gust inputs.
- [x] VLM `solve_generic` input validation: done (primal check of alpha, speed, density in `solve_core`). Still open: `solve_state` at |beta| = 90 deg normalises a zero vector; check reference area, chord and span **Partly done:** primal checks on alpha, speed and density. Still open: `solve_state` at |β| = 90° and reference geometry checks. (Validation, the |β| = 90° lift direction, and the reference checks are all in.)
- [x] VLM `from_planform_compressible`: `with_reference` double-counts beta (beta^2 S); add `from_surfaces_compressible` (`with_reference` stretches a physical reference; `from_surfaces_compressible` added.)
- [x] VLM `from_lattice`: span counts only the first surface, so chord = area/span is wrong for multi-surface lattices (Reference area, span and chord come from the first surface; the span was n_chord times too large.)
- [x] VLM/core `refine`: uniform refinement drops interior section nodes (kinks fall inside a strip, area differs from `planform.area()`) (Every original section is kept as a station; the strip count must cover the planform's segments.)
- [x] VLM `with_incidence` silently ignores a `from_sections` failure (now returns `Result`; rejects non-finite incidence)
- [ ] VLM: no vortex core, absolute 1e-12 singularity guards; coplanar tail-in-wing-plane conditioning **Partly done:** segment guards now scale with segment length (test `segment_velocity_is_length_unit_invariant`). Still open: vortex core and coplanar conditioning; a vortex core would move the AVL-matched results, so it waits for an AVL re-run
- [ ] VLM: `solve()` and `prepare().solve()` use different singularity criteria; `lu_solve_generic` misses NaN pivots **Partly done:** `lu_solve_generic` now rejects NaN pivots and non-finite input. Still open: the two solve paths use different singularity criteria. **Partly done:** both paths now use one factorization and one NaN rule. **Done:** the two paths share one pivot criterion (test `cached_and_generic_lu_share_singularity_threshold`). Vortex core and coplanar conditioning are still open.
- [x] VLM: ground plane intersecting the geometry is accepted; add a minimum clearance check (Clearance check added: bound and collocation points must lie above z = −h.)
- [x] VLM NACA camber: `p = 1` divides by zero (VLM part: `LatticeError::InvalidCamber`); core NACA parser dropped camber for `p = 0` silently (now `InvalidNacaCode` for any nonzero camber with `p = 0`; test `camber_with_zero_position_is_rejected_not_flattened`)
- [ ] Rotor BEMT: windmill and descent floored to zero induction though the docs claim them; no swirl (a'); bisection assumes monotone `g`; `RotorError::NotConverged` unreachable; dead `_cl_slope`/`_relaxation` parameters; wrong error variant for bad operating conditions **Partly done:** descent bracket now keeps λ_tot ≥ 0; NotConverged is reachable; dead parameters removed; operating-condition error variant added. Still open: windmill-brake state, swirl (a'). **Partly done:** swirl (a') is solved when built `with_swirl(true)`, opt-in so the validated hover numbers do not move (tests in `bemt_validation.rs`). Still open: the windmill-brake state, which needs an empirical turbulent-wake closure that cannot be validated here.
- [x] Rotor `SkewedWake`: factor is `1 + tan(chi)*x*cos(psi)`, Coleman/Glauert first harmonic is about `(4/3)tan(chi/2)`; can go negative (Coleman `(4/3)·tan(χ/2)`, clamped at zero.)
- [ ] Rotor `CylinderWake`: NaN at r = R, 4001 elliptic evaluations per query; `Slipstream` NaN at T = 0, V = 0 **Partly done:** NaN on the wake edge and `Slipstream` NaN at T = V = 0 fixed. Still open: the 4001-point quadrature per query. **Done:** the quadrature is adaptive (doubling Simpson to 1e-9 relative agreement, capped at 4096 panels; test `adaptive_integral_matches_the_fixed_fine_grid`).
- [x] Rotor `RotorGeometry::rectangular` with hub radius 0 produces an invalid first station; `validate()` misses NaN and `r_frac >= hub_radius/tip_radius`
- [ ] Atmo turbulence: low-altitude scale lengths use feet coefficients on metres; `vertical` argument discarded; Dryden v/w need the second-order filter; von Karman is only a tau scaling; Euler step unstable for dt >= 2 tau (use the exact discrete update); validate sigma and length **Partly done:** feet and exponent fixed, vertical flag used, exact first-order update, second-order Dryden v/w. Still open: von Kármán is a first-order approximation, not the spectrum.
- [x] Atmo ISA: 84.852-86 km mixes geopotential and geometric altitude; `altitude_from_pressure` range; `dynamic_viscosity(T <= 0)`; unbounded `at_isa_offset` (Range is 0 to 84 852 m geopotential; ΔISA bounded; `altitude_from_pressure` bounded; viscosity NaN for T ≤ 0.)
- [x] Atmo `shear.rs` doc: peak gradient is at x = h/2, not h
- [ ] Panel: `doublet` accepts NaN, zero chord or zero `v_inf`; centroid vortex can fall outside thin or cambered sections; `source.rs` comment says 8-point Gauss, table has 4 **Partly done:** NaN, zero chord and zero `v_inf` are rejected and the `source.rs` Gauss comment is corrected. Still open: centroid vortex placement on thin or cambered sections. **Done:** the bound circulation is a uniform vortex sheet on the contour, not a centroid vortex (see the panel CHANGELOG), so it stays on the surface for thin and cambered sections.
- [x] Panel `compress::karman_tsien` doc claims weaker suction, code and tests give stronger (doc corrected to match the formula)
- [x] Stability: Durand-Kerner has no convergence flag and burns 2000 iterations on repeated roots (replace with Hessenberg QR); mode classification by exact `im == 0.0`; phugoid/short-period return `None` when overdamped (return the eigenvalues) (Hessenberg reduction and Francis double-shift QR; a failed iteration returns NaN; repeated and defective roots tested.)
- [ ] Stability: `LongitudinalPlant::from_derivatives` has no validation; two longitudinal plants with overlapping names; `with_integrator`/`with_feedback` panic on bad indices; trim Newton has no damping or CL_max check; turn trim uses absolute tolerances and fixed continuation **Partly done:** `try_from_derivatives` validates inputs; `from_derivatives` is unchanged and documented. Panics on bad indices are documented under `# Panics`. Overlapping plant names: not applicable, since `LongitudinalPlant` carries no name field in the current code. **Partly done:** `try_from_derivatives` validates. The trim Newton items: the level trim has a damped step, an envelope check and a relative tolerance; the turn trim has damping and adaptive continuation and now tests its force residuals relative to W/(qS).
- [ ] Stability: `cl_delta`/`cm_delta` use a uniform incidence change, not the named control **Partly done:** the functions are uniform by design and now say so in their docs; `control_derivatives` gives a named surface.
- [x] Gravity defined four times (`turn.rs`, `plant.rs`, `phugoid.rs`, `lateral.rs`); keep one (One `STANDARD_GRAVITY` at the crate root; `turn::G` re-exports it.)
- [ ] Core: `AirfoilPolar::new` accepts NaN; `StallModel` raises lift toward the plateau for already-stalled tabulated polars; `PolarSet` blends linearly in Re and allocates per call; `MachNumber::new`/`ReynoldsNumber::new` panic and no solver uses them; `viscous` has no validation; `cholesky_spd` passes NaN **Partly done:** NaN and invalid Mach or Reynolds are rejected; viscous functions return NaN for invalid input; `cholesky_spd` rejects NaN. Still open: `StallModel` plateau, `PolarSet` blending, `MachNumber::new` panic, and the core NACA `p = 0` parser. **Done:** the NACA `p = 0` parser rejects camber; the `StallModel` plateau is clamped to the raw lift; `PolarSet` blends in ln(Re) and does not allocate per lookup. Decision: `MachNumber::new` keeps its documented panic; `try_new` is the checked constructor and new code should use it.

### 2. CLI usability
- [x] Adopt `clap`: per-subcommand `--help`, `--version`, `--out`/`-o` on every command, duplicate-flag errors, flags in any position (Derive API; `--help` and `--version` per command; `-o/--out` on every command.)
- [x] Display, not `{e:?}`, for errors; include file path and valid range (ISA 0-86 km) (CLI errors use Display and name the file; ISA range reported.)
- [x] Alpha sweep prints `0.30000000000000004`; make the point-count tolerance proportional to the step (Sweep alphas rounded to 1e-9°.)
- [x] Aircraft schema: CG, inertia, reference point; serde enums for `kind` and `symmetry` (CG, inertias, measured or estimated; serde enums for kind and symmetry.)
- [x] Rotor: negative speed (descent), hub loss, polar input, per-station chord and twist, `rotor sweep` over rpm, J or collective; fix column alignment **Partly done:** descent, hub loss, polar input and a taper ratio are in, plus an RPM sweep. **Update:** the J (axial-speed) and collective sweeps are in (`--speed-sweep`, `--collective-sweep`). **Done:** per-station chord and twist arrays (`chord_profile`, `twist_profile_deg`).
- [ ] Output: units column or header notes (per q-hat, per pb/2V), CL, CD, L/D, static margin, CL_max margin in the trim report; spanwise loading from `wing`; SVG polar; CSV for `trim` **Partly done:** trim report now has CL, CD, L/D, static margin and CL_max margin, and the wing has an SVG polar. **Update:** the trim report has a units column and `trim --csv`; `wing --spanwise` gives the strip loading. **Done:** units for every CLI output are in the table in `docs/CONVENTIONS.md`.
- [x] Mach correction and `cd0` option in `wing`; warn on large alpha and on aspect ratio, taper and sweep extremes; bound camber `m` (Mach selects the Göthert model; `cd0` adds columns; warnings; camber bounded.)
- [x] Imperial unit layer and `isa_offset` (decision: no full imperial layer; SI throughout, with the units table in `docs/CONVENTIONS.md`; `isa --imperial` and `at_isa_offset` remain) **Partly done:** ISA only (`isa --imperial`). Trim and wing are still SI.
- [x] Binary-level CLI tests (`assert_cmd`) for argument parsing and exit codes (`tests/binary.rs`: version, help, exit codes, `-o`, negative values.)

### 3. Doc and spec drift
- [ ] spec.txt: license (says MIT only), phases and crate list stale, `tpt-control` coupling claimed but absent, typos **Partly done:** licence line corrected to MIT OR Apache-2.0; `tpt-control` now marked planned. Still open: the phases and crate list. **Partly done:** licence and `tpt-control` are fixed. Done since: phases note and crate list updated in `spec.txt` §3 and §6. Still open: a typo audit of `spec.txt`.
- [ ] todo.md: duplicate v0.1.0 tag items, `[x]` items under "Open:", "differentiable" claim covers the VLM only, Kani proofs not run in CI **Partly done:** not edited in this pass (the tag and publish items are on hold). **Partly done:** not edited in this pass; the tag and publish items are on hold.
- [x] Reconcile `docs/VALIDATION.md` with the "XFLR5 not run" note; say what Cholesky is for (XFLR5 note is consistent; `cholesky_spd` is documented as a utility that no solver calls.)
- [x] "Eight examples" count in `crates/tpt-aero/README.md`; CLI README omits `rotor`, `trim`, `--json` (Nine examples; CLI README covers `rotor`, `trim` and `--json`.)
- [x] Move `todo.md` and `spec.txt` to `docs/dev/` (superseded: both files stay at the repository root by decision, and `docs/dev/` was removed)

### 4. Features
- [ ] Stall, post-stall and CL_max in the VLM; CD profile drag with Mach and Re *Not in this pass:* the VLM is inviscid and has no stall model yet.
- [ ] Swirl-capable, windmill-capable BEMT; forward-flight BEMT using `SkewedWake` and `DynamicInflow`; differentiable rotor path *Not in this pass:* not started.
- [ ] Atmosphere, wind and gust inputs in `SixDofConfig`; quaternion attitude option **Partly done:** steady wind and body-fixed thrust are in. **Update:** gusts are in (`simulate_partial_with_wind`, `one_minus_cosine_gust`). Still open: atmosphere profiles, and quaternions. Decision: quaternions are deferred; the Euler formulation is already guarded near gimbal lock.
- [x] Prelude gaps (`trim_steady_turn`, `sixdof`, `StateSpace`, modes, gusts, panel and rotor results); rename the `core` alias to `aero_core` (Prelude now exports the turn trim, six-DOF, state-space, mode, gust, panel and rotor result types; `core` renamed `aero_core`.)
- [x] `tpt-aero-panel` `serde` feature; polar and stability/trim strategies in `tpt-aero-verify` **Partly done:** the panel serde feature is in. **Update:** polar, flight-condition and control strategies are in `tpt-aero-verify`, with property tests. **Done:** a trim-request strategy (`trim_request`, `TrimRequestParams`), with property tests that generated requests are accepted and returned trims satisfy level flight.
- [ ] Tests: atmo `tests/`, compressible `from_surfaces`, ground effect with dihedral/fin, NaN and zero inputs, windmill, stalled polar, ISA 84-86 km **Partly done:** NaN and zero-input tests added for the VLM, rotor, core, panel, atmosphere and CLI. Still open: `tests/` for atmo, compressible `from_surfaces`, windmill, stalled polar. **Partly done:** atmo `tests/` directory added, ground effect with dihedral and a fin, NaN and zero inputs. **Update:** the stalled-polar behaviour is covered by `stall_blend_never_raises_lift_magnitude` and the polar tests. Still open: windmill, which waits on the windmill-brake closure.

### 5. Innovation
- [x] `tpt-aero-wasm` crate and a browser demo (polar, ISA, trim) (Builds for wasm32-unknown-unknown; the demo needs `wasm-bindgen` run once, see its README. Not run in a browser here.)
- [x] `tpt-aero-py` (PyO3/maturin) and a pandas recipe on `--json` (Wheel built and exercised in a scratch venv. The pandas recipe reads the CSV sweep; there is no JSON path.)
- [x] `tpt-aero report aircraft.toml` to Markdown/HTML (neutral point, static margin, modes, trim, derivatives, plots) (`tpt-aero report`: trim, static margin, derivatives, modes; no plots yet.)
- [x] AVL `.avl` importer and exporter (decks exist in `docs/validation/avl/`) (Subset: surfaces, sections, controls, reference geometry, vertical fins; YDUPLICATE and transforms are refused.)
- [ ] Optimisation commands on the differentiable VLM: min induced drag twist, CG-range trim, speed/CG/mesh sweeps **Partly done:** one-parameter linear washout minimising induced drag (`optimize_linear_twist`, example `twist_optimization`). **Update:** `trim --speed-sweep`, `wing --mesh-sweep N,N,...`, and a cubic spanwise twist law (`optimize_twist_shape`) are in. Still open: CG-range trim (CG enters only through the static margin, so a CG sweep needs a model change first).
- [x] Mesh-convergence estimate in results, with a suggested `n_span` (`wing --check-mesh`, and `wing_mesh_check`.)
- [ ] Validity-envelope flags (alpha, Mach, Re, beta) on every result **Partly done:** warnings cover aspect ratio, taper, sweep, |α| and Mach, and trim refuses |α| > 20°. **Update:** `validity_flags` checks α, β, Mach and Re against documented limits. Decision: the flags stay a function, not fields on `VlmSolution`, to avoid breaking that struct.
- [ ] Validation: swept, cambered and control-surface cases against AVL; XFLR5; measured data for 6-DOF and compressible **Partly done:** the T-tail and rectangular decks now import and match AVL's lift slope within 5 %. Still open: swept, cambered and control-surface cases, XFLR5, and measured data.

### 6. Adoption
- [ ] Publish `tpt-math-*`, then `tpt-aero-*` (unblocks sibling-repo clone, cargo-generate template, `cargo install`); until then the README quick-start says to clone both repos side by side
- [x] Pin the `tpt-math` ref in `.github/actions/setup-tpt-math`: pinned to `2b86a62` (on `origin/master`), fetched by SHA.
- [ ] Release automation: `cargo-dist` binaries, `[package.metadata.binstall]`, `release-plz` *Not in this pass:* release tooling, which is part of publishing.
- [ ] CI: Windows and macOS, MSRV 1.84, beta; rustdoc `-D warnings`; `mdbook build` and `test`, Pages deploy; build the cargo-generate template; auto-discover examples (`slipstream_tail` is skipped); coverage; criterion tracking; `cargo deny` advisories on a schedule; lychee and typos; Kani; thumbv6m builds without the `|| fallback` mask **Partly done:** matrix job, MSRV raised to 1.85 (clap and proptest need it), beta job, rustdoc with warnings denied, mdBook build and test, and a Kani job were added and are not yet run on a runner. **Update:** the thumbv6m fallback mask is removed; the strict serde build exposed and fixed three no_std regressions (core `ln`, VLM `round`/`sqrt`, stability `hypot`). **Update:** typos, lychee, a scheduled advisory check and a scheduled coverage job are added. Still open: criterion benchmark tracking. Decision: defer it until a baseline store exists, since a benchmark without stored history gives no signal.
- [x] Track `Cargo.lock`; dependabot (`Cargo.lock` no longer ignored; `.github/dependabot.yml` added.)
- [ ] README: badges, runnable `fn main()` quick-start tested as a doctest, sample CLI output and a plot, book link *Not in this pass:* the runnable quick-start as a doctest and badges are not done.
- [ ] `examples/README.md` ladder; CONTRIBUTING, issue templates (including "validation discrepancy"), SECURITY; docs.rs metadata for the umbrella and CLI; keywords (uav, drone, propeller, glider) **Partly done:** CONTRIBUTING, SECURITY, issue and PR templates added. Still open: the examples ladder and docs.rs metadata. **Partly done:** ladder, CONTRIBUTING, SECURITY, issue and PR templates, docs.rs metadata and keywords are done.
- [x] Examples: UAV wing sizing, drone prop selection with motor curve, glider polar, trim and stability report (Four examples: UAV sizing, drone propeller selection, glider polar, trim and stability report.)
- [ ] Templates: `cases/` directory of ready TOML files with expected outputs; `template uav|glider|quadrotor`; cargo-generate variants and prompts; reusable GitHub Action that runs a case file and uploads CSV **Partly done:** `cases/` with expected outputs and `template uav|glider|quadrotor` are done. **Update:** the reusable action is in (`.github/actions/run-case`). Still open: cargo-generate variants.
- [x] Docs: coming from AVL, XFLR5, OpenVSP and JSBSim; FAQ; cookbook; mesh-selection guide; error and units guide (Coming-from page, FAQ, cookbook, mesh-selection and errors-and-units pages in the book.)

### Order
1. Correctness fixes. 2. Publish prep and doc drift. 3. CLI hardening. 4. Physics gaps. 5. Showcase (wasm, Python, cases, AVL importer).
