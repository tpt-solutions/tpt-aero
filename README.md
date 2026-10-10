# tpt-aero

A pure-Rust **aerodynamics and flight-mechanics** library. From a planform to a
handling-qualities analysis in one workspace: vortex-lattice wing and
whole-aircraft aerodynamics, panel methods, propeller/rotor BEMT, the standard
atmosphere with gusts and turbulence, stability and control derivatives, trim,
flight modes, time-domain simulation and simple control design.

No C/C++ FFI and no external solver binaries: every influence matrix and linear
solve is plain Rust. The core crates are `no_std` + `alloc`.

```rust
use tpt_aero::prelude::*;

// Wing + tailplane with an elevator; trim for level flight at 50 m/s.
let aircraft = VlmModel::from_surfaces(&[
    Surface::horizontal(WingPlanform::rectangular(10.0, 1.2), 20, 3),
    Surface::horizontal(WingPlanform::rectangular(4.0, 0.6), 10, 3)
        .at(4.0, 0.0, 0.0)
        .with_control("elevator", 0.6, 0.0, 2.0, ControlSymmetry::Symmetric),
])?;
let trim = trim_level_flight(&aircraft, &TrimRequest {
    weight_n: 600.0 * 9.80665, v_inf: 50.0, density: 1.225,
    elevator: "elevator", cd0: 0.03,
})?;
let d = longitudinal_fd(&aircraft, trim.alpha_rad, 50.0, 1.225, 1e-4);
println!("alpha {:.1} deg, elevator {:.1} deg, Cm_alpha {:.2}",
    trim.alpha_rad.to_degrees(), trim.elevator_rad.to_degrees(), d.cm_alpha);
```

## What you can do

| Task | Where |
|---|---|
| Wing aerodynamics: CL, induced drag, moments, spanwise loading, alpha sweeps (about 20x faster with the prepared solver) | `VlmModel`, `solve_sweep` |
| Whole aircraft: wing, tailplane, fin with downwash interference; incidence, camber, twist | `Surface`, `VlmModel::from_surfaces` |
| Control surfaces: elevator, aileron, rudder, and their derivatives | `Surface::with_control`, `control_derivatives` |
| Stability derivatives: `CLa, Cma, CLq, Cmq` and `Cyb, Clb, Cnb, Clp, Cnp, Clr, Cnr` | `longitudinal_fd`, `lateral_fd` |
| Level-flight trim (alpha, elevator, thrust required) | `trim_level_flight` |
| Flight modes: short period, phugoid, Dutch roll, roll subsidence, spiral | `FullLongitudinalPlant`, `LateralPlant` |
| Time-domain responses to controls and gusts; closed-loop state feedback | `StateSpace`, `simulate`, `with_feedback` |
| 2-D airfoils: source/vortex panels, NACA 4-digit generator, Prandtl-Glauert / Karman-Tsien | `tpt-aero-panel`, `naca4_closed` |
| Propellers and rotors: BEMT, tip/hub loss, dynamic inflow, performance maps | `tpt-aero-rotor` |
| Atmosphere: ISA, wind shear, 1-cos gusts, Dryden / von Karman turbulence | `tpt-aero-atmo` |
| Gradients through the solver (forward-mode autodiff) | `VlmModel::solve_generic` |
| Visualize in ParaView: lattice with circulation | `export::to_vtk` |

## Try it in two minutes

```text
cargo run --example aircraft_stability    # neutral point, derivatives, trim, modes
cargo run --example gust_response         # 1-cos gust vs FAR 25, Dryden turbulence
cargo run --example control_design        # yaw damper and pitch-attitude hold
cargo run --example wing_vlm              # wing polar and spanwise loading
cargo run --example airfoil_cp            # NACA 0012 pressure distribution
cargo run --example propeller_sweep       # BEMT propeller map
cargo run --example stability_derivatives # CG sweep and short-period mode
cargo run --example isa_table             # standard atmosphere
```

Or without writing Rust, with the `tpt-aero` command (`cargo install --path crates/tpt-aero-cli`):

```text
tpt-aero template wing > wing.toml          # a commented case file
tpt-aero wing wing.toml --out polar.csv     # alpha sweep to CSV
tpt-aero derivatives wing.toml --alpha 4    # stability derivatives
tpt-aero vtk wing.toml --alpha 5 --out lattice.vtk
tpt-aero isa 0 5000 11000
```

## Crates

| Crate | Purpose |
|---|---|
| [`tpt-aero`](crates/tpt-aero) | Feature-gated umbrella, prelude and examples |
| [`tpt-aero-core`](crates/tpt-aero-core) | Unit-safe quantities, airfoil polars and geometry, planforms |
| [`tpt-aero-vlm`](crates/tpt-aero-vlm) | Vortex lattice method, multi-surface aircraft, controls, VTK export |
| [`tpt-aero-panel`](crates/tpt-aero-panel) | 2-D panel methods, slender bodies, compressibility corrections |
| [`tpt-aero-rotor`](crates/tpt-aero-rotor) | BEMT propellers and rotors |
| [`tpt-aero-atmo`](crates/tpt-aero-atmo) | ISA, gusts, turbulence |
| [`tpt-aero-stability`](crates/tpt-aero-stability) | Derivatives, trim, modes, simulation, feedback |
| [`tpt-aero-verify`](crates/tpt-aero-verify) | proptest strategies and Kani proof harnesses |
| [`tpt-aero-cli`](crates/tpt-aero-cli) | The `tpt-aero` command |

## Read this before trusting a number

- **Conventions matter.** Frames, signs and the control-command sense are
  documented in [docs/CONVENTIONS.md](docs/CONVENTIONS.md). A positive aileron
  command rolls *left* and a positive rudder command yaws *right*.
- **What has been checked, and against what** is listed honestly in
  [docs/VALIDATION.md](docs/VALIDATION.md), including known weaknesses. The
  checks are against analytic theory (lifting line, thin/thick airfoil,
  Lanchester, FAR 25), not flight-test or wind-tunnel data.
- **Scope.** Potential-flow, low-speed, attached-flow aerodynamics. No stall,
  no viscous drag build-up (you supply `cd0`), no fuselage or propeller
  interference, no ground effect, no compressibility in the VLM. The
  aircraft is wing + tail + fin; the dynamics are linear about a trim point.

## Design principles

- **Pure Rust, small dependency surface.** The solver crates depend only on the
  in-house [`tpt-math`](https://github.com/tpt-solutions/tpt-math) substrate;
  the CLI adds `serde` and `toml`.
- **Strict licensing.** MIT OR Apache-2.0; no Apache-2.0-only dependencies
  (enforced by `cargo-deny`, see `deny.toml`).
- **Safe by construction.** `unsafe` is forbidden. Library entry points return
  `Result` for bad input rather than panicking, with documented exceptions.
- **`no_std` + `alloc`** for the solver crates; the `libm` feature covers
  hosts without hardware floats.
- **Verification.** Proptest invariants (symmetry gives zero roll and yaw,
  freestream reversal reverses lift) and Kani harnesses for the BEMT induction
  bounds and linear-solve soundness.

## Status

Pre-release (v0.1.0 pending): see [CHANGELOG.md](CHANGELOG.md) for what changed
and [todo.md](todo.md) for the open items. MSRV 1.84.

## Contributing

`cargo fmt --all`, `cargo clippy --workspace --all-targets -- -D warnings` and
`cargo test --workspace` must pass; CI also runs `cargo-deny`, a `no_std`
build, and every example plus a CLI smoke test. The workspace path-depends on
a sibling checkout of `tpt-math` (CI clones it automatically).

## License

Dual-licensed under MIT or Apache-2.0: see [LICENSE-MIT](LICENSE-MIT) and
[LICENSE-APACHE](LICENSE-APACHE). (c) 2026 TPT Solutions.
