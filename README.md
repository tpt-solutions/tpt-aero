# tpt-aero

A pure-Rust, AI-native specialized **aerodynamics and flight mechanics** library.
Fast, analytical, semi-empirical solvers — Vortex Lattice Method (VLM), panel
methods, Blade Element Momentum Theory (BEMT) — plus atmospheric modeling and
stability derivative extraction. No C/C++ FFI, no external solver binaries.

## Crate map

| Crate | Purpose |
|---|---|
| [`tpt-aero-core`](crates/tpt-aero-core) | Unit-safe quantities, airfoil polars, lifting-surface geometry |
| [`tpt-aero-vlm`](crates/tpt-aero-vlm) | Horseshoe-vortex VLM with differentiable AIC assembly |
| [`tpt-aero-panel`](crates/tpt-aero-panel) | Source/doublet panel methods, slender body theory, compressibility |
| [`tpt-aero-rotor`](crates/tpt-aero-rotor) | BEMT propellers/rotors, tip losses, dynamic inflow, performance maps |
| [`tpt-aero-atmo`](crates/tpt-aero-atmo) | ISA atmosphere, wind shear, gusts, Dryden/von Kármán turbulence |
| [`tpt-aero-stability`](crates/tpt-aero-stability) | Stability & control derivatives, plant matrices |
| [`tpt-aero-verify`](crates/tpt-aero-verify) | proptest strategies, invariant checks, Kani proof harnesses |
| [`tpt-aero`](crates/tpt-aero) | Feature-gated umbrella re-exporting the whole stack |

## Design principles

- **Zero external solver dependencies** — all aerodynamic influence coefficient
  (AIC) assembly and linear solves are from-scratch pure Rust on the in-house
  [`tpt-math`](https://github.com/tpt-solutions/tpt-math) substrate.
- **Strict licensing** — MIT OR Apache-2.0; no Apache-2.0-only runtime
  dependencies (enforced by `cargo-deny`, see `deny.toml`).
- **Type-state & dimensional safety** — phantom types separate freestream
  conditions from local flow; non-dimensional coefficients from dimensional
  forces.
- **no_std + alloc** at the core; `libm` feature for hosts without hardware
  floats.
- **Native formal verification** — Kani bounded model checking (induction
  factor bounds, AIC non-singularity) and proptest physical invariants
  (symmetry ⇒ zero rolling/yawing moment; freestream reversal reverses lift).

## Status

Phase 0–4 build-out tracked in [todo.md](todo.md) per [spec.txt](spec.txt).
v0.1.0 target: VLM + BEMT + ISA + stability derivatives with full verification.

## License

Dual-licensed under MIT or Apache-2.0 — see [LICENSE-MIT](LICENSE-MIT) and
[LICENSE-APACHE](LICENSE-APACHE). © 2026 TPT Solutions.
