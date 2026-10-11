# tpt-aero-core

The foundation crate of [tpt-aero](https://github.com/tpt-solutions/tpt-aero).
It defines the types every solver exchanges: unit-safe quantities, force and
moment coefficients, lifting-surface geometry, airfoil section data and a
native dense linear solver. It does not solve any flow itself, so it is cheap
to depend on from embedded or `no_std` code.

## What it provides

| Module | Contents |
|---|---|
| `quantities` | Newtypes for `MachNumber`, `ReynoldsNumber`, `AngleOfAttack`, `AngleOfSideslip`, `DynamicPressure`, `Freestream`, `FreestreamVelocity`, `LocalVelocity`. Constructors such as `MachNumber::try_new` return `Option` for out-of-range input. |
| `coefficients` | `LiftCoefficient`, `DragCoefficient`, `SideForceCoefficient`, pitching, rolling and yawing moment coefficients, plus `ForceVector` and `MomentVector`. |
| `geometry` | `Vec3`, `WingPlanform` (built with `rectangular`, `tapered` or `from_sections` for arbitrary section lists), `WingSection`, and the lifting-line and surface panel types used by the VLM. |
| `polar` | `AirfoilPolar`, `PolarSet` and `StallModel` for tabulated section data. |
| `airfoil` | NACA 4-digit generators: `naca4` (open trailing edge) and `naca4_closed` (sharp trailing edge, for panel methods). |
| `linalg` | `LuFactors`, `lu_solve`, `lu_solve_generic`, `cholesky_spd`. Partial-pivot LU written from scratch; no external solver or C/C++ FFI. `cholesky_spd` is for symmetric positive-definite systems (normal equations, mass matrices); no solver in this release calls it, so it is a utility kept for callers. |
| `viscous` | Flat-plate skin friction (laminar and Schlichting mixed), Hoerner form factor and a `C_D0` component build-up. |
| `num` | The `Real` trait, so that solvers run unchanged over `f64` or forward-mode dual numbers for gradients. |

## Quick start

```rust
use tpt_aero_core::{airfoil::naca4_closed, MachNumber, WingPlanform};

// A 10 m span, 1 m chord rectangular wing.
let wing = WingPlanform::rectangular(10.0, 1.0);
println!("area {:.1} m^2, mean aerodynamic chord {:.2} m", wing.area(), wing.mac());

// Closed-trailing-edge NACA 2412, the input the panel crate expects.
let foil = naca4_closed("2412", 80).expect("valid NACA 4-digit code");

// Mach numbers are checked on construction.
assert!(MachNumber::try_new(0.45).is_some());
assert!(MachNumber::try_new(-0.1).is_none());
```

## Examples

```text
cargo run -p tpt-aero-core --example quickstart
```

`examples/quickstart.rs` tours the crate: a wing's geometry, an airfoil, a
checked Mach number, a skin-friction estimate, a zero-lift drag build-up, and
a small linear solve whose answer you can check by hand.

## Feature flags

| Feature | Default | Effect |
|---|---|---|
| `std` | yes | Standard library support in the numeric dependencies. |
| `alloc` | no | `no_std` + `alloc` build (the geometry and polar types use `Vec`). |
| `libm` | no | Pure-software math for targets without hardware floats. Implies `alloc`. |
| `serde` | no | `Serialize` / `Deserialize` on result and parameter types. Implies `alloc`. |

The crate is `#![no_std]` and `#![forbid(unsafe_code)]`. Constructors that can
fail return `Result` or `Option`; a small number of `new` helpers document
their panic conditions and are intended for inputs known to be valid.

## Conventions

All quantities are SI (m, kg, s, rad, N, Pa) unless a name says otherwise.
Coefficients are normalised by the reference area, chord and span of the
surface they describe. Frames, signs and the control-command sense are defined
once in [docs/CONVENTIONS.md](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md).
Read it before comparing a coefficient with a published value.

## Validation

The types have unit tests in each module. The aerodynamic checks that depend on
these types (lifting-line agreement, panel-method pressure distributions, the
viscous build-up against its own formulas) live in the solver crates. What is
checked against what, and what is not, is listed in
[docs/VALIDATION.md](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md).

## Limits

* `WingPlanform` describes straight-line segments between sections. It is not a
  CAD kernel and carries no airfoil thickness distribution beyond what
  `WingSection` stores.
* The linear solver is dense. It is appropriate for the few-thousand-unknown
  systems the VLM and panel methods produce, not for large sparse problems.
* The viscous helpers are first-order component estimates for early design,
  not a substitute for a boundary-layer solution.

## Minimum supported Rust version

1.85, inherited from the workspace.

Part of the [tpt-aero](https://github.com/tpt-solutions/tpt-aero) workspace.
See [conventions](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)
(frames and signs) and [validation](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md)
(what is checked, and what is not). Release notes are in the
[changelog](https://github.com/tpt-solutions/tpt-aero/blob/master/crates/tpt-aero-core/CHANGELOG.md).
Licensed MIT OR Apache-2.0, (c) 2026 TPT Solutions.
