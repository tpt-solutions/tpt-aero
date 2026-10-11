# tpt-aero-rotor

Blade Element Momentum Theory (BEMT) for propellers and rotors. The crate
divides a blade into radial stations, balances each annulus's blade-element
forces against the momentum-theory induced velocity, and integrates thrust,
torque and power. It includes the standard tip and hub loss corrections,
dynamic inflow for time-domain work, and advance-ratio performance maps.

## What it provides

| Item | Purpose |
|---|---|
| `RotorGeometry`, `BladeSection` | Radial blade description: `rectangular(...)` builds a constant-chord blade, `with_linear_twist(...)` sets root-to-tip twist. |
| `SectionModel` | Section aerodynamics. `Linear` is the classic lift-slope and quadratic-drag closure with `cl_slope`, `cl_max`, `cd0`, `cd_cl2`. |
| `BemtSolver`, `OperatingCondition` | Bind geometry and closures once, then solve many operating points. `OperatingCondition` holds axial speed, shaft speed, density and collective pitch. |
| `BemtSolution`, `RotorLoads`, `SectionResult` | Thrust, torque, power and figure of merit, plus per-station results and convergence information. |
| `TipLoss`, `HubLoss` | Prandtl and Glauert tip loss; Prandtl hub loss with root cutout fraction. |
| `performance_map` | Thrust, torque and efficiency across an advance-ratio sweep. |
| `DynamicInflow`, `SkewedWake` | First-order inflow dynamics and skewed-wake modulation for rotor-on-airframe work. |
| `Slipstream`, `CylinderWake`, `SlipstreamModel` | Momentum-theory slipstream, a cylinder-wake induced-velocity field, and a choice between them. |

## Quick start

A small UAV propeller: 10 inch diameter, two blades, 20 mm chord, root twist
25 degrees tapering to 5 degrees at the tip.

```rust
use tpt_aero_rotor::{
    BemtSolver, HubLoss, OperatingCondition, RotorGeometry, SectionModel, TipLoss,
};

let geometry = RotorGeometry::rectangular(0.127, 0.0127, 2, 0.020, 30)
    .with_linear_twist(25.0_f64.to_radians(), 5.0_f64.to_radians());

let solver = BemtSolver::new(
    geometry,
    SectionModel::Linear {
        cl_slope: 2.0 * core::f64::consts::PI,
        cl_max: 1.2,
        cd0: 0.012,
        cd_cl2: 0.02,
    },
    TipLoss::Prandtl,
    HubLoss::Prandtl { x_hub: 0.1 },
)
.expect("valid geometry");

// Hover at about 7600 rpm, sea-level density.
let hover = OperatingCondition {
    axial_speed: 0.0,
    omega: 800.0,
    density: 1.225,
    collective: 0.0,
};
let solution = solver.solve(&hover).expect("converges");
println!("converged in {} iterations", solution.iterations);
```

The umbrella crate's `propeller_sweep` example runs the full advance-ratio map,
and `propeller_flight` feeds these loads into a six-degree-of-freedom run.

## Examples

```text
cargo run -p tpt-aero-rotor --example hover_bemt
cargo run -p tpt-aero-rotor --example slipstream_tail
```

* `hover_bemt` solves one hover point of a 10 inch propeller, converts `C_T` to
  newtons, and prints an advance-ratio map. Thrust reverses between `J` of about
  0.31 and 0.46, and efficiency is reported as zero from there on.
* `slipstream_tail` shows how a propeller slipstream raises the dynamic
  pressure, and so the lift, of a tailplane behind it.

## Feature flags

| Feature | Default | Effect |
|---|---|---|
| `std` | yes | Standard library support in dependencies. |
| `alloc` | no | `no_std` + `alloc` build. Confirmed to build for `thumbv6m`. |
| `libm` | no | Software math for targets without hardware floats. Implies `alloc`. |
| `serde` | no | `Serialize` / `Deserialize` on parameter and result types. Implies `alloc`. |

The crate is `#![no_std]` and `#![forbid(unsafe_code)]`.

## Conventions

* `axial_speed` is positive for climb or propeller inflow, zero in hover, and
  small negative for incipient descent.
* Efficiency and figure of merit are reported as zero once thrust reverses.
* No thrust is integrated inside the hub cutout.
* Advance ratio is `J = V / (n D)`, the convention used by `performance_map`.

## Validation

The solver is checked against momentum-theory limits and mass-flow
conservation. The Caradonna-Tung hover case is compared with published
measurements: BEMT over-predicts hover thrust by about 23 %, a known
consequence of the simple section model. Full results are in
[docs/VALIDATION.md](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md).
The induction bound is checked with Kani in `tpt-aero-verify`.

## Limits

* The default section model is linear with a quadratic drag polar. It has no
  stall delay, no Reynolds or Mach dependence and no three-dimensional lift
  correction. Use `SectionModel` with your own polar data for better accuracy.
* Hover thrust is over-predicted with the simple closures (see Validation).
* Slipstream and cylinder-wake models are inviscid and linearised.
* Single isolated rotor. No blade-to-blade or rotor-to-fuselage interference
  beyond what the slipstream models carry.

## Benchmarks

`cargo bench -p tpt-aero-rotor` runs the BEMT criterion benchmark
(`benches/bemt.rs`).

## Minimum supported Rust version

1.85, inherited from the workspace.

Part of the [tpt-aero](https://github.com/tpt-solutions/tpt-aero) workspace.
See [conventions](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)
(frames and signs) and [validation](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md)
(what is checked, and what is not). Release notes are in the
[changelog](https://github.com/tpt-solutions/tpt-aero/blob/master/crates/tpt-aero-rotor/CHANGELOG.md).
Licensed MIT OR Apache-2.0, (c) 2026 TPT Solutions.
