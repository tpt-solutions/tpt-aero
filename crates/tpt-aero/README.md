# tpt-aero

Pure-Rust aerodynamics and flight mechanics: vortex-lattice wings and whole
aircraft, panel methods, propeller and rotor BEMT, the standard atmosphere with
gusts and turbulence, stability and control derivatives, trim, flight modes,
time-domain simulation and control design.

This is the umbrella crate. It re-exports the sub-crates behind Cargo features,
provides a `prelude` with the types you use most, and carries the runnable
examples. Depend on it when you want the whole stack; depend on a single
sub-crate when you want less.

## Quick start

A wing with a tailplane and elevator, trimmed for level flight at 50 m/s,
with its static stability read off:

```rust
use tpt_aero::prelude::*;

let aircraft = VlmModel::from_surfaces(&[
    Surface::horizontal(WingPlanform::rectangular(10.0, 1.2), 20, 3),
    Surface::horizontal(WingPlanform::rectangular(4.0, 0.6), 10, 3)
        .at(4.0, 0.0, 0.0)
        .with_control("elevator", 0.6, 0.0, 2.0, ControlSymmetry::Symmetric),
]).expect("valid configuration");

let trim = trim_level_flight(&aircraft, &TrimRequest {
    weight_n: 600.0 * 9.80665,
    v_inf: 50.0,
    density: 1.225,
    elevator: "elevator",
    cd0: 0.03,
}).expect("trims");

let d = longitudinal_fd(&aircraft, trim.alpha_rad, 50.0, 1.225, 1e-4).expect("solves");
println!("alpha {:.1} deg, elevator {:.1} deg, Cm_alpha {:.2}",
    trim.alpha_rad.to_degrees(), trim.elevator_rad.to_degrees(), d.cm_alpha);
```

A single wing needs even less:

```rust
use tpt_aero::{CamberModel, VlmModel, WingPlanform};

let wing = WingPlanform::rectangular(10.0, 1.0);
let model = VlmModel::from_planform(&wing, 24, 4, CamberModel::None).expect("valid");
let sol = model.solve(5.0_f64.to_radians(), 50.0, 1.225, &[]).expect("solvable");
println!("CL = {:.3}", sol.cl);
```

## Features

Every solver sub-crate is on by default. Turn off `default-features` and pick
what you need.

| Feature | Default | Enables |
|---|---|---|
| `core` | yes | `tpt-aero-core`: quantities, coefficients, airfoil polars, geometry. Available as `tpt_aero::aero_core`. |
| `vlm` | yes | `tpt-aero-vlm`: vortex lattice method. Available as `tpt_aero::vlm`. |
| `panel` | yes | `tpt-aero-panel`: source and doublet panels, slender bodies, compressibility. Available as `tpt_aero::panel`. |
| `rotor` | yes | `tpt-aero-rotor`: BEMT, inflow models, performance maps. Available as `tpt_aero::rotor`. |
| `atmo` | yes | `tpt-aero-atmo`: ISA, shear, gusts, turbulence. Available as `tpt_aero::atmo`. |
| `stability` | yes | `tpt-aero-stability`: derivatives, trim, modes, simulation. Available as `tpt_aero::stability`. |
| `verify` | no | `tpt-aero-verify`: proptest strategies and Kani harnesses. Requires `std`. |
| `serde` | no | `Serialize` / `Deserialize` on result and parameter types across the enabled crates. |
| `std` | yes | Standard library support across the enabled crates. |
| `alloc` | no | `no_std` + `alloc` build. |
| `libm` | no | Software math for targets without hardware floats. Implies `alloc`. |

The solver crates are `no_std` + `alloc`. Disable `std` for embedded targets.
The full feature powerset (52 combinations) was confirmed to compile when the
feature wiring was last fixed; re-run the check after changing features.

```toml
# Stability analysis only, for a no_std + alloc target with software math:
tpt-aero = { version = "0.1", default-features = false, features = ["stability", "alloc", "libm"] }
```

## Prelude

`use tpt_aero::prelude::*;` brings in the common types: `VlmModel`, `Surface`,
`ControlSymmetry`, `WingPlanform`, `WingSection`, `AngleOfAttack`,
`DynamicPressure`, `MachNumber`, `ReynoldsNumber`, `naca4`, `naca4_closed`, `isa`,
`DiscreteGust`, `GustShape`, `TrimRequest`, `trim_level_flight` and
`longitudinal_fd`, among others. See the module rustdoc for the full list.

## Examples

Run from the repository root with `cargo run --example <name>`:

| Example | Shows |
|---|---|
| `aircraft_stability` | Neutral point, derivatives, trim and short-period mode for a whole aircraft. |
| `gust_response` | 1-cos gust against the FAR 25 alleviation factor, and Dryden turbulence load factors. |
| `control_design` | A yaw damper that cuts the ringing yaw rate about nine-fold, and a pitch-attitude hold with about a 1 s rise and no overshoot. |
| `wing_vlm` | Wing polar and spanwise loading. |
| `airfoil_cp` | NACA 0012 pressure distribution from the panel method. |
| `propeller_sweep` | BEMT propeller performance map over advance ratio. |
| `propeller_flight` | BEMT propeller thrust and torque in a six-degree-of-freedom run, with the slipstream raising the tail's dynamic pressure. Use `--release`. |
| `stability_derivatives` | Centre-of-gravity sweep and short-period mode. |
| `isa_table` | Standard atmosphere table. |

## Serialisation

With the `serde` feature, results and parameters serialise to JSON or any other
serde format. Use `serde_json` with its `float_roundtrip` feature for bit-exact
round trips. `WingPlanform` is serialise-only: deserialise it through the
validating constructors so invalid geometry is rejected.

## Conventions and validation

Frames, signs and the meaning of a positive control command are fixed in
[docs/CONVENTIONS.md](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md).
Read it before comparing a coefficient with a published value.

What each part has been checked against, and what it has not, is in
[docs/VALIDATION.md](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md).
The checks are against analytic theory (lifting line, thin and thick airfoil,
Lanchester, FAR 25), the AVL 3.52 reference and a small set of published data.
None is a comparison with wind-tunnel or flight-test data for a specific
aircraft.

## Scope

Potential-flow, low-speed, attached-flow aerodynamics. There is no stall, no
viscous drag build-up beyond what you supply (`cd0`), no fuselage or propeller
interference in the lattice, and no compressibility in the VLM beyond the
Göthert option. The aircraft model is a wing plus tail plus fin, and the
dynamics are linear about a trim point unless you use the nonlinear six-DOF
model.

## Minimum supported Rust version

1.85.

Part of the [tpt-aero](https://github.com/tpt-solutions/tpt-aero) workspace.
See [conventions](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)
(frames and signs) and [validation](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md)
(what is checked, and what is not). Release notes are in the
[changelog](https://github.com/tpt-solutions/tpt-aero/blob/master/crates/tpt-aero/CHANGELOG.md).
Licensed MIT OR Apache-2.0, (c) 2026 TPT Solutions.
