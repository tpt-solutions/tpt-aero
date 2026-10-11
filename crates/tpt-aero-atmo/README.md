# tpt-aero-atmo

Atmospheric models for flight mechanics and gust work: the International
Standard Atmosphere, wind-shear profiles, discrete gusts, and continuous
turbulence from Dryden and von Kármán filters. Everything is `no_std` and runs
with a reproducible internal noise source, so simulations repeat exactly.

## What it provides

| Module | Contents |
|---|---|
| `isa` | `at(altitude)` returns an `IsaState` with temperature, pressure, density and speed of sound from 0 to 84 852 m geopotential. `at_isa_offset` gives a non-standard day that stays hydrostatically consistent. Sutherland `dynamic_viscosity` and `kinematic_viscosity`. |
| `shear` | `power_law` and `log_law` wind profiles. `DiscreteGust` with `GustShape` for 1-cos and sharp-edged gusts, with `speed_at` and `gradient_at`. |
| `turbulence` | `TurbulenceModel` (Dryden or von Kármán), `TurbulenceAxis`, `TurbulenceField` with low- and medium-altitude presets, `WhiteNoise` (uniform and Gaussian from a seeded generator), and `simulate_axis`. |

## Quick start

```rust
use tpt_aero_atmo::isa;

// Sea-level standard atmosphere.
let sea_level = isa::at(0.0).expect("within the standard atmosphere range");
assert!((sea_level.density - 1.225).abs() < 1e-3);

// The same call at 5000 m.
let fl164 = isa::at(5000.0).expect("within range");
println!("density at 5 km: {:.3} kg/m^3", fl164.density);
```

The umbrella crate's `gust_response` example combines a 1-cos gust with a
Dryden turbulence field and reports load factors.

## Examples

```text
cargo run -p tpt-aero-atmo --example atmosphere_gusts
```

`examples/atmosphere_gusts.rs` prints the standard atmosphere at four altitudes,
a non-standard day, a power-law wind profile and a 1-cos gust's profile over
its 2H length. Each value can be checked by hand: the gust peaks at H (15 m in
the example) and returns to zero at 2H.

## Feature flags

| Feature | Default | Effect |
|---|---|---|
| `std` | yes | Standard library support in dependencies. |
| `alloc` | no | `no_std` + `alloc` build. |
| `libm` | no | Software math for targets without hardware floats. Implies `alloc`. |
| `serde` | no | `Serialize` / `Deserialize` on state and parameter types. Implies `alloc`. |

The crate is `#![no_std]` and `#![forbid(unsafe_code)]`.

## Conventions

* `isa::at` takes geopotential altitude in metres.
* Gust `speed_at` is the vertical gust speed; positive is an updraft.
* `WhiteNoise` is seeded, so a given seed gives the same sequence on every run.

## Validation

The gust response is checked against the FAR 25 alleviation factor
`Kg = 0.88 mu / (5.3 + mu)` at 12.5 chords, through the umbrella's
`gust_response` example (0.59 against 0.62). The atmosphere and turbulence
models are not compared against external tables in the validation suite. Details are in
[docs/VALIDATION.md](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md).

## Limits

* The ISA covers 0 to 84 852 m geopotential. Requests above that return an
  error rather than an extrapolated value.
* Turbulence filters are linear, shaping-filter models. They are not a fitted
  description of any particular atmosphere, and they carry no intermittency or
  altitude-dependent scale beyond the two presets.
* Shear profiles are single-variable, steady and horizontally uniform.

## Minimum supported Rust version

1.85, inherited from the workspace.

Part of the [tpt-aero](https://github.com/tpt-solutions/tpt-aero) workspace.
See [conventions](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)
(frames and signs) and [validation](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md)
(what is checked, and what is not). Release notes are in the
[changelog](https://github.com/tpt-solutions/tpt-aero/blob/master/crates/tpt-aero-atmo/CHANGELOG.md).
Licensed MIT OR Apache-2.0, (c) 2026 TPT Solutions.
