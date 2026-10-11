# Examples

Run any example with `cargo run --release -p tpt-aero --example <name>`. They
are listed from the simplest to the most complete. Each one prints its
assumptions at the top.

## 1. First steps

| Example | What it shows |
|---|---|
| `isa_table` | The standard atmosphere at a list of altitudes. |
| `airfoil_cp` | Surface pressure on a NACA section from the panel method. |
| `wing_vlm` | A rectangular wing swept through angle of attack. |

## 2. Wings, sizing and polars

| Example | What it shows |
|---|---|
| `uav_wing_sizing` | Wing area for a cruise lift coefficient, traded against span. |
| `glider_polar` | A glider polar with the best L/D and the minimum-sink angle. |

## 3. Propellers and rotors

| Example | What it shows |
|---|---|
| `propeller_sweep` | Thrust, torque and power over the advance ratio. |
| `drone_prop_selection` | Propeller choice against a motor torque–speed curve. |
| `propeller_flight` | A propeller on an aircraft: thrust, torque reaction and slipstream at trim. |

## 4. Stability, trim and control

| Example | What it shows |
|---|---|
| `stability_derivatives` | The longitudinal and lateral derivatives of a wing. |
| `aircraft_stability` | Derivatives and modes for a wing and tail. |
| `trim_stability_report` | Trim, static margin, short-period and Dutch-roll modes, and the derivatives. |
| `control_design` | Closed-loop state feedback: a yaw damper and a pitch-attitude hold. |
| `gust_response` | A discrete gust and a Dryden turbulence record on a complete aircraft. |

The rotor crate has one more: `cargo run -p tpt-aero-rotor --example slipstream_tail`.

## Adding an example

Put a `.rs` file in this directory. Cargo finds it without a `Cargo.toml`
entry, and CI runs every example in both crates, so a broken one fails the
build.
