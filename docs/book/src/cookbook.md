# Cookbook

Each recipe is a command or a few lines of Rust. The case files are in
`cases/` and the examples are in `crates/tpt-aero/examples/`.

## Size a UAV wing for a cruise lift coefficient

```text
tpt-aero template uav -o uav.toml
tpt-aero wing uav.toml --svg polar.svg -o polar.csv
```

or run `cargo run --release -p tpt-aero --example uav_wing_sizing`, which
trades span against induced drag at a fixed lift coefficient.

## Pick a propeller for a drone

`cargo run --release -p tpt-aero --example drone_prop_selection` matches each
candidate's torque to a motor's linear torque–speed curve and ranks the ones
that meet the thrust. To see the whole thrust and power curve for one
propeller, sweep its speed:

```text
tpt-aero rotor quadrotor.toml --rpm-sweep 4000:9000:500
```

## Find the best glide ratio

`cargo run --release -p tpt-aero --example glider_polar` prints the polar with
the best L/D and the minimum-sink angle.

## Trim and check stability

```text
tpt-aero template aircraft -o aircraft.toml
tpt-aero trim aircraft.toml --speed 50
tpt-aero report aircraft.toml --speed 50 -o report.md
```

`report` gives the trim, the static margin, the derivatives and the modes in
Markdown. Set `mass.cg_x` and the measured inertias first.

## Check the mesh

```text
tpt-aero wing uav.toml --check-mesh
```

Prints the lift at 1×, 2×, 4× and 8× the strip count and the first level that
changes it by less than 1 %. See [mesh selection](mesh-selection.md).

## Standard atmosphere at an offset

```text
tpt-aero isa 0 5000 --delta-isa 15          # hot day, SI
tpt-aero isa 0 16000 --imperial             # feet in, °R, psf, slug/ft³, ft/s out
```

## Turbulence and gusts in code

`cargo run --release -p tpt-aero --example gust_response` drives a complete
aircraft with a discrete gust and a Dryden record, and prints the load-factor
statistics.
