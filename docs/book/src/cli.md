# Command line

The `tpt-aero` binary wraps the library. Each command reads a TOML case file
(see [Case files](cases.md)) and writes text, CSV, JSON or a VTK file.

| Command | What it does |
|---|---|
| `tpt-aero wing <case.toml> [--out file.csv]` | Alpha sweep of a single wing; CSV of CL, CD, CM and dimensional forces. |
| `tpt-aero derivatives <case.toml> [--alpha <deg>] [--json]` | Longitudinal and lateral stability derivatives at one angle of attack. |
| `tpt-aero vtk <case.toml> [--alpha <deg>] [--out file.vtk]` | Lattice and circulation for ParaView. |
| `tpt-aero rotor <rotor.toml> [--json]` | One BEMT operating point: thrust, torque, power, figure of merit. |
| `tpt-aero trim <aircraft.toml> --speed <m/s> [--altitude <m>] [--bank <deg>] [--json]` | Steady level trim, or a coordinated turn when `--bank` is given. |
| `tpt-aero isa <altitude_m>...` | ISA properties at the given altitudes, as CSV. |
| `tpt-aero template wing`, `template rotor`, `template aircraft` | Print a commented starting case file. |

## Worked example: a turn trim

```sh
tpt-aero template aircraft > aircraft.toml
tpt-aero trim aircraft.toml --speed 50 --bank 20
```

The output gives the angle of attack, sideslip, the three control deflections,
the load factor and the turn radius. A coordinated turn has sideslip near zero
and a load factor of `1/cos φ`. The rudder and aileron that hold the turn
appear as `rudder_deg` and `aileron_deg`.

Errors exit with a non-zero status and a message that names the file, the
field or the control at fault.
