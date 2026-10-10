# tpt-aero-cli

The `tpt-aero` command: aerodynamics from a TOML case file, no Rust required.

```text
tpt-aero template wing > wing.toml          # a commented starting point
tpt-aero wing wing.toml --out polar.csv     # alpha sweep to CSV
tpt-aero derivatives wing.toml --alpha 4    # stability derivatives
tpt-aero vtk wing.toml --alpha 5 --out lattice.vtk   # open in ParaView
tpt-aero isa 0 5000 11000                   # standard atmosphere
```

`wing` columns: `alpha_deg,cl,cdi,cm,cl_roll,cn,lift_n` (the last roll/yaw
columns are body-axis; `derivatives` uses stability axes). Unknown or
misspelled case-file keys are rejected.

Part of the [tpt-aero](https://github.com/tpt-solutions/tpt-aero) workspace.
See [conventions](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)
(frames and signs) and [validation](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md)
(what is checked, and what is not). Licensed MIT OR Apache-2.0, (c) 2026 TPT Solutions.
