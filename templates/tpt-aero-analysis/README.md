# {{project-name}}

An aerodynamic analysis generated from the `tpt-aero-analysis` template.

```sh
cargo run --release
```

The program sweeps a tapered, cambered wing through angle of attack and prints
the lift coefficient, induced drag and pitching moment. Edit `src/main.rs` to
change the planform, the mesh, the flight condition or the sweep.

For the command line and the full library, see the tpt-aero documentation at
https://github.com/tpt-solutions/tpt-aero.
