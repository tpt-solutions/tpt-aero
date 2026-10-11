# Ready case files

Each case is a runnable TOML file with its expected output next to it. The
test `crates/tpt-aero-cli/tests/cases.rs` runs every case and compares the
numbers with the expected file to a relative tolerance of 1e-6, so a change
to the solver that moves a published number fails the test.

| Case | Command | Expected output |
|---|---|---|
| `uav.toml` | `tpt-aero wing uav.toml` | `uav.expected.csv` |
| `glider.toml` | `tpt-aero wing glider.toml` | `glider.expected.csv` |
| `quadrotor.toml` | `tpt-aero rotor quadrotor.toml --json` | `quadrotor.expected.json` |
| `aircraft.toml` | `tpt-aero trim aircraft.toml --speed 50` | `aircraft.trim-50ms.expected.txt` |

To refresh an expected file after an intended change, rerun the command with
`-o` and review the diff before committing it.
