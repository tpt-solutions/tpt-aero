# tpt-aero-panel

Two-dimensional panel methods for airfoils and bodies, plus the compressibility
corrections that move low-speed pressures to subsonic Mach numbers.

## What it provides

| Item | Purpose |
|---|---|
| `source::solve_source_panels` | Constant-strength source panels for a non-lifting body. |
| `doublet::solve_doublet_panels` | Hess-Smith style lifting solver with a Kutta condition. Returns `DoubletSolution` with surface pressure coefficients and lift. |
| `slender::SlenderBody` | Slender-body theory for bodies of revolution: Munk lift and moment in incidence, side force and yaw at sideslip. |
| `compress::{prandtl_glauert, karman_tsien, beta}` | Subsonic compressibility corrections. Return `Result`, and reject `M` outside `0 <= M < 1`. |

Geometry comes from `tpt-aero-core`: `airfoil::naca4_closed` produces the
closed-trailing-edge contour the solvers expect.

## Quick start

```rust
use tpt_aero_core::airfoil::naca4_closed;
use tpt_aero_panel::{compress, solve_doublet_panels};

// NACA 0012 with a closed trailing edge.
let foil = naca4_closed("0012", 100).expect("valid NACA code");
let points: Vec<(f64, f64)> = foil.x.iter().zip(&foil.y).map(|(&x, &y)| (x, y)).collect();

// Lifting solve at 4 degrees, 50 m/s, 1 m chord.
let sol = solve_doublet_panels(&points, 4.0_f64.to_radians(), 50.0, 1.0).expect("solvable");
println!("cl = {:.3}", sol.cl);

// Correct a low-speed pressure coefficient to Mach 0.6.
let cp_corrected = compress::prandtl_glauert(-0.5, 0.6).expect("subsonic");
```

## Examples

```text
cargo run -p tpt-aero-panel --example airfoil_pressure
```

`examples/airfoil_pressure.rs` computes the lift of a NACA 0012 at 4 degrees
from the doublet solver, the pressure distribution at zero incidence from the
source solver (peak suction near x/c 0.1), and a Prandtl-Glauert correction to
Mach 0.6. The umbrella crate's `airfoil_cp` example shows the zero-incidence
pressure distribution.

## Feature flags

| Feature | Default | Effect |
|---|---|---|
| `std` | yes | Standard library support in dependencies. |
| `alloc` | no | `no_std` + `alloc` build. |
| `libm` | no | Software math for targets without hardware floats. Implies `alloc`. |

The crate is `#![no_std]` and `#![forbid(unsafe_code)]`. It has no `serde`
feature.

## Conventions

* Contours are **counter-clockwise**. Clockwise input is rejected with
  `PanelError::ClockwiseContour`.
* Use a **closed trailing edge** (`naca4_closed`). Across the gap of an open
  trailing edge, the pressures on the bridging panels are not meaningful.
* Angle of attack is positive nose-up, as in the rest of the workspace (see
  [docs/CONVENTIONS.md](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)).

## Validation

The source solver's NACA 0012 pressure at zero incidence is checked against
published inviscid pressure: the minimum pressure coefficient is within 0.04 of
the reference (-0.43 near x/c 0.1). The lifting solver's lift at 5 degrees is
checked against thick-airfoil theory, within 0.04. At 4 degrees its aft-chord
pressure jump (lower minus upper) matches the flat-plate thin-airfoil value
`4 alpha sqrt((1 - x)/x)` within 15 % for x/c between 0.6 and 0.9, and its
surface pressure is smooth with the trailing-edge Kutta condition satisfied. A
symmetry property test confirms a symmetric section at zero incidence has zero
lift.
Slender-body and compressibility formulas are checked against their closed
forms, and have no external reference data in the repository. See [docs/VALIDATION.md](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md).

## Limits

* Two-dimensional, inviscid. No boundary layer, no separation, no stall.
* The lifting solver is a constant-strength doublet panel method. It is not a
  three-dimensional panel code.
* Compressibility corrections are linear and subsonic only. They do not model
  shock formation; `karman_tsien` extends the usable range somewhat but is still
  not a transonic method.
* Slender-body results are first-order in thickness and incidence.

## Benchmarks

`cargo bench -p tpt-aero-panel` runs a criterion benchmark (`benches/panel.rs`).

## Minimum supported Rust version

1.85, inherited from the workspace.

Part of the [tpt-aero](https://github.com/tpt-solutions/tpt-aero) workspace.
See [conventions](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)
(frames and signs) and [validation](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md)
(what is checked, and what is not). Release notes are in the
[changelog](https://github.com/tpt-solutions/tpt-aero/blob/master/crates/tpt-aero-panel/CHANGELOG.md).
Licensed MIT OR Apache-2.0, (c) 2026 TPT Solutions.
