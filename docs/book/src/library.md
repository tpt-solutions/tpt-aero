# Library tour

The umbrella crate `tpt-aero` re-exports the solver crates. Import the common
items from the prelude:

```rust
use tpt_aero::prelude::*;
```

| Crate | Purpose |
|---|---|
| `tpt-aero-core` | Unit-safe quantities, polars, planform geometry, native linear solves, viscous drag build-up. |
| `tpt-aero-vlm` | Vortex lattice method for wings and whole aircraft. |
| `tpt-aero-panel` | Source and doublet panel methods, slender bodies, compressibility corrections. |
| `tpt-aero-rotor` | BEMT, loss models, dynamic inflow, performance maps. |
| `tpt-aero-atmo` | ISA, wind shear, gusts, turbulence. |
| `tpt-aero-stability` | Derivatives, trim, turns, modes, simulation, control design, six-DOF. |
| `tpt-aero-verify` | Property tests and Kani harnesses. |
| `tpt-aero-cli` | The `tpt-aero` binary. |

Every solver is generic over the `Real` trait in `tpt-aero-core`, so the same
code runs on `f64` and on forward-mode dual numbers.

See the next chapters for the main entry points in each area.
