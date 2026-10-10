# tpt-aero

Pure-Rust aerodynamics and flight mechanics: vortex-lattice wings and whole
aircraft, panel methods, propeller/rotor BEMT, the standard atmosphere with
gusts and turbulence, stability derivatives, trim, flight modes, time-domain
simulation and control design. This is the umbrella crate: it re-exports the
sub-crates behind features and provides `tpt_aero::prelude`.

```rust
use tpt_aero::prelude::*;

let wing = WingPlanform::rectangular(10.0, 1.0);
let model = VlmModel::from_surfaces(&[Surface::horizontal(wing, 24, 4)]).unwrap();
let sol = model.solve(5.0_f64.to_radians(), 50.0, 1.225, &[]).unwrap();
println!("CL = {:.3}", sol.cl);
```

Features: `core`, `vlm`, `panel`, `rotor`, `atmo`, `stability` (all default),
`verify` (opt-in), `std`, `alloc`, `libm`. Eight runnable examples live in
`examples/` (`cargo run --example aircraft_stability`).

Part of the [tpt-aero](https://github.com/tpt-solutions/tpt-aero) workspace.
See [conventions](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)
(frames and signs) and [validation](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md)
(what is checked, and what is not). Licensed MIT OR Apache-2.0, (c) 2026 TPT Solutions.
