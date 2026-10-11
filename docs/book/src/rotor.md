# Rotors and propellers (BEMT)

```rust
use tpt_aero::rotor::{BemtSolver, HubLoss, OperatingCondition, RotorGeometry, SectionModel, TipLoss};

let geometry = RotorGeometry::rectangular(1.0, 0.1, 2, 0.08, 20).with_linear_twist(0.4, 0.1);
let section = SectionModel::Linear { cl_slope: 6.283, cl_max: 1.4, cd0: 0.01, cd_cl2: 0.01 };
let solver = BemtSolver::new(geometry, section, TipLoss::Prandtl, HubLoss::None)?;
let sol = solver.solve(&OperatingCondition::hover(300.0, 1.225))?;
println!("C_T = {:.5}, figure of merit = {:.3}", sol.loads.ct, sol.loads.figure_of_merit);
```

Use `performance_map` to sweep axial speed and get thrust, power and efficiency
against advance ratio. The propeller example,
`crates/tpt-aero/examples/propeller_flight.rs`, feeds BEMT thrust into the
six-degree-of-freedom simulation.
