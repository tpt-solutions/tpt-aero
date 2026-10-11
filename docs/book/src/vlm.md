# Wings and aircraft (VLM)

Build a model from a planform, then solve it:

```rust
use tpt_aero::prelude::*;

let wing = WingPlanform::rectangular(10.0, 1.0);
let model = VlmModel::from_planform(&wing, 24, 4, CamberModel::None)?;
let sol = model.solve(5.0_f64.to_radians(), 50.0, 1.225, &[])?;
println!("CL = {:.3}", sol.cl);
```

## Options

Each option is opt-in and leaves the default behaviour unchanged.

* **Spanwise spacing.** `VlmModel::from_planform_spaced` with
  `SpanSpacing::Cosine` clusters strips toward the tips.
* **Trailing wake.** `with_wake(WakeModel::Freestream)` aligns the trailing
  legs with the flow. The default body-axis wake keeps the solve linear in
  speed and enables the cached `prepare()` path.
* **Ground effect.** `with_ground_height(h)` adds image vortices below a plane
  `h` metres under the model.
* **Compressibility.** `VlmModel::from_planform_compressible(…, mach)` applies
  Göthert's rule: the planform is stretched by `β = √(1 − M²)` and force
  coefficients are referenced accordingly. This reproduces the Helmbold–
  Diederich lift slope for a straight wing.

## Derivatives and trim

`control_derivatives` gives the force and moment derivatives of a named
control. `lateral_fd` and `longitudinal_fd` give the stability derivatives,
and `trim_level_flight` finds the level-flight angle of attack and elevator.
