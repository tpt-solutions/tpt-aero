# tpt-aero-stability

Everything between a VLM model and handling qualities.

* **Derivatives**: `longitudinal_fd` (`CLa, Cma, CLq, Cmq`, flap effectiveness),
  `lateral_fd` (`Cyb, Clb, Cnb, Clp, Cnp, Clr, Cnr, ...`),
  `control_derivatives` for any named control.
* **Trim**: `trim_level_flight` (angle of attack and elevator, thrust required).
* **Modes**: `FullLongitudinalPlant` (short period and phugoid) and
  `LateralPlant` (Dutch roll, roll subsidence, spiral), with a dependency-free
  4x4 eigenvalue routine.
* **Simulation**: RK4 `StateSpace` with control, gust and sidewash inputs, and
  `with_feedback` for closed-loop design.

Linear, small-perturbation, level-flight models. See the examples in the
umbrella crate: `aircraft_stability`, `gust_response`, `control_design`.

Part of the [tpt-aero](https://github.com/tpt-solutions/tpt-aero) workspace.
See [conventions](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)
(frames and signs) and [validation](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md)
(what is checked, and what is not). Licensed MIT OR Apache-2.0, (c) 2026 TPT Solutions.
