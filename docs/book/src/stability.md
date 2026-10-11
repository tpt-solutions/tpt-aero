# Stability, trim and simulation

* **Trim.** `trim_level_flight` for level flight. `trim_steady_turn` for a
  steady turn or steady sideslip, from a bank angle and turn rate. The turn
  solve continues from level flight, so the root it finds stays on the
  physical branch.
* **Linear plants.** `FullLongitudinalPlant` and `LateralPlant` give the state
  matrices, and `modes` gives their eigenmodes.
* **Time domain.** `simulate` integrates a linear `StateSpace` with RK4.
  `with_feedback` closes a state-feedback loop. `with_integrator` adds
  integral action, the autotrim of a pitch hold.
* **Six degrees of freedom.** `sixdof::simulate` integrates the nonlinear
  rigid-body equations (Stevens and Lewis, Euler angles) with VLM aerodynamics
  and a thrust closure. Thrust acts along the airspeed vector, which matches
  the trim solvers. Parasite drag enters as a constant `C_D0`. Further loads
  (fuselage, propeller torque reaction) go through `simulate_with_loads`.
