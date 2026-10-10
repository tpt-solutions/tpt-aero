# tpt-aero-rotor

Blade element momentum theory for propellers and rotors: bounded induction
solve (the bound is formally checked with Kani in `tpt-aero-verify`), Prandtl
and Glauert tip loss, hub loss, polar or linear section models, dynamic
inflow, skewed-wake modulation and advance-ratio performance maps.

No thrust is integrated inside the hub cutout. Efficiency and figure of merit
are reported as zero once thrust reverses. Checked against momentum-theory
limits; not yet against measured rotor data.

Part of the [tpt-aero](https://github.com/tpt-solutions/tpt-aero) workspace.
See [conventions](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)
(frames and signs) and [validation](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md)
(what is checked, and what is not). Licensed MIT OR Apache-2.0, (c) 2026 TPT Solutions.
