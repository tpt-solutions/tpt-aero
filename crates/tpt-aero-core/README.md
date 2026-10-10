# tpt-aero-core

Foundations shared by every other crate: unit-safe quantities (`MachNumber`,
`ReynoldsNumber`, `AngleOfAttack`, dynamic pressure), force and moment
coefficients, 3-D vectors generic over `Real` (so solvers can run on dual
numbers for gradients), lifting-surface geometry (`WingPlanform`,
`WingSection`), airfoil polars and the NACA 4-digit generators (`naca4`,
`naca4_closed`).

`no_std` + `alloc`. Constructors that can fail return `Result`; `new`
variants with a documented panic are provided where the input is known good.

Part of the [tpt-aero](https://github.com/tpt-solutions/tpt-aero) workspace.
See [conventions](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)
(frames and signs) and [validation](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md)
(what is checked, and what is not). Licensed MIT OR Apache-2.0, (c) 2026 TPT Solutions.
