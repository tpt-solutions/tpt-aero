# tpt-aero-panel

Two-dimensional panel methods for airfoils and bodies: constant-strength
source panels (non-lifting), a Hess-Smith lifting solver with a Kutta condition
(`solve_doublet_panels`), slender-body theory for bodies of revolution, and
the Prandtl-Glauert and Karman-Tsien compressibility corrections (which return
`Result` outside `0 <= M < 1`).

Contours must be counter-clockwise (clockwise input is rejected). Use a
closed trailing edge (`tpt_aero_core::naca4_closed`): across the gap of an
open trailing edge the pressures on the bridging panels are not meaningful.

Part of the [tpt-aero](https://github.com/tpt-solutions/tpt-aero) workspace.
See [conventions](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)
(frames and signs) and [validation](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md)
(what is checked, and what is not). Licensed MIT OR Apache-2.0, (c) 2026 TPT Solutions.
