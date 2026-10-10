# tpt-aero-verify

Verification support: proptest strategies for wings, angles of attack, speeds
and lattice resolutions, property tests (symmetry gives zero roll and yaw,
freestream reversal reverses lift) and Kani proof harnesses for the BEMT
induction bounds and linear-solve soundness (run with `cargo kani`).

`std` only; not part of the default feature set of the umbrella crate.

Part of the [tpt-aero](https://github.com/tpt-solutions/tpt-aero) workspace.
See [conventions](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)
(frames and signs) and [validation](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md)
(what is checked, and what is not). Licensed MIT OR Apache-2.0, (c) 2026 TPT Solutions.
