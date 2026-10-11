# tpt-aero-verify

Verification support for the tpt-aero solvers. It holds three kinds of check:

* **proptest strategies** that generate valid wings, angles of attack, speeds
  and lattice resolutions for property-based tests;
* **physical-invariant property tests**: symmetry gives zero roll and yaw,
  freestream reversal reverses lift, and dynamic-pressure scaling holds;
* **Kani proof harnesses** (compiled only under `cfg(kani)`) for the BEMT
  induction bounds and the soundness of the linear solve.

This crate is for the test suite and for verifiers. Application code should
not depend on it.

## What it provides

| Item | Purpose |
|---|---|
| `symmetric_wing()` | Strategy producing `SymmetricWingParams`: span, chord, sweep and taper within valid ranges. |
| `angle_of_attack()` | Strategy over angles of attack. |
| `freestream_speed()` | Strategy over flight speeds. |
| `lattice_resolution()` | Strategy over spanwise and chordwise panel counts. |
| `tests/properties.rs` | The property tests above, run by `cargo test -p tpt-aero-verify`. |
| `kani_proofs` | Kani harnesses (`cfg(kani)`) for BEMT induction bounds and LU soundness. |

## Usage

Run the property tests:

```text
cargo test -p tpt-aero-verify
```

Run the Kani proofs (requires the [Kani verifier](https://model-checking.github.io/kani/)):

```text
cargo kani -p tpt-aero-verify
```

Kani harnesses are compiled only when `cfg(kani)` is set, so ordinary builds
and tests do not need the verifier installed.

## Examples

```text
cargo run -p tpt-aero-verify --example property_sweep
```

`examples/property_sweep.rs` shows the strategies used outside a test suite. It
draws random symmetric wings and angles of attack, solves each one, and checks
that the roll and yaw moments vanish in symmetric flight. Copy its pattern to
write your own property checks against the solvers.

## Feature flags

| Feature | Default | Effect |
|---|---|---|
| `std` | yes | Required: `proptest` needs the standard library. |

This crate is `std`-only. It is not part of the default feature set of the
umbrella crate; enable it there with `tpt-aero/verify`.

## Validation scope

This crate is itself the validation layer for the other crates. The checks it
makes are invariants and formal bounds, not comparisons with flight or
wind-tunnel data. What each solver is checked against is in
[docs/VALIDATION.md](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md).

## Minimum supported Rust version

1.85, inherited from the workspace. The Kani harnesses need a recent Kani
release.

Part of the [tpt-aero](https://github.com/tpt-solutions/tpt-aero) workspace.
See [conventions](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)
(frames and signs) and [validation](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md)
(what is checked, and what is not). Release notes are in the
[changelog](https://github.com/tpt-solutions/tpt-aero/blob/master/crates/tpt-aero-verify/CHANGELOG.md).
Licensed MIT OR Apache-2.0, (c) 2026 TPT Solutions.
