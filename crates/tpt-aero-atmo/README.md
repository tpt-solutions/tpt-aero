# tpt-aero-atmo

International Standard Atmosphere (0-86 km), wind shear profiles (power and
log law), discrete gusts (1-cos and sharp-edged) and continuous turbulence
(Dryden and von Karman filters with a reproducible `no_std` noise source).

```rust
let sea_level = tpt_aero_atmo::isa::at(0.0).unwrap();
assert!((sea_level.density - 1.225).abs() < 1e-3);
```

`no_std` + `alloc`.

Part of the [tpt-aero](https://github.com/tpt-solutions/tpt-aero) workspace.
See [conventions](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)
(frames and signs) and [validation](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md)
(what is checked, and what is not). Licensed MIT OR Apache-2.0, (c) 2026 TPT Solutions.
