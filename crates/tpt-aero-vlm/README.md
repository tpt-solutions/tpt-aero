# tpt-aero-vlm

A horseshoe-vortex lattice method (VLM) for wings and whole aircraft. It
builds a lattice of bound and trailing vortices on each lifting surface,
assembles the aerodynamic influence coefficient (AIC) matrix, solves for the
circulation and integrates forces and moments. Everything runs in plain Rust
on a native dense LU solver.

## What it provides

* **Wings and whole aircraft.** `Surface::horizontal` and `Surface::vertical`
  describe a planform with its spanwise and chordwise resolution. Compose
  several with `VlmModel::from_surfaces` (wing, tailplane, fin) and the
  downwash and sidewash between them come from the shared influence matrix.
* **Fast alpha sweeps.** `VlmModel::solve_sweep` factors the AIC once and caches
  the induced-velocity influences, then costs O(N²) per angle. It is about
  20x faster than solving each angle separately for a 17-angle sweep at 500
  panels. `VlmModel::prepare` returns a `PreparedVlm` for the same reuse when
  you need more than one kind of solve.
* **Sideslip and body rates.** `PreparedVlm::solve_state(FlightState)` returns
  wind-axis forces and body moments in aircraft-standard signs, with moments
  about stability axes.
* **Control surfaces.** `Surface::with_control` defines an elevator, aileron
  or rudder as a change of incidence aft of a hinge. `prepare_with_twists`
  solves with a deflection applied.
* **Planform twist and camber.** `WingSection::twist` is applied on top of any
  per-solve twist slice.
* **Options, all off by default.** Cosine spanwise spacing
  (`SpanSpacing::Cosine`), a freestream-aligned wake (`WakeModel::Freestream`),
  ground effect by image horseshoes (`with_ground_height`), and subsonic
  compressibility by Göthert's rule (`from_planform_compressible`).
* **Gradients.** `solve_generic` runs on dual numbers from `tpt-aero-core`.
* **Export.** `export::to_vtk` writes the lattice with circulation and surface
  data as legacy ASCII VTK for ParaView or PyVista.

## Quick start

```rust
use tpt_aero_core::WingPlanform;
use tpt_aero_vlm::{Surface, VlmModel};

// 10 m span, 1 m chord, 24 spanwise by 4 chordwise panels.
let planform = WingPlanform::rectangular(10.0, 1.0);
let model = VlmModel::from_surfaces(&[Surface::horizontal(planform, 24, 4)])
    .expect("valid configuration");

// One solve at 5 degrees, 50 m/s, sea-level density.
let sol = model.solve(5.0_f64.to_radians(), 50.0, 1.225, &[]).expect("solvable");
println!("CL = {:.3}, CDi = {:.4}", sol.cl, sol.cdi);

// A sweep factors the influence matrix once and reuses it for every angle.
let alphas: Vec<f64> = (0..=8).map(|i| (i as f64).to_radians()).collect();
let sweep = model.solve_sweep(&alphas, 50.0, 1.225).expect("solvable");
assert_eq!(sweep.len(), alphas.len());
```

## Examples

```text
cargo run -p tpt-aero-vlm --example wing_polar
```

`examples/wing_polar.rs` runs an alpha sweep of a 10:1 rectangular wing, prints
the polar, and compares the fitted lift slope with the elliptic-wing value
`2 pi AR / (AR + 2)`. At this mesh the VLM comes out about 6 % below that
value, so treat the comparison as a sanity check, not a precision claim.

The umbrella crate has more: `wing_vlm` for the polar and spanwise loading, and
`aircraft_stability` for a whole aircraft.

## Feature flags

| Feature | Default | Effect |
|---|---|---|
| `std` | yes | Standard library support in the core and numeric dependencies. |
| `alloc` | no | `no_std` + `alloc` build. |
| `libm` | no | Software math for targets without hardware floats. Implies `alloc`. |
| `serde` | no | `Serialize` / `Deserialize` on solution and parameter types. Implies `alloc`. |

The crate is `#![no_std]` and `#![forbid(unsafe_code)]`.

## Conventions

* Bound vortices are on the quarter chord; control points are at three
  quarters chord.
* Every trailing leg sheds along **body +x**, regardless of alpha and beta.
  This keeps the influence matrix independent of the flow state, so a sweep can
  factor it once. It is accurate for small angles and degrades at large alpha
  or sideslip (see the `WakeModel::Freestream` option).
* `VlmSolution::cl_roll` and `cn` are body-axis moments. Use
  `VlmSolution::stability_roll_yaw(alpha)` when forming stability derivatives.
* Coefficients are normalised by the first surface's reference area, mean
  aerodynamic chord and span. Override with `VlmModel::with_reference`.

The full set of frame and sign rules is in
[docs/CONVENTIONS.md](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md).

## Validation

Checked against Prandtl lifting-line theory for rectangular and tapered wings
(CL, section loading and Oswald efficiency), against the elliptic-wing limit,
for aileron roll authority and adverse yaw, and for prepared-versus-one-shot
agreement. Details and known weaknesses, including the near-field induced-drag
residual and the zero-lift offset that decays with aspect ratio, are in
[docs/VALIDATION.md](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md).

## Limits

* Potential flow, small angles, attached flow. There is no stall, no
  viscosity and no thickness.
* Induced drag is a near-field residual, not a far-field Trefftz-plane
  integral. Treat `cdi` as indicative and read the validation notes before
  relying on it, especially near zero lift where the Oswald ratio is
  ill-conditioned.
* Compressibility is available only through the Göthert option, for subsonic
  flow.
* No fuselage, nacelle or propeller interference.
* Dense AIC: the cost grows with the square of the panel count. The CLI caps a
  surface at 2 500 panels.

## Minimum supported Rust version

1.85, inherited from the workspace.

Part of the [tpt-aero](https://github.com/tpt-solutions/tpt-aero) workspace.
See [conventions](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)
(frames and signs) and [validation](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md)
(what is checked, and what is not). Release notes are in the
[changelog](https://github.com/tpt-solutions/tpt-aero/blob/master/crates/tpt-aero-vlm/CHANGELOG.md).
Licensed MIT OR Apache-2.0, (c) 2026 TPT Solutions.
