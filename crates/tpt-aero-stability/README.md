# tpt-aero-stability

Flight mechanics on top of a VLM model: stability and control derivatives,
level-flight and coordinated-turn trim, the linearised longitudinal and lateral
plants with their modes, time-domain simulation, closed-loop state feedback,
and a nonlinear six-degree-of-freedom integrator for cases the linear models
do not cover.

## What it provides

| Module | Contents |
|---|---|
| `longitudinal` | `longitudinal_fd` returns `LongitudinalDerivatives`: `C_Lα`, `C_mα`, `C_Lq`, `C_mq`, `C_Lδ`, `C_mδ` by finite difference. `cl_alpha_dual` gives the lift slope from dual numbers. |
| `lateral` | `lateral_fd` returns `LateralDerivatives`: `C_Yβ`, `C_lβ`, `C_nβ`, the roll, yaw and sideforce rate derivatives. `LateralPlant` with `from_derivatives`, `eigenvalues`, `modes` (Dutch roll, roll subsidence, spiral) and `state_space`. |
| `controls` | `control_derivatives` for any named control: lift, drag, side force and the three moments per radian of command. |
| `trim` | `trim_level_flight` with `TrimRequest`: Newton solve for angle of attack and elevator, reporting the thrust required. |
| `turn` | `trim_steady_turn` and `coordinated_turn_rate`: steady turns and steady sideslip by continuation from level flight. |
| `plant` | `LongitudinalPlant` and `PlantParameters`: the linearised longitudinal `A` and `B` matrices and the short-period characteristics. |
| `phugoid` | `FullLongitudinalPlant`: four-state longitudinal model with speed and drag couplings, giving short-period and phugoid modes. |
| `modes` | `eigenvalues4`: a dependency-free 4x4 eigenvalue routine. |
| `sim` | `StateSpace`, RK4 `simulate`, `with_feedback` for state-feedback design, `eigenvalues`, and `with_integrator` for integral action (autotrim). |
| `sixdof` | Nonlinear six-degree-of-freedom simulation in aircraft axes with Euler angles, VLM aerodynamics, thrust along the airspeed vector, and parasite drag. `simulate_with_loads` adds external force and moment at the centre of gravity. |

## Quick start

Longitudinal derivatives of a wing. The moment reference is set ahead of the
neutral point so that the static stability sign is meaningful.

```rust
use tpt_aero_core::{Vec3, WingPlanform};
use tpt_aero_stability::longitudinal::longitudinal_fd;
use tpt_aero_vlm::{CamberModel, ReferenceGeometry, VlmModel};

let wing = WingPlanform::rectangular(10.0, 1.0);
let model = VlmModel::from_planform(&wing, 16, 4, CamberModel::None)
    .expect("valid")
    .with_reference(ReferenceGeometry {
        area: wing.area(),
        chord: wing.mac(),
        span: wing.span(),
        moment_center: Vec3::new(0.0, 0.0, 0.0),
    });

let d = longitudinal_fd(&model, 2.0_f64.to_radians(), 50.0, 1.225, 1e-4).expect("solves");
println!("C_Lα = {:.2} /rad, C_mα = {:.3} /rad", d.cl_alpha, d.cm_alpha);
```

Trim a whole aircraft (wing, tailplane with an elevator, fin) in level flight
and read off its static stability:

```rust
use tpt_aero_core::WingPlanform;
use tpt_aero_stability::{longitudinal::longitudinal_fd, trim::{trim_level_flight, TrimRequest}};
use tpt_aero_vlm::{ControlSymmetry, Surface, VlmModel};

let aircraft = VlmModel::from_surfaces(&[
    Surface::horizontal(WingPlanform::rectangular(10.0, 1.2), 20, 3),
    Surface::horizontal(WingPlanform::rectangular(4.0, 0.6), 10, 3)
        .at(4.0, 0.0, 0.0)
        .with_control("elevator", 0.6, 0.0, 2.0, ControlSymmetry::Symmetric),
])
.expect("valid configuration");

let trim = trim_level_flight(&aircraft, &TrimRequest {
    weight_n: 600.0 * 9.80665,
    v_inf: 50.0,
    density: 1.225,
    elevator: "elevator",
    cd0: 0.03,
})
.expect("trims");
println!("alpha {:.1} deg, elevator {:.1} deg",
    trim.alpha_rad.to_degrees(), trim.elevator_rad.to_degrees());

let d = longitudinal_fd(&aircraft, trim.alpha_rad, 50.0, 1.225, 1e-4).expect("solves");
assert!(d.cm_alpha < 0.0, "statically stable");
```

The `aircraft_stability`, `control_design` and `gust_response` examples in the
umbrella crate show the full workflow, including modes and time responses.

## Examples

```text
cargo run -p tpt-aero-stability --example short_period
```

`examples/short_period.rs` takes a wing from the VLM to its longitudinal
derivatives, builds the linearised plant for a 500 kg aircraft, and reports the
short-period natural frequency and damping. It is the smallest complete path
from geometry to handling qualities.

For trim, lateral modes, time-domain responses and control design, see the
umbrella crate's `aircraft_stability`, `control_design` and `gust_response`
examples.

## Feature flags

| Feature | Default | Effect |
|---|---|---|
| `std` | yes | Standard library support, forwarded to the VLM and numeric crates. |
| `alloc` | no | `no_std` + `alloc` build. |
| `libm` | no | Software math for targets without hardware floats. Implies `alloc`. |
| `serde` | no | `Serialize` / `Deserialize` on derivative sets, trim results and related types. Implies `alloc`. |

The crate is `#![no_std]` and `#![forbid(unsafe_code)]`. The `sixdof` and
`turn` modules use the same feature set.

## Conventions

* Models are small-perturbation, **level flight** (`gamma0 = 0`), in stability
  axes. The longitudinal state is `[u, alpha, q, theta]` with inputs
  `[elevator, vertical gust]`. The lateral state is `[beta, p, r, phi]` with
  inputs `[aileron, rudder, sidewash]`.
* Rate derivatives are per non-dimensional rate: `p̂ = p b / 2V`,
  `q̂ = q c / 2V`, `r̂ = r b / 2V`.
* Control derivatives are per radian of command. A positive aileron command
  rolls left and a positive rudder command yaws right. See
  [docs/CONVENTIONS.md](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md).
* Thrust is held constant in the phugoid model.

## Validation

Derivatives are checked against the AVL 3.52 reference run (decks and script in
`docs/validation/avl/`) for `C_Lα`, `C_Dᵢ`, `C_Yβ`, `C_lβ`, `C_nβ`, `C_lp`,
`C_lr` and `C_nr`. The phugoid is checked against Lanchester's approximations,
the gust response against FAR 25, and the simulation against the eigenvalues it
should reproduce. The Navion wing-tail lift slope and pitch stiffness are
compared with NASA TN D-6643 flight data, with the differences recorded in
[docs/VALIDATION.md](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md).

## Limits

* Linear, small-perturbation, level-flight models. The linear plants are about a
  single trim point and are not valid at large angles or in a turn (use `turn`
  or `sixdof` for those).
* The aerodynamics come from the VLM: no stall, no viscous drag beyond what you
  supply as `cd0`, no fuselage or propeller interference. The `sixdof` model
  carries a parasite drag term; any other external load must be passed through
  `simulate_with_loads`.
* The phugoid model holds thrust constant.
* No actuator dynamics or sensor models in the simulation.

## Minimum supported Rust version

1.85, inherited from the workspace.

Part of the [tpt-aero](https://github.com/tpt-solutions/tpt-aero) workspace.
See [conventions](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)
(frames and signs) and [validation](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md)
(what is checked, and what is not). Release notes are in the
[changelog](https://github.com/tpt-solutions/tpt-aero/blob/master/crates/tpt-aero-stability/CHANGELOG.md).
Licensed MIT OR Apache-2.0, (c) 2026 TPT Solutions.
