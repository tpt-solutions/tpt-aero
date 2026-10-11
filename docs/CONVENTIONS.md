# Conventions

Mistakes in aerodynamics code are usually sign and frame mistakes. This page
is the single source of truth for tpt-aero. All quantities are SI (m, kg, s,
rad, N, Pa) unless a name says otherwise.

## Body frame

The model frame is **x aft, y out the right wing, z up**: the nose points to
negative x. (This is a right-handed frame, but not the aircraft-standard
x-forward, z-down one.)

* Leading edges sit at smaller `x` than trailing edges.
* A horizontal surface spans `y`; a fin is a planform rolled 90 degrees about
  `x`, so its span runs along `z`.
* Chords always run along `x`, for every surface.

## Flow angles

* **Angle of attack `alpha`**: positive nose-up. The relative wind in the body
  frame is `V (cos(alpha) cos(beta), -sin(beta), sin(alpha) cos(beta))`: it
  comes from below the nose and moves up through the body.
* **Sideslip `beta`**: positive when the wind comes from the right of the
  nose (the air moves toward `-y`).
* **Rates `p, q, r`**: aircraft-standard: `p` right wing down, `q` nose up,
  `r` nose right, in rad/s, about the stability axes (below).

## Wake model

Every horseshoe vortex sheds its trailing legs **along body +x**, regardless
of `alpha` and `beta`. This keeps the influence matrix independent of the flow
state (so a sweep factors it once) and is accurate for small angles; it
degrades at large `alpha` or sideslip.

## Stability axes and moment signs

Stability derivatives are defined about the **stability axes**: body axes
rotated by `alpha` about the pitch axis, so the roll axis lies along the
relative wind projected on the plane of symmetry. Using body axes instead
mixes roll and yaw by `alpha` and flips the sign of induced (adverse) yaw.

`FlightState` / `StateLoads` and `VlmSolution::stability_roll_yaw` use
aircraft-standard signs:

| Quantity | Positive means |
|---|---|
| lift `C_L` | normal to the relative wind, up |
| drag `C_D` | along the relative wind, aft |
| side force `C_Y` | to the right |
| pitching moment `C_m` | nose up |
| rolling moment `C_l` | right wing down |
| yawing moment `C_n` | nose right |

The raw `VlmSolution::cl_roll` and `cn` fields are **body-axis, model-frame**
moments about `x` (aft) and `z` (up) and have opposite signs: use
`stability_roll_yaw(alpha)` when forming derivatives.

## Reference geometry

Coefficients are normalized by the **first surface's** area, mean aerodynamic
chord and span. The moment reference point defaults to the quarter chord of
the mean chord of that surface (shifted by its `offset`) and is the centre of
gravity for trim and stability. Override it with `VlmModel::with_reference`.

## Twist and incidence

`WingSection::twist` is a rotation about the span axis, **positive nose-up**,
applied by every solver on top of any per-solve `twists` slice. A tail with
negative incidence carries download.

## Controls

A control surface is a change of local incidence over the panels aft of its
hinge (a plain-flap approximation). A **positive command is trailing edge
down on the +y half** of a horizontal surface:

| Control | Positive command | Aircraft-standard positive |
|---|---|---|
| elevator | nose-down moment, more lift | trailing edge **up** (nose up): negative command |
| aileron (antisymmetric) | right half more lift: rolls **left** | right wing down: negative command |
| rudder on an up-pointing fin | side force to the left: yaws **right** | trailing edge left (nose left): negative command |

`control_derivatives` returns derivatives per radian of this *command*.

## Derivative normalization

* `beta` derivatives: per radian.
* Rate derivatives: per non-dimensional rate, `p_hat = p b / 2V`,
  `q_hat = q c / 2V`, `r_hat = r b / 2V`.
* Control derivatives: per radian of command.

## Dynamics

The plants are small-perturbation, level-flight (`gamma0 = 0`), stability-axes
models. Longitudinal state `[u, alpha, q, theta]` with inputs
`[elevator, vertical gust]` (gust positive = updraft); lateral state
`[beta, p, r, phi]` with inputs `[aileron, rudder, sidewash]` (sidewash
positive = air arriving from the right). Thrust is held constant in the
phugoid model.

## CLI output units

Every column the command-line tool writes has a fixed unit. The text and CSV
headers carry the names; this table is the reference.

| Output | Column | Unit |
|---|---|---|
| `wing` CSV | `alpha_deg` | deg |
| `wing` CSV | `cl`, `cdi`, `cm`, `cl_roll`, `cn`, `cd` | dimensionless (`cl_roll`, `cn` body-axis, see above) |
| `wing` CSV | `lift_n` | N |
| `wing` CSV | `l_over_d` | dimensionless, blank at zero drag |
| `wing --spanwise` | `y_m`, `chord_m` | m |
| `wing --spanwise` | `cl_section` | dimensionless |
| `trim` (text, CSV) | `speed_ms` | m/s |
| `trim` (text, CSV) | `density` | kg/m³ |
| `trim` (text, CSV) | `bank_deg`, `alpha_deg`, `beta_deg`, `elevator_deg`, `aileron_deg`, `rudder_deg` | deg |
| `trim` (text, CSV) | `thrust_n` | N |
| `trim` (text, CSV) | `load_factor` | g |
| `trim` (text, CSV) | `static_margin_mac` | fraction of MAC, positive = stable |
| `trim` (text, CSV) | `cl_max_margin` | C_Lmax/C_L − 1 |
| `trim` (text, CSV) | `turn_rate_deg_s` | deg/s |
| `trim` (text, CSV) | `turn_radius_m` | m |
| `rotor` sweeps | `rpm` | rpm |
| `rotor` sweeps | `collective_deg` | deg |
| `rotor` sweeps | `axial_speed_ms` | m/s |
| `rotor` sweeps | `advance_ratio` | J = π λ_c, with D = 2R |
| `rotor` sweeps | `thrust_n`, `torque_nm`, `power_w` | N, N·m, W |
| `rotor` sweeps | `ct`, `cq`, `cp` | dimensionless (ΩR and ρ A (ΩR)² normalised) |

The only imperial option is `isa --imperial`, which reads altitudes in feet.
Trim and wing inputs and outputs are SI. A full imperial layer is deliberately
not added: converting every case field doubles the surface for little gain,
since the case files are TOML with explicit SI names.
