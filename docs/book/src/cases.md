# Case files

Case files are TOML and reject unknown fields, so a typo is an error rather
than a silently ignored setting. Print a full, commented template with
`tpt-aero template <kind>`.

## `wing`

A single tapered wing swept through angle of attack:

```toml
[wing]
span = 10.0
root_chord = 1.0
tip_chord = 1.0
dihedral_deg = 0.0
# camber = { m = 0.02, p = 0.4 }

[mesh]
n_span = 24
n_chord = 4

[conditions]
speed = 50.0
altitude = 0.0
alpha = { start = -4.0, stop = 10.0, step = 2.0 }
```

## `rotor`

One blade-element operating point, with a linear section closure
(`c_l = cl_slope · α`, clamped to `±cl_max`). The rotor's `hub_radius`, twist
and tip loss are optional with defaults. `speed` is the axial inflow in m/s:
zero for hover, positive for a propeller.

## `aircraft`

A list of `[[surface]]` tables, each with optional `[[surface.control]]`
tables. Surfaces are `horizontal` (wings, tails) or `vertical` (fins). A
control is named (`elevator`, `aileron`, `rudder`), starts at a chord
fraction `hinge` and runs from `inboard` to `outboard` metres from the
centreline. `symmetry` is `symmetric` (elevator, rudder) or `antisymmetric`
(aileron).

Level trim needs an `elevator`. A turn trim needs an `elevator`, an `aileron`
and a `rudder`, all with these names.
