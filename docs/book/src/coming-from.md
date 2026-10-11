# Coming from other tools

tpt-aero covers the same ground as AVL, XFLR5, OpenVSP and JSBSim, but it is a
library and a command line, not a GUI, and it does not import their geometry
files yet. This page maps the concepts so you can move a model across by hand.

## From AVL

AVL is the closest match: a vortex lattice with stability derivatives.

| AVL | tpt-aero |
|---|---|
| `Sref`, `Cref`, `Bref` in the header | `ReferenceGeometry` (`area`, `chord`, `span`); set by the planform, or override with `with_reference` |
| `Xref Yref Zref` (moment reference) | `moment_center`; for an aircraft case, `mass.cg_x` |
| `SURFACE` with `NACA`, `SECTION` (`Xle Yle Zle Chord Ainc`) | `WingPlanform` sections (`le`, `chord`, `twist`); `CamberModel::Naca4` |
| `CONTROL` on a section | `with_control(name, hinge, y_in, y_out, symmetry)` on a `Surface` |
| Trefftz plane, `alpha` sweep | `VlmModel::solve`, `solve_sweep`; `wing` case `[conditions] alpha` |
| `ELEVATOR`, `AILERON`, `RUDDER` derivatives | `control_derivatives` for any named control |

Sign conventions follow [`docs/CONVENTIONS.md`](../../CONVENTIONS.md): alpha
nose-up positive, roll positive right-wing-down, yaw positive nose-right. The
AVL comparisons in `crates/tpt-aero-stability/tests/avl_reference.rs` agree to
1 to 4 % on a converged grid (see [validation](validation.md)). An `.avl`
importer is not written yet; the decks in `docs/validation/avl/` show the
geometry to copy.

## From XFLR5

XFLR5 polars are the wing sweep. The `wing` command writes the same kind of
table: `alpha_deg,cl,cdi,cm,cl_roll,cn,lift_n,cd,l_over_d`. The `cd` column
adds your parasite drag (`conditions.cd0`) to the induced drag, which is what
XFLR5's total drag column gives you.

XFLR5 uses a panel method with viscous airfoil polars. tpt-aero's wing solver
is inviscid, so profile drag comes from `cd0` or the `viscous` module
(`flat_plate_cf` and `cd0_buildup`). Use XFLR5's polars for airfoils that need
stall or viscous detail; tpt-aero's `AirfoilPolar` reads that table directly.

## From OpenVSP

OpenVSP builds geometry; tpt-aero solves it. Export a planform by hand: the
span, the root and tip chords, the sweep at the leading edge, and the
dihedral. Each lifting surface becomes a `Surface` with its own offset (`x`,
`y`, `z`) in body axes. Fuselage and nacelle shapes are not modelled as
geometry; `SlenderBody` and the fuselage drag in `sixdof` cover their first
order effects.

## From JSBSim

JSBSim is a six-degree-of-freedom simulator driven by aerodynamic coefficient
tables. tpt-aero's `sixdof` module integrates the same equations (aircraft
axes, Euler angles, RK4) with the VLM supplying the aerodynamics, so the
coefficient tables are not needed. `sixdof::simulate_partial` returns the
history up to any failure, such as `GimbalLock` near a 90° pitch.

For a linear model of an aircraft in JSBSim's trim, `trim_level_flight` and
`longitudinal_fd` give the same trim point and derivatives.
