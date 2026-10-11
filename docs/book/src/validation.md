# Validation and limits

Every number here comes from an automated test in the repository. The full
table, with the reference for each row, is in `docs/VALIDATION.md`.

The checks compare against analytic theory, independent formulations and
published formulas. None is a comparison with flight-test or wind-tunnel data.
The AVL 3.52 checks compare two lattice codes on matched geometry (straight and
dihedral wings, a T-tail and fin); they are in `docs/VALIDATION.md`.

## Known limits

* Induced drag comes from near-field forces. Treat `cdi` as approximate.
* The VLM trailing wake is along body `+x` by default. Large incidences and
  sideslip need `WakeModel::Freestream`.
* A wing with only a fin as its lateral surface cannot separate sideslip from
  rudder deflection: both act through the same fin moment arm. Add dihedral or
  another lateral surface to trim in a turn.
* There is no fuselage, nacelle or viscous lift model in the VLM. Parasite drag
  is supplied through `C_D0`. The propeller slipstream is available as a linear
  vortex-cylinder field (`CylinderWake`, the default in the six-DOF example) or
  as a momentum-theory contracting tube (`Slipstream`). Neither feeds back into
  the VLM.
* The VLM converges first order in spanwise strips. Default grids sit 3–5 % from
  converged values, and a fin and tail that meet at the fin root give unreliable
  sideslip derivatives.
