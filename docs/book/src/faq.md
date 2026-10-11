# FAQ

**Why are my lift and drag numbers different from XFLR5 or AVL?**
Inviscid VLM and XFLR5's viscous panel method disagree on profile drag and on
the lift near stall; the VLM has no stall at all. Compare lift slopes and
induced drag, which should agree within a few per cent on a converged mesh.

**The static margin is negative. Is the solver wrong?**
Probably not. Check the CG: `mass.cg_x` is the moment reference, and a CG
behind the neutral point is unstable. The `report` command prints the neutral
point distance, so you can see how far to move the CG.

**Why is `l_over_d` blank in the wing output?**
At zero drag the ratio is undefined, which happens at α = 0 for a symmetric
wing with no `cd0`. The cell is left empty on purpose.

**Why does trim fail with `OutsideEnvelope`?**
The trimmed angle of attack is beyond ±20°, where the attached-flow linear VLM
is no longer meaningful. Move the CG, add tail area, or lower the design lift.

**Did compressible results change?**
Yes, for one case. Supplying a physical reference to a compressible model used
to count β twice. The reference is now stretched internally, so a physical
area gives the same coefficients as the default. Rerun any compressible result
made with `with_reference` before 0.1.0.

**Can I use a fin that is one-sided?**
The `vertical` surface kind is documented as one-sided, but the CLI currently
builds every surface symmetric about the body axis. This is an open item; see
`todo.md` before comparing fin-heavy results with AVL.

**Are the mode frequencies right?**
They depend on the inertias. Without `mass.inertia_x/y/z` the `report` command
uses estimates from the span and says so. Use measured inertias for anything
that matters.

**Is the library `no_std`?**
The core crates are `no_std` with `alloc`, and CI builds them for an embedded
target. The CLI and the examples need `std`.
