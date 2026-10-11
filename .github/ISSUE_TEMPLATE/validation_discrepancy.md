---
name: Validation discrepancy
about: A result disagrees with a reference (AVL, XFLR5, a textbook formula, measured data)
labels: validation
---

**Reference and its version**
Name the source, the version, and the exact figure or table used. Attach the input deck if it is an AVL or XFLR5 comparison.

**tpt-aero input**
The case file or code that produced the number, with the flight condition (speed, density, alpha, beta).

**Numbers**
| Quantity | tpt-aero | Reference | Difference |
|---|---|---|---|

**Mesh**
`n_span`, `n_chord`, and whether a convergence check was run. See [docs/VALIDATION.md](../../docs/VALIDATION.md) for the bands already recorded.
