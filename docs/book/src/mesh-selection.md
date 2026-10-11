# Mesh selection

A vortex lattice converges with the number of spanwise strips (`mesh.n_span`)
and chordwise panels (`mesh.n_chord`). The lift converges first, then the
induced drag, then the pitching moment.

## Rules of thumb

- Start with `n_span` at 20 to 30 for a wing and 4 chordwise panels. This
  resolves the trailing-vortex system for most planforms.
- Check convergence with `tpt-aero wing case.toml --check-mesh` before you
  quote a number. It doubles `n_span` up to three times and reports the first
  level that moves the lift by less than 1 %.
- Give the case an explicit `n_span`. The check edits that line.
- Cosine spanwise spacing (`SpanSpacing::Cosine`) resolves the wing tips better
  at the same count. It is available from the library.
- The lattice needs at least one strip per planform segment, so a wing with a
  kink needs at least as many strips as it has sections minus one.
- The strip count is capped at 500 per direction and 2500 panels per surface.
  The dense influence matrix grows with the square of the panel count, so the
  cap keeps a run inside memory.

## Chordwise panels

Four panels give the lift slope to a few per cent. Eight is the point where
the zero-lift angle stops moving; the zero-lift study in
[validation](validation.md) shows the convergence.

## When the mesh is not enough

If `--check-mesh` does not converge, the problem is usually the geometry, not
the count: a sharp tip, a very thin tip chord, or a control hinge at a panel
boundary. The taper and aspect-ratio warnings flag the first two cases. A
hinge on a panel boundary is not yet warned about; it is a known limit, listed
in `todo.md`.
