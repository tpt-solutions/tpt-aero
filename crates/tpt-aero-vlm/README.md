# tpt-aero-vlm

Horseshoe-vortex lattice method.

* Wings and **whole aircraft**: `Surface::horizontal` / `Surface::vertical`
  composed with `VlmModel::from_surfaces` (downwash and sidewash interference
  come from the shared influence matrix).
* **Fast sweeps**: `VlmModel::prepare` factors the influence matrix once;
  `solve_sweep` then costs O(N^2) per angle (about 20x faster for a 17-angle
  sweep at 500 panels).
* **Sideslip and body rates**: `PreparedVlm::solve_state` with aircraft-standard
  signs and stability-axis moments.
* **Control surfaces** (`with_control`, `control_twists`) and planform twist.
* **Gradients**: `solve_generic` runs on dual numbers.
* **Export**: `export::to_vtk` writes the lattice with circulation for ParaView.

Potential flow, small angles, trailing vortices along body +x. Induced drag is
a near-field residual: see the validation notes before relying on it.

Part of the [tpt-aero](https://github.com/tpt-solutions/tpt-aero) workspace.
See [conventions](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/CONVENTIONS.md)
(frames and signs) and [validation](https://github.com/tpt-solutions/tpt-aero/blob/master/docs/VALIDATION.md)
(what is checked, and what is not). Licensed MIT OR Apache-2.0, (c) 2026 TPT Solutions.
