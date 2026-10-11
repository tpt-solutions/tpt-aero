# Introduction

tpt-aero is a pure-Rust library and command-line tool for aerodynamics and
flight mechanics. It takes a planform through to handling qualities in one
workspace:

* **Wings and whole aircraft** with the vortex lattice method (VLM): lift,
  induced drag, moments, spanwise loading, control surfaces, and the
  stability derivatives they produce.
* **2-D airfoils** with source and doublet panel methods and compressibility
  corrections.
* **Propellers and rotors** with blade-element momentum theory (BEMT).
* **The atmosphere**: the standard atmosphere (with temperature offsets and
  viscosity), wind shear, discrete gusts and continuous turbulence.
* **Flight mechanics**: trim, steady turns, linear and nonlinear (six-degree-of-
  freedom) simulation, flight modes and closed-loop control.

There are no C or C++ dependencies and no external solver binaries. Every
influence matrix and linear solve is plain Rust. The solver crates are
`no_std` with `alloc`.

## Where to go next

* New to the tool? Start with [installing](install.md), then the
  [command line](cli.md).
* Writing Rust? The [library tour](library.md) lists the modules and the entry
  points for each analysis.
* Want to know how far to trust a number? Read [validation and limits](validation.md)
  first.
