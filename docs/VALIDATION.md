# Validation

What each part of tpt-aero has been checked against, and where it is weak.
Every row below is an automated test in the repository. **None of these are
comparisons with flight-test or wind-tunnel data**; they check the numerics
against analytic theory, independent formulations and regulatory formulas.

## Checked

| Area | Reference | Result |
|---|---|---|
| VLM, rectangular wing AR 10 | Prandtl lifting line (Fourier) | CL within 0.02 |
| VLM, tapered wings (taper 0.5 / AR 8, 0.3 / AR 12) | Lifting line, generalized to any chord | CL within 4 %, section loading within 7 %, Oswald e within 5 % |
| Lifting-line reference itself | Elliptic-wing limit | e = 1, CL = 2 pi AR/(AR+2) alpha within 1 % |
| VLM prepared sweep vs one-shot solve | Same model, two code paths | agree to 1e-9 |
| Rate derivatives `CLq`, `Cmq` | Independent incidence-perturbation method | agree to 1 % |
| Aileron roll authority | Lifting line, antisymmetric incidence | within 15 % (VLM ~10 % lower near the tips) |
| Aileron induced (adverse) yaw | Lifting line | correct sense, magnitude within 50 % |
| Wing-alone neutral point | Theory (quarter chord) | within 5 % of a chord |
| Roll damping `Clp`, AR 10 rectangular | Textbook range | -0.63 (inside -0.75 .. -0.45) |
| `Clr`, `Cnp` (observed, sign tested) | `CL/4`, `-CL/8` | 0.071 vs 0.07; -0.032 vs -0.035 |
| Dihedral effect | Theory | `Clb` negative, linear in dihedral |
| Panel method, NACA 0012 | Published inviscid Cp and thick-airfoil theory | `Cp_min` within 0.04 of -0.43; lift slope within 0.04 |
| BEMT | Momentum theory limits, mass-flow conservation | tests in `tpt-aero-rotor` |
| Phugoid | Lanchester (`wn = sqrt2 g/V`, `zeta = 1/(sqrt2 L/D)`) | ~13 % low in frequency, damping within ~25 %; `wn V` constant across speed |
| Gust response | FAR 25 alleviation factor `Kg = 0.88 mu/(5.3+mu)` at 12.5 chords | 0.59 vs 0.62 |
| Time integration | Analytic solutions; step halving | RK4 error < 1e-6 |
| Lateral modes | Spiral/dihedral trade, eigenvalue periods | simulation matches eigenvalues to 2-3 % |
| Symmetry invariants | Property tests, Kani harnesses | zero roll/yaw for symmetric flight; freestream reversal reverses lift |

## Known weaknesses

* **Induced drag is optimistic.** The VLM computes drag from near-field
  forces. The Oswald efficiency converges to within a few percent of lifting
  line but can exceed 1 (impossible for a planar wing) by up to ~6 % at low
  resolution. Do not read three-digit accuracy into `cdi`.
* **VLM lift is a few percent below lifting line** (about 4.5 % at AR 8), the
  finite-chord effect; aileron authority is ~10 % low for the same reason.
* **Small-angle wake.** Trailing vortices stay along body +x. Large angles of
  attack or sideslip are outside the model's accuracy.
* **Wing-alone `Cnr`, `Cnb`** are near zero because there is no profile drag
  or fuselage; with a fin the values are plausible but unvalidated.
* **Control surfaces** are an incidence-change approximation: no flap gap
  losses, no hinge moments, no stall.
* **The phugoid model** assumes constant thrust, `C_Lu = C_Du = C_mu = 0` and
  you supply the zero-lift drag `cd0`.
* **Aircraft scope.** Wing + tail + fin only: no fuselage, nacelles, propeller
  slipstream or ground effect.
* **The 2-D panel solver** uses constant-strength panels; the pressure on the
  panels next to a sharp trailing edge is not accurate, and the open-
  trailing-edge NACA generator (`naca4`) gives a spurious gap panel, so use
  `naca4_closed` for panel work.

## Not yet validated

* AVL or XFLR5 comparisons for swept, cambered or multi-surface configurations.
* BEMT against measured propeller or rotor data (Caradonna-Tung, UIUC).
* Any published light-aircraft data set (for example Navion or Cessna 172
  derivatives).

If you have such data and are willing to share it, a comparison test would be
the most valuable contribution.
