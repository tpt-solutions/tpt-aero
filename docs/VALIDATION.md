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
| Rudder yawing moment `C_nδr` | Fin side force times its moment arm | within 1 % |
| VLM compressibility (Göthert rule) | Helmbold-Diederich lift slope, Mach dependence as a ratio | within 2 % at M 0.5 and 0.7 |
| VLM cosine spacing | Strip widths, agreement with uniform spacing | tip strip < 25 % of centre strip; CL within 1 % |
| VLM freestream wake | Body-axis wake at zero incidence; small change at 5° | identical at 0°; CL within 3 % at 5° |
| VLM ground effect | Trends: far-field limit, CL rise, monotone in height | CL +13 % at h/b 0.1; induced drag per unit CL² roughly halves at h/b 0.1 (rule of thumb) |
| ISA temperature offset | Sea-level and standard-day agreement; hydrostatic balance | pressure within 1e-4 of ISA at ΔISA = 0; `dp/dh = −ρg` within 1e-3 |
| Sutherland viscosity | Sea-level reference values | μ within 0.1 %, ν within 0.7 % |
| Viscous drag | Blasius laminar; Schlichting mixed law; Hoerner form factor | Re = 10⁶ laminar 0.001328 exact; turbulent 0.00297 ± 5e-5 |
| Steady turn trim | Coordinated-turn load factor 1/cos φ; mirror symmetry in bank | load factor within 1e-4; aileron and rudder flip, α and elevator unchanged |
| Six-DOF simulation | Level trim is an equilibrium; pitch-rate kick damps; steady turn holds under integration | Δu < 0.05 m/s and Δθ < 2e-3 rad over 5 s; φ within 3 %, r within 5 %, Δh < 1 m over 3 s |
| Autotrim (integral action) | Steady pitch error of a proportional hold | integral error < 0.5 % of the command, under a tenth of the proportional error |
| Rotor CLI | Dimensional forms follow the coefficients | T = C_T ρ A (ΩR)² to 1e-6 |
| Munk fuselage lift, side force and moments | Definition: integral of `ρVα dS/dx`, sign of weathercock yaw | integrals exact to 2 %; signs tested |
| Additive loads in six-DOF | Constant roll moment gives `M / I_xx` | within 2 % |
| AVL 3.52, wing alone (rectangular, AR 8.33, 20 × 4 panels) | AVL stability derivatives on the same geometry (`docs/validation/avl/`) | C_Lα 4.76 vs 4.63 (+3 %, band 5 %); C_Di 0.001037 vs 0.001030 (+1 %, band 3 %) |
| AVL, wing with 5° dihedral | same | C_Yβ −0.0257 vs −0.0262; C_lβ −0.081 vs −0.077 (bands 10 %) |
| AVL, wing + T-tail + fin, refined grid (40 / 20 / 24 strips) | same | C_Lα +2 %, C_mα +3 %, C_Yβ −2 %, C_lβ +2 %, C_nβ −1 %, C_lp +4 %, C_lr +4 %, C_nr +2 % (bands 3–5 %). C_np not banded: AVL itself has not converged on it |
| Propeller slipstream, momentum tube (`Slipstream`) | Momentum theory: `T = 2ρA v_d (V + v_d)`; far wake `V + 2 v_d`; mass flow through the contracting tube | thrust balance and mass flow exact to 1e-9; on-axis increment equals `v_d` at the disk and `2 v_d` far downstream |
| Propeller slipstream, vortex cylinder (`CylinderWake`) | On-axis closed form `v_d(1 + z/√(z²+R²))` (exact for the cylinder of ring strength `2 v_d`); divergence-free; axisymmetric | on-axis to 1e-9 for z ≤ 10R; disk-plane `u_z = v_d` for r < R to 1e-6; ∂u_z/∂z + div_r(u_r) < 1e-4; quadrature converged to 1e-6 off-axis and 1e-4 at the wake edge |

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
* **Fin-only lateral control.** A configuration whose only lateral surface is a
  fin with a full-span rudder cannot separate sideslip from rudder: both act
  through the same fin moment arm, so the turn and sideslip trim equations are
  singular. The test aircraft adds dihedral to the wing to break the tie. A
  real trim needs a second lateral source at a different station.
* **Ground effect** is checked for trends and the rule-of-thumb induced-drag
  halving at h/b 0.1, not against a published formula. The Cheeseman–Bennett
  expression `(16h/b)²/(1 + (16h/b)²)` does not match this VLM for either a
  rectangular or an elliptic planform, so it was not used as a reference.
* **Six-DOF scope.** Thrust acts along the airspeed vector (the assumption of
  the trim solvers). The fuselage enters through parasite drag and a first-order
  Munk lift and pitching moment (`ρ V α dS/dx`), which is small for a typical
  fuselage. Propeller torque is applied as a reaction moment.
* **Slipstream models are linear or first-order.** `CylinderWake` is the
  linearized vortex-cylinder wake of a uniformly loaded disk (Jewel & Heyson,
  NASA Memo 4-15-59L, "a single vortex cylinder"). It has no wake contraction,
  so it is not Conway's exact nonlinear theory. Its on-axis values are exact to
  1e-9 up to z = 10R; the quadrature loses about 0.5 % by z = 40R. `Slipstream`
  adds the momentum-theory contraction on the axis and has no crossflow. The
  `propeller_flight` example uses the cylinder by default (`SLIPSTREAM`). At
  its trim point the two models give tail lift within 1.2 % of each other
  (36.6 vs 36.2 N). The tail's angle of attack is the wing's (no downwash),
  which overstates the tail share at trim, by an estimated 1.5× for that example.
  The VLM itself does not see the slipstream.
* **VLM convergence is first order in the spanwise strip count.** The wing's
  C_Lα is 4.88, 4.76 and 4.70 at 10, 20 and 40 strips, against AVL's 4.63 at every
  count. Default grids in the examples and CLI are 3–5 % from converged values.
  Refine before comparing against another code or data.
* **Fin and tail that meet at the fin root give unreliable C_Yβ and C_nβ.** The
  `tpt-aero trim` aircraft template puts the horizontal tail at z = 0, on the fin
  root. Both codes see the overlap, but on that layout tpt-aero's C_Yβ is 31 %
  and C_nβ 37 % away from AVL, while a T-tail clear of the fin agrees within
  about 2 %. Keep the tail clear of the fin.
* **Propeller sensitivity.** Near the trim point the example propeller's thrust
  is steep in shaft speed (a 1 % rpm change gives roughly 16 % more thrust), so
  throttle steps must be small. This is a property of the BEMT closure with a
  linear lift model, not checked against measured propeller data.
* **The phugoid model** assumes constant thrust, `C_Lu = C_Du = C_mu = 0` and
  you supply the zero-lift drag `cd0`.
* **Aircraft scope.** Wing + tail + fin only: no fuselage, nacelles or ground
  effect beyond the image model. The propeller slipstream is a tail-only estimate (above). The Munk fuselage terms in six-DOF are a separate
  first-order estimate and are not part of the VLM.
* **The 2-D panel solver** uses constant-strength panels; the pressure on the
  panels next to a sharp trailing edge is not accurate, and the open-
  trailing-edge NACA generator (`naca4`) gives a spurious gap panel, so use
  `naca4_closed` for panel work.

## Compared with published data (results recorded, not agreement)

| Case | Source | Model | Measurement | Result |
|---|---|---|---|---|
| Caradonna-Tung rotor, hover, collective 8°, tip Mach 0.877, `C_T` | Caradonna & Tung, NASA TM-81232 (1981); value as quoted in [arXiv:2107.13895](https://arxiv.org/pdf/2107.13895) (the NASA table itself was not checked) | BEMT, linear lift closure, no Mach or stall correction | 0.00473 | 0.00584: +23 %. Tested within +30 % (`tpt-aero-rotor/tests/caradonna_tung.rs`) |
| Navion wing and tail, C_Lα at 73.2 m/s | Seckel & Morris, NASA TN D-6643, Table VII (clean), Table I and Figure 2 | VLM wing + tail, no fuselage; tail arm 4.7 m read from Figure 2 | C_Zα = −4.33 | C_Lα 4.98: +15 %. Band +5 % to +25 % (`tpt-aero-stability/tests/navion.rs`) |
| Navion, C_mα about 25 % c | same | same | −0.77 | −1.37: 1.8× the flight value. Band 1.4× to 2.2× |

The Navion C_mα gap is expected: the model has no fuselage, which is normally destabilising. The C_Lα gap is not explained by the missing fuselage, which would add lift. The flight value is a variable-stability extraction and the cause is not resolved here.

## Not yet validated

* Comparisons with AVL for swept, cambered, or control-surface configurations. The AVL checks above use straight, flat or dihedral wings with a tail and fin. XFLR5 has not been run.
* UIUC or other measured rotor data beyond the Caradonna-Tung hover point above.
* Cessna 172 or other light-aircraft derivative sets with geometry.
* The compressible, ground-effect and six-DOF results have not been compared with published wind-tunnel or flight-test data for any configuration.

If you have such data and are willing to share it, a comparison test would be
the most valuable contribution.
