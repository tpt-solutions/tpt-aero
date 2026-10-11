//! Stability derivatives against AVL 3.52 (Drela & Youngren) on matched geometry.
//!
//! The AVL decks are in `docs/validation/avl/` and the run script is
//! `derivatives.in` (`avl <deck> < derivatives.in`, which runs `oper`, `a a 2`,
//! `x`, `st`). Both codes are horseshoe-vortex lattices with the bound vortex at
//! the quarter chord and the control point at three quarters, and both run at
//! alpha = 2 deg, so they should agree closely once the grids are fine enough.
//!
//! The lattice is first-order in the spanwise panel count here: the wing's
//! CL_alpha is 4.88 / 4.76 / 4.70 at 10 / 20 / 40 strips, against AVL's 4.63 at
//! every count. The bands are therefore tight on the refined T-tail model and
//! looser on the coarse wing-alone model, which is the one the documentation
//! quotes. The AVL yaw-rate-from-roll derivative `Cnp` is not banded: AVL itself
//! has not converged on it between 20 and 40 strips (-0.0068 and -0.0113).

use tpt_aero_core::{Vec3, WingPlanform, WingSection};
use tpt_aero_stability::lateral::lateral_fd;
use tpt_aero_stability::longitudinal::longitudinal_fd;
use tpt_aero_vlm::{Surface, VlmModel};

const ALPHA: f64 = 2.0 * core::f64::consts::PI / 180.0;
const V: f64 = 50.0;
const RHO: f64 = 1.225;

fn assert_rel(name: &str, got: f64, want: f64, tol: f64) {
    let rel = (got - want).abs() / want.abs();
    assert!(
        rel <= tol,
        "{name}: tpt-aero {got:.5} vs AVL {want:.5} ({:+.1} %, band ±{:.0} %)",
        100.0 * (got - want) / want,
        100.0 * tol
    );
}

/// Rectangular wing, span 10 m, chord 1.2 m, flat, 20 x 4 panels.
fn rectangular_wing() -> VlmModel {
    let wing = WingPlanform::rectangular(10.0, 1.2);
    VlmModel::from_surfaces(&[Surface::horizontal(wing, 20, 4)]).expect("valid wing")
}

/// The same wing with 5° of dihedral (tips raised).
fn dihedral_wing() -> VlmModel {
    let wing = WingPlanform::tapered(10.0, 1.2, 1.2, 0.0, 5.0_f64.to_radians());
    VlmModel::from_surfaces(&[Surface::horizontal(wing, 20, 4)]).expect("valid wing")
}

/// Wing with 5° dihedral, a horizontal tail raised to z = 1.8 m (a T-tail, clear
/// of the fin), and a fin from z = 0 to 1.5 m. `n` is the spanwise strip count
/// of (wing, tail, fin); the AVL deck uses the same counts.
fn t_tail(n: (usize, usize, usize)) -> VlmModel {
    let wing = WingPlanform::tapered(10.0, 1.2, 1.2, 0.0, 5.0_f64.to_radians());
    let tail = WingPlanform::tapered(4.0, 0.6, 0.6, 0.0, 0.0);
    let fin = WingPlanform::from_sections(vec![
        WingSection {
            le: Vec3::new(0.0, 0.0, 0.0),
            chord: 0.9,
            twist: 0.0,
        },
        WingSection {
            le: Vec3::new(0.0, 1.5, 0.0),
            chord: 0.5,
            twist: 0.0,
        },
    ])
    .expect("valid fin");
    VlmModel::from_surfaces(&[
        Surface::horizontal(wing, n.0, 3),
        Surface::horizontal(tail, n.1, 3).at(4.0, 0.0, 1.8),
        Surface::vertical(fin, n.2, 3).at(3.8, 0.0, 0.0),
    ])
    .expect("valid aircraft")
}

#[test]
fn rectangular_wing_lift_and_induced_drag_match_avl() {
    // AVL 3.52, wing_rectangular.avl, 20 x 4 panels: CLa 4.6279, CDind 0.0010297.
    let m = rectangular_wing();
    let d = longitudinal_fd(&m, ALPHA, V, RHO, 1e-4).expect("solve");
    let s = m.solve(ALPHA, V, RHO, &[]).expect("solve");
    assert_rel("wing C_La", d.cl_alpha, 4.627_851, 0.05);
    assert_rel("wing C_Di", s.cdi, 0.001_029_7, 0.03);
}

#[test]
fn dihedral_wing_side_force_and_dihedral_effect_match_avl() {
    // AVL 3.52, wing_dihedral.avl, 20 x 4 panels: CYb -0.026217, Clb -0.076758, Cnb -0.001637.
    let m = dihedral_wing();
    let l = lateral_fd(&m, ALPHA, V, RHO).expect("lateral");
    assert_rel("dihedral wing C_Yb", l.cy_beta, -0.026_217, 0.10);
    assert_rel("dihedral wing C_lb", l.cl_beta, -0.076_758, 0.10);
    assert!(
        (l.cn_beta - (-0.001_637)).abs() < 5.0e-4,
        "dihedral wing C_nb = {}",
        l.cn_beta
    );
}

#[test]
fn refined_t_tail_derivatives_match_avl() {
    // AVL 3.52, t_tail_fine.avl (wing 40, tail 20, fin 24 strips):
    // CLa 5.302540, Cma -2.126698, CYb -0.243056, Clb -0.092155, Cnb 0.078154,
    // Clp -0.544663, Clr 0.065340, Cnr -0.067860. Cmα is about the wing
    // quarter chord, the model's reference point.
    let m = t_tail((40, 20, 24));
    let d = longitudinal_fd(&m, ALPHA, V, RHO, 1e-4).expect("solve");
    let l = lateral_fd(&m, ALPHA, V, RHO).expect("lateral");
    assert_rel("T-tail C_La", d.cl_alpha, 5.302_540, 0.03);
    assert_rel("T-tail C_ma", d.cm_alpha, -2.126_698, 0.05);
    assert_rel("T-tail C_Yb", l.cy_beta, -0.243_056, 0.05);
    assert_rel("T-tail C_lb", l.cl_beta, -0.092_155, 0.05);
    assert_rel("T-tail C_nb", l.cn_beta, 0.078_154, 0.05);
    assert_rel("T-tail C_lp", l.cl_p, -0.544_663, 0.05);
    assert_rel("T-tail C_lr", l.cl_r, 0.065_340, 0.05);
    assert_rel("T-tail C_nr", l.cn_r, -0.067_860, 0.05);
}
