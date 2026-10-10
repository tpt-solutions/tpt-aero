//! Multi-surface configurations: wing + tailplane + fin.

use tpt_aero_core::{Vec3, WingPlanform, WingSection};
use tpt_aero_vlm::{FlightState, Lattice, LatticeError, Surface, VlmError, VlmModel};

const V: f64 = 50.0;
const RHO: f64 = 1.225;
const DEG: f64 = core::f64::consts::PI / 180.0;

fn wing() -> Surface {
    Surface::horizontal(WingPlanform::rectangular(10.0, 1.2), 20, 3)
}

fn htail() -> Surface {
    Surface::horizontal(WingPlanform::rectangular(4.0, 0.6), 10, 3).at(4.0, 0.0, 0.0)
}

fn fin() -> Surface {
    let s = |y: f64, chord: f64| WingSection {
        le: Vec3::new(0.0, y, 0.0),
        chord,
        twist: 0.0,
    };
    let planform = WingPlanform::from_sections(vec![s(0.0, 0.9), s(1.5, 0.5)]).unwrap();
    Surface::vertical(planform, 6, 3).at(3.8, 0.0, 0.0)
}

/// Lift-curve and pitch-moment slopes by central difference on `solve`.
fn longitudinal_slopes(m: &VlmModel, alpha: f64) -> (f64, f64) {
    let h = 1.0e-4;
    let at = |a: f64| m.solve(a, V, RHO, &[]).unwrap();
    let (hi, lo) = (at(alpha + h), at(alpha - h));
    ((hi.cl - lo.cl) / (2.0 * h), (hi.cm - lo.cm) / (2.0 * h))
}

fn beta_slopes(m: &VlmModel, alpha: f64) -> (f64, f64, f64) {
    let prepared = m.prepare().unwrap();
    let h = 1.0e-3;
    let at = |beta: f64| {
        prepared
            .solve_state(
                &FlightState {
                    alpha,
                    beta,
                    ..FlightState::default()
                },
                V,
                RHO,
            )
            .unwrap()
    };
    let (hi, lo) = (at(h), at(-h));
    let d = |a: f64, b: f64| (a - b) / (2.0 * h);
    (
        d(hi.c_side, lo.c_side),
        d(hi.c_roll, lo.c_roll),
        d(hi.c_yaw, lo.c_yaw),
    )
}

#[test]
fn wing_alone_neutral_point_is_at_the_quarter_chord() {
    let m = VlmModel::from_surfaces(&[wing()]).unwrap();
    let r = m.reference();
    assert!((r.moment_center.x - 0.3).abs() < 1e-12, "reference at c/4");
    let (cla, cma) = longitudinal_slopes(&m, 3.0 * DEG);
    let x_np_offset = -cma / cla * r.chord;
    assert!(
        x_np_offset.abs() < 0.05 * r.chord,
        "NP offset {x_np_offset}"
    );
}

#[test]
fn tailplane_adds_stability_and_is_reduced_by_the_wing_downwash() {
    let alpha = 3.0 * DEG;
    let wing_only = VlmModel::from_surfaces(&[wing()]).unwrap();
    let with_tail = VlmModel::from_surfaces(&[wing(), htail()]).unwrap();
    let (cla_w, _) = longitudinal_slopes(&wing_only, alpha);
    let (cla_wt, cma_wt) = longitudinal_slopes(&with_tail, alpha);

    assert!(
        cma_wt < -1.0,
        "statically stable with a tail: Cma = {cma_wt}"
    );
    let r = with_tail.reference();
    let np_aft = -cma_wt / cla_wt * r.chord;
    assert!(
        np_aft > 0.2,
        "neutral point moves aft of the CG: {np_aft} m"
    );

    // The tail's own lift-curve slope (in the wing-sized coefficient) and the
    // fraction of it that survives the wing's downwash (1 - d eps/d alpha).
    let tail_alone = VlmModel::from_surfaces(&[htail()]).unwrap();
    let (cla_t, _) = longitudinal_slopes(&tail_alone, alpha);
    let tail_contribution = cla_t * tail_alone.reference().area / r.area;
    let gained = cla_wt - cla_w;
    let survives = gained / tail_contribution;
    assert!(
        (0.45..=0.95).contains(&survives),
        "downwash leaves {survives:.2} of the isolated tail lift slope"
    );
}

#[test]
fn tail_incidence_changes_the_pitching_moment_in_the_expected_direction() {
    let level = VlmModel::from_surfaces(&[wing(), htail()]).unwrap();
    // Trailing-edge-down tail (positive incidence) lifts; negative incidence
    // carries download, which pitches the nose up about a CG ahead of it.
    let download = VlmModel::from_surfaces(&[wing(), htail().with_incidence(-3.0 * DEG)]).unwrap();
    let a = 0.0;
    let c_level = level.solve(a, V, RHO, &[]).unwrap();
    let c_down = download.solve(a, V, RHO, &[]).unwrap();
    assert!(c_down.cl < c_level.cl, "download reduces total lift");
    assert!(c_down.cm > c_level.cm, "download pitches the nose up");
}

#[test]
fn fin_provides_weathercock_stability_and_side_force() {
    let alpha = 3.0 * DEG;
    let no_fin = VlmModel::from_surfaces(&[wing(), htail()]).unwrap();
    let with_fin = VlmModel::from_surfaces(&[wing(), htail(), fin()]).unwrap();
    let (_, _, cn_beta_none) = beta_slopes(&no_fin, alpha);
    let (cy_beta, _, cn_beta) = beta_slopes(&with_fin, alpha);
    assert!(
        cn_beta_none.abs() < 1e-6,
        "no fin, no Cn_beta: {cn_beta_none}"
    );
    assert!(cn_beta > 0.05, "Cn_beta = {cn_beta}");
    assert!(cy_beta < -0.15, "Cy_beta = {cy_beta}");
    // The fin does not disturb the longitudinal behaviour in symmetric flight.
    let (cla_a, cma_a) = longitudinal_slopes(&no_fin, alpha);
    let (cla_b, cma_b) = longitudinal_slopes(&with_fin, alpha);
    assert!((cla_a - cla_b).abs() < 0.02 * cla_a);
    assert!((cma_a - cma_b).abs() < 0.02 * cma_a.abs());
}

#[test]
fn symmetric_flight_produces_no_roll_or_yaw_for_a_full_aircraft() {
    let m = VlmModel::from_surfaces(&[wing(), htail(), fin()]).unwrap();
    let s = m.solve(4.0 * DEG, V, RHO, &[]).unwrap();
    assert!(s.cl > 0.2);
    assert!(s.cl_roll.abs() < 1e-9 && s.cn.abs() < 1e-9);
    // And the prepared path agrees with the one-shot solve.
    let prepared = m.prepare().unwrap();
    let p = prepared.solve(4.0 * DEG, V, RHO).unwrap();
    assert!((p.cl - s.cl).abs() < 1e-9 && (p.cm - s.cm).abs() < 1e-9);
}

#[test]
fn rolled_surface_geometry_is_a_rigid_motion() {
    let lat = Lattice::build(
        &WingPlanform::from_sections(vec![
            WingSection {
                le: Vec3::new(0.0, 0.0, 0.0),
                chord: 1.0,
                twist: 0.0,
            },
            WingSection {
                le: Vec3::new(0.0, 2.0, 0.0),
                chord: 1.0,
                twist: 0.0,
            },
        ])
        .unwrap(),
        4,
        2,
        Default::default(),
    )
    .unwrap();
    let up = lat.transformed(core::f64::consts::FRAC_PI_2, Vec3::new(1.0, 0.0, 0.5));
    for (before, after) in lat.elements().iter().zip(up.elements()) {
        // y -> z, translated; lengths preserved.
        assert!((after.bound_a.y).abs() < 1e-12 && (after.bound_b.y).abs() < 1e-12);
        assert!((after.bound_a.z - (0.5 + before.bound_a.y)).abs() < 1e-12);
        assert!((after.bound_a.x - (1.0 + before.bound_a.x)).abs() < 1e-12);
        let len = |a: Vec3<f64>, b: Vec3<f64>| (b - a).norm();
        assert!(
            (len(before.bound_a, before.bound_b) - len(after.bound_a, after.bound_b)).abs() < 1e-12
        );
        // The wing normal (+z) becomes -y; the span axis (+y) becomes +z.
        assert!((after.flat_normal.y + 1.0).abs() < 1e-12);
        assert!((after.span_axis.z - 1.0).abs() < 1e-12);
    }
}

#[test]
fn combine_keeps_strips_separate_and_rejects_empty_input() {
    let a = Lattice::build(
        &WingPlanform::rectangular(4.0, 1.0),
        3,
        2,
        Default::default(),
    )
    .unwrap();
    let b = Lattice::build(
        &WingPlanform::rectangular(2.0, 0.5),
        5,
        2,
        Default::default(),
    )
    .unwrap();
    let both = Lattice::combine(&[a.clone(), b.clone()]).unwrap();
    assert_eq!(both.len(), a.len() + b.len());
    assert_eq!(both.n_span(), 8);
    let max_strip = both.elements().iter().map(|e| e.strip).max().unwrap();
    assert_eq!(max_strip, 7);
    // The second surface's strips start after the first's.
    assert!(both.elements()[a.len()..].iter().all(|e| e.strip >= 3));
    assert_eq!(
        Lattice::combine(&[]).unwrap_err(),
        LatticeError::InvalidResolution
    );
    assert!(matches!(
        VlmModel::from_surfaces(&[]),
        Err(VlmError::Lattice(LatticeError::InvalidResolution))
    ));
}
