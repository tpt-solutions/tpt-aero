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
    let download =
        VlmModel::from_surfaces(&[wing(), htail().with_incidence(-3.0 * DEG).unwrap()]).unwrap();
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

use tpt_aero_vlm::ControlSymmetry;

fn controlled_aircraft() -> VlmModel {
    let s = |y: f64, chord: f64| WingSection {
        le: Vec3::new(0.0, y, 0.0),
        chord,
        twist: 0.0,
    };
    let fin_planform = WingPlanform::from_sections(vec![s(0.0, 0.9), s(1.5, 0.5)]).unwrap();
    VlmModel::from_surfaces(&[
        wing().with_control("aileron", 0.75, 2.5, 4.8, ControlSymmetry::Antisymmetric),
        htail().with_control("elevator", 0.6, 0.0, 2.0, ControlSymmetry::Symmetric),
        Surface::vertical(fin_planform, 6, 3)
            .at(3.8, 0.0, 0.0)
            .with_control("rudder", 0.6, 0.0, 1.5, ControlSymmetry::Symmetric),
    ])
    .unwrap()
}

/// d(coefficients)/d(deflection) for one control by central difference:
/// (CL, Cm, Cy, Cl_roll_std, Cn_std).
fn control_derivative(m: &VlmModel, name: &str) -> [f64; 5] {
    let alpha = 3.0 * DEG;
    let index = m.control_index(name).expect("control exists");
    let h = 1.0e-3;
    let q_s = 0.5 * RHO * V * V * m.reference().area;
    let at = |d: f64| {
        let mut defl = vec![0.0; m.controls().len()];
        defl[index] = d;
        let sol = m
            .solve(alpha, V, RHO, &m.control_twists(&defl).unwrap())
            .unwrap();
        let (roll, yaw) = sol.stability_roll_yaw(alpha);
        [sol.cl, sol.cm, sol.force.y / q_s, roll, yaw]
    };
    let (hi, lo) = (at(h), at(-h));
    core::array::from_fn(|i| (hi[i] - lo[i]) / (2.0 * h))
}

#[test]
fn elevator_lifts_and_pitches_nose_down_for_trailing_edge_down() {
    let m = controlled_aircraft();
    let [cl, cm, cy, croll, cn] = control_derivative(&m, "elevator");
    assert!(cl > 0.3, "CL_de = {cl}");
    assert!(cm < -1.0, "Cm_de = {cm}");
    assert!(cy.abs() < 1e-6 && croll.abs() < 1e-6 && cn.abs() < 1e-6);
}

#[test]
fn aileron_rolls_with_adverse_yaw_and_leaves_the_longitudinal_axis_alone() {
    let m = controlled_aircraft();
    let [cl, cm, _cy, croll, cn] = control_derivative(&m, "aileron");
    assert!(cl.abs() < 1e-6 && cm.abs() < 1e-6, "CL {cl} Cm {cm}");
    // Right trailing edge down (more lift on the right) rolls left: negative
    // right-wing-down roll moment, and adverse yaw is nose-right (positive).
    assert!(croll < -0.2, "Cl_da = {croll}");
    assert!(cn > 0.0 && cn < 0.1 * croll.abs(), "Cn_da = {cn}");
}

#[test]
fn rudder_side_force_and_yaw_are_consistent() {
    let m = controlled_aircraft();
    let [cl, cm, cy, _croll, cn] = control_derivative(&m, "rudder");
    assert!(cl.abs() < 1e-6 && cm.abs() < 1e-6);
    // The fin sits behind the CG: a side force to the left yaws the nose right.
    assert!(cy < -0.1, "Cy_dr = {cy}");
    assert!(cn > 0.05, "Cn_dr = {cn}");
}

#[test]
fn control_twists_select_panels_and_respect_symmetry() {
    let m = controlled_aircraft();
    assert_eq!(m.controls().len(), 3);
    assert_eq!(m.control_index("aileron"), Some(0));
    assert_eq!(m.control_index("elevator"), Some(1));
    assert_eq!(m.control_index("rudder"), Some(2));
    assert_eq!(m.control_index("flap"), None);

    let aileron = m.control_twists(&[0.1, 0.0, 0.0]).unwrap();
    let elements = m.lattice().elements();
    let moved: Vec<_> = elements
        .iter()
        .zip(&aileron)
        .filter(|(_, t)| **t != 0.0)
        .collect();
    assert!(!moved.is_empty());
    for (e, t) in &moved {
        assert_eq!(e.surface, 0, "only wing panels move");
        assert!(e.chord_frac >= 0.75, "only panels aft of the hinge");
        assert!((2.5..=4.8).contains(&e.eta.abs()));
        // Antisymmetric: right half +, left half -.
        assert_eq!(t.signum(), e.eta.signum());
    }
    assert!(moved.iter().any(|(e, _)| e.eta > 0.0) && moved.iter().any(|(e, _)| e.eta < 0.0));

    let elevator = m.control_twists(&[0.0, 0.1, 0.0]).unwrap();
    assert!(elevator
        .iter()
        .zip(elements)
        .all(|(t, e)| (*t != 0.0) == (e.surface == 1 && e.chord_frac >= 0.6)));
    assert!(elevator.iter().all(|t| *t == 0.0 || *t == 0.1));

    assert!(matches!(
        m.control_twists(&[0.0]),
        Err(VlmError::DimensionMismatch)
    ));
}

#[test]
fn naca_camber_domain_is_checked_before_building() {
    use tpt_aero_vlm::{CamberModel, Lattice, LatticeError};
    let wing = WingPlanform::rectangular(10.0, 1.0);
    // p = 1 would divide by zero inside the camber line.
    let p_one = CamberModel::Naca4 { m: 0.02, p: 1.0 };
    assert_eq!(
        Lattice::build(&wing, 4, 2, p_one).err(),
        Some(LatticeError::InvalidCamber)
    );
    // p = 0 with non-zero camber would silently be treated as flat.
    let p_zero = CamberModel::Naca4 { m: 0.02, p: 0.0 };
    assert!(Lattice::build(&wing, 4, 2, p_zero).is_err());
    let ok = CamberModel::Naca4 { m: 0.02, p: 0.4 };
    assert!(Lattice::build(&wing, 4, 2, ok).is_ok());
    // m = 0 means a flat plate, so p is not used and any value is fine.
    let flat = CamberModel::Naca4 { m: 0.0, p: 1.0 };
    assert!(Lattice::build(&wing, 4, 2, flat).is_ok());
}

#[test]
fn incidence_shift_rejects_non_finite_angle() {
    let wing = WingPlanform::rectangular(10.0, 1.0);
    assert!(Surface::horizontal(wing, 4, 2)
        .with_incidence(f64::NAN)
        .is_err());
}

#[test]
fn lattice_reference_span_is_the_surface_span_not_a_multiple_of_chord_panels() {
    use tpt_aero_vlm::{CamberModel, Lattice};
    // The span used to be the first strip's width summed over every chordwise
    // panel, i.e. n_chord times too large.
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let lattice = Lattice::build(&wing, 8, 4, CamberModel::None).expect("valid");
    let model = VlmModel::from_lattice(lattice);
    let r = model.reference();
    assert!((r.span - 10.0).abs() < 1e-9, "span {}", r.span);
    assert!((r.area - 10.0).abs() < 1e-9, "area {}", r.area);
    assert!((r.chord - 1.0).abs() < 1e-9, "chord {}", r.chord);
}

#[test]
fn ground_below_the_geometry_is_rejected() {
    // A wing at z = 0 with the ground at z = -h solves for h > 0. A non-positive
    // ground height is rejected as invalid.
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let model =
        VlmModel::from_planform(&wing, 8, 2, tpt_aero_vlm::CamberModel::None).expect("valid");
    assert!(model
        .clone()
        .with_ground_height(2.0)
        .solve(0.05, 30.0, 1.225, &[])
        .is_ok());
    assert!(model
        .with_ground_height(0.0)
        .solve(0.05, 30.0, 1.225, &[])
        .is_err());
}

#[test]
fn sideslip_at_ninety_degrees_gives_finite_loads() {
    use tpt_aero_vlm::FlightState;
    // At β = 90° the wind is along ŷ; the lift direction used to be NaN.
    let model = VlmModel::from_planform(
        &WingPlanform::rectangular(10.0, 1.0),
        8,
        2,
        tpt_aero_vlm::CamberModel::None,
    )
    .expect("valid");
    let prepared = model.prepare().expect("prepares");
    let st = FlightState {
        alpha: 0.05,
        beta: core::f64::consts::FRAC_PI_2,
        p: 0.0,
        q: 0.0,
        r: 0.0,
    };
    let loads = prepared.solve_state(&st, 30.0, 1.225).expect("solves");
    assert!(loads.c_lift.is_finite() && loads.c_drag.is_finite() && loads.c_side.is_finite());
}

#[test]
fn zero_reference_area_is_rejected() {
    let model = VlmModel::from_planform(
        &WingPlanform::rectangular(10.0, 1.0),
        8,
        2,
        tpt_aero_vlm::CamberModel::None,
    )
    .expect("valid")
    .with_reference(tpt_aero_vlm::ReferenceGeometry {
        area: 0.0,
        chord: 1.0,
        span: 10.0,
        moment_center: tpt_aero_core::Vec3::zero(),
    });
    assert_eq!(
        model.solve(0.05, 30.0, 1.225, &[]).unwrap_err(),
        tpt_aero_vlm::VlmError::InvalidReference
    );
}

#[test]
fn too_few_strips_for_a_kinked_planform_is_an_error() {
    use tpt_aero_core::{Vec3, WingSection};
    use tpt_aero_vlm::{CamberModel, Lattice, LatticeError};
    let kinked = WingPlanform::from_sections(vec![
        WingSection {
            le: Vec3::new(0.0, -5.0, 0.0),
            chord: 1.0,
            twist: 0.0,
        },
        WingSection {
            le: Vec3::new(0.0, -2.0, 0.0),
            chord: 1.0,
            twist: 0.0,
        },
        WingSection {
            le: Vec3::new(0.5, 5.0, 0.0),
            chord: 0.4,
            twist: 0.0,
        },
    ])
    .expect("valid");
    // Two segments need at least two strips.
    assert_eq!(
        Lattice::build(&kinked, 1, 2, CamberModel::None).err(),
        Some(LatticeError::InvalidResolution)
    );
    assert!(Lattice::build(&kinked, 2, 2, CamberModel::None).is_ok());
}

#[test]
fn ground_effect_with_dihedral_and_a_fin_raises_lift_near_the_ground() {
    use tpt_aero_vlm::Surface;
    // A dihedral wing with a fin, over a ground plane: lift grows as the
    // ground comes closer, and stays finite at every height tried.
    let wing = WingPlanform::tapered(10.0, 1.0, 0.7, 0.2, 5.0_f64.to_radians());
    let fin = WingPlanform::rectangular(1.2, 0.8);
    let surfaces = [
        Surface::horizontal(wing, 16, 4).at(0.0, 0.0, 0.0),
        Surface::vertical(fin, 4, 2).at(4.0, 0.0, 0.0),
    ];
    let model = VlmModel::from_surfaces(&surfaces).expect("valid aircraft");
    let a = 3.0_f64.to_radians();
    let cl_at = |h: f64| {
        model
            .clone()
            .with_ground_height(h)
            .solve(a, 40.0, 1.225, &[])
            .expect("solves")
            .cl
    };
    let (far, mid, near) = (cl_at(50.0), cl_at(4.0), cl_at(1.5));
    assert!(far.is_finite() && mid.is_finite() && near.is_finite());
    assert!(
        near > mid && mid > far,
        "CL {far} (50 m) < {mid} (4 m) < {near} (1.5 m) expected"
    );
    // The fin is built symmetric about the body axis (its span is ±0.6 m in z),
    // so its lower half reaches z = −0.6 m. A ground plane at z = −0.5 m therefore
    // intersects the geometry, and the clearance check must refuse it. This is the
    // symmetric-fin item in todo.md; a one-sided fin would pass.
    assert_eq!(
        model
            .clone()
            .with_ground_height(0.5)
            .solve(a, 40.0, 1.225, &[])
            .unwrap_err(),
        tpt_aero_vlm::VlmError::InvalidGroundHeight
    );
}

#[test]
fn kernel_path_matches_direct_prepare_for_any_twist() {
    // A kernel built once must reproduce `prepare_with_twists` exactly for
    // every deflection, so the six-DOF loop can reuse it across control changes.
    let m = controlled_aircraft();
    let kernel = m.influence_kernel().expect("kernel");
    assert_eq!(kernel.len(), m.lattice().len());
    let state = FlightState {
        alpha: 3.0 * DEG,
        beta: 0.0,
        p: 0.0,
        q: 0.0,
        r: 0.0,
    };
    for (name, d) in [
        ("elevator", 0.0),
        ("elevator", 0.1),
        ("aileron", -0.05),
        ("rudder", 0.08),
    ] {
        let index = m.control_index(name).expect("control");
        let mut defl = vec![0.0; m.controls().len()];
        defl[index] = d;
        let twists = m.control_twists(&defl).unwrap();
        let direct = m
            .prepare_with_twists(&twists)
            .unwrap()
            .solve_state(&state, V, RHO)
            .unwrap();
        let via_kernel = m
            .prepare_with_kernel(&kernel, &twists)
            .unwrap()
            .solve_state(&state, V, RHO)
            .unwrap();
        assert!(
            (direct.c_lift - via_kernel.c_lift).abs() < 1e-12
                && (direct.c_pitch - via_kernel.c_pitch).abs() < 1e-12
                && (direct.c_side - via_kernel.c_side).abs() < 1e-12,
            "{name} = {d}: {direct:?} vs {via_kernel:?}"
        );
    }
}

#[test]
fn kernel_rejects_a_mismatched_lattice() {
    let m = controlled_aircraft();
    let other = VlmModel::from_surfaces(&[wing()]).unwrap();
    let kernel = other.influence_kernel().unwrap();
    assert_eq!(
        m.prepare_with_kernel(&kernel, &[]).unwrap_err(),
        VlmError::DimensionMismatch
    );
}
