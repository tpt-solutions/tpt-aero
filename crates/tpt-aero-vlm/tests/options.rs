//! Opt-in solver options: cosine spanwise spacing and the freestream wake.

use tpt_aero_core::WingPlanform;
use tpt_aero_vlm::{CamberModel, SpanSpacing, VlmError, VlmModel, WakeModel};

const RHO: f64 = 1.225;
const V: f64 = 50.0;

fn rectangular(n_span: usize, spacing: SpanSpacing) -> VlmModel {
    let wing = WingPlanform::rectangular(10.0, 1.0);
    VlmModel::from_planform_spaced(&wing, n_span, 4, CamberModel::None, spacing).expect("valid")
}

#[test]
fn cosine_spacing_clusters_strips_toward_the_tips() {
    // A full-span wing: strips run tip to tip, so both end strips are tip
    // strips. Cosine spacing makes them much narrower than the middle ones.
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let uniform = VlmModel::from_planform(&wing, 16, 4, CamberModel::None).expect("valid");
    let cosine = rectangular(16, SpanSpacing::Cosine);

    let widths = |m: &VlmModel| -> Vec<f64> {
        (0..m.lattice().n_span())
            .map(|k| {
                m.lattice()
                    .elements()
                    .iter()
                    .find(|e| e.strip == k)
                    .expect("strip")
                    .strip_dy
            })
            .collect()
    };
    let spread = |w: &[f64]| {
        let max = w.iter().copied().fold(f64::MIN, f64::max);
        let min = w.iter().copied().fold(f64::MAX, f64::min);
        (min, max)
    };
    let (u_min, u_max) = spread(&widths(&uniform));
    assert!(
        (u_max - u_min).abs() < 1e-12,
        "uniform strips must be equal"
    );

    let c = widths(&cosine);
    let (c_min, c_max) = spread(&c);
    assert!(
        c_min < 0.25 * c_max,
        "tip strip {c_min} vs centre strip {c_max}"
    );
    // Symmetric about the centre of the span.
    assert!((c[0] - c[c.len() - 1]).abs() < 1e-12);
}

#[test]
fn cosine_spacing_agrees_with_uniform_on_a_rectangular_wing() {
    let a = 5.0_f64.to_radians();
    let uniform = rectangular(24, SpanSpacing::Uniform)
        .solve(a, V, RHO, &[])
        .unwrap();
    let cosine = rectangular(24, SpanSpacing::Cosine)
        .solve(a, V, RHO, &[])
        .unwrap();
    assert!(
        (cosine.cl / uniform.cl - 1.0).abs() < 0.01,
        "{} vs {}",
        cosine.cl,
        uniform.cl
    );
    assert!((cosine.cdi / uniform.cdi - 1.0).abs() < 0.05);
}

#[test]
fn default_spacing_is_uniform() {
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let default = VlmModel::from_planform(&wing, 12, 4, CamberModel::None).unwrap();
    let explicit = rectangular(12, SpanSpacing::Uniform);
    let a = 4.0_f64.to_radians();
    let s1 = default.solve(a, V, RHO, &[]).unwrap();
    let s2 = explicit.solve(a, V, RHO, &[]).unwrap();
    assert_eq!(s1.cl, s2.cl);
}

#[test]
fn freestream_wake_matches_body_axis_at_zero_incidence() {
    let body = rectangular(12, SpanSpacing::Uniform);
    let free = body.clone().with_wake(WakeModel::Freestream);
    assert_eq!(free.wake(), WakeModel::Freestream);
    let s1 = body.solve(0.0, V, RHO, &[]).unwrap();
    let s2 = free.solve(0.0, V, RHO, &[]).unwrap();
    assert!((s1.cl - s2.cl).abs() < 1e-12, "{} vs {}", s1.cl, s2.cl);
}

#[test]
fn freestream_wake_is_a_small_correction_at_moderate_incidence() {
    let body = rectangular(16, SpanSpacing::Uniform);
    let free = body.clone().with_wake(WakeModel::Freestream);
    let a = 5.0_f64.to_radians();
    let s1 = body.solve(a, V, RHO, &[]).unwrap();
    let s2 = free.solve(a, V, RHO, &[]).unwrap();
    assert!(s2.cl > 0.0);
    assert!((s2.cl / s1.cl - 1.0).abs() < 0.03, "{} vs {}", s2.cl, s1.cl);
}

#[test]
fn freestream_wake_rejects_cached_paths_but_sweeps_per_angle() {
    let free = rectangular(10, SpanSpacing::Uniform).with_wake(WakeModel::Freestream);
    assert!(matches!(
        free.prepare(),
        Err(VlmError::FreestreamWakeNotPrepared)
    ));

    let angles = [0.0, 2.0_f64.to_radians(), 4.0_f64.to_radians()];
    let sweep = free.solve_sweep(&angles, V, RHO).unwrap();
    for (sol, &a) in sweep.iter().zip(&angles) {
        let direct = free.solve(a, V, RHO, &[]).unwrap();
        assert_eq!(sol.cl, direct.cl);
    }
}

#[test]
fn ground_effect_vanishes_far_from_the_ground() {
    let free = rectangular(16, SpanSpacing::Uniform);
    let far = free.clone().with_ground_height(1.0e3);
    let a = 5.0_f64.to_radians();
    let s1 = free.solve(a, V, RHO, &[]).unwrap();
    let s2 = far.solve(a, V, RHO, &[]).unwrap();
    assert!((s2.cl / s1.cl - 1.0).abs() < 1e-3);
    assert!((s2.cdi / s1.cdi - 1.0).abs() < 1e-2);
}

#[test]
fn ground_effect_raises_lift_and_cuts_induced_drag_near_the_ground() {
    // Ten-metre span: h/b = 0.1 is a classic ground-effect distance.
    let free = rectangular(16, SpanSpacing::Uniform);
    let near = free.clone().with_ground_height(1.0);
    let a = 5.0_f64.to_radians();
    let s1 = free.solve(a, V, RHO, &[]).unwrap();
    let s2 = near.solve(a, V, RHO, &[]).unwrap();
    assert!(s2.cl > 1.05 * s1.cl, "CL {} vs free {}", s2.cl, s1.cl);
    // Induced drag at equal lift, k = CDi / CL²: it roughly halves at h/b = 0.1
    // (the rule of thumb for ground effect); the bound is deliberately loose.
    let k_ratio = (s2.cdi / (s2.cl * s2.cl)) / (s1.cdi / (s1.cl * s1.cl));
    assert!((0.35..0.6).contains(&k_ratio), "k ratio {k_ratio}");
}

#[test]
fn ground_effect_is_monotone_in_height() {
    let free = rectangular(16, SpanSpacing::Uniform);
    let a = 5.0_f64.to_radians();
    let cl_at = |h: f64| {
        free.clone()
            .with_ground_height(h)
            .solve(a, V, RHO, &[])
            .unwrap()
            .cl
    };
    let (c05, c1, c2, c5) = (cl_at(0.5), cl_at(1.0), cl_at(2.0), cl_at(5.0));
    assert!(c05 > c1 && c1 > c2 && c2 > c5, "{c05} {c1} {c2} {c5}");
}

#[test]
fn ground_effect_is_consistent_between_solve_and_prepared_paths() {
    let model = rectangular(12, SpanSpacing::Uniform).with_ground_height(1.5);
    let a = 4.0_f64.to_radians();
    let direct = model.solve(a, V, RHO, &[]).unwrap();
    let prepared = model.prepare().unwrap().solve(a, V, RHO).unwrap();
    assert!(
        (direct.cl - prepared.cl).abs() < 1e-9,
        "{} vs {}",
        direct.cl,
        prepared.cl
    );
}

#[test]
fn invalid_ground_height_is_rejected() {
    let bad = rectangular(8, SpanSpacing::Uniform).with_ground_height(-1.0);
    assert_eq!(
        bad.solve(0.05, V, RHO, &[]).unwrap_err(),
        VlmError::InvalidGroundHeight
    );
    assert_eq!(bad.prepare().unwrap_err(), VlmError::InvalidGroundHeight);
}

fn helmbold_lift_slope(aspect_ratio: f64, mach: f64) -> f64 {
    let beta = (1.0 - mach * mach).sqrt();
    let ar_beta = beta * aspect_ratio;
    2.0 * core::f64::consts::PI * aspect_ratio / (2.0 + (4.0 + ar_beta * ar_beta).sqrt())
}

#[test]
fn compressible_lift_slope_follows_helmbold_diederich() {
    // AR = 10 rectangular wing, linear regime (1°). The lattice has a small
    // discretization offset from the closed form at M = 0, so compare the
    // Mach-dependence as a ratio to the incompressible result, which cancels it.
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let alpha = 1.0_f64.to_radians();
    let cl_at = |mach: f64| {
        VlmModel::from_planform_compressible(&wing, 24, 6, CamberModel::None, mach)
            .unwrap()
            .solve(alpha, V, RHO, &[])
            .unwrap()
            .cl
    };
    let cl0 = cl_at(0.0);
    for &mach in &[0.5, 0.7] {
        let ratio = cl_at(mach) / cl0;
        let expected = helmbold_lift_slope(10.0, mach) / helmbold_lift_slope(10.0, 0.0);
        assert!(
            (ratio / expected - 1.0).abs() < 0.02,
            "M {mach}: CL ratio {ratio} vs Helmbold {expected}"
        );
    }
}

#[test]
fn compressibility_raises_the_lift_slope_and_vanishes_at_zero_mach() {
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let alpha = 2.0_f64.to_radians();
    let at = |mach: f64| {
        VlmModel::from_planform_compressible(&wing, 16, 4, CamberModel::None, mach)
            .unwrap()
            .solve(alpha, V, RHO, &[])
            .unwrap()
            .cl
    };
    let incompressible = VlmModel::from_planform(&wing, 16, 4, CamberModel::None)
        .unwrap()
        .solve(alpha, V, RHO, &[])
        .unwrap()
        .cl;
    assert!((at(0.0) - incompressible).abs() < 1e-12);
    assert!(at(0.3) > incompressible && at(0.6) > at(0.3) && at(0.8) > at(0.6));
    let model = VlmModel::from_planform_compressible(&wing, 8, 4, CamberModel::None, 0.6).unwrap();
    assert!((model.compressibility_beta() - 0.8).abs() < 1e-12);
}

#[test]
fn compressible_model_rejects_supersonic_and_invalid_mach() {
    let wing = WingPlanform::rectangular(10.0, 1.0);
    for &mach in &[1.0, 1.5, -0.1, f64::NAN] {
        assert_eq!(
            VlmModel::from_planform_compressible(&wing, 8, 4, CamberModel::None, mach).unwrap_err(),
            VlmError::InvalidFlightCondition,
            "mach {mach}"
        );
    }
}

#[test]
fn compressible_reference_is_physical_and_not_double_counted() {
    use tpt_aero_vlm::ReferenceGeometry;
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let mach = 0.6;
    let default =
        VlmModel::from_planform_compressible(&wing, 16, 4, CamberModel::None, mach).expect("valid");
    // Supplying the true (physical) reference must give the same coefficients.
    let physical = ReferenceGeometry {
        area: wing.area(),
        chord: wing.mac(),
        span: wing.span(),
        moment_center: tpt_aero_core::Vec3::zero(),
    };
    let overridden = default.clone().with_reference(physical);
    let a = 3.0_f64.to_radians();
    let s1 = default.solve(a, 50.0, 1.225, &[]).expect("solves");
    let s2 = overridden.solve(a, 50.0, 1.225, &[]).expect("solves");
    assert!((s1.cl - s2.cl).abs() < 1e-9, "{} vs {}", s1.cl, s2.cl);
}

#[test]
fn surfaces_compressible_matches_planform_compressible_for_a_wing() {
    use tpt_aero_vlm::Surface;
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let single =
        VlmModel::from_planform_compressible(&wing, 16, 4, CamberModel::None, 0.5).expect("valid");
    let multi =
        VlmModel::from_surfaces_compressible(&[Surface::horizontal(wing.clone(), 16, 4)], 0.5)
            .expect("valid");
    let a = 2.0_f64.to_radians();
    let s1 = single.solve(a, 50.0, 1.225, &[]).expect("solves");
    let s2 = multi.solve(a, 50.0, 1.225, &[]).expect("solves");
    assert!((s1.cl - s2.cl).abs() < 1e-9);
    assert!(VlmModel::from_surfaces_compressible(&[Surface::horizontal(wing, 4, 2)], 1.0).is_err());
}

#[test]
fn validity_flags_mark_each_limit_independently() {
    use tpt_aero_vlm::{validity_flags, VALIDITY_ALPHA_LIMIT};
    let ok = validity_flags(0.05, 0.02, 0.2, 1.0e6);
    assert!(ok.all_ok());
    // Each limit fails on its own; the others stay valid.
    let high_alpha = validity_flags(VALIDITY_ALPHA_LIMIT + 0.01, 0.0, 0.2, 1.0e6);
    assert!(
        !high_alpha.alpha_ok && high_alpha.beta_ok && high_alpha.mach_ok && high_alpha.reynolds_ok
    );
    assert!(!validity_flags(0.0, 0.3, 0.2, 1.0e6).beta_ok);
    assert!(!validity_flags(0.0, 0.0, 0.8, 1.0e6).mach_ok);
    assert!(!validity_flags(0.0, 0.0, 0.2, 5.0e4).reynolds_ok);
    // NaN never reads as valid.
    assert!(!validity_flags(0.0, 0.0, f64::NAN, 1.0e6).mach_ok);
}
