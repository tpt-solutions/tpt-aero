//! Numerical acceptance tests for the VLM against analytic references.

use tpt_aero_core::WingPlanform;
use tpt_aero_vlm::{CamberModel, VlmModel};

use tpt_math_autodiff_fwd::Dual;

const RHO: f64 = 1.225;
const V_INF: f64 = 50.0;
const DEG: f64 = core::f64::consts::PI / 180.0;

/// Exact Prandtl lifting-line solution (Fourier collocation) for a
/// rectangular wing — the ground truth for the VLM acceptance band.
fn fourier_ll_rectangular(span: f64, chord: f64, alpha: f64, n: usize) -> f64 {
    let pi = core::f64::consts::PI;
    let mut a = vec![vec![0.0; n]; n];
    let rhs = vec![alpha; n];
    for (i, row) in a.iter_mut().enumerate() {
        let theta = pi * (i + 1) as f64 / (n + 1) as f64;
        let st = theta.sin();
        for (nn, aij) in row.iter_mut().enumerate() {
            let k = nn as i32 + 1;
            *aij = k as f64 * (k as f64 * theta).sin() / st
                + 2.0 * span * (k as f64 * theta).sin() / (pi * chord);
        }
    }
    let sol = tpt_aero_vlm::lu_solve(&a, &rhs).expect("LL solve");
    let ar = span * span / (span * chord);
    pi * ar * sol[0]
}

#[test]
fn rectangular_wing_matches_lifting_line() {
    // AR = 10 rectangular wing: the exact Fourier lifting-line solution gives
    // C_L ~= 0.4404 at alpha = 5 deg; the classic horseshoe VLM converges
    // ~2-3% below it (in-plane wake induction).
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let model = VlmModel::from_planform(&wing, 24, 4, CamberModel::None).expect("valid");
    let alpha = 5.0 * DEG;
    let sol = model.solve(alpha, V_INF, RHO, &[]).expect("solvable");

    let cl_expected = fourier_ll_rectangular(10.0, 1.0, alpha, 48);
    assert!(
        (sol.cl - cl_expected).abs() < 0.02,
        "cl = {} vs lifting-line {cl_expected}",
        sol.cl
    );
    assert!(sol.cl > 0.0, "positive alpha must give positive lift");

    // Induced drag: e ~= 1 for a rectangular planform.
    let cdi_expected = sol.cl * sol.cl / (core::f64::consts::PI * 10.0);
    assert!(
        (sol.cdi - cdi_expected).abs() < 0.10 * cdi_expected,
        "cdi = {} vs {cdi_expected}",
        sol.cdi
    );

    // Spanwise distribution: symmetric and tip-unloaded.
    assert_eq!(sol.strip_cl.len(), 24);
    let mid = sol.strip_cl[12];
    let tip = sol.strip_cl[23];
    assert!(
        tip < mid,
        "tip loading {tip} must be below mid {mid} (tip loss)"
    );
    assert!(tip > 0.0, "tip must keep positive loading, got {tip}");
    for i in 0..12 {
        assert!(
            (sol.strip_cl[i] - sol.strip_cl[23 - i]).abs() < 1e-9,
            "spanwise loading must be symmetric"
        );
    }
}

#[test]
fn symmetric_wing_zero_alpha_zero_roll_yaw() {
    // Symmetric geometry at zero AoA: C_L ≈ 0, C_l ≈ 0, C_n ≈ 0.
    let wing = WingPlanform::tapered(8.0, 1.5, 0.8, 0.3, 0.12);
    let model = VlmModel::from_planform(&wing, 30, 5, CamberModel::None).expect("valid");
    let sol = model.solve(0.0, V_INF, RHO, &[]).expect("solvable");
    assert!(sol.cl.abs() < 1e-10, "cl = {}", sol.cl);
    assert!(sol.cl_roll.abs() < 1e-10, "cl_roll = {}", sol.cl_roll);
    assert!(sol.cn.abs() < 1e-10, "cn = {}", sol.cn);
    assert!(sol.cdi.abs() < 1e-12, "cdi = {}", sol.cdi);
}

#[test]
fn positive_camber_gives_positive_lift_at_zero_alpha() {
    // NACA 24-series-like camber (2% at 40%): thin-airfoil C_L(0) ≈ +0.22.
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let model = VlmModel::from_planform(&wing, 24, 4, CamberModel::Naca4 { m: 0.02, p: 0.4 })
        .expect("valid");
    let sol = model.solve(0.0, V_INF, RHO, &[]).expect("solvable");
    assert!(
        sol.cl > 0.1 && sol.cl < 0.35,
        "cambered cl(0) = {} (expected ≈ 0.2)",
        sol.cl
    );
}

#[test]
fn washout_unloads_tips() {
    // Negative (washout) twist toward the tips reduces tip loading relative
    // to the untwisted wing at the same alpha.
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let refined = wing.refine(25);
    let n = refined.sections().len();
    let secs: Vec<tpt_aero_core::WingSection> = refined
        .sections()
        .iter()
        .enumerate()
        .map(|(i, s)| tpt_aero_core::WingSection {
            le: s.le,
            chord: s.chord,
            twist: -4.0 * DEG * (i as f64 / (n - 1) as f64),
        })
        .collect();
    let twisted = WingPlanform::from_sections(secs).expect("valid");

    let model_plain = VlmModel::from_planform(&wing, 24, 4, CamberModel::None).expect("valid");
    let model_twist = VlmModel::from_planform(&twisted, 24, 4, CamberModel::None).expect("valid");

    let alpha = 5.0 * DEG;
    let sol_plain = model_plain.solve(alpha, V_INF, RHO, &[]).expect("ok");
    // Twist vector: one entry per element, taken from the lattice sections.
    let twists: Vec<f64> = model_twist
        .lattice()
        .elements()
        .iter()
        .map(|e| e.twist)
        .collect();
    let sol_twist = model_twist.solve(alpha, V_INF, RHO, &twists).expect("ok");

    assert!(
        sol_twist.strip_cl[23] < sol_plain.strip_cl[23],
        "washout must unload the tip: {} vs {}",
        sol_twist.strip_cl[23],
        sol_plain.strip_cl[23]
    );
    assert!(sol_twist.cl > 0.0 && sol_twist.cl < sol_plain.cl * 1.05);
}

#[test]
fn freestream_reversal_flips_lift() {
    // spec.txt §2 property: reversing the freestream reverses the lift vector.
    let wing = WingPlanform::tapered(8.0, 1.5, 0.8, 0.3, 0.12);
    let model = VlmModel::from_planform(&wing, 20, 4, CamberModel::None).expect("valid");
    let alpha = 4.0 * DEG;
    let fwd = model.solve(alpha, V_INF, RHO, &[]).expect("ok");
    let rev = model
        .solve(alpha + core::f64::consts::PI, V_INF, RHO, &[])
        .expect("ok");
    assert!(
        (rev.cl + fwd.cl).abs() < 1e-9 * fwd.cl.abs().max(1.0),
        "cl fwd {} rev {}",
        fwd.cl,
        rev.cl
    );
}

#[test]
fn spanwise_loading_shape() {
    // The inboard half of the spanwise distribution must rise monotonically
    // toward midspan (elliptic-like shape), symmetric about the root.
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let model = VlmModel::from_planform(&wing, 40, 4, CamberModel::None).expect("valid");
    let sol = model.solve(5.0 * DEG, V_INF, RHO, &[]).expect("ok");
    let n = sol.strip_cl.len();
    for i in 0..n / 2 - 1 {
        assert!(
            sol.strip_cl[i] < sol.strip_cl[i + 1],
            "loading must rise inboard: strip {i} = {} vs {}",
            sol.strip_cl[i],
            sol.strip_cl[i + 1]
        );
    }
    // Span efficiency: the near-field drag must match C_L^2/(pi*AR*e) with
    // e close to 1 for the rectangular planform.
    let e = sol.cl * sol.cl / (core::f64::consts::PI * 10.0 * sol.cdi);
    assert!((e - 1.0).abs() < 0.08, "span efficiency e = {e}");
}

#[test]
fn differentiable_alpha_gradient_matches_thin_wing() {
    // dC_L/dα from a dual-number solve must match both the finite difference
    // and the lifting-line estimate 2πAR/(AR+2) ≈ 5.236 /rad.
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let model = VlmModel::from_planform(&wing, 16, 4, CamberModel::None).expect("valid");

    let alpha0 = 3.0 * DEG;
    let a = Dual::<f64>::variable(alpha0, 0);
    let fields = model
        .solve_generic::<Dual<f64, 1>>(a, Dual::constant(V_INF), Dual::constant(RHO), &[])
        .expect("ok");
    let ad_grad = fields.cl.du(0);

    let eps = 1e-6;
    let lo = model.solve(alpha0 - eps, V_INF, RHO, &[]).expect("ok").cl;
    let hi = model.solve(alpha0 + eps, V_INF, RHO, &[]).expect("ok").cl;
    let fd_grad = (hi - lo) / (2.0 * eps);

    assert!(
        (ad_grad - fd_grad).abs() < 1e-4,
        "ad {ad_grad} vs fd {fd_grad}"
    );
    // The VLM lift-curve slope for AR=10 sits slightly below the
    // 2*pi*AR/(AR+2)=5.24 /rad elliptic estimate (exact LL ~ 5.05 /rad here).
    assert!((4.6..=5.1).contains(&ad_grad), "dC_L/dalpha = {ad_grad}");
}

#[test]
fn differentiable_twist_gradient_matches_finite_difference() {
    // Seeding one element's twist as a dual variable must give the same
    // dC_Di/dtwist as a central finite difference — this is the gradient the
    // induced-drag-minimization loop consumes.
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let model = VlmModel::from_planform(&wing, 12, 4, CamberModel::None).expect("valid");
    let alpha = 3.0 * DEG;

    // Pick an inboard chordwise panel of a mid strip.
    let seed_idx = model
        .lattice()
        .elements()
        .iter()
        .position(|e| e.strip == 6)
        .expect("strip 6 exists");
    let base: Vec<f64> = model.lattice().elements().iter().map(|e| e.twist).collect();

    let twist_val = 1.0 * DEG;
    let dual_twists: Vec<Dual<f64, 1>> = base
        .iter()
        .enumerate()
        .map(|(i, &t)| {
            if i == seed_idx {
                Dual::<f64>::variable(t + twist_val, 0)
            } else {
                Dual::constant(t)
            }
        })
        .collect();
    let fields = model
        .solve_generic::<Dual<f64, 1>>(
            Dual::<f64, 1>::constant(alpha),
            Dual::constant(V_INF),
            Dual::constant(RHO),
            &dual_twists,
        )
        .expect("ok");
    let ad_grad = fields.cdi.du(0);

    // Finite difference w.r.t. that single element's twist.
    let mut lo_t = base.clone();
    let mut hi_t = base.clone();
    lo_t[seed_idx] += twist_val - 1e-6;
    hi_t[seed_idx] += twist_val + 1e-6;
    let lo = model.solve(alpha, V_INF, RHO, &lo_t).expect("ok").cdi;
    let hi = model.solve(alpha, V_INF, RHO, &hi_t).expect("ok").cdi;
    let fd_grad = (hi - lo) / 2e-6;

    assert!(
        (ad_grad - fd_grad).abs() < 1e-6,
        "ad {ad_grad} vs fd {fd_grad}"
    );
}

#[test]
fn flat_wing_induced_drag_scales_with_cl_squared() {
    // For a flat rectangular wing the induced drag is proportional to C_L²
    // away from zero lift, so C_Di / C_L² is constant across the sweep
    // (measured spread ≈ 0.2 % over −4°..+8°). A regression here would mean
    // the induced-drag projection picked up a non-quadratic residual.
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let model = VlmModel::from_planform(&wing, 32, 4, CamberModel::None).expect("valid");
    let ratios: Vec<f64> = [-4.0, -2.0, -1.0, 2.0, 4.0, 8.0]
        .iter()
        .map(|deg: &f64| {
            let sol = model.solve(deg * DEG, V_INF, RHO, &[]).expect("solvable");
            sol.cdi / (sol.cl * sol.cl)
        })
        .collect();
    let lo = ratios.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = ratios.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert!(
        (hi - lo) / lo < 0.01,
        "C_Di/C_L² spread {:.4}% over sweep: {ratios:?}",
        100.0 * (hi - lo) / lo
    );
}

#[test]
fn rejects_invalid_flight_conditions() {
    let wing = WingPlanform::rectangular(10.0, 1.0);
    let model = VlmModel::from_planform(&wing, 8, 2, CamberModel::None).expect("valid");
    for (alpha, v, rho) in [
        (0.05, 0.0, RHO),
        (0.05, -1.0, RHO),
        (0.05, V_INF, 0.0),
        (f64::NAN, V_INF, RHO),
        (0.05, f64::NAN, RHO),
        (0.05, V_INF, f64::INFINITY),
    ] {
        assert!(
            matches!(
                model.solve(alpha, v, rho, &[]),
                Err(tpt_aero_vlm::VlmError::InvalidFlightCondition)
            ),
            "({alpha}, {v}, {rho}) must be rejected"
        );
    }
    assert!(model.solve(0.05, V_INF, RHO, &[]).is_ok());
}

#[test]
fn prepared_sweep_matches_one_shot_solve() {
    // Tapered, swept, dihedral, cambered: nothing symmetric to hide a bug.
    let wing = WingPlanform::tapered(9.0, 1.6, 0.7, 0.6, 4.0 * DEG);
    let model =
        VlmModel::from_planform(&wing, 10, 4, CamberModel::Naca4 { m: 0.03, p: 0.45 }).unwrap();
    let alphas: Vec<f64> = [-3.0, 0.0, 2.5, 7.0].iter().map(|a| a * DEG).collect();
    let sweep = model.solve_sweep(&alphas, V_INF, RHO).unwrap();
    assert_eq!(sweep.len(), alphas.len());
    for (&alpha, fast) in alphas.iter().zip(&sweep) {
        let slow = model.solve(alpha, V_INF, RHO, &[]).unwrap();
        let close = |a: f64, b: f64, what: &str| {
            assert!(
                (a - b).abs() <= 1e-9 * (1.0 + b.abs()),
                "{what} at {alpha}: {a} vs {b}"
            );
        };
        close(fast.cl, slow.cl, "cl");
        close(fast.cdi, slow.cdi, "cdi");
        close(fast.cm, slow.cm, "cm");
        close(fast.cl_roll, slow.cl_roll, "cl_roll");
        close(fast.cn, slow.cn, "cn");
        close(fast.force.x, slow.force.x, "fx");
        close(fast.force.z, slow.force.z, "fz");
        assert_eq!(fast.strip_cl.len(), slow.strip_cl.len());
        for (a, b) in fast.strip_cl.iter().zip(&slow.strip_cl) {
            close(*a, *b, "strip_cl");
        }
        for (a, b) in fast.gamma.iter().zip(&slow.gamma) {
            close(*a, *b, "gamma");
        }
    }
}

#[test]
fn prepared_solve_validates_flight_condition() {
    let wing = WingPlanform::rectangular(8.0, 1.0);
    let model = VlmModel::from_planform(&wing, 6, 2, CamberModel::None).unwrap();
    let prepared = model.prepare().unwrap();
    assert!(matches!(
        prepared.solve(0.1, 0.0, RHO),
        Err(tpt_aero_vlm::VlmError::InvalidFlightCondition)
    ));
    assert!(matches!(
        model.solve_sweep(&[0.1, f64::NAN], V_INF, RHO),
        Err(tpt_aero_vlm::VlmError::InvalidFlightCondition)
    ));
}
