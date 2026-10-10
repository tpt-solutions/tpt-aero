//! VLM against an independent Prandtl lifting-line solution for straight
//! tapered wings (Fourier collocation, any chord distribution).

use tpt_aero_core::WingPlanform;
use tpt_aero_vlm::{CamberModel, VlmModel};

const V_INF: f64 = 50.0;
const RHO: f64 = 1.225;
const DEG: f64 = core::f64::consts::PI / 180.0;
const PI: f64 = core::f64::consts::PI;

/// Prandtl lifting-line solution for a flat, untwisted, unswept wing.
struct LiftingLine {
    span: f64,
    area: f64,
    a: Vec<f64>,
    chord: Box<dyn Fn(f64) -> f64>,
}

impl LiftingLine {
    /// Solve for Fourier coefficients `A_k` at angle of attack `alpha` using
    /// `n` collocation stations; `chord(y)` is the chord at spanwise `y`.
    fn solve(span: f64, chord: impl Fn(f64) -> f64 + 'static, alpha: f64, n: usize) -> Self {
        // Area by midpoint quadrature of the chord distribution.
        let m = 4000;
        let area: f64 = (0..m)
            .map(|i| chord(-0.5 * span + span * (i as f64 + 0.5) / m as f64) * span / m as f64)
            .sum();
        let mut a = vec![vec![0.0; n]; n];
        for (i, row) in a.iter_mut().enumerate() {
            let theta = PI * (i + 1) as f64 / (n + 1) as f64;
            let y = 0.5 * span * theta.cos();
            let st = theta.sin();
            for (k0, aij) in row.iter_mut().enumerate() {
                let k = (k0 + 1) as f64;
                // alpha = (4b / (a0 c)) sum A_k sin(k theta) + sum k A_k sin(k theta)/sin(theta)
                // with the thin-airfoil section slope a0 = 2 pi.
                *aij =
                    k * (k * theta).sin() / st + 2.0 * span * (k * theta).sin() / (PI * chord(y));
            }
        }
        let rhs = vec![alpha; n];
        let coeffs = tpt_aero_vlm::lu_solve(&a, &rhs).expect("lifting-line system");
        Self {
            span,
            area,
            a: coeffs,
            chord: Box::new(chord),
        }
    }

    fn aspect_ratio(&self) -> f64 {
        self.span * self.span / self.area
    }

    fn cl(&self) -> f64 {
        PI * self.aspect_ratio() * self.a[0]
    }

    fn cdi(&self) -> f64 {
        PI * self.aspect_ratio()
            * self
                .a
                .iter()
                .enumerate()
                .map(|(k0, a)| (k0 + 1) as f64 * a * a)
                .sum::<f64>()
    }

    /// Section lift coefficient (based on the local chord) at spanwise `y`.
    fn section_cl(&self, y: f64) -> f64 {
        let theta = (2.0 * y / self.span).acos();
        let gamma_over_2bv: f64 = self
            .a
            .iter()
            .enumerate()
            .map(|(k0, a)| a * ((k0 + 1) as f64 * theta).sin())
            .sum();
        4.0 * self.span * gamma_over_2bv / (self.chord)(y)
    }
}

fn tapered_chord(span: f64, root: f64, tip: f64) -> impl Fn(f64) -> f64 {
    move |y: f64| root - (root - tip) * y.abs() / (0.5 * span)
}

fn check_tapered(span: f64, root: f64, tip: f64, alpha_deg: f64, n_span: usize) {
    let alpha = alpha_deg * DEG;
    let ll = LiftingLine::solve(span, tapered_chord(span, root, tip), alpha, 80);
    let wing = WingPlanform::tapered(span, root, tip, 0.0, 0.0);
    assert!(
        (wing.area() - ll.area).abs() < 1e-6 * ll.area,
        "planform area {} vs {}",
        wing.area(),
        ll.area
    );
    let model = VlmModel::from_planform(&wing, n_span, 4, CamberModel::None).expect("valid");
    let sol = model.solve(alpha, V_INF, RHO, &[]).expect("solvable");

    // Lift: horseshoe VLM sits within a few percent of lifting line at AR 8-12.
    let cl_err = (sol.cl / ll.cl() - 1.0).abs();
    assert!(
        cl_err < 0.04,
        "AR {:.1} taper {:.2}: CL {} vs LL {} ({:.1}%)",
        ll.aspect_ratio(),
        tip / root,
        sol.cl,
        ll.cl(),
        100.0 * cl_err
    );
    // Induced drag, compared through the Oswald efficiency e = CL^2/(pi AR CDi)
    // so the (mesh-independent) few-percent CL offset does not leak into it.
    // The mesh-converged VLM e is within about 3% of lifting line.
    let ar = ll.aspect_ratio();
    let e_vlm = sol.cl * sol.cl / (PI * ar * sol.cdi);
    let e_ll = ll.cl() * ll.cl() / (PI * ar * ll.cdi());
    assert!(
        (e_vlm / e_ll - 1.0).abs() < 0.05,
        "Oswald e: VLM {e_vlm:.4} vs LL {e_ll:.4}"
    );

    // Spanwise loading: section cl at strip centres, away from the tips.
    let strips = sol.strip_cl.len();
    for (i, &vlm_cl) in sol.strip_cl.iter().enumerate() {
        let y = -0.5 * span + span * (i as f64 + 0.5) / strips as f64;
        if y.abs() > 0.9 * 0.5 * span {
            continue;
        }
        let want = ll.section_cl(y);
        assert!(
            (vlm_cl / want - 1.0).abs() < 0.07,
            "strip {i} (y = {y:.2}): c_l {vlm_cl} vs LL {want}"
        );
    }
}

#[test]
fn tapered_ar8_matches_lifting_line() {
    // Span 8, root 1.333, tip 0.667: taper 0.5, AR 8.
    check_tapered(8.0, 4.0 / 3.0, 2.0 / 3.0, 5.0, 32);
}

#[test]
fn strongly_tapered_ar12_matches_lifting_line() {
    // Taper 0.3, AR 12.
    check_tapered(12.0, 1.0 / 0.65 * 1.0, 0.3 / 0.65, 4.0, 40);
}

#[test]
fn lifting_line_reference_reproduces_the_elliptic_limit() {
    // Sanity-check the reference itself: an elliptic chord gives e = 1,
    // CDi = CL^2 / (pi AR), and CL = a AR/(AR+2) for a0 = 2 pi.
    let span = 10.0;
    let c0 = 1.2;
    let ll = LiftingLine::solve(
        span,
        move |y: f64| c0 * (1.0 - (2.0 * y / span).powi(2)).max(1e-12).sqrt(),
        5.0 * DEG,
        120,
    );
    let ar = ll.aspect_ratio();
    let e = ll.cl() * ll.cl() / (PI * ar * ll.cdi());
    assert!((e - 1.0).abs() < 0.01, "elliptic Oswald e = {e}");
    let cl_theory = 2.0 * PI * ar / (ar + 2.0) * 5.0 * DEG;
    assert!(
        (ll.cl() / cl_theory - 1.0).abs() < 0.01,
        "{} vs {cl_theory}",
        ll.cl()
    );
}
