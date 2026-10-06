//! Finite-difference and dual-number extraction of aerodynamic derivatives
//! from a VLM plant model.
//!
//! Longitudinal derivatives extracted here (per spec.txt §4):
//!
//! * `C_{Lα}`, `C_{mα}` — angle-of-attack sensitivity,
//! * `C_{Lq}`, `C_{mq}` — pitch-rate sensitivity, modeled by the classic
//!   quasi-steady approach: a pitch rate q about the reference point changes
//!   each lattice element's local incidence by `−q·Δx/V∞` where Δx is the
//!   element's x-offset from the reference point,
//! * `C_{Lδe}`, `C_{mδe}` — elevator (flap) effectiveness via a uniform
//!   camber-line incidence change δ (a camber flap spanning the full chord
//!   has δ_e ≡ effective incidence).

use alloc::{vec, vec::Vec};

#[allow(unused_imports)]
use tpt_aero_core::Real as _;
use tpt_aero_vlm::VlmModel;

/// A set of extracted stability derivatives (1/rad units unless noted).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LongitudinalDerivatives {
    /// Lift-curve slope `C_{Lα}`.
    pub cl_alpha: f64,
    /// Pitch-moment slope `C_{mα}` (negative = statically stable).
    pub cm_alpha: f64,
    /// Pitch-rate lift derivative `C_{Lq}` (per rad/s, normalized: the
    /// classic `C_{Lq}` uses q·c/2V; here the raw dC_L/dq̂ with q̂ = q·c/2V∞).
    pub cl_q: f64,
    /// Pitch-rate moment derivative `C_{mq}`.
    pub cm_q: f64,
    /// Flap-effectiveness lift derivative `C_{Lδ}`.
    pub cl_delta: f64,
    /// Flap-effectiveness moment derivative `C_{mδ}`.
    pub cm_delta: f64,
}

/// A small central-difference helper.
fn central_difference(f: impl Fn(f64) -> f64, x0: f64, eps: f64) -> f64 {
    (f(x0 + eps) - f(x0 - eps)) / (2.0 * eps)
}

/// Extract longitudinal derivatives from a VLM model by central finite
/// differences at the given trim point.
#[must_use]
pub fn longitudinal_fd(
    model: &VlmModel,
    alpha_rad: f64,
    v_inf: f64,
    density: f64,
    eps_alpha: f64,
) -> LongitudinalDerivatives {
    let cl_at = |a: f64| model.solve(a, v_inf, density, &[]).map(|s| s.cl).unwrap_or(f64::NAN);
    let cm_at = |a: f64| model.solve(a, v_inf, density, &[]).map(|s| s.cm).unwrap_or(f64::NAN);
    let cl_alpha = central_difference(cl_at, alpha_rad, eps_alpha);
    let cm_alpha = central_difference(cm_at, alpha_rad, eps_alpha);

    // Pitch rate: local incidence change −q̂·(x_i − x_ref)·2/c_ref... with
    // the normalized rate q̂ = q·c/2V, the incidence change per unit q̂ is
    // Δα_i = −2·q̂·(x_i − x_ref)/c_ref. Enforce via the per-element twist
    // field (twist adds incidence linearly in this VLM).
    let ref_x = model.reference().moment_center.x;
    let c_ref = model.reference().chord;
    let pitch_twists = |qhat: f64| -> Vec<f64> {
        model
            .lattice()
            .elements()
            .iter()
            .map(|e| 2.0 * qhat * (e.collocation.x - ref_x) / c_ref)
            .collect()
    };
    let cl_q_at = |qhat: f64| {
        let t = pitch_twists(qhat);
        model
            .solve(alpha_rad, v_inf, density, &t)
            .map(|s| s.cl)
            .unwrap_or(f64::NAN)
    };
    let cm_q_at = |qhat: f64| {
        let t = pitch_twists(qhat);
        model
            .solve(alpha_rad, v_inf, density, &t)
            .map(|s| s.cm)
            .unwrap_or(f64::NAN)
    };
    let cl_q = central_difference(cl_q_at, 0.0, 1.0e-3);
    let cm_q = central_difference(cm_q_at, 0.0, 1.0e-3);

    // Flap effectiveness: uniform incidence increment (full-chord camber
    // flap: δ_e shifts every element's zero-lift line by δ_e).
    let cl_d_at = |d: f64| {
        let t = vec![d; model.lattice().len()];
        model
            .solve(alpha_rad, v_inf, density, &t)
            .map(|s| s.cl)
            .unwrap_or(f64::NAN)
    };
    let cm_d_at = |d: f64| {
        let t = vec![d; model.lattice().len()];
        model
            .solve(alpha_rad, v_inf, density, &t)
            .map(|s| s.cm)
            .unwrap_or(f64::NAN)
    };
    let cl_delta = central_difference(cl_d_at, 0.0, 1.0e-3);
    let cm_delta = central_difference(cm_d_at, 0.0, 1.0e-3);

    LongitudinalDerivatives {
        cl_alpha,
        cm_alpha,
        cl_q,
        cm_q,
        cl_delta,
        cm_delta,
    }
}

/// Extract `C_{Lα}` exactly (to machine precision) with the forward-mode
/// dual path of the VLM — differentiable assembly, no step-size tuning.
#[must_use]
pub fn cl_alpha_dual(model: &VlmModel, alpha_rad: f64, v_inf: f64, density: f64) -> f64 {
    use tpt_math_autodiff_fwd::Dual;
    let a = Dual::<f64>::variable(alpha_rad, 0);
    model
        .solve_generic::<Dual<f64, 1>>(
            a,
            Dual::constant(v_inf),
            Dual::constant(density),
            &[],
        )
        .map(|f| f.cl.du(0))
        .unwrap_or(f64::NAN)
}

/// Classic analytical estimates for a rectangular (or high-aspect tapered)
/// wing — the sanity references for the extracted values.
pub mod analytical {
    #[allow(unused_imports)]
    use tpt_aero_core::Real as _;

    /// Lifting-line estimate of the lift-curve slope (1/rad):
    /// `C_{Lα} = a_0·AR/(AR + 2)` with `a_0` the section lift-curve slope.
    #[must_use]
    pub fn cl_alpha_wing(ar: f64, a0: f64) -> f64 {
        a0 * ar / (ar + 2.0)
    }

    /// Thin-airfoil flap effectiveness `dC_l/dδ_e` for a plain flap hinged at
    /// `x_h/c`: `2(π − θ_h + sin θ_h)` with `cos θ_h = 2·x_h/c − 1`
    /// (Glauert). Full-chord "flap" (x_h → 0) gives 2π; hinge at the TE
    /// gives 0.
    #[must_use]
    pub fn cl_delta_flap(x_h_over_c: f64) -> f64 {
        let theta_h = (1.0 - 2.0 * x_h_over_c).acos();
        2.0 * (core::f64::consts::PI - theta_h + theta_h.sin())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_aero_core::WingPlanform;
    use tpt_aero_vlm::CamberModel;

    #[test]
    fn finite_difference_matches_dual_and_thin_wing() {
        let wing = WingPlanform::rectangular(10.0, 1.0);
        let model = VlmModel::from_planform(&wing, 20, 4, CamberModel::None).expect("valid");
        let alpha = 2.0_f64.to_radians();
        let fd = longitudinal_fd(&model, alpha, 50.0, 1.225, 1.0e-4);
        let dual = cl_alpha_dual(&model, alpha, 50.0, 1.225);
        assert!(
            (fd.cl_alpha - dual).abs() < 1e-3,
            "fd {} vs dual {}",
            fd.cl_alpha,
            dual
        );
        // Compare against the exact Fourier lifting-line value ≈ 5.05 /rad.
        assert!(
            (4.7..=5.2).contains(&fd.cl_alpha),
            "C_Lα = {}",
            fd.cl_alpha
        );
        // A camberless wing referenced at its aerodynamic center (MAC
        // quarter chord) has C_mα ≈ 0; static stability arises once the
        // reference point is ahead of the AC.
        assert!(
            fd.cm_alpha.abs() < 0.1,
            "C_mα = {} (should vanish at the AC)",
            fd.cm_alpha
        );
    }

    #[test]
    fn pitch_rate_damping() {
        // C_mq must be negative (pitch damping) and C_Lq positive.
        let wing = WingPlanform::rectangular(10.0, 1.0);
        let model = VlmModel::from_planform(&wing, 16, 4, CamberModel::None).expect("valid");
        let d = longitudinal_fd(&model, 2.0_f64.to_radians(), 50.0, 1.225, 1e-4);
        assert!(d.cm_q < 0.0, "C_mq = {}", d.cm_q);
        assert!(d.cl_q > 0.0, "C_Lq = {}", d.cl_q);
    }

    #[test]
    fn flap_effectiveness_sign_and_scale() {
        let wing = WingPlanform::rectangular(10.0, 1.0);
        let model = VlmModel::from_planform(&wing, 16, 4, CamberModel::None).expect("valid");
        let d = longitudinal_fd(&model, 2.0_f64.to_radians(), 50.0, 1.225, 1e-4);
        // A full-chord camber flap shifts every element's zero-lift line:
        // its effectiveness equals the wing lift-curve slope.
        assert!(
            (d.cl_delta - d.cl_alpha).abs() < 0.1,
            "C_Lδ = {} vs C_Lα = {}",
            d.cl_delta,
            d.cl_alpha
        );
        // Analytical reference: monotone decreasing in hinge position,
        // 2π at the LE and 0 at the TE.
        assert!(analytical::cl_delta_flap(0.75) > analytical::cl_delta_flap(0.9));
        assert!((analytical::cl_delta_flap(0.0) - core::f64::consts::TAU).abs() < 1e-9);
    }
}
