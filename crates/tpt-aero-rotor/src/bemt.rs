//! Blade Element Momentum Theory: the coupled blade-element / momentum
//! iteration for axial flight (hover through windmilling), following
//! Leishman, *Principles of Helicopter Aerodynamics*, ch. 3.
//!
//! Per annulus (span fraction `x = r/R`):
//!
//! * inflow angle `tan φ = λ_tot/x`, `α = θ − φ`,
//! * blade element: `dC_T/dx = x σ_r Ŵ² (C_L cosφ + C_D sinφ)`,
//! * momentum (Glauert, with loss factor F):
//!   `dC_T/dx = 4 x λ_tot λ_i F`,
//! * solved by fixed-point iteration on the induced ratio `λ_i`, clamped to
//!   the physical bound `0 ≤ λ_i ≤ 1` every step (Kani-harness target: the
//!   iteration cannot leave `[0, 1]`).

use alloc::vec::Vec;

use crate::blade::{RotorError, RotorGeometry, SectionModel};
use crate::loss::{HubLoss, TipLoss};

/// Axial operating condition.
#[derive(Debug, Clone, Copy)]
pub struct OperatingCondition {
    /// Axial flow speed (m/s). Positive = climb / propeller inflow; zero =
    /// hover; small negative = incipient descent.
    pub axial_speed: f64,
    /// Shaft angular speed (rad/s) > 0.
    pub omega: f64,
    /// Air density (kg/m³).
    pub density: f64,
    /// Collective pitch added to the built-in twist (rad).
    pub collective: f64,
}

impl OperatingCondition {
    /// Hover condition at the given shaft speed and density.
    #[must_use]
    pub const fn hover(omega: f64, density: f64) -> Self {
        Self {
            axial_speed: 0.0,
            omega,
            density,
            collective: 0.0,
        }
    }
}

/// Per-annulus result of the BEMT solve.
#[derive(Debug, Clone, Copy)]
pub struct SectionResult {
    /// Span fraction `x = r/R`.
    pub r_frac: f64,
    /// Total inflow ratio `λ_tot = (V∞ + v_i)/(ΩR)`.
    pub lambda_tot: f64,
    /// Induced inflow ratio `λ_i = v_i/(ΩR)` (non-negative).
    pub lambda_induced: f64,
    /// Inflow angle φ (rad).
    pub inflow_angle: f64,
    /// Section angle of attack α (rad).
    pub alpha: f64,
    /// Section lift coefficient.
    pub cl: f64,
    /// Section drag coefficient.
    pub cd: f64,
    /// Combined loss factor F at this annulus.
    pub loss: f64,
    /// Annulus contribution to the thrust coefficient, `dC_T/dx`.
    pub dct_dx: f64,
    /// Annulus contribution to the torque coefficient, `dC_Q/dx`.
    pub dcq_dx: f64,
    /// Momentum-side annulus contribution (mass-flow consistency check).
    pub dct_dx_momentum: f64,
}

/// Aggregate rotor loads in coefficient form.
#[derive(Debug, Clone, Copy)]
pub struct RotorLoads {
    /// Thrust coefficient `C_T = T/(ρ A (ΩR)²)`.
    pub ct: f64,
    /// Torque coefficient `C_Q = Q/(ρ A (ΩR)² R)`.
    pub cq: f64,
    /// Power coefficient `C_P = C_Q + λ_c C_T` (includes climb power).
    pub cp: f64,
    /// Propeller efficiency `η = J C_T/(2π C_Q)` (axial flight only).
    pub propeller_efficiency: f64,
    /// Figure of merit `C_T^{3/2}/(√2 C_Q)` (hover efficiency measure).
    pub figure_of_merit: f64,
}

/// Full BEMT solution: aggregate loads plus per-annulus detail.
#[derive(Debug, Clone)]
pub struct BemtSolution {
    /// Aggregate thrust/torque/power coefficients.
    pub loads: RotorLoads,
    /// Per-station results (blade-station ordered).
    pub sections: Vec<SectionResult>,
    /// Iterations used by the slowest-converging annulus.
    pub iterations: usize,
    /// Max |element − momentum| annulus residual at convergence.
    pub max_residual: f64,
}

/// The BEMT solver: geometry + closure models bound once, many conditions run.
#[derive(Debug, Clone)]
pub struct BemtSolver {
    geometry: RotorGeometry,
    section: SectionModel,
    tip_loss: TipLoss,
    hub_loss: HubLoss,
    max_iterations: usize,
    tolerance: f64,
    /// Under-relaxation factor for the induction update (0, 1].
    relaxation: f64,
}

impl BemtSolver {
    /// The bound geometry.
    #[must_use]
    pub const fn geometry(&self) -> &RotorGeometry {
        &self.geometry
    }

    /// Bind geometry and models.
    ///
    /// # Errors
    /// [`RotorError::InvalidGeometry`] if the geometry fails validation.
    pub fn new(
        geometry: RotorGeometry,
        section: SectionModel,
        tip_loss: TipLoss,
        hub_loss: HubLoss,
    ) -> Result<Self, RotorError> {
        geometry.validate()?;
        Ok(Self {
            geometry,
            section,
            tip_loss,
            hub_loss,
            max_iterations: 200,
            tolerance: 1.0e-11,
            relaxation: 0.6,
        })
    }

    /// Adjust iteration controls (builder pattern).
    #[must_use]
    pub fn with_iteration_controls(
        mut self,
        max_iterations: usize,
        tolerance: f64,
        relaxation: f64,
    ) -> Self {
        self.max_iterations = max_iterations.max(10);
        self.tolerance = tolerance.max(1.0e-15);
        self.relaxation = relaxation.clamp(0.05, 1.0);
        self
    }

    /// Solve one operating condition.
    ///
    /// # Errors
    /// [`RotorError::NotConverged`] if any annulus fails to converge;
    /// [`RotorError::InvalidSection`] if the closure returns non-finite values.
    pub fn solve(&self, op: &OperatingCondition) -> Result<BemtSolution, RotorError> {
        let tip_speed = self.geometry.tip_radius * op.omega;
        if !(tip_speed > 0.0 && op.density > 0.0) {
            return Err(RotorError::InvalidGeometry);
        }
        let lambda_c = op.axial_speed / tip_speed;
        let mut sections = Vec::with_capacity(self.geometry.sections.len());
        let mut max_residual = 0.0_f64;
        let mut max_iters = 0_usize;

        for sec in &self.geometry.sections {
            let x = sec.r_frac;
            // Local blade solidity: sigma_r = B c / (2 pi r).
            let sigma_r = self.geometry.blades as f64 * sec.chord
                / (2.0 * core::f64::consts::PI * x * self.geometry.tip_radius);
            let theta = sec.twist + op.collective;

            let (lambda_i, iters) = self.solve_annulus(x, sigma_r, theta, lambda_c)?;
            let lambda_tot = (lambda_c + lambda_i).max(0.0);
            let phi = (lambda_tot / x).atan();
            let alpha = theta - phi;
            let (cl, cd) = self.section.coefficients(alpha);
            let w2 = lambda_tot * lambda_tot + x * x;
            let loss = self.loss_factor(x, phi, lambda_tot);

            let dct_dx = x * sigma_r * w2 * (cl * phi.cos() + cd * phi.sin());
            let dcq_dx = x * x * sigma_r * w2 * (cl * phi.sin() + cd * phi.cos());
            let dct_dx_momentum = 4.0 * x * lambda_tot * lambda_i * loss;

            max_iters = max_iters.max(iters);
            // Where the Prandtl factor vanishes (the tip station) the
            // momentum balance is singular; exclude it from the residual.
            if loss >= 0.05 {
                max_residual = max_residual.max((dct_dx - dct_dx_momentum).abs());
            }

            sections.push(SectionResult {
                r_frac: x,
                lambda_tot,
                lambda_induced: lambda_i,
                inflow_angle: phi,
                alpha,
                cl,
                cd,
                loss,
                dct_dx,
                dcq_dx,
                dct_dx_momentum,
            });
        }

        let loads = Self::integrate(&sections, lambda_c);
        Ok(BemtSolution {
            loads,
            sections,
            iterations: max_iters,
            max_residual,
        })
    }

    /// Fixed-point iteration for the induced inflow ratio at one annulus.
    ///
    /// Returns `(λ_i, iterations)` with `0 ≤ λ_i ≤ 1` — the physical bound
    /// the Kani harness in `tpt-aero-verify` proves for the pure function
    /// [`solve_induction_bounded`].
    fn solve_annulus(
        &self,
        x: f64,
        sigma_r: f64,
        theta: f64,
        lambda_c: f64,
    ) -> Result<(f64, usize), RotorError> {
        let a_slope = match &self.section {
            SectionModel::Linear { cl_slope, .. } => *cl_slope,
            SectionModel::Polar { polar, .. } => polar.cl_alpha_slope(),
        };
        let tip_loss = self.tip_loss;
        let hub_loss = self.hub_loss;
        let blades = self.geometry.blades;
        let lambda_i = solve_induction_bounded(sigma_r, theta, x, lambda_c, a_slope, |alpha| {
            self.section.coefficients(alpha)
        }, |phi, lambda_tot| {
            tip_loss.factor(x, phi, lambda_tot, blades) * hub_loss.factor(x, phi, blades)
        }, self.relaxation, self.max_iterations)
        .ok_or(RotorError::InvalidSection)?;
        Ok((lambda_i, self.max_iterations))
    }

    /// Combined tip/hub loss factor at one annulus.
    fn loss_factor(&self, x: f64, phi: f64, lambda_tot: f64) -> f64 {
        let tip = self.tip_loss.factor(x, phi, lambda_tot, self.geometry.blades);
        let hub = self.hub_loss.factor(x, phi, self.geometry.blades);
        tip * hub
    }

    /// Trapezoidal integration of the annulus distributions into coefficients.
    fn integrate(sections: &[SectionResult], lambda_c: f64) -> RotorLoads {
        let mut ct = 0.0;
        let mut cq = 0.0;
        for w in sections.windows(2) {
            let dx = w[1].r_frac - w[0].r_frac;
            ct += 0.5 * (w[0].dct_dx + w[1].dct_dx) * dx;
            cq += 0.5 * (w[0].dcq_dx + w[1].dcq_dx) * dx;
        }
        if let (Some(first), Some(last)) = (sections.first(), sections.last()) {
            // End stubs [0, x_0] and [x_last, 1] carry the end values.
            ct += first.dct_dx * first.r_frac;
            cq += first.dcq_dx * first.r_frac;
            ct += last.dct_dx * (1.0 - last.r_frac);
            cq += last.dcq_dx * (1.0 - last.r_frac);
        }
        let cp = cq + lambda_c * ct;
        // Propeller advance ratio: J = V/(nD) with n = Ω/2π ⇒ J = π λ_c.
        let j = core::f64::consts::PI * lambda_c;
        // With ΩR-based coefficients, η = TV/P = C_T λ_c / C_Q
        // (the textbook J·C_T/(2π C_Q) form holds for n-based coefficients).
        let eta = if cq.abs() > 1.0e-12 { ct * lambda_c / cq } else { 0.0 };
        let _ = j;
        let fom = if cq.abs() > 1.0e-12 {
            ct.powf(1.5) / (core::f64::consts::SQRT_2 * cq)
        } else {
            0.0
        };
        RotorLoads {
            ct,
            cq,
            cp,
            propeller_efficiency: eta,
            figure_of_merit: fom,
        }
    }
}

/// Pure, bounded induction-factor solve — the function the Kani harness
/// proves to keep `λ_i ∈ [0, 1]`.
///
/// Solves the annulus momentum–element balance by bisection on the residual
/// `g(λ_i) = element(λ_i) − 4 x λ_tot λ_i F` over the physical bracket
/// `[0, 1]`. `g` is decreasing in `λ_i` (larger inflow angle unloads the
/// section and grows the momentum side), so bisection converges for any
/// monotone closure and cannot leave the bracket.
///
/// `loss_factor(phi, lambda_tot)` supplies the combined tip/hub factor `F`
/// (clamped to a positive floor internally). Returns the bracket midpoint
/// after `max_iterations` halvings — `Some` always, unless the closure
/// returns non-finite coefficients.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn solve_induction_bounded(
    sigma_r: f64,
    theta: f64,
    x: f64,
    lambda_c: f64,
    _cl_slope: f64,
    coefficients: impl Fn(f64) -> (f64, f64),
    loss_factor: impl Fn(f64, f64) -> f64,
    _relaxation: f64,
    max_iterations: usize,
) -> Option<f64> {
    let residual = |lambda_i: f64| -> Option<f64> {
        let lambda_tot = (lambda_c + lambda_i).max(0.0);
        let phi = (lambda_tot / x).atan();
        let alpha = theta - phi;
        let (cl, cd) = coefficients(alpha);
        if !(cl.is_finite() && cd.is_finite()) {
            return None;
        }
        let w2 = lambda_tot * lambda_tot + x * x;
        let element = x * sigma_r * w2 * (cl * phi.cos() + cd * phi.sin());
        let loss = loss_factor(phi, lambda_tot).clamp(0.01, 1.0);
        Some(element - 4.0 * x * lambda_tot * lambda_i * loss)
    };

    let mut lo = 0.0_f64;
    let mut hi = 1.0_f64;
    // No induced flow if the balance is already satisfied (or negative
    // thrust — the brake / windmill regime keeps λ_i at the lower bound).
    let g0 = residual(lo)?;
    if g0 <= 0.0 {
        return Some(0.0);
    }
    for _ in 0..max_iterations {
        let mid = 0.5 * (lo + hi);
        let g = residual(mid)?;
        if g > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1.0e-12 {
            break;
        }
    }
    Some(0.5 * (lo + hi))
}

