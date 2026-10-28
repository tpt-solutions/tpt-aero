//! Rotor loss models: tip and hub corrections to the momentum theory.

// Brings the `no_std` float methods (`atan`, `powf`, ...) into scope.
#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// Tip-loss model applied to the momentum (induced) term of BEMT.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum TipLoss {
    /// No tip loss (actuator disk).
    #[default]
    None,
    /// Prandtl's finite-number-of-blades tip loss (Leishman eq. 3.171).
    Prandtl,
    /// Glauert's series approximation of the same effect (Leishman eq. 3.170),
    /// using the local inflow ratio in place of sin φ.
    Glauert,
}

impl TipLoss {
    /// Tip-loss factor `F_tip ∈ (0, 1]` at span fraction `x = r/R`.
    ///
    /// `phi` is the local inflow angle (rad), `lambda_tot` the total inflow
    /// ratio, `blades` the blade count.
    #[must_use]
    pub fn factor(&self, x: f64, phi: f64, lambda_tot: f64, blades: usize) -> f64 {
        let b = blades as f64;
        match *self {
            TipLoss::None => 1.0,
            TipLoss::Prandtl => {
                let sin_phi = phi.sin().abs().max(1.0e-6);
                let f = 0.5 * b * (1.0 - x) / (x * sin_phi);
                prandtl_factor(f)
            }
            TipLoss::Glauert => {
                let sin_phi = lambda_tot.abs().max(1.0e-6);
                let f = 0.5 * b * (1.0 - x) / (x * sin_phi);
                prandtl_factor(f)
            }
        }
    }
}

/// Prandtl's closed form `(2/π)·acos(exp(−f))`.
fn prandtl_factor(f: f64) -> f64 {
    core::f64::consts::FRAC_2_PI * (-f.max(0.0)).exp().acos()
}

/// Hub-loss model (same functional form applied to the inboard station).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum HubLoss {
    /// No hub loss.
    #[default]
    None,
    /// Prandtl-type hub loss with root cutout fraction `x_hub = r_hub/R`.
    Prandtl {
        /// Root cutout fraction (0 .. 0.3).
        x_hub: f64,
    },
}

impl HubLoss {
    /// Hub-loss factor `F_hub ∈ (0, 1]` at span fraction `x`.
    #[must_use]
    pub fn factor(&self, x: f64, phi: f64, blades: usize) -> f64 {
        match *self {
            HubLoss::None => 1.0,
            HubLoss::Prandtl { x_hub } => {
                let b = blades as f64;
                let sin_phi = phi.sin().abs().max(1.0e-6);
                let f = 0.5 * b * (x - x_hub) / (x * sin_phi);
                prandtl_factor(f)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prandtl_tip_loss_limits() {
        // Far from the tip: factor → 1.
        let f = TipLoss::Prandtl.factor(0.5, 0.1, 0.05, 2);
        assert!(f > 0.99, "F = {f}");
        // Near the tip the factor drops steeply...
        let f = TipLoss::Prandtl.factor(0.95, 0.1, 0.05, 2);
        assert!((0.3..1.0).contains(&f), "F = {f}");
        // ...and vanishes exactly at the tip station (f → 0 ⇒ acos(1) = 0).
        let f = TipLoss::Prandtl.factor(1.0, 0.1, 0.05, 2);
        assert!(f.abs() < 1e-9, "F = {f}");
        // More blades ⇒ less loss near the tip.
        let b2 = TipLoss::Prandtl.factor(0.95, 0.1, 0.05, 2);
        let b4 = TipLoss::Prandtl.factor(0.95, 0.1, 0.05, 4);
        assert!(b4 > b2);
    }

    #[test]
    fn hub_loss_limits() {
        let hl = HubLoss::Prandtl { x_hub: 0.2 };
        // At the root cutout the factor vanishes (same functional form).
        let f = hl.factor(0.2, 0.1, 2);
        assert!(f.abs() < 1e-9, "F = {f}");
        // Mid-span: close to 1.
        let f = hl.factor(0.6, 0.1, 2);
        assert!(f > 0.99, "F = {f}");
    }
}

/// Momentum-side thrust coefficient `C_T(a)` with the tip/hub loss factor `F`,
/// including the turbulent windmill state (Buhl, NREL/TP-500-36834, eqs. 3, 18).
///
/// Up to `a = 0.4` this is the classical `4Fa(1 − a)`. Above that, Buhl's
/// parabola `C_T = 8/9 + (4F − 40/9)a + (50/9 − 4F)a²` takes over. It touches the
/// classical curve at `a = 0.4` (value `0.96F`, slope `0.8F`) and reaches 2 at
/// `a = 1`, so the momentum balance has no gap for the root search to step
/// across.
///
/// Chapman et al. (2013, *Renewable Energy* 60, 472–480) apply the same
/// correction to tidal turbines with tip/hub losses. Their abstract reports
/// good power prediction and overpredicted thrust against laboratory data. The
/// paper's equations and numbers were not available for this check.
///
/// This is the closure only. The BEMT solver does not use it yet. Its
/// induction variable `λ_i` maps to `a` only in the windmill regime, and that
/// mapping needs a validation case before it is wired in.
///
/// Returns NaN for a non-finite `a` or `loss`.
#[must_use]
pub fn buhl_momentum_ct(a: f64, loss: f64) -> f64 {
    if !(a.is_finite() && loss.is_finite()) {
        return f64::NAN;
    }
    if a <= 0.4 {
        4.0 * loss * a * (1.0 - a)
    } else {
        8.0 / 9.0 + (4.0 * loss - 40.0 / 9.0) * a + (50.0 / 9.0 - 4.0 * loss) * a * a
    }
}

#[cfg(test)]
mod buhl_tests {
    use super::buhl_momentum_ct;

    #[test]
    fn touches_the_classical_curve_at_forty_per_cent() {
        for loss in [0.5_f64, 0.9, 1.0] {
            let classical = 4.0 * loss * 0.4 * 0.6;
            let buhl =
                8.0 / 9.0 + (4.0 * loss - 40.0 / 9.0) * 0.4 + (50.0 / 9.0 - 4.0 * loss) * 0.16;
            assert!((buhl - classical).abs() < 1e-12, "F={loss}");
            assert!((buhl_momentum_ct(0.4, loss) - classical).abs() < 1e-12);
        }
    }

    #[test]
    fn slope_matches_the_classical_curve_at_forty_per_cent() {
        // The classical slope is 4F(1 − 2a) = 0.8F at a = 0.4.
        let h = 1e-6;
        for loss in [0.7, 0.9] {
            let slope =
                (buhl_momentum_ct(0.4 + h, loss) - buhl_momentum_ct(0.4 - h, loss)) / (2.0 * h);
            assert!((slope - 0.8 * loss).abs() < 1e-5, "F={loss}: slope {slope}");
        }
    }

    #[test]
    fn reaches_two_at_full_induction() {
        for loss in [0.9, 1.0] {
            assert!((buhl_momentum_ct(1.0, loss) - 2.0).abs() < 1e-12);
        }
    }

    #[test]
    fn classical_branch_below_forty_per_cent() {
        assert!((buhl_momentum_ct(0.2, 1.0) - 4.0 * 0.2 * 0.8).abs() < 1e-12);
        assert!(buhl_momentum_ct(f64::NAN, 1.0).is_nan());
    }
}
