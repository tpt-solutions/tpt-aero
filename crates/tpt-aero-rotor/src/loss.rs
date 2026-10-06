//! Rotor loss models: tip and hub corrections to the momentum theory.

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
