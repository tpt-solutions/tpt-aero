//! Full longitudinal plant with speed and drag couplings: short period *and*
//! phugoid.
//!
//! [`crate::plant::LongitudinalPlant`] keeps only the `(α, q)` short-period
//! block coupled to a decoupled speed row. This module adds the terms that
//! create the long-period mode: speed sensitivity of lift and drag (`2C_L`,
//! `2C_D`), the drag-due-to-α slope, and the `C_D`/`C_L` pieces of `X_α` and
//! `Z_α`. Low-Mach assumptions: `C_{Lu} = C_{Du} = C_{mu} = 0`, thrust held
//! constant, level trim (`γ₀ = 0`), stability axes.
//!
//! The classical Lanchester approximation (`ω_p = √2·g/V`,
//! `ζ_p = 1/(√2·L/D)`) is the analytic cross-check used in the tests.

#[allow(unused_imports)]
use tpt_aero_core::Real as _;
use tpt_aero_vlm::VlmModel;

use crate::controls::ControlDerivatives;
use crate::longitudinal::LongitudinalDerivatives;
use crate::modes::{eigenvalues4, Complex};
use crate::sim::StateSpace;

/// Standard gravity (m/s²).
const GRAVITY: f64 = 9.80665;

/// Flight condition and mass properties for [`FullLongitudinalPlant`].
#[derive(Debug, Clone, Copy)]
pub struct FullLongitudinalParameters {
    /// Trim airspeed (m/s).
    pub v_inf: f64,
    /// Reference chord (m).
    pub chord: f64,
    /// Reference area (m²).
    pub area: f64,
    /// Air density (kg/m³).
    pub density: f64,
    /// Mass (kg).
    pub mass: f64,
    /// Pitch inertia `I_y` (kg m²).
    pub inertia_y: f64,
    /// Trim lift coefficient (`W/(qS)` in level flight).
    pub cl: f64,
    /// Total trim drag coefficient (induced plus zero-lift).
    pub cd: f64,
    /// Drag-due-to-α slope `∂C_D/∂α` (per rad); see [`drag_alpha_fd`].
    pub cd_alpha: f64,
}

/// `∂C_{Di}/∂α` by central difference of the VLM near-field induced drag.
/// Zero-lift drag does not depend on α in this model, so this is the whole
/// slope unless the caller adds profile-drag growth.
///
/// Returns `None` if either solve fails.
#[must_use]
pub fn drag_alpha_fd(model: &VlmModel, alpha_rad: f64, v_inf: f64, density: f64) -> Option<f64> {
    let h = 1.0e-4;
    let cdi = |a: f64| model.solve(a, v_inf, density, &[]).ok().map(|s| s.cdi);
    Some((cdi(alpha_rad + h)? - cdi(alpha_rad - h)?) / (2.0 * h))
}

/// Linearized longitudinal plant `x' = A x`, state `[u, α, q, θ]`
/// (m/s, rad, rad/s, rad), about level flight:
///
/// ```text
/// u'  = -g θ - (qS/mV) 2 C_D u - (qS/m)(C_Dα - C_L) α
/// α'  = q - (qS/mV²) 2 C_L u - (qS/mV)(C_Lα + C_D) α - (qS/mV)(c/2V) C_Lq q
/// q'  = (qS c/I_y)(C_mα α + C_mq (c/2V) q)
/// θ'  = q
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FullLongitudinalPlant {
    /// Dynamics matrix (4×4, row-major).
    pub a: [[f64; 4]; 4],
}

/// The two classical longitudinal oscillations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LongitudinalModes {
    /// Short-period natural frequency (rad/s).
    pub short_period_wn: f64,
    /// Short-period damping ratio.
    pub short_period_zeta: f64,
    /// Phugoid natural frequency (rad/s).
    pub phugoid_wn: f64,
    /// Phugoid damping ratio (negative = divergent).
    pub phugoid_zeta: f64,
}

impl FullLongitudinalPlant {
    /// Assemble the plant. Returns `None` for non-positive speed, chord,
    /// area, density, mass or inertia, or non-finite coefficients.
    #[must_use]
    pub fn from_derivatives(
        d: &LongitudinalDerivatives,
        p: &FullLongitudinalParameters,
    ) -> Option<Self> {
        let positive = |x: f64| x.is_finite() && x > 0.0;
        if !(positive(p.v_inf)
            && positive(p.chord)
            && positive(p.area)
            && positive(p.density)
            && positive(p.mass)
            && positive(p.inertia_y)
            && p.cl.is_finite()
            && p.cd.is_finite()
            && p.cd_alpha.is_finite())
        {
            return None;
        }
        let q_s = 0.5 * p.density * p.v_inf * p.v_inf * p.area;
        let v = p.v_inf;
        let half_c_over_v = p.chord / (2.0 * v);
        let per_mass = q_s / p.mass;
        let m_scale = q_s * p.chord / p.inertia_y;
        let a = [
            [
                -per_mass * 2.0 * p.cd / v,
                -per_mass * (p.cd_alpha - p.cl),
                0.0,
                -GRAVITY,
            ],
            [
                -per_mass * 2.0 * p.cl / (v * v),
                -per_mass * (d.cl_alpha + p.cd) / v,
                1.0 - per_mass * d.cl_q * half_c_over_v / v,
                0.0,
            ],
            [
                0.0,
                d.cm_alpha * m_scale,
                d.cm_q * m_scale * half_c_over_v,
                0.0,
            ],
            [0.0, 0.0, 1.0, 0.0],
        ];
        Some(Self { a })
    }

    /// Eigenvalues of the plant.
    #[must_use]
    pub fn eigenvalues(&self) -> [Complex; 4] {
        eigenvalues4(&self.a)
    }

    /// The plant as a [`StateSpace`] with inputs `[elevator, w_gust]`:
    /// `elevator` is the control command (rad, VLM sign convention, see
    /// [`crate::controls`]) with effectiveness `elevator`; `w_gust` is the
    /// vertical air velocity (m/s, positive = updraft), which raises the
    /// aerodynamic angle of attack by `w_gust / V` without moving the
    /// airframe's own attitude.
    #[must_use]
    pub fn state_space(
        &self,
        p: &FullLongitudinalParameters,
        elevator: &ControlDerivatives,
    ) -> StateSpace<4, 2> {
        let q_s = 0.5 * p.density * p.v_inf * p.v_inf * p.area;
        let per_mass = q_s / p.mass;
        let m_scale = q_s * p.chord / p.inertia_y;
        let mut b = [[0.0; 2]; 4];
        b[0][0] = -per_mass * elevator.cd;
        b[1][0] = -per_mass * elevator.cl / p.v_inf;
        b[2][0] = m_scale * elevator.cm;
        // Every aerodynamic dependence on alpha also sees the gust.
        for i in 0..3 {
            b[i][1] = self.a[i][1] / p.v_inf;
        }
        StateSpace { a: self.a, b }
    }

    /// Split the spectrum into short period and phugoid (the two complex
    /// pairs, ordered by frequency). `None` if the spectrum is not two
    /// complex pairs, e.g. a statically unstable airframe or an overdamped
    /// mode.
    #[must_use]
    pub fn modes(&self) -> Option<LongitudinalModes> {
        let eig = self.eigenvalues();
        let mut pairs: [Option<Complex>; 2] = [None, None];
        for z in eig.iter().filter(|z| z.im > 0.0) {
            if pairs[0].is_none() {
                pairs[0] = Some(*z);
            } else if pairs[1].is_none() {
                pairs[1] = Some(*z);
            } else {
                return None;
            }
        }
        let (a, b) = (pairs[0]?, pairs[1]?);
        let (phugoid, short) = if a.abs() < b.abs() { (a, b) } else { (b, a) };
        let describe = |z: Complex| (z.abs(), -z.re / z.abs());
        let (short_period_wn, short_period_zeta) = describe(short);
        let (phugoid_wn, phugoid_zeta) = describe(phugoid);
        Some(LongitudinalModes {
            short_period_wn,
            short_period_zeta,
            phugoid_wn,
            phugoid_zeta,
        })
    }
}
