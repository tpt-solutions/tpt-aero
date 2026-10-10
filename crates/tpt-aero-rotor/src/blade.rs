//! Rotor/propeller blade geometry and section aerodynamics.

// Brings the `no_std` float methods (`atan`, `powf`, ...) into scope.
use alloc::vec::Vec;
#[allow(unused_imports)]
use tpt_aero_core::Real as _;

use tpt_aero_core::{AirfoilPolar, StallModel};

/// Errors from rotor model construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RotorError {
    /// Radius ratio bounds violated (need 0 ≤ r_hub < r_tip, chord > 0).
    InvalidGeometry,
    /// The BEMT iteration did not converge within the iteration budget.
    NotConverged,
    /// The section model returned a non-physical lift coefficient.
    InvalidSection,
}

impl core::fmt::Display for RotorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            RotorError::InvalidGeometry => {
                write!(
                    f,
                    "invalid rotor geometry (need 0 ≤ r_hub < r_tip, chord > 0)"
                )
            }
            RotorError::NotConverged => write!(f, "BEMT induction iteration did not converge"),
            RotorError::InvalidSection => write!(f, "section lift coefficient out of range"),
        }
    }
}

/// One blade station: chord and geometric twist at a span fraction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BladeSection {
    /// Span fraction `r/R ∈ (0, 1]`.
    pub r_frac: f64,
    /// Chord (m).
    pub chord: f64,
    /// Geometric twist (rad, positive = local pitch increase).
    pub twist: f64,
}

/// Section aerodynamic model: 2-D lift/drag vs local angle of attack.
#[derive(Debug, Clone)]
pub enum SectionModel {
    /// Linear lift with quadratic profile drag — the classic BEMT closure.
    ///
    /// `c_l = cl_slope·α` clamped to ±`cl_max`; `c_d = cd0 + k·c_l²`.
    Linear {
        /// Lift-curve slope (1/rad), typically 2π.
        cl_slope: f64,
        /// Stall lift magnitude.
        cl_max: f64,
        /// Profile (form) drag coefficient at zero lift.
        cd0: f64,
        /// Induced-drag-like quadratic coefficient.
        cd_cl2: f64,
    },
    /// Tabulated airfoil polar (from `tpt-aero-core`), with optional stall
    /// blending beyond the table.
    Polar {
        /// The polar table (at the rotor's operating Mach/Reynolds).
        polar: AirfoilPolar,
        /// Optional stall model applied on top of the table.
        stall: Option<StallModel>,
    },
}

impl SectionModel {
    /// Lift and drag coefficients at local angle of attack `alpha` (rad).
    #[must_use]
    pub fn coefficients(&self, alpha: f64) -> (f64, f64) {
        match self {
            SectionModel::Linear {
                cl_slope,
                cl_max,
                cd0,
                cd_cl2,
            } => {
                let cl = (cl_slope * alpha).clamp(-*cl_max, *cl_max);
                (cl, cd0 + cd_cl2 * cl * cl)
            }
            SectionModel::Polar { polar, stall } => {
                let cl_raw = polar.cl_at(alpha);
                let cl = stall.map_or(cl_raw, |m| m.apply(alpha, cl_raw));
                (cl, polar.cd_at(alpha))
            }
        }
    }
}

/// Rotor geometry: blade count, radii, and spanwise stations.
#[derive(Debug, Clone)]
pub struct RotorGeometry {
    /// Tip radius (m).
    pub tip_radius: f64,
    /// Hub (root cutout) radius (m), 0 ≤ r_hub < r_tip.
    pub hub_radius: f64,
    /// Number of blades.
    pub blades: usize,
    /// Blade stations, ascending in span fraction.
    pub sections: Vec<BladeSection>,
}

impl RotorGeometry {
    /// Validate the geometry.
    ///
    /// # Errors
    /// [`RotorError::InvalidGeometry`] on out-of-range radii, chords, or
    /// non-ascending sections.
    pub fn validate(&self) -> Result<(), RotorError> {
        if !(self.tip_radius > 0.0
            && self.hub_radius >= 0.0
            && self.hub_radius < self.tip_radius
            && self.blades >= 1)
        {
            return Err(RotorError::InvalidGeometry);
        }
        if self.sections.len() < 2 {
            return Err(RotorError::InvalidGeometry);
        }
        if let Some(w) = self
            .sections
            .windows(2)
            .find(|w| w[1].r_frac <= w[0].r_frac)
        {
            let _ = w;
            return Err(RotorError::InvalidGeometry);
        }
        if self.sections.iter().any(|s| {
            !(s.chord > 0.0 && s.chord.is_finite()) || !(s.r_frac > 0.0 && s.r_frac <= 1.0)
        }) {
            return Err(RotorError::InvalidGeometry);
        }
        Ok(())
    }

    /// Rectangular, untwisted blade (constant chord), `n` evenly spaced
    /// stations between the hub and the tip.
    #[must_use]
    pub fn rectangular(
        tip_radius: f64,
        hub_radius: f64,
        blades: usize,
        chord: f64,
        n: usize,
    ) -> Self {
        let n = n.max(2);
        let r0 = hub_radius / tip_radius;
        let sections = (0..n)
            .map(|i| {
                let t = i as f64 / (n - 1) as f64;
                BladeSection {
                    r_frac: r0 + t * (1.0 - r0),
                    chord,
                    twist: 0.0,
                }
            })
            .collect();
        Self {
            tip_radius,
            hub_radius,
            blades,
            sections,
        }
    }

    /// Linear (ideal) twist from `root_twist` at the first station to
    /// `tip_twist` at the last, interpolated on the existing stations.
    #[must_use]
    pub fn with_linear_twist(mut self, root_twist: f64, tip_twist: f64) -> Self {
        let n = self.sections.len();
        for (i, s) in self.sections.iter_mut().enumerate() {
            let t = i as f64 / (n - 1).max(1) as f64;
            s.twist = root_twist + t * (tip_twist - root_twist);
        }
        self
    }
}
