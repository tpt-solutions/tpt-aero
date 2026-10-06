//! Airfoil polar tables: `C_l`/`C_d`/`C_m` interpolation over α, Mach and
//! Reynolds number, with smooth stall modeling.

use alloc::vec::Vec;

/// Errors from polar table construction.
#[derive(Debug, Clone, PartialEq)]
pub enum PolarError {
    /// Table arrays differ in length or are too short.
    InvalidTable,
    /// The α grid is not strictly increasing.
    NonMonotonicAlpha,
}

impl core::fmt::Display for PolarError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            PolarError::InvalidTable => write!(f, "polar arrays must be equal length ≥ 2"),
            PolarError::NonMonotonicAlpha => write!(f, "alpha grid must be strictly increasing"),
        }
    }
}

/// Linear interpolation on a strictly increasing grid, clamped at the ends.
pub(crate) fn interp_linear(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    debug_assert_eq!(xs.len(), ys.len());
    if x <= xs[0] {
        return ys[0];
    }
    let last = xs.len() - 1;
    if x >= xs[last] {
        return ys[last];
    }
    let mut lo = 0usize;
    let mut hi = last;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if xs[mid] <= x {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let t = (x - xs[lo]) / (xs[hi] - xs[lo]);
    ys[lo] * (1.0 - t) + ys[hi] * t
}

/// A tabulated airfoil polar at fixed Mach and Reynolds number.
///
/// The α grid must be strictly increasing. Outside the tabulated range the end
/// values are held (constant extrapolation), which the caller is expected to
/// combine with a [`StallModel`] before trusting at high α.
#[derive(Debug, Clone, PartialEq)]
pub struct AirfoilPolar {
    alphas: Vec<f64>,
    cl: Vec<f64>,
    cd: Vec<f64>,
    cm: Vec<f64>,
    mach: f64,
    reynolds: f64,
}

impl AirfoilPolar {
    /// Build a polar table.
    ///
    /// # Errors
    /// [`PolarError::InvalidTable`] for mismatched/short arrays,
    /// [`PolarError::NonMonotonicAlpha`] for a non-increasing α grid.
    pub fn new(
        alphas: Vec<f64>,
        cl: Vec<f64>,
        cd: Vec<f64>,
        cm: Vec<f64>,
        mach: f64,
        reynolds: f64,
    ) -> Result<Self, PolarError> {
        let n = alphas.len();
        if n < 2 || cl.len() != n || cd.len() != n || cm.len() != n {
            return Err(PolarError::InvalidTable);
        }
        if alphas.windows(2).any(|w| w[1] <= w[0]) {
            return Err(PolarError::NonMonotonicAlpha);
        }
        Ok(Self {
            alphas,
            cl,
            cd,
            cm,
            mach,
            reynolds,
        })
    }

    /// The tabulated Mach number.
    #[must_use]
    pub const fn mach(&self) -> f64 {
        self.mach
    }

    /// The tabulated Reynolds number.
    #[must_use]
    pub const fn reynolds(&self) -> f64 {
        self.reynolds
    }

    /// The α grid in radians.
    #[must_use]
    pub fn alphas(&self) -> &[f64] {
        &self.alphas
    }

    /// Section lift coefficient at α (linear, clamped).
    #[must_use]
    pub fn cl_at(&self, alpha: f64) -> f64 {
        interp_linear(&self.alphas, &self.cl, alpha)
    }

    /// Section drag coefficient at α (linear, clamped).
    #[must_use]
    pub fn cd_at(&self, alpha: f64) -> f64 {
        interp_linear(&self.alphas, &self.cd, alpha)
    }

    /// Section pitching-moment coefficient at α (linear, clamped).
    #[must_use]
    pub fn cm_at(&self, alpha: f64) -> f64 {
        interp_linear(&self.alphas, &self.cm, alpha)
    }

    /// Lift-curve slope dC_l/dα from a central difference at the middle of
    /// the α grid (1/rad).
    #[must_use]
    pub fn cl_alpha_slope(&self) -> f64 {
        let n = self.alphas.len();
        if n < 3 {
            return (self.cl[n - 1] - self.cl[0]) / (self.alphas[n - 1] - self.alphas[0]);
        }
        let i = n / 2;
        (self.cl[i + 1] - self.cl[i - 1]) / (self.alphas[i + 1] - self.alphas[i - 1])
    }

    /// Apply a smooth stall saturation to a raw (linear) lift value at α.
    ///
    /// See [`StallModel`].
    #[must_use]
    pub fn cl_with_stall(&self, alpha: f64, model: &StallModel) -> f64 {
        model.apply(alpha, self.cl_at(alpha))
    }
}

/// Smooth lift saturation beyond the stall angles.
///
/// Between −α_stall and +α_stall the table (or linear) value passes through
/// unchanged; outside, the lift is driven smoothly toward the plateau `c_l_max`
/// (resp. `c_l_min`) with a smoothstep over `width` radians. This is a
/// engineering approximation adequate for stability-and-control work — it is
/// *not* a separated-flow model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StallModel {
    /// Absolute α beyond which saturation starts (rad).
    pub alpha_stall: f64,
    /// Saturation band width (rad).
    pub width: f64,
    /// Upper lift plateau.
    pub cl_max: f64,
    /// Lower lift plateau.
    pub cl_min: f64,
}

impl StallModel {
    /// A symmetric model with equal and opposite plateaus.
    #[must_use]
    pub const fn symmetric(alpha_stall: f64, width: f64, cl_max: f64) -> Self {
        Self {
            alpha_stall,
            width,
            cl_max,
            cl_min: -cl_max,
        }
    }

    /// Effective section lift at angle of attack `alpha` for the raw
    /// (unstalled) value `cl_raw`.
    ///
    /// For |α| ≤ α_stall the raw value passes through unchanged. Beyond, the
    /// value is blended smoothly (smoothstep over `width` radians) toward the
    /// plateau `cl_max`/`cl_min`, preserving the sign of the raw lift.
    #[must_use]
    pub fn apply(&self, alpha: f64, cl_raw: f64) -> f64 {
        let plateau = if cl_raw >= 0.0 { self.cl_max } else { self.cl_min };
        let t = ((alpha.abs() - self.alpha_stall) / self.width.max(1.0e-9)).clamp(0.0, 1.0);
        let s = t * t * (3.0 - 2.0 * t);
        cl_raw + (plateau - cl_raw) * s
    }
}

/// A family of polars over a (Mach, Reynolds) grid.
///
/// Lookup interpolates linearly between bracketing tables; tables may have
/// different α grids because values are evaluated per-table before blending.
#[derive(Debug, Clone, Default)]
pub struct PolarSet {
    polars: Vec<AirfoilPolar>,
}

impl PolarSet {
    /// Build from a list of polars (need not be pre-sorted).
    ///
    /// # Errors
    /// [`PolarError::InvalidTable`] if the set is empty or contains duplicate
    /// (Mach, Re) stations.
    pub fn new(polars: Vec<AirfoilPolar>) -> Result<Self, PolarError> {
        if polars.is_empty() {
            return Err(PolarError::InvalidTable);
        }
        let mut polars = polars;
        polars.sort_by(|a, b| {
            a.mach()
                .total_cmp(&b.mach())
                .then(a.reynolds().total_cmp(&b.reynolds()))
        });
        for w in polars.windows(2) {
            if w[0].mach() == w[1].mach() && w[0].reynolds() == w[1].reynolds() {
                return Err(PolarError::InvalidTable);
            }
        }
        Ok(Self { polars })
    }

    /// Number of tabulated (Mach, Re) stations.
    #[must_use]
    pub fn len(&self) -> usize {
        self.polars.len()
    }

    /// True if empty (cannot happen via `new`).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.polars.is_empty()
    }

    /// Weighted blend over the up-to-four bracketing (Mach, Re) tables.
    fn blend<F>(&self, mach: f64, reynolds: f64, eval: F) -> f64
    where
        F: Fn(&AirfoilPolar) -> f64,
    {
        // Bracket mach among distinct mach values.
        let machs: Vec<f64> = {
            let mut m = Vec::new();
            for p in &self.polars {
                if m.last() != Some(&p.mach()) {
                    m.push(p.mach());
                }
            }
            m
        };
        let (m0, m1, wm) = bracket(&machs, mach);

        let section = |mach_idx: usize| -> f64 {
            let at_mach: Vec<&AirfoilPolar> =
                self.polars.iter().filter(|p| p.mach() == machs[mach_idx]).collect();
            let res: Vec<f64> = at_mach.iter().map(|p| p.reynolds()).collect();
            let (r0, r1, wr) = bracket(&res, reynolds);
            eval(at_mach[r0]) * (1.0 - wr) + eval(at_mach[r1]) * wr
        };

        if m0 == m1 {
            section(m0)
        } else {
            section(m0) * (1.0 - wm) + section(m1) * wm
        }
    }

    /// Blended section lift coefficient at (M, Re, α).
    #[must_use]
    pub fn cl_at(&self, mach: f64, reynolds: f64, alpha: f64) -> f64 {
        self.blend(mach, reynolds, |p| p.cl_at(alpha))
    }

    /// Blended section drag coefficient at (M, Re, α).
    #[must_use]
    pub fn cd_at(&self, mach: f64, reynolds: f64, alpha: f64) -> f64 {
        self.blend(mach, reynolds, |p| p.cd_at(alpha))
    }

    /// Blended section pitching-moment coefficient at (M, Re, α).
    #[must_use]
    pub fn cm_at(&self, mach: f64, reynolds: f64, alpha: f64) -> f64 {
        self.blend(mach, reynolds, |p| p.cm_at(alpha))
    }
}

/// Bracket `x` in an ascending list: returns (lo_index, hi_index, weight).
/// Outside the range, the nearest segment is used with weight 0 or 1.
fn bracket(xs: &[f64], x: f64) -> (usize, usize, f64) {
    debug_assert!(!xs.is_empty());
    if xs.len() == 1 || x <= xs[0] {
        return (0, 0, 0.0);
    }
    let last = xs.len() - 1;
    if x >= xs[last] {
        return (last, last, 0.0);
    }
    let mut lo = 0usize;
    let mut hi = last;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if xs[mid] <= x {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let w = (x - xs[lo]) / (xs[hi] - xs[lo]);
    (lo, hi, w)
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn flat_polar(cl_slope: f64, mach: f64, re: f64) -> AirfoilPolar {
        let alphas: Vec<f64> = (-10..=15).map(|d| f64::from(d).to_radians()).collect();
        let cl: Vec<f64> = alphas.iter().map(|&a| cl_slope * a).collect();
        let cd: Vec<f64> = alphas.iter().map(|&a| 0.01 + 0.05 * a * a).collect();
        let cm = vec![-0.05; alphas.len()];
        AirfoilPolar::new(alphas, cl, cd, cm, mach, re).expect("valid")
    }

    #[test]
    fn linear_interpolation_midpoint() {
        let p = flat_polar(TWO_PI, 0.0, 1e6);
        let mid = 2.5_f64.to_radians();
        assert!((p.cl_at(mid) - TWO_PI * mid).abs() < 1e-12);
    }

    const TWO_PI: f64 = 2.0 * core::f64::consts::PI;

    #[test]
    fn clamped_extrapolation() {
        let p = flat_polar(2.0, 0.0, 1e6);
        let a_max = *p.alphas().last().expect("non-empty");
        assert_eq!(p.cl_at(a_max + 0.5), p.cl_at(a_max));
        assert_eq!(p.cl_at(-10.0), p.cl_at(p.alphas()[0]));
    }

    #[test]
    fn stall_saturates() {
        let p = flat_polar(TWO_PI, 0.0, 1e6);
        let model = StallModel::symmetric(15.0_f64.to_radians(), 3.0_f64.to_radians(), 1.2);
        let unstalled = p.cl_at(10.0_f64.to_radians());
        let stalled = model.apply(25.0_f64.to_radians(), p.cl_at(25.0_f64.to_radians()));
        assert!(stalled <= 1.2 + 1e-9);
        assert!((model.apply(10.0_f64.to_radians(), unstalled) - unstalled).abs() < 1e-12);
    }

    #[test]
    fn polar_set_blends_mach() {
        let set = PolarSet::new(vec![flat_polar(4.0, 0.0, 1e6), flat_polar(6.0, 0.4, 1e6)])
            .expect("valid");
        let alpha = 3.0_f64.to_radians();
        let cl_mid = set.cl_at(0.2, 1e6, alpha);
        assert!((cl_mid - 5.0 * alpha).abs() < 1e-9, "{cl_mid}");
    }

    #[test]
    fn polar_set_rejects_duplicates() {
        assert_eq!(
            PolarSet::new(vec![flat_polar(4.0, 0.0, 1e6), flat_polar(6.0, 0.0, 1e6)])
                .unwrap_err(),
            PolarError::InvalidTable
        );
    }
}
