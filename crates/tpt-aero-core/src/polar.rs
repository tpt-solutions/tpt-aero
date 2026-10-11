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

impl core::error::Error for PolarError {}

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
        // NaN must not reach the interpolants; Mach and Reynolds key the table.
        let finite = alphas
            .iter()
            .chain(&cl)
            .chain(&cd)
            .chain(&cm)
            .all(|v| v.is_finite());
        if !(finite && mach.is_finite() && mach >= 0.0 && reynolds.is_finite() && reynolds > 0.0) {
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
        // The blend targets the plateau, but never a value larger in magnitude
        // than the raw lift: a tabulated polar that is already below the plateau
        // (or past it) is left alone rather than pushed up toward the plateau.
        let target = if cl_raw >= 0.0 {
            cl_raw.min(self.cl_max)
        } else {
            cl_raw.max(self.cl_min)
        };
        let t = ((alpha.abs() - self.alpha_stall) / self.width.max(1.0e-9)).clamp(0.0, 1.0);
        let s = t * t * (3.0 - 2.0 * t);
        cl_raw + (target - cl_raw) * s
    }
}

/// A family of polars over a (Mach, Reynolds) grid.
///
/// Lookup interpolates linearly between bracketing tables; tables may have
/// different α grids because values are evaluated per-table before blending.
#[derive(Debug, Clone, Default)]
pub struct PolarSet {
    polars: Vec<AirfoilPolar>,
    /// Distinct Mach numbers, ascending.
    machs: Vec<f64>,
    /// For each Mach value, the Reynolds numbers of its tables (ascending).
    reynolds_at: Vec<Vec<f64>>,
    /// For each Mach value, indices into `polars` in the same Reynolds order.
    index_at: Vec<Vec<usize>>,
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
        // Precompute the bracketing grid so lookups do not allocate.
        let mut machs: Vec<f64> = Vec::new();
        let mut reynolds_at: Vec<Vec<f64>> = Vec::new();
        let mut index_at: Vec<Vec<usize>> = Vec::new();
        for (i, p) in polars.iter().enumerate() {
            if machs.last() != Some(&p.mach()) {
                machs.push(p.mach());
                reynolds_at.push(Vec::new());
                index_at.push(Vec::new());
            }
            let g = machs.len() - 1;
            reynolds_at[g].push(p.reynolds());
            index_at[g].push(i);
        }
        Ok(Self {
            polars,
            machs,
            reynolds_at,
            index_at,
        })
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

    /// Weighted blend over the up-to-four bracketing (Mach, Re) tables. NaN
    /// for an empty set (the constructor refuses one).
    fn blend<F>(&self, mach: f64, reynolds: f64, eval: F) -> f64
    where
        F: Fn(&AirfoilPolar) -> f64,
    {
        if self.machs.is_empty() {
            return f64::NAN;
        }
        let (m0, m1, wm) = bracket(&self.machs, mach);
        let section = |mach_idx: usize| -> f64 {
            let (r0, r1, wr) = bracket(&self.reynolds_at[mach_idx], reynolds);
            let idx = &self.index_at[mach_idx];
            eval(&self.polars[idx[r0]]) * (1.0 - wr) + eval(&self.polars[idx[r1]]) * wr
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
            PolarSet::new(vec![flat_polar(4.0, 0.0, 1e6), flat_polar(6.0, 0.0, 1e6)]).unwrap_err(),
            PolarError::InvalidTable
        );
    }
}

#[cfg(test)]
mod validation_tests {
    use super::*;
    use alloc::vec;

    #[test]
    fn polar_rejects_nan_and_bad_flight_condition() {
        let ok = || vec![-1.0, 0.0, 1.0];
        let v = || vec![0.1, 0.0, 0.2];
        assert!(AirfoilPolar::new(ok(), v(), v(), v(), 0.2, 1.0e6).is_ok());
        let nan = vec![-1.0, f64::NAN, 1.0];
        assert!(AirfoilPolar::new(ok(), nan, v(), v(), 0.2, 1.0e6).is_err());
        assert!(AirfoilPolar::new(ok(), v(), v(), v(), f64::NAN, 1.0e6).is_err());
        assert!(AirfoilPolar::new(ok(), v(), v(), v(), 0.2, 0.0).is_err());
    }
}

#[cfg(test)]
mod stall_and_set_tests {
    use super::*;
    use alloc::vec;

    #[test]
    fn stall_blend_never_raises_lift_magnitude() {
        let m = StallModel::symmetric(0.2, 0.1, 1.2);
        // Below the plateau and past the stall angle: must not be pushed up.
        assert!(m.apply(0.5, 0.4) <= 0.4 + 1e-12);
        // Above the plateau: moves toward it, downward.
        assert!(m.apply(0.5, 1.5) < 1.5 && m.apply(0.5, 1.5) >= 1.2 - 1e-12);
        // Negative side mirrors it.
        assert!(m.apply(-0.5, -0.4) >= -0.4 - 1e-12);
    }

    #[test]
    fn polar_set_blend_matches_direct_lookup_on_the_grid() {
        let mk = |mach: f64, re: f64, cl: f64| {
            AirfoilPolar::new(
                vec![-1.0, 0.0, 1.0],
                vec![-cl, 0.0, cl],
                vec![0.01, 0.01, 0.01],
                vec![0.0, 0.0, 0.0],
                mach,
                re,
            )
            .expect("valid")
        };
        let set = PolarSet::new(vec![
            mk(0.0, 1e6, 1.0),
            mk(0.0, 2e6, 1.2),
            mk(0.5, 1e6, 1.4),
        ])
        .expect("valid");
        assert!((set.cl_at(0.0, 1e6, 0.5) - 0.5).abs() < 1e-12);
        assert!((set.cl_at(0.0, 1.5e6, 0.5) - 0.55).abs() < 1e-9);
        assert!((set.cl_at(0.25, 1e6, 0.5) - 0.6).abs() < 1e-9);
        assert!(PolarSet::default().cl_at(0.0, 1e6, 0.5).is_nan());
    }
}
