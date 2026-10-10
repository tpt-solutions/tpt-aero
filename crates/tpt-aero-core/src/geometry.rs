//! Lifting-surface geometry primitives.
//!
//! Body frame convention (used across tpt-aero): *x* aft along the chord,
//! *y* out the right wing, *z* up. Leading edges sit at smaller `x` than
//! trailing edges. The freestream at angle of attack α flows as
//! `V·(cos α, 0, +sin α)` — positive α pitches the nose up, so the relative
//! wind arrives from below and moves up through the body frame.

use alloc::{vec, vec::Vec};

use crate::num::Real;

/// 3-component vector (also used for points, as is customary in aero codes).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Vec3<R: Real = f64> {
    /// x component (aft).
    pub x: R,
    /// y component (right span).
    pub y: R,
    /// z component (up).
    pub z: R,
}

impl<R: Real> Vec3<R> {
    /// New vector from components.
    #[must_use]
    pub const fn new(x: R, y: R, z: R) -> Self {
        Self { x, y, z }
    }

    /// Zero vector.
    #[must_use]
    pub fn zero() -> Self {
        Self::new(R::zero(), R::zero(), R::zero())
    }

    /// Euclidean dot product.
    #[must_use]
    pub fn dot(&self, o: &Self) -> R {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Euclidean cross product `self × o`.
    #[must_use]
    pub fn cross(&self, o: &Self) -> Self {
        Self::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }

    /// Euclidean norm.
    #[must_use]
    pub fn norm(&self) -> R {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Unit vector in the same direction (undefined at zero norm).
    #[must_use]
    pub fn normalized(&self) -> Self {
        let n = self.norm();
        Self::new(self.x / n, self.y / n, self.z / n)
    }

    /// Component-wise scaling.
    #[must_use]
    pub fn scaled(&self, s: R) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }
}

impl<R: Real> core::ops::Add for Vec3<R> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}

impl<R: Real> core::ops::Sub for Vec3<R> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}

impl<R: Real> core::ops::Neg for Vec3<R> {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(self.x.neg(), self.y.neg(), self.z.neg())
    }
}

/// A trapezoidal strip of a lifting surface between two span stations.
///
/// The strip is defined by the quarter-chord line endpoints and the chords;
/// `twist` is the geometric (washout/washin) angle about each section's
/// quarter-chord, positive nose-up.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LiftingLinePanel {
    /// Quarter-chord point at the root-side station.
    pub root_qc: Vec3<f64>,
    /// Quarter-chord point at the tip-side station.
    pub tip_qc: Vec3<f64>,
    /// Root-side chord.
    pub root_chord: f64,
    /// Tip-side chord.
    pub tip_chord: f64,
    /// Root-side twist (rad, positive nose-up).
    pub root_twist: f64,
    /// Tip-side twist (rad, positive nose-up).
    pub tip_twist: f64,
}

impl LiftingLinePanel {
    /// Planform reference area of this strip (projected on the x–y plane).
    #[must_use]
    pub fn area(&self) -> f64 {
        let dy = (self.tip_qc.y - self.root_qc.y).abs();
        0.5 * (self.root_chord + self.tip_chord) * dy
    }

    /// Quarter-chord point at fractional span `t ∈ [0, 1]`.
    #[must_use]
    pub fn qc_at(&self, t: f64) -> Vec3<f64> {
        Vec3::new(
            self.root_qc.x + t * (self.tip_qc.x - self.root_qc.x),
            self.root_qc.y + t * (self.tip_qc.y - self.root_qc.y),
            self.root_qc.z + t * (self.tip_qc.z - self.root_qc.z),
        )
    }

    /// Local chord at fractional span `t ∈ [0, 1]`.
    #[must_use]
    pub fn chord_at(&self, t: f64) -> f64 {
        self.root_chord + t * (self.tip_chord - self.root_chord)
    }

    /// Local twist (rad) at fractional span `t ∈ [0, 1]`.
    #[must_use]
    pub fn twist_at(&self, t: f64) -> f64 {
        self.root_twist + t * (self.tip_twist - self.root_twist)
    }
}

/// A quadrilateral surface panel with corners ordered consistently (e.g.
/// counter-clockwise when viewed from the outside of the body).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfacePanel {
    /// Corner points `p0 → p1 → p2 → p3`.
    pub corners: [Vec3<f64>; 4],
}

impl SurfacePanel {
    /// Build from four corners in order.
    #[must_use]
    pub const fn new(p0: Vec3<f64>, p1: Vec3<f64>, p2: Vec3<f64>, p3: Vec3<f64>) -> Self {
        Self {
            corners: [p0, p1, p2, p3],
        }
    }

    /// Area-weighted normal vector (magnitude ≈ 2× area for a planar quad).
    #[must_use]
    pub fn normal(&self) -> Vec3<f64> {
        let d1 = self.corners[2] - self.corners[0];
        let d2 = self.corners[3] - self.corners[1];
        d1.cross(&d2)
    }

    /// Planar-quad area (half the magnitude of the area-weighted normal).
    #[must_use]
    pub fn area(&self) -> f64 {
        self.normal().norm() * 0.5
    }

    /// Corner average.
    #[must_use]
    pub fn centroid(&self) -> Vec3<f64> {
        let mut c = Vec3::zero();
        for p in &self.corners {
            c = c + *p;
        }
        c.scaled(0.25)
    }
}

/// One span station of a wing planform.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct WingSection {
    /// Leading-edge point (body frame).
    pub le: Vec3<f64>,
    /// Chord length (> 0).
    pub chord: f64,
    /// Geometric twist (rad, positive nose-up). Washout is negative twist.
    pub twist: f64,
}

impl WingSection {
    /// Trailing-edge point of the section.
    #[must_use]
    pub fn te(&self) -> Vec3<f64> {
        Vec3::new(self.le.x + self.chord, self.le.y, self.le.z)
    }

    /// Quarter-chord point of the section.
    #[must_use]
    pub fn quarter_chord(&self) -> Vec3<f64> {
        Vec3::new(self.le.x + 0.25 * self.chord, self.le.y, self.le.z)
    }
}

/// A wing planform: an ordered list of span stations from one tip to the
/// other (monotone in `y`).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct WingPlanform {
    sections: Vec<WingSection>,
}

impl WingPlanform {
    /// Build from ordered sections; the span stations must be strictly
    /// increasing in `y` and every chord > 0.
    ///
    /// # Errors
    /// Returns [`GeometryError::InvalidPlanform`] otherwise.
    pub fn from_sections(sections: Vec<WingSection>) -> Result<Self, GeometryError> {
        if sections.len() < 2 {
            return Err(GeometryError::InvalidPlanform);
        }
        if let Some(w) = sections.windows(2).find(|w| w[1].le.y <= w[0].le.y) {
            let _ = w;
            return Err(GeometryError::InvalidPlanform);
        }
        if sections
            .iter()
            .any(|s| !(s.chord > 0.0 && s.chord.is_finite()))
        {
            return Err(GeometryError::InvalidPlanform);
        }
        Ok(Self { sections })
    }

    /// Symmetric rectangular wing centred on `y = 0`.
    #[must_use]
    pub fn rectangular(span: f64, chord: f64) -> Self {
        let half = span * 0.5;
        Self {
            sections: vec![
                WingSection {
                    le: Vec3::new(0.0, -half, 0.0),
                    chord,
                    twist: 0.0,
                },
                WingSection {
                    le: Vec3::new(0.0, half, 0.0),
                    chord,
                    twist: 0.0,
                },
            ],
        }
    }

    /// Symmetric trapezoidal wing centred on `y = 0`: three sections (left
    /// tip, root, right tip), so both semi-spans taper, sweep and rise
    /// identically.
    ///
    /// `sweep_le` is the aft offset of each tip leading edge relative to the
    /// root leading edge (m); `dihedral` (rad) raises (+) both tips in z (a
    /// V-dihedral).
    #[must_use]
    pub fn tapered(
        span: f64,
        root_chord: f64,
        tip_chord: f64,
        sweep_le: f64,
        dihedral: f64,
    ) -> Self {
        let half = span * 0.5;
        let dz = half * dihedral.tan();
        Self {
            sections: vec![
                WingSection {
                    le: Vec3::new(sweep_le, -half, dz),
                    chord: tip_chord,
                    twist: 0.0,
                },
                WingSection {
                    le: Vec3::new(0.0, 0.0, 0.0),
                    chord: root_chord,
                    twist: 0.0,
                },
                WingSection {
                    le: Vec3::new(sweep_le, half, dz),
                    chord: tip_chord,
                    twist: 0.0,
                },
            ],
        }
    }

    /// The span stations, root to tip.
    #[must_use]
    pub fn sections(&self) -> &[WingSection] {
        &self.sections
    }

    /// Total projected planform area (sum over station pairs).
    #[must_use]
    pub fn area(&self) -> f64 {
        self.sections
            .windows(2)
            .map(|w| {
                let dy = w[1].le.y - w[0].le.y;
                0.5 * (w[0].chord + w[1].chord) * dy
            })
            .sum()
    }

    /// Span (tip-to-tip extent in y).
    #[must_use]
    pub fn span(&self) -> f64 {
        self.sections.last().map_or(0.0, |s| s.le.y) - self.sections.first().map_or(0.0, |s| s.le.y)
    }

    /// Mean aerodynamic chord (panel area-weighted standard formula).
    #[must_use]
    pub fn mac(&self) -> f64 {
        let mut num = 0.0;
        let mut den = 0.0;
        for w in self.sections.windows(2) {
            let dy = w[1].le.y - w[0].le.y;
            let (c0, c1) = (w[0].chord, w[1].chord);
            let mac_pair = (2.0 / 3.0) * (c0 + c1 - c0 * c1 / (c0 + c1));
            let area_pair = 0.5 * (c0 + c1) * dy;
            num += mac_pair * area_pair;
            den += area_pair;
        }
        if den.abs() < 1.0e-12 {
            0.0
        } else {
            num / den
        }
    }

    /// Insert `n` evenly *y*-spaced stations by linear interpolation between
    /// the given sections (a lattice generator's first step).
    ///
    /// # Panics
    /// If `n < 2` (a planform always has at least two sections).
    #[must_use]
    pub fn refine(&self, n: usize) -> Self {
        assert!(n >= 2, "refine needs at least 2 stations");
        let mut out = Vec::with_capacity(n);
        let last = self.sections.len() - 1;
        for i in 0..n {
            let t = i as f64 / (n - 1) as f64;
            let s = t * (last as f64);
            let i0 = (s as usize).min(last - 1);
            let f = s - i0 as f64;
            let (a, b) = (&self.sections[i0], &self.sections[i0 + 1]);
            out.push(WingSection {
                le: Vec3::new(
                    a.le.x + f * (b.le.x - a.le.x),
                    a.le.y + f * (b.le.y - a.le.y),
                    a.le.z + f * (b.le.z - a.le.z),
                ),
                chord: a.chord + f * (b.chord - a.chord),
                twist: a.twist + f * (b.twist - a.twist),
            });
        }
        Self { sections: out }
    }

    /// One [`LiftingLinePanel`] per adjacent section pair.
    #[must_use]
    pub fn lifting_line_panels(&self) -> Vec<LiftingLinePanel> {
        self.sections
            .windows(2)
            .map(|w| LiftingLinePanel {
                root_qc: w[0].quarter_chord(),
                tip_qc: w[1].quarter_chord(),
                root_chord: w[0].chord,
                tip_chord: w[1].chord,
                root_twist: w[0].twist,
                tip_twist: w[1].twist,
            })
            .collect()
    }

    /// Quadrilateral surface panels covering the planform (one per strip).
    #[must_use]
    pub fn surface_panels(&self) -> Vec<SurfacePanel> {
        self.sections
            .windows(2)
            .map(|w| SurfacePanel::new(w[0].le, w[0].te(), w[1].te(), w[1].le))
            .collect()
    }
}

/// Errors from geometry construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeometryError {
    /// Sections must be ≥2, strictly increasing in `y`, with positive chords.
    InvalidPlanform,
}

impl core::fmt::Display for GeometryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            GeometryError::InvalidPlanform => {
                write!(
                    f,
                    "planform needs ≥2 sections, strictly increasing y, chord > 0"
                )
            }
        }
    }
}

/// Linear interpolation over `t ∈ [0,1]` between two scalars.
#[must_use]
pub fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + t * (b - a)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vec3_basics() {
        let a = Vec3::new(1.0, 0.0, 0.0);
        let b = Vec3::new(0.0, 1.0, 0.0);
        assert_eq!(a.cross(&b), Vec3::new(0.0, 0.0, 1.0));
        assert_eq!(a.dot(&b), 0.0);
        assert_eq!(Vec3::new(3.0, 4.0, 0.0).norm(), 5.0);
    }

    #[test]
    fn rectangular_planform() {
        let w = WingPlanform::rectangular(10.0, 1.0);
        assert!((w.area() - 10.0).abs() < 1e-12);
        assert!((w.span() - 10.0).abs() < 1e-12);
        assert!((w.mac() - 1.0).abs() < 1e-12);
        let panels = w.lifting_line_panels();
        assert_eq!(panels.len(), 1);
        assert!((panels[0].area() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn tapered_planform_area() {
        // Projected area is exact regardless of dihedral.
        let w = WingPlanform::tapered(10.0, 2.0, 1.0, 0.5, 0.1);
        assert!((w.area() - 0.5 * (2.0 + 1.0) * 10.0).abs() < 1e-12);
        // Surface panels are 3D quads: their area matches the projected area
        // only when there is no dihedral.
        let flat = WingPlanform::tapered(10.0, 2.0, 1.0, 0.5, 0.0);
        let panels = flat.surface_panels();
        assert_eq!(panels.len(), 2, "one panel per semi-span");
        let total: f64 = panels.iter().map(|p| p.area()).sum();
        assert!((total - flat.area()).abs() < 1e-9);
    }

    #[test]
    fn tapered_wing_is_mirror_symmetric() {
        let w = WingPlanform::tapered(10.0, 2.0, 1.0, 0.5, 0.1);
        let s = w.sections();
        assert_eq!(s.len(), 3);
        // Root at y = 0 with the root chord; both tips mirror each other.
        assert_eq!(s[1].le.y, 0.0);
        assert_eq!(s[1].chord, 2.0);
        assert_eq!(s[0].chord, 1.0);
        assert_eq!(s[2].chord, 1.0);
        assert_eq!(s[0].le.y, -s[2].le.y);
        assert_eq!(s[0].le.x, s[2].le.x);
        assert_eq!(s[0].le.z, s[2].le.z);
        assert!((s[2].le.z - 5.0 * 0.1_f64.tan()).abs() < 1e-12);
        // The planform mean chord and the MAC lie between tip and root chord.
        assert!(w.mac() > 1.0 && w.mac() < 2.0);
    }

    #[test]
    fn refine_preserves_endpoints() {
        let w = WingPlanform::tapered(10.0, 2.0, 1.0, 0.5, 0.1).refine(9);
        let s = w.sections();
        assert_eq!(s.len(), 9);
        assert!((s[0].le.y - -5.0).abs() < 1e-12);
        assert!((s[8].le.y - 5.0).abs() < 1e-12);
        // The centre station is the root; the tips carry the tip chord.
        assert!((s[4].chord - 2.0).abs() < 1e-12);
        assert!((s[0].chord - 1.0).abs() < 1e-12 && (s[8].chord - 1.0).abs() < 1e-12);
        assert!((s[0].le.z - 5.0 * 0.1_f64.tan()).abs() < 1e-12);
        assert!((s[8].le.z - 5.0 * 0.1_f64.tan()).abs() < 1e-12);
    }

    #[test]
    fn rejects_bad_planform() {
        let bad = WingPlanform::from_sections(vec![WingSection {
            le: Vec3::zero(),
            chord: 1.0,
            twist: 0.0,
        }]);
        assert_eq!(bad.unwrap_err(), GeometryError::InvalidPlanform);
    }

    #[test]
    fn surface_panel_normal_and_centroid() {
        let p = SurfacePanel::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(1.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        );
        assert!((p.area() - 1.0).abs() < 1e-12);
        let n = p.normal();
        assert!((n.z - 2.0).abs() < 1e-12);
        let c = p.centroid();
        assert!((c.x - 0.5).abs() < 1e-12 && (c.y - 0.5).abs() < 1e-12);
    }
}
