//! Lattice generation: discretize a [`WingPlanform`] into horseshoe elements.
//!
//! Convention (Katz & Plotkin): the bound vortex sits on the panel's
//! quarter-chord line, the collocation point at the panel's three-quarter
//! chord, mid-span. Trailing legs extend to infinity along the freestream.
//! Camber and twist enter through the surface normal at the collocation
//! point, so the AIC assembly stays linear in the unknowns.

use alloc::vec::Vec;

use tpt_aero_core::geometry::WingPlanform;
use tpt_aero_core::{Real, Vec3};

/// Errors from lattice construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LatticeError {
    /// Fewer than 1 spanwise or chordwise panel requested.
    InvalidResolution,
    /// Planform validation failed.
    InvalidPlanform,
    /// NACA camber outside its domain: `m` must be finite, and when it is
    /// non-zero `p` must lie strictly between 0 and 1 (`p = 1` divides by zero).
    InvalidCamber,
}

impl core::fmt::Display for LatticeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            LatticeError::InvalidResolution => {
                write!(f, "need at least 1 spanwise and 1 chordwise panel")
            }
            LatticeError::InvalidPlanform => write!(f, "invalid wing planform"),
            LatticeError::InvalidCamber => {
                write!(f, "NACA camber needs finite m and 0 < p < 1 when m != 0")
            }
        }
    }
}

impl core::error::Error for LatticeError {}

/// Camber model applied to the lattice surface normals.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum CamberModel {
    /// Flat plate (zero camber).
    #[default]
    None,
    /// NACA 4-digit-style circular-arc-free camber line: maximum camber `m`
    /// (chord fraction) at chordwise position `p`.
    Naca4 {
        /// Maximum camber as a chord fraction (e.g. 0.02).
        m: f64,
        /// Chordwise position of maximum camber (e.g. 0.4).
        p: f64,
    },
}

impl CamberModel {
    /// Camber-line slope dz/dx at chordwise fraction `x`.
    #[must_use]
    pub fn slope_at(&self, x: f64) -> f64 {
        match *self {
            CamberModel::None => 0.0,
            CamberModel::Naca4 { m, p } => {
                if m == 0.0 || p == 0.0 {
                    0.0
                } else if x < p {
                    2.0 * m / (p * p) * (p - x)
                } else {
                    2.0 * m / ((1.0 - p) * (1.0 - p)) * (p - x)
                }
            }
        }
    }

    /// Camber-line ordinate z/c at chordwise fraction `x`.
    #[must_use]
    pub fn z_at(&self, x: f64) -> f64 {
        match *self {
            CamberModel::None => 0.0,
            CamberModel::Naca4 { m, p } => {
                if m == 0.0 || p == 0.0 {
                    0.0
                } else if x < p {
                    m / (p * p) * (2.0 * p * x - x * x)
                } else {
                    m / ((1.0 - p) * (1.0 - p)) * (1.0 - 2.0 * p + 2.0 * p * x - x * x)
                }
            }
        }
    }
}

/// Spanwise distribution of the lattice strips.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SpanSpacing {
    /// Equal strip widths (the default; matches [`Lattice::build`]).
    #[default]
    Uniform,
    /// Cosine-clustered toward the tips: narrow strips where the loading
    /// gradient is steepest. Opt-in, see [`WingPlanform::refine_cosine`].
    Cosine,
}

/// One horseshoe element of the lattice (fixed `f64` geometry).
///
/// Twist is *not* baked into the geometry: it is applied when the AIC is
/// assembled, so a differentiable assembly can seed twist with dual numbers.
#[derive(Debug, Clone)]
pub struct LatticeElement {
    /// Left (smaller `y`) end of the bound vortex (quarter chord).
    pub bound_a: Vec3<f64>,
    /// Right (larger `y`) end of the bound vortex (quarter chord).
    pub bound_b: Vec3<f64>,
    /// Collocation point (three-quarter chord, mid-span).
    pub collocation: Vec3<f64>,
    /// Surface normal at the collocation point ignoring camber and twist
    /// (i.e. from the flat, dihedral-only surface), unit length.
    pub flat_normal: Vec3<f64>,
    /// Unit vector along the bound vortex (`a → b`), the twist/camber
    /// rotation axis at this panel.
    pub span_axis: Vec3<f64>,
    /// Camber slope dz/dx at the collocation chordwise station.
    pub camber_slope: f64,
    /// Panel chord.
    pub chord: f64,
    /// Panel planform (projected) area.
    pub area: f64,
    /// Spanwise strip index (0-based, from the −y tip to the +y tip).
    pub strip: usize,
    /// Trailing-edge `x` of the panel at the collocation station (moment arm).
    pub x_te: f64,
    /// Local span of the strip containing this panel.
    pub strip_dy: f64,
    /// Local strip chord (mid-span).
    pub strip_chord: f64,
    /// Geometric twist at the collocation station (rad, positive nose-up).
    pub twist: f64,
    /// Index of the surface this element belongs to (see [`Lattice::combine`]).
    pub surface: usize,
    /// Spanwise coordinate of the strip centre in the surface's own planform
    /// frame (m, before any placement transform).
    pub eta: f64,
    /// Chordwise fraction (0 = leading edge, 1 = trailing edge) of the
    /// collocation point on the local chord.
    pub chord_frac: f64,
}

impl LatticeElement {
    /// Rigidly roll this element about the body *x* axis by `roll_rad`
    /// (right-hand rule: `+y` toward `+z`), then translate by `offset`.
    ///
    /// Directions (normal, span axis) rotate only; points translate too. The
    /// strip index is untouched (see [`Lattice::combine`]).
    #[must_use]
    pub fn transformed(&self, roll_rad: f64, offset: Vec3<f64>) -> Self {
        let (s, c) = roll_rad.sin_cos();
        let rot = |v: Vec3<f64>| Vec3::new(v.x, v.y * c - v.z * s, v.y * s + v.z * c);
        let place = |p: Vec3<f64>| rot(p) + offset;
        Self {
            bound_a: place(self.bound_a),
            bound_b: place(self.bound_b),
            collocation: place(self.collocation),
            flat_normal: rot(self.flat_normal),
            span_axis: rot(self.span_axis),
            x_te: self.x_te + offset.x,
            ..self.clone()
        }
    }

    /// Effective surface slope angle (rad) at the collocation point for a
    /// given local twist: geometric twist minus camber-line incidence
    /// (a surface sloping up toward +x, dz/dx > 0, tilts its unit normal to
    /// `n ∝ (−dz/dx, 0, 1)`, i.e. rotation `ε = twist − atan(dz/dx)` about
    /// the span axis; nose-up twist is positive and adds incidence).
    #[must_use]
    pub fn effective_twist<R: Real>(&self, twist: R) -> R {
        twist - R::from_f64(self.camber_slope.atan2(1.0))
    }

    /// Surface normal at the collocation point for a given local twist
    /// (Rodrigues rotation of the flat normal about the span axis).
    #[must_use]
    pub fn normal_for_twist<R: Real>(&self, twist: R) -> Vec3<R> {
        let n0 = Vec3::new(
            R::from_f64(self.flat_normal.x),
            R::from_f64(self.flat_normal.y),
            R::from_f64(self.flat_normal.z),
        );
        let s = Vec3::new(
            R::from_f64(self.span_axis.x),
            R::from_f64(self.span_axis.y),
            R::from_f64(self.span_axis.z),
        );
        let eps = self.effective_twist(twist);
        let (sin_e, cos_e) = eps.sin_cos();
        let sn = s.cross(&n0);
        let s_dot_n = s.dot(&n0);
        // n = n0 cos ε + (s × n0) sin ε + s (s·n0)(1 − cos ε)
        n0.scaled(cos_e) + sn.scaled(sin_e) + s.scaled(s_dot_n * (R::one() - cos_e))
    }

    /// Bound-vortex endpoints in `R` (constant geometry).
    #[must_use]
    pub fn bounds<R: Real>(&self) -> (Vec3<R>, Vec3<R>) {
        (
            Vec3::new(
                R::from_f64(self.bound_a.x),
                R::from_f64(self.bound_a.y),
                R::from_f64(self.bound_a.z),
            ),
            Vec3::new(
                R::from_f64(self.bound_b.x),
                R::from_f64(self.bound_b.y),
                R::from_f64(self.bound_b.z),
            ),
        )
    }

    /// Collocation point in `R` (constant geometry).
    #[must_use]
    pub fn cp<R: Real>(&self) -> Vec3<R> {
        Vec3::new(
            R::from_f64(self.collocation.x),
            R::from_f64(self.collocation.y),
            R::from_f64(self.collocation.z),
        )
    }
}

/// A horseshoe lattice for one lifting surface.
#[derive(Debug, Clone)]
pub struct Lattice {
    elements: Vec<LatticeElement>,
    n_span: usize,
    n_chord: usize,
}

impl Lattice {
    /// Build a lattice from a planform.
    ///
    /// `n_span` strips spanwise (uniform in *y*), `n_chord` panels chordwise
    /// (cosine clustered at leading/trailing edges, which keeps the
    /// quarter-chord line dense where the load gradient is largest).
    ///
    /// # Errors
    /// [`LatticeError::InvalidResolution`] if either count is zero,
    /// [`LatticeError::InvalidPlanform`] if the planform is invalid.
    pub fn build(
        planform: &WingPlanform,
        n_span: usize,
        n_chord: usize,
        camber: CamberModel,
    ) -> Result<Self, LatticeError> {
        Self::build_spaced(planform, n_span, n_chord, camber, SpanSpacing::Uniform)
    }

    /// Like [`Lattice::build`], with an explicit spanwise strip distribution.
    ///
    /// # Errors
    /// As [`Lattice::build`].
    pub fn build_spaced(
        planform: &WingPlanform,
        n_span: usize,
        n_chord: usize,
        camber: CamberModel,
        spacing: SpanSpacing,
    ) -> Result<Self, LatticeError> {
        if n_span == 0 || n_chord == 0 {
            return Err(LatticeError::InvalidResolution);
        }
        // Every planform segment needs at least one spanwise strip, so that
        // kinks stay on strip edges.
        if n_span < planform.sections().len() - 1 {
            return Err(LatticeError::InvalidResolution);
        }
        if let CamberModel::Naca4 { m, p } = camber {
            let p_ok = m == 0.0 || (p > 0.0 && p < 1.0);
            if !(m.is_finite() && p.is_finite() && p_ok) {
                return Err(LatticeError::InvalidCamber);
            }
        }
        let refined = match spacing {
            SpanSpacing::Uniform => planform.refine(n_span + 1),
            SpanSpacing::Cosine => planform.refine_cosine(n_span + 1),
        };
        let sections = refined.sections();

        let mut elements = Vec::with_capacity(n_span * n_chord);
        for i in 0..n_span {
            let (a, b) = (&sections[i], &sections[i + 1]);
            let strip_dy = b.le.y - a.le.y;
            for k in 0..n_chord {
                // Cosine chordwise stations across the section chord (LE→TE).
                let (t_le, t_te) = (chordwise(k, n_chord), chordwise(k + 1, n_chord));
                let dt = t_te - t_le;
                let t_qb = t_le + 0.25 * dt; // panel quarter chord
                let t_cp = t_le + 0.75 * dt; // collocation station

                // Spanwise-mid LE point, chord, twist.
                let le_mid = Vec3::new(
                    a.le.x + 0.5 * (b.le.x - a.le.x),
                    a.le.y + 0.5 * (b.le.y - a.le.y),
                    a.le.z + 0.5 * (b.le.z - a.le.z),
                );
                let chord_mid = 0.5 * (a.chord + b.chord);
                let twist = a.twist + 0.5 * (b.twist - a.twist);

                // Bound vortex at the panel's quarter chord, spanning the
                // strip in y; x follows each end section's own geometry.
                let x_qb_root = a.le.x + t_qb * a.chord;
                let x_qb_tip = b.le.x + t_qb * b.chord;
                let bound_a = Vec3::new(x_qb_root, a.le.y, a.le.z);
                let bound_b = Vec3::new(x_qb_tip, b.le.y, b.le.z);

                // Collocation point at the panel's three-quarter chord,
                // mid-span, on the camber line.
                let collocation = Vec3::new(
                    le_mid.x + t_cp * chord_mid,
                    le_mid.y,
                    le_mid.z + camber.z_at(t_cp) * chord_mid,
                );

                // Flat normal: from the strip's dihedral tilt only.
                let flat_normal = Vec3::new(0.0, -dz_tilt(a, b), 1.0).normalized();
                let span_axis = (bound_b - bound_a).normalized();
                let panel_chord = dt * chord_mid;
                elements.push(LatticeElement {
                    bound_a,
                    bound_b,
                    collocation,
                    flat_normal,
                    span_axis,
                    camber_slope: camber.slope_at(t_cp),
                    chord: panel_chord,
                    area: panel_chord * strip_dy,
                    strip: i,
                    x_te: le_mid.x + t_te * chord_mid,
                    strip_dy,
                    strip_chord: chord_mid,
                    twist,
                    surface: 0,
                    eta: le_mid.y,
                    chord_frac: t_cp,
                });
            }
        }
        Ok(Self {
            elements,
            n_span,
            n_chord,
        })
    }

    /// A copy of this lattice rolled about the body *x* axis by `roll_rad` and
    /// translated by `offset` (see [`LatticeElement::transformed`]).
    #[must_use]
    pub fn transformed(&self, roll_rad: f64, offset: Vec3<f64>) -> Self {
        Self {
            elements: self
                .elements
                .iter()
                .map(|e| e.transformed(roll_rad, offset))
                .collect(),
            n_span: self.n_span,
            n_chord: self.n_chord,
        }
    }

    /// Merge several lattices into one model of a multi-surface configuration.
    ///
    /// Strip indices are offset so every surface keeps its own strips, and
    /// `n_span` becomes the total strip count. `n_chord` is informational
    /// only (taken from the first part).
    ///
    /// # Errors
    /// [`LatticeError::InvalidResolution`] if `parts` is empty.
    pub fn combine(parts: &[Lattice]) -> Result<Self, LatticeError> {
        let first = parts.first().ok_or(LatticeError::InvalidResolution)?;
        let mut elements = Vec::with_capacity(parts.iter().map(Lattice::len).sum());
        let mut strips = 0;
        let mut surfaces = 0;
        for part in parts {
            elements.extend(part.elements.iter().cloned().map(|mut e| {
                e.strip += strips;
                e.surface += surfaces;
                e
            }));
            strips += part.n_span;
            surfaces += part
                .elements
                .iter()
                .map(|e| e.surface)
                .max()
                .map_or(0, |m| m + 1);
        }
        Ok(Self {
            elements,
            n_span: strips,
            n_chord: first.n_chord,
        })
    }

    /// All elements, spanwise strips outer loop, chordwise inner.
    #[must_use]
    pub fn elements(&self) -> &[LatticeElement] {
        &self.elements
    }

    /// Number of spanwise strips.
    #[must_use]
    pub const fn n_span(&self) -> usize {
        self.n_span
    }

    /// Number of chordwise panels per strip.
    #[must_use]
    pub const fn n_chord(&self) -> usize {
        self.n_chord
    }

    /// Total number of horseshoe elements.
    #[must_use]
    pub fn len(&self) -> usize {
        self.elements.len()
    }

    /// True if empty (cannot happen via `build`).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// Distinct strip count helper for downstream reduction.
    #[must_use]
    pub fn strip_count(&self) -> usize {
        self.n_span
    }
}

fn chordwise(k: usize, n: usize) -> f64 {
    // Cosine spacing: dense at LE (k=0) and TE (k=n).
    let theta = core::f64::consts::PI * k as f64 / n as f64;
    0.5 * (1.0 - theta.cos())
}

fn dz_tilt(a: &tpt_aero_core::WingSection, b: &tpt_aero_core::WingSection) -> f64 {
    let dy = b.le.y - a.le.y;
    if dy.abs() < 1.0e-12 {
        0.0
    } else {
        (b.le.z - a.le.z) / dy
    }
}
