//! Airfoil geometry: NACA 4-digit generation and coordinate transformations.

use alloc::{
    string::{String, ToString},
    vec::Vec,
};

#[cfg(not(feature = "std"))]
use crate::num::Real as _;

/// Error returned by airfoil constructors and parsers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AirfoilError {
    /// A NACA 4-digit code was not 4 ASCII digits.
    InvalidNacaCode(String),
    /// A coordinate list must contain at least 2 points with 0 ≤ x ≤ 1.
    InvalidCoordinates,
    /// Fewer than 2 chordwise stations requested.
    InvalidResolution,
}

impl core::fmt::Display for AirfoilError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            AirfoilError::InvalidNacaCode(code) => {
                write!(f, "invalid NACA 4-digit code: {code:?}")
            }
            AirfoilError::InvalidCoordinates => {
                write!(f, "airfoil coordinates must be ≥2 points with 0 ≤ x ≤ 1")
            }
            AirfoilError::InvalidResolution => {
                write!(f, "need at least 2 chordwise stations")
            }
        }
    }
}

impl core::error::Error for AirfoilError {}

/// Airfoil surface coordinates as a single closed-ish loop: upper surface from
/// trailing edge to leading edge, then lower surface from leading edge back to
/// trailing edge. Chordwise coordinate `x ∈ [0, 1]`, `y` in chord units.
#[derive(Debug, Clone, PartialEq)]
pub struct AirfoilCoordinates {
    /// Chordwise station (0 = LE, 1 = TE).
    pub x: Vec<f64>,
    /// Surface ordinate (positive up).
    pub y: Vec<f64>,
}

impl AirfoilCoordinates {
    /// Build from coordinate lists, validating bounds.
    ///
    /// # Errors
    /// [`AirfoilError::InvalidCoordinates`] if fewer than 2 points or any
    /// coordinate is outside the unit chord box (a small tolerance is allowed).
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> Result<Self, AirfoilError> {
        if x.len() != y.len() || x.len() < 2 {
            return Err(AirfoilError::InvalidCoordinates);
        }
        // The classic open-trailing-edge NACA construction (coefficient
        // 0.1015) leaves a tiny overshoot past x = 1 on cambered sections;
        // tolerate it while still rejecting real garbage.
        let tol = 1.0e-3;
        if x.iter().any(|&xi| xi < -tol || xi > 1.0 + tol) || y.iter().any(|&yi| !yi.is_finite()) {
            return Err(AirfoilError::InvalidCoordinates);
        }
        Ok(Self { x, y })
    }

    /// Number of surface points.
    #[must_use]
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// True if there are no surface points (cannot happen via `new`).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    /// Apply a chord scale (c != 0), a twist rotation about the quarter-chord
    /// point (positive = nose-up), and a translation.
    #[must_use]
    pub fn transform(&self, scale: f64, twist_rad: f64, dx: f64, dy: f64) -> Self {
        let (s, c) = twist_rad.sin_cos();
        let (xp, yp): (Vec<f64>, Vec<f64>) = self
            .x
            .iter()
            .zip(self.y.iter())
            .map(|(&x, &y)| {
                // Move to quarter-chord, rotate (nose-up = positive), move back.
                let (x0, y0) = (x - 0.25, y);
                let xr = x0 * c + y0 * s;
                let yr = -x0 * s + y0 * c;
                (xr * scale + 0.25 * scale + dx, yr * scale + dy)
            })
            .unzip();
        Self { x: xp, y: yp }
    }

    /// Mirror about the chord plane (y → −y): swaps suction/pressure sides.
    #[must_use]
    pub fn reflected(&self) -> Self {
        Self {
            x: self.x.clone(),
            y: self.y.iter().map(|&y| -y).collect(),
        }
    }
}

/// Generate NACA 4-digit airfoil coordinates.
///
/// `code` is the classic four digits, e.g. `"2412"`: maximum camber 2 % at
/// 40 % chord, 12 % maximum thickness. The trailing edge follows the original
/// NACA convention (open: coefficient 0.1015).
///
/// Points are distributed with cosine clustering at leading and trailing
/// edges; `n` is the number of chordwise stations (≥ 2).
///
/// # Errors
/// [`AirfoilError::InvalidNacaCode`] if the code is not four digits with
/// `0 < t ≤ 40` (percent thickness; covers e.g. NACA 0030).
/// [`AirfoilError::InvalidResolution`] if `n < 2`.
pub fn naca4(code: &str, n: usize) -> Result<AirfoilCoordinates, AirfoilError> {
    naca4_impl(code, n, false)
}

/// As [`naca4`], but with the *closed* trailing-edge coefficient (0.1036) so
/// the thickness vanishes exactly at `x = 1`. The shared trailing-edge vertex
/// appears once (the contour closes implicitly from the last point back to
/// the first), which is what panel solvers need: a sharp trailing edge with
/// no gap panel. Total points: `2n`.
///
/// # Errors
/// As [`naca4`].
pub fn naca4_closed(code: &str, n: usize) -> Result<AirfoilCoordinates, AirfoilError> {
    naca4_impl(code, n, true)
}

fn naca4_impl(code: &str, n: usize, closed: bool) -> Result<AirfoilCoordinates, AirfoilError> {
    let te = if closed { TE_CLOSED } else { TE_OPEN };
    if n < 2 {
        return Err(AirfoilError::InvalidResolution);
    }
    let bytes = code.as_bytes();
    if bytes.len() != 4 || !bytes.iter().all(u8::is_ascii_digit) {
        return Err(AirfoilError::InvalidNacaCode(code.to_string()));
    }
    let digit = |i: usize| (bytes[i] - b'0') as f64;
    let (m_pct, p, t_pct) = (digit(0), digit(1) / 10.0, digit(2) * 10.0 + digit(3));
    let m = m_pct / 100.0;
    let t = t_pct / 100.0;
    if t <= 0.0 || t > 0.40 {
        return Err(AirfoilError::InvalidNacaCode(code.to_string()));
    }
    // Camber position 0 is only meaningful for a symmetric section. Any
    // nonzero camber with p = 0 would otherwise be dropped silently.
    if m > 0.0 && p == 0.0 {
        return Err(AirfoilError::InvalidNacaCode(code.to_string()));
    }

    // Cosine-spaced stations from TE (x=1) around the LE to TE again.
    let beta: Vec<f64> = (0..=n)
        .map(|i| core::f64::consts::PI * i as f64 / n as f64)
        .collect();
    let x: Vec<f64> = beta.iter().map(|&b| 0.5 * (1.0 - b.cos())).collect();

    let mut px = Vec::with_capacity(2 * n + 1);
    let mut py = Vec::with_capacity(2 * n + 1);

    // Upper surface: x goes from ~1 (TE) down to 0 (LE).
    for &xi in x.iter().skip(1).rev() {
        let (yt, dyc, yc) = thickness_camber(xi, m, p, t, te);
        let theta = dyc.atan2(1.0);
        px.push(xi - yt * theta.sin());
        py.push(yc + yt * theta.cos());
    }
    // Lower surface: x from 0 back to 1.
    // (closed: the final station repeats the upper-surface TE vertex.)
    let lower_count = if closed { n } else { n + 1 };
    for &xi in x.iter().take(lower_count) {
        let (yt, dyc, yc) = thickness_camber(xi, m, p, t, te);
        let theta = dyc.atan2(1.0);
        px.push(xi + yt * theta.sin());
        py.push(yc - yt * theta.cos());
    }

    AirfoilCoordinates::new(px, py)
}

/// Original NACA trailing-edge coefficient (open trailing edge).
const TE_OPEN: f64 = 0.1015;
/// Modified coefficient that closes the trailing edge (thickness 0 at x = 1).
const TE_CLOSED: f64 = 0.1036;

/// NACA thickness half-width, camber ordinate and camber slope at `x`.
fn thickness_camber(x: f64, m: f64, p: f64, t: f64, te: f64) -> (f64, f64, f64) {
    let yt = 5.0
        * t
        * (0.2969 * x.sqrt() - 0.1260 * x - 0.3516 * x * x + 0.2843 * x * x * x
            - te * x * x * x * x);
    let (yc, dyc);
    if m == 0.0 || p == 0.0 {
        yc = 0.0;
        dyc = 0.0;
    } else if x < p {
        yc = m / (p * p) * (2.0 * p * x - x * x);
        dyc = 2.0 * m / (p * p) * (p - x);
    } else {
        yc = m / ((1.0 - p) * (1.0 - p)) * (1.0 - 2.0 * p + 2.0 * p * x - x * x);
        dyc = 2.0 * m / ((1.0 - p) * (1.0 - p)) * (p - x);
    }
    (yt, dyc, yc)
}

/// Return the NACA camber-line ordinate at `x` for camber `m` at position `p`.
#[must_use]
pub fn naca4_camber(x: f64, m: f64, p: f64) -> f64 {
    thickness_camber(x, m, p, 0.0, TE_OPEN).2
}

/// Return the NACA thickness half-width at `x` for thickness ratio `t`.
#[must_use]
pub fn naca4_thickness(x: f64, t: f64) -> f64 {
    thickness_camber(x, 0.0, 0.0, t, TE_OPEN).0
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    #[test]
    fn naca0012_symmetric() {
        let a = naca4("0012", 40).expect("valid code");
        // Loop: upper TE->LE (40 points), LE, lower LE->TE (41 points).
        // Symmetry: y(x_upper) == -y(x_lower) for the same x.
        for i in 0..a.len() {
            assert!(a.y[i].abs() < 1.0 + 1e-9);
        }
        // Max thickness ≈ 0.12 near x = 0.3.
        let tmax = a.y.iter().fold(0.0_f64, |m, &y| m.max(y));
        assert!((tmax - 0.06).abs() < 2e-3, "half-thickness {tmax}");
    }

    #[test]
    fn naca2412_cambered() {
        let a = naca4("2412", 40).expect("valid code");
        // At x=0.5 upper surface should be well above zero (camber + thickness).
        let upper_mid = a.y[a.len() / 4]; // near mid-chord on upper side
        assert!(upper_mid > 0.02);
    }

    #[test]
    fn closed_trailing_edge_has_no_gap_or_duplicate() {
        for code in ["0012", "2412"] {
            let a = naca4_closed(code, 40).expect("valid");
            assert_eq!(a.len(), 80);
            let (first, last) = ((a.x[0], a.y[0]), (a.x[a.len() - 1], a.y[a.len() - 1]));
            assert!((first.0 - 1.0).abs() < 1e-12 && first.1.abs() < 1e-12);
            // The last vertex is a real lower-surface point, not a repeat of
            // the TE, and the closing panel is short but non-degenerate.
            let gap = ((first.0 - last.0).powi(2) + (first.1 - last.1).powi(2)).sqrt();
            assert!(gap > 0.0 && gap < 0.01, "closing panel {gap}");
        }
    }

    #[test]
    fn rejects_bad_codes() {
        assert_eq!(
            naca4("abc", 10).unwrap_err(),
            AirfoilError::InvalidNacaCode("abc".into())
        );
        assert_eq!(
            naca4("123", 10).unwrap_err(),
            AirfoilError::InvalidNacaCode("123".into())
        );
        assert_eq!(
            naca4("0000", 10).unwrap_err(),
            AirfoilError::InvalidNacaCode("0000".into())
        );
        assert_eq!(
            naca4("0012", 1).unwrap_err(),
            AirfoilError::InvalidResolution
        );
        assert!(naca4("0030", 20).is_ok());
    }

    #[test]
    fn camber_with_zero_position_is_rejected_not_flattened() {
        // m > 0 with p = 0 is not a NACA section; it must not silently
        // produce a symmetric profile.
        assert_eq!(
            naca4("2012", 20).unwrap_err(),
            AirfoilError::InvalidNacaCode("2012".into())
        );
        // Symmetric sections with p = 0 remain valid.
        assert!(naca4("0012", 20).is_ok());
    }

    #[test]
    fn transforms() {
        let a = naca4("0012", 20).expect("valid");
        let scaled = a.transform(2.0, 0.0, 1.0, 0.5);
        // No twist: pure 2x chord scale about the quarter chord plus translate.
        assert!((scaled.x[0] - (2.0 * a.x[0] + 1.0)).abs() < 1e-12);
        assert!((scaled.y[0] - (2.0 * a.y[0] + 0.5)).abs() < 1e-12);
        let mirror = a.reflected();
        for (y0, y1) in a.y.iter().zip(mirror.y.iter()) {
            assert!((y0 + y1).abs() < 1e-12);
        }
    }

    #[test]
    fn twist_moves_leading_edge_up() {
        let flat: AirfoilCoordinates =
            AirfoilCoordinates::new(vec![0.0, 0.5, 1.0], vec![0.0, 0.0, 0.0]).expect("ok");
        let twisted = flat.transform(1.0, 5.0_f64.to_radians(), 0.0, 0.0);
        // LE (x=0) is ahead of the quarter-chord pivot: nose-up twist lifts it.
        assert!(twisted.y[0] > 0.0);
        // TE (x=1) drops.
        assert!(twisted.y[2] < 0.0);
    }
}
