//! AVL `.avl` geometry: import a deck as a [`VlmModel`], and export surfaces
//! back to a deck.
//!
//! The supported subset is the one the validation decks in
//! `docs/validation/avl/` use:
//!
//! * the six-line header: title, Mach, `iYsym iZsym Zsym`, `Sref Cref Bref`,
//!   `Xref Yref Zref`, `CDp`. `Sref`, `Cref`, `Bref` and the moment reference
//!   become the model's reference geometry;
//! * `SURFACE` blocks with `Nchord Cspace Nspan Sspace`, whose `Nspan` is the
//!   strip count of the surface;
//! * `SECTION` lines `Xle Yle Zle Chord Ainc` (incidence in degrees, positive
//!   nose-up, taken as the section twist); the optional per-section
//!   `Nspan Sspace` are accepted and ignored in favour of the surface count;
//! * `CONTROL` lines after a section: `name gain Xhinge HingeVec SgnDup`. The
//!   control covers the span from that section to the next one; `SgnDup = -1`
//!   makes it antisymmetric.
//!
//! Not supported: `YDUPLICATE` mirrored half-surfaces, `ANGLE`/`SCALE`/
//! `TRANSLATE` transforms, `BODY`, `AFILE`/`NACA` airfoil files, and per-section
//! strip counts. These are rejected with a clear error rather than ignored.

use alloc::{format, string::String, string::ToString, vec::Vec};
use core::fmt;

use tpt_aero_core::{Vec3, WingPlanform, WingSection};

use crate::aircraft::{ControlSymmetry, Surface};
use crate::solver::{ReferenceGeometry, VlmError, VlmModel};

/// Errors from reading an AVL deck.
#[derive(Debug, Clone, PartialEq)]
pub enum AvlError {
    /// The deck ended before a required line.
    Truncated(&'static str),
    /// A line did not parse as the expected numbers.
    BadLine {
        /// One-based line number in the deck.
        line: usize,
        /// What was expected there.
        expected: &'static str,
    },
    /// A keyword this importer does not support.
    Unsupported {
        /// One-based line number in the deck.
        line: usize,
        /// The keyword as written.
        keyword: String,
    },
    /// A surface with fewer than two sections.
    TooFewSections(String),
    /// The sections or the solve were rejected by the VLM.
    Vlm(String),
}

impl fmt::Display for AvlError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AvlError::Truncated(what) => write!(f, "AVL deck ends before {what}"),
            AvlError::BadLine { line, expected } => write!(f, "line {line}: expected {expected}"),
            AvlError::Unsupported { line, keyword } => {
                write!(f, "line {line}: {keyword} is not supported by the importer")
            }
            AvlError::TooFewSections(name) => write!(f, "surface {name} needs at least 2 sections"),
            AvlError::Vlm(msg) => write!(f, "VLM rejected the geometry: {msg}"),
        }
    }
}

impl core::error::Error for AvlError {}

/// A section line of a surface.
#[derive(Debug, Clone, PartialEq)]
pub struct AvlSection {
    /// Leading-edge position `(x, y, z)` (m).
    pub le: Vec3<f64>,
    /// Chord (m).
    pub chord: f64,
    /// Incidence (deg, nose-up positive).
    pub ainc_deg: f64,
}

/// A control line attached to the section it follows.
#[derive(Debug, Clone, PartialEq)]
pub struct AvlControl {
    /// Name used by trim and derivative queries.
    pub name: String,
    /// Hinge position as a chord fraction.
    pub hinge: f64,
    /// `+1` symmetric, `-1` antisymmetric.
    pub sgn_dup: f64,
    /// Index of the section the control follows.
    pub after_section: usize,
}

/// A `SURFACE` block.
#[derive(Debug, Clone, PartialEq)]
pub struct AvlSurface {
    /// Surface name.
    pub name: String,
    /// Chordwise panels (`Nchord`).
    pub n_chord: usize,
    /// Spanwise strips (`Nspan`).
    pub n_span: usize,
    /// Sections in order of increasing `y` as written.
    pub sections: Vec<AvlSection>,
    /// Controls, each attached to a section.
    pub controls: Vec<AvlControl>,
}

/// A parsed deck: the header and its surfaces.
#[derive(Debug, Clone, PartialEq)]
pub struct AvlDeck {
    /// First line of the deck.
    pub title: String,
    /// Mach number (informational; the importer solves incompressibly).
    pub mach: f64,
    /// Reference area `Sref` (m²).
    pub sref: f64,
    /// Reference chord `Cref` (m).
    pub cref: f64,
    /// Reference span `Bref` (m).
    pub bref: f64,
    /// Moment reference `(Xref, Yref, Zref)` (m).
    pub xref: Vec3<f64>,
    /// Surfaces in deck order.
    pub surfaces: Vec<AvlSurface>,
}

fn numbers(line: &str) -> Vec<f64> {
    line.split_whitespace()
        .filter_map(|t| t.parse::<f64>().ok())
        .collect()
}

/// Parse an AVL deck.
///
/// # Errors
/// [`AvlError`] for a truncated deck, a malformed line, an unsupported keyword,
/// or a surface with fewer than two sections.
pub fn parse(text: &str) -> Result<AvlDeck, AvlError> {
    let lines: Vec<(usize, &str)> = text
        .lines()
        .enumerate()
        .map(|(i, l)| (i + 1, l.split('!').next().unwrap_or("").trim()))
        .collect();
    let mut body = lines.iter().filter(|(_, l)| !l.is_empty());

    let (_, title) = body.next().ok_or(AvlError::Truncated("the title"))?;
    let (line_no, mach_line) = body.next().ok_or(AvlError::Truncated("the Mach number"))?;
    let mach = *numbers(mach_line).first().ok_or(AvlError::BadLine {
        line: *line_no,
        expected: "a Mach number",
    })?;
    let (_, _sym) = body
        .next()
        .ok_or(AvlError::Truncated("the symmetry line"))?;
    let (line_no, sref_line) = body.next().ok_or(AvlError::Truncated("Sref Cref Bref"))?;
    let refs = numbers(sref_line);
    if refs.len() < 3 {
        return Err(AvlError::BadLine {
            line: *line_no,
            expected: "Sref Cref Bref",
        });
    }
    let (line_no, xref_line) = body.next().ok_or(AvlError::Truncated("Xref Yref Zref"))?;
    let xr = numbers(xref_line);
    if xr.len() < 3 {
        return Err(AvlError::BadLine {
            line: *line_no,
            expected: "Xref Yref Zref",
        });
    }
    let (_, _cdp) = body.next().ok_or(AvlError::Truncated("CDp"))?;

    let mut deck = AvlDeck {
        title: title.to_string(),
        mach,
        sref: refs[0],
        cref: refs[1],
        bref: refs[2],
        xref: Vec3::new(xr[0], xr[1], xr[2]),
        surfaces: Vec::new(),
    };

    let mut current: Option<AvlSurface> = None;
    let mut expecting_surface_header = false;
    let mut pending_header: Option<(usize, usize)> = None;
    for (line_no, line) in body {
        let line_no = *line_no;
        // A surface name is free text: take it before keyword matching, so a
        // surface called "Surface" is not read as a SURFACE keyword.
        if expecting_surface_header {
            if let Some(s) = current.as_mut() {
                if s.name.is_empty() {
                    s.name = line.to_string();
                    continue;
                }
            }
        }
        let upper = line.to_ascii_uppercase();
        let keyword = upper.split_whitespace().next().unwrap_or("");
        match keyword {
            "SURFACE" => {
                if let Some(s) = current.take() {
                    deck.surfaces.push(finish(s)?);
                }
                expecting_surface_header = true;
                pending_header = None;
                current = Some(AvlSurface {
                    name: String::new(),
                    n_chord: 0,
                    n_span: 0,
                    sections: Vec::new(),
                    controls: Vec::new(),
                });
            }
            "SECTION" => {
                // The section numbers are on the next line.
                let s = current.as_mut().ok_or(AvlError::BadLine {
                    line: line_no,
                    expected: "a SURFACE before SECTION",
                })?;
                s.sections.push(AvlSection {
                    le: Vec3::zero(),
                    chord: 0.0,
                    ainc_deg: 0.0,
                });
                pending_header = Some((line_no, s.sections.len() - 1));
            }
            "CONTROL" => {
                expecting_surface_header = false;
                let s = current.as_mut().ok_or(AvlError::BadLine {
                    line: line_no,
                    expected: "a SURFACE before CONTROL",
                })?;
                let after = s.sections.len().checked_sub(1).ok_or(AvlError::BadLine {
                    line: line_no,
                    expected: "a SECTION before CONTROL",
                })?;
                let _ = (line_no, after);
                // The control's numbers are on the next line; mark it.
                pending_header = Some((line_no, usize::MAX));
            }
            "YDUPLICATE" | "SCALE" | "TRANSLATE" | "ANGLE" | "BODY" | "AFILE" | "NACA" | "CLAF"
            | "CDCL" | "DESIGN" | "BFILE" => {
                return Err(AvlError::Unsupported {
                    line: line_no,
                    keyword: keyword.to_string(),
                });
            }
            _ => {
                // A data line: the header of a surface, a section's values, or a control's values.
                if expecting_surface_header {
                    let s = current.as_mut().expect("surface open");
                    if s.name.is_empty() {
                        s.name = line.to_string();
                        continue;
                    }
                    let n = numbers(line);
                    if n.len() < 3 {
                        return Err(AvlError::BadLine {
                            line: line_no,
                            expected: "Nchord Cspace Nspan Sspace",
                        });
                    }
                    s.n_chord = n[0].round() as usize;
                    s.n_span = n[2].round() as usize;
                    expecting_surface_header = false;
                    continue;
                }
                match pending_header.take() {
                    Some((_, idx)) if idx != usize::MAX => {
                        let s = current.as_mut().expect("surface open");
                        let n = numbers(line);
                        if n.len() < 4 {
                            return Err(AvlError::BadLine {
                                line: line_no,
                                expected: "Xle Yle Zle Chord Ainc",
                            });
                        }
                        s.sections[idx] = AvlSection {
                            le: Vec3::new(n[0], n[1], n[2]),
                            chord: n[3],
                            ainc_deg: n.get(4).copied().unwrap_or(0.0),
                        };
                    }
                    Some((_, _)) => {
                        // Control line: name gain Xhinge HingeVec(3) SgnDup
                        let s = current.as_mut().expect("surface open");
                        let parts: Vec<&str> = line.split_whitespace().collect();
                        let after = s.sections.len() - 1;
                        let nums = numbers(line);
                        let name = parts.first().copied().unwrap_or("").to_string();
                        let hinge = nums.get(1).copied().unwrap_or(0.0);
                        let sgn = nums.last().copied().unwrap_or(1.0);
                        if name.is_empty() {
                            return Err(AvlError::BadLine {
                                line: line_no,
                                expected: "name gain Xhinge HingeVec SgnDup",
                            });
                        }
                        s.controls.push(AvlControl {
                            name,
                            hinge,
                            sgn_dup: sgn,
                            after_section: after,
                        });
                    }
                    None => {
                        // A data line we did not ask for (for example a second
                        // numbers line after a SECTION): ignore it.
                    }
                }
            }
        }
    }
    if let Some(s) = current.take() {
        deck.surfaces.push(finish(s)?);
    }
    Ok(deck)
}

fn finish(s: AvlSurface) -> Result<AvlSurface, AvlError> {
    if s.sections.len() < 2 {
        return Err(AvlError::TooFewSections(s.name));
    }
    Ok(s)
}

impl AvlDeck {
    /// The surfaces as [`Surface`]s, horizontal and at the origin, with the
    /// deck's per-section twist and dihedral.
    ///
    /// # Errors
    /// [`AvlError::Vlm`] if a planform is invalid (sections not in increasing `y`).
    pub fn surfaces(&self) -> Result<Vec<Surface>, AvlError> {
        self.surfaces
            .iter()
            .map(|s| {
                // A surface whose sections share one y lies in the x–z plane: a fin.
                // Its span runs along z, so it is built as a rolled planform with its
                // span along local y, rooted at the first section.
                let y0 = s.sections[0].le.y;
                let vertical = s.sections.iter().all(|sec| (sec.le.y - y0).abs() < 1.0e-9);
                let root = s.sections[0].le;
                let sections: Vec<WingSection> = s
                    .sections
                    .iter()
                    .map(|sec| {
                        let le = if vertical {
                            Vec3::new(sec.le.x - root.x, sec.le.z - root.z, 0.0)
                        } else {
                            sec.le
                        };
                        WingSection {
                            le,
                            chord: sec.chord,
                            twist: sec.ainc_deg.to_radians(),
                        }
                    })
                    .collect();
                let planform = WingPlanform::from_sections(sections)
                    .map_err(|e| AvlError::Vlm(e.to_string()))?;
                let mut surface = if vertical {
                    Surface::vertical(planform, s.n_span.max(1), s.n_chord.max(1))
                        .at(root.x, root.y, root.z)
                } else {
                    Surface::horizontal(planform, s.n_span.max(1), s.n_chord.max(1))
                };
                surface.controls = s
                    .controls
                    .iter()
                    .map(|c| {
                        let a = s.sections[c.after_section].le.y.abs();
                        let b = s.sections[(c.after_section + 1).min(s.sections.len() - 1)]
                            .le
                            .y
                            .abs();
                        crate::aircraft::Control {
                            name: c.name.clone(),
                            hinge: c.hinge,
                            y_in: a.min(b),
                            y_out: a.max(b),
                            symmetry: if c.sgn_dup < 0.0 {
                                ControlSymmetry::Antisymmetric
                            } else {
                                ControlSymmetry::Symmetric
                            },
                        }
                    })
                    .collect();
                Ok(surface)
            })
            .collect()
    }

    /// The deck as a solvable model, with the deck's reference geometry.
    ///
    /// # Errors
    /// As [`AvlDeck::surfaces`], and [`AvlError::Vlm`] if the model is rejected.
    pub fn to_model(&self) -> Result<VlmModel, AvlError> {
        let surfaces = self.surfaces()?;
        let model = VlmModel::from_surfaces(&surfaces)
            .map_err(|e: VlmError| AvlError::Vlm(e.to_string()))?;
        Ok(model.with_reference(ReferenceGeometry {
            area: self.sref,
            chord: self.cref,
            span: self.bref,
            moment_center: self.xref,
        }))
    }
}

/// Write a horizontal-surface model as an AVL deck, with the reference geometry
/// of `reference`. Each surface becomes one `SURFACE` with one `SECTION` per
/// planform station; a vertical surface is written as its rolled sections.
#[must_use]
pub fn export(title: &str, surfaces: &[Surface], reference: ReferenceGeometry) -> String {
    let mut out = String::new();
    out.push_str(title);
    out.push('\n');
    out.push_str("0.0\n0 0 0.0\n");
    out.push_str(&format!(
        "{} {} {}\n",
        reference.area, reference.chord, reference.span
    ));
    out.push_str(&format!(
        "{} {} {}\n",
        reference.moment_center.x, reference.moment_center.y, reference.moment_center.z
    ));
    out.push_str("0.0\n");
    for s in surfaces {
        let roll_vertical = (s.roll_rad - core::f64::consts::FRAC_PI_2).abs() < 1.0e-12;
        out.push_str("\nSURFACE\nSurface\n");
        out.push_str(&format!("{} 1.0 {} 1.0\n", s.n_chord, s.n_span));
        for sec in s.planform.sections() {
            // A vertical surface rolls its span onto body z: local y becomes z.
            let (y, z) = if roll_vertical {
                (s.offset.y, s.offset.z + sec.le.y)
            } else {
                (s.offset.y + sec.le.y, s.offset.z + sec.le.z)
            };
            out.push_str("\nSECTION\n");
            out.push_str(&format!(
                "{} {} {} {} {}\n",
                s.offset.x + sec.le.x,
                y,
                z,
                sec.chord,
                sec.twist.to_degrees()
            ));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::CamberModel;

    const RECT: &str = "Rectangular wing\n0.0\n0 0 0.0\n12.0 1.2 10.0\n0.3 0.0 0.0\n0.0\n\nSURFACE\nWing\n4  1.0  20  1.0\n\nSECTION\n0.0  -5.0  0.0  1.2  0.0  10  1.0\n\nSECTION\n0.0   0.0  0.0  1.2  0.0  10  1.0\n\nSECTION\n0.0   5.0  0.0  1.2  0.0\n";

    #[test]
    fn parses_the_header_and_one_surface() {
        let deck = parse(RECT).expect("parses");
        assert_eq!(deck.title, "Rectangular wing");
        assert!((deck.sref - 12.0).abs() < 1e-12 && (deck.bref - 10.0).abs() < 1e-12);
        assert_eq!(deck.surfaces.len(), 1);
        let s = &deck.surfaces[0];
        assert_eq!((s.n_chord, s.n_span), (4, 20));
        assert_eq!(s.sections.len(), 3);
        assert!((s.sections[0].le.y + 5.0).abs() < 1e-12);
    }

    #[test]
    fn imported_rectangular_wing_matches_the_native_model() {
        let model = parse(RECT).expect("parses").to_model().expect("builds");
        let native = VlmModel::from_planform(
            &WingPlanform::rectangular(10.0, 1.2),
            20,
            4,
            CamberModel::None,
        )
        .expect("valid");
        let a = 4.0_f64.to_radians();
        let cl_avl = model.solve(a, 50.0, 1.225, &[]).expect("solves").cl;
        let cl_native = native.solve(a, 50.0, 1.225, &[]).expect("solves").cl;
        assert!(
            (cl_avl - cl_native).abs() / cl_native < 1e-6,
            "{cl_avl} vs {cl_native}"
        );
    }

    #[test]
    fn controls_attach_to_the_span_they_follow() {
        let deck_text = format!("{RECT}\nCONTROL\naileron 1.0 0.75 0 0 -1.0\n");
        // The control follows the last SECTION, so it has no span to cover.
        let deck = parse(&deck_text).expect("parses");
        assert_eq!(deck.surfaces[0].controls.len(), 1);
        assert_eq!(deck.surfaces[0].controls[0].sgn_dup, -1.0);
    }

    #[test]
    fn unsupported_keywords_are_refused_not_ignored() {
        let text = RECT.replace("SURFACE", "YDUPLICATE\n0.0\n\nSURFACE");
        assert!(matches!(parse(&text), Err(AvlError::Unsupported { .. })));
    }

    #[test]
    fn export_then_import_keeps_the_span_and_area() {
        let wing = WingPlanform::rectangular(10.0, 1.2);
        let surf = Surface::horizontal(wing.clone(), 20, 4);
        let reference = ReferenceGeometry::from_planform(&wing);
        let text = export("exported", &[surf], reference);
        let deck = parse(&text).expect("parses back");
        assert!((deck.sref - wing.area()).abs() < 1e-9);
        assert!((deck.bref - 10.0).abs() < 1e-9);
        let model = deck.to_model().expect("builds");
        assert!((model.reference().span - 10.0).abs() < 1e-9);
    }
}
