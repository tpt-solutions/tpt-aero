//! Multi-surface configurations: wing + tail + fin as one lattice.
//!
//! Each [`Surface`] is a planform placed in body axes (*x* aft, *y* right,
//! *z* up). Horizontal surfaces span *y*; a vertical fin is the same kind of
//! planform rolled 90° about *x*, so its span runs along *z* (build it
//! one-sided, with section `y` from 0 to the fin height). All surfaces shed
//! their trailing vortices along body +x, so downwash and sidewash
//! interference between them is captured by the shared AIC matrix.

use alloc::vec::Vec;

use tpt_aero_core::{Vec3, WingPlanform, WingSection};

use crate::lattice::{CamberModel, Lattice, LatticeError};
use crate::solver::{ReferenceGeometry, VlmError, VlmModel};

/// One lifting surface of a configuration.
#[derive(Debug, Clone)]
pub struct Surface {
    /// Planform in its own frame (span along *y*, zero dihedral roll).
    pub planform: WingPlanform,
    /// Spanwise strips.
    pub n_span: usize,
    /// Chordwise panels.
    pub n_chord: usize,
    /// Camber line.
    pub camber: CamberModel,
    /// Position of the planform origin (root leading edge) in body axes (m).
    pub offset: Vec3<f64>,
    /// Roll about body *x* (rad): 0 for a horizontal surface, `π/2` for a
    /// fin pointing up.
    pub roll_rad: f64,
}

impl Surface {
    /// A horizontal surface (wing, tailplane, canard) with the origin at the body origin.
    #[must_use]
    pub fn horizontal(planform: WingPlanform, n_span: usize, n_chord: usize) -> Self {
        Self {
            planform,
            n_span,
            n_chord,
            camber: CamberModel::None,
            offset: Vec3::zero(),
            roll_rad: 0.0,
        }
    }

    /// A fin pointing up: the planform's span (*y*, from 0 to the fin height)
    /// becomes body *z*.
    #[must_use]
    pub fn vertical(planform: WingPlanform, n_span: usize, n_chord: usize) -> Self {
        Self {
            roll_rad: core::f64::consts::FRAC_PI_2,
            ..Self::horizontal(planform, n_span, n_chord)
        }
    }

    /// Place the surface's root leading edge at `(x, y, z)` in body axes.
    #[must_use]
    pub fn at(mut self, x: f64, y: f64, z: f64) -> Self {
        self.offset = Vec3::new(x, y, z);
        self
    }

    /// Use a camber line.
    #[must_use]
    pub fn with_camber(mut self, camber: CamberModel) -> Self {
        self.camber = camber;
        self
    }

    /// Set the incidence (rad, nose-up positive) by adding it to every
    /// section's twist; a horizontal tail with negative incidence carries
    /// download.
    #[must_use]
    pub fn with_incidence(mut self, incidence_rad: f64) -> Self {
        let sections: Vec<WingSection> = self
            .planform
            .sections()
            .iter()
            .map(|s| WingSection {
                twist: s.twist + incidence_rad,
                ..*s
            })
            .collect();
        // The sections came from a valid planform and keep their y order.
        if let Ok(p) = WingPlanform::from_sections(sections) {
            self.planform = p;
        }
        self
    }
}

impl VlmModel {
    /// Build a multi-surface model. The reference geometry (area, chord,
    /// span, moment centre at the quarter chord of the mean chord) comes from
    /// the **first** surface — conventionally the wing — and can be replaced
    /// with [`VlmModel::with_reference`], e.g. to move the centre of gravity.
    ///
    /// # Errors
    /// [`VlmError::Lattice`] if `surfaces` is empty or a surface has an
    /// invalid planform or resolution.
    pub fn from_surfaces(surfaces: &[Surface]) -> Result<Self, VlmError> {
        let first = surfaces
            .first()
            .ok_or(VlmError::Lattice(LatticeError::InvalidResolution))?;
        let parts: Vec<Lattice> = surfaces
            .iter()
            .map(|s| {
                Lattice::build(&s.planform, s.n_span, s.n_chord, s.camber)
                    .map(|l| l.transformed(s.roll_rad, s.offset))
            })
            .collect::<Result<_, _>>()
            .map_err(VlmError::Lattice)?;
        let lattice = Lattice::combine(&parts).map_err(VlmError::Lattice)?;

        let mut reference = ReferenceGeometry::from_planform(&first.planform);
        reference.moment_center = reference.moment_center + first.offset;
        Ok(VlmModel::from_lattice(lattice).with_reference(reference))
    }
}
