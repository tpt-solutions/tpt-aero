//! # tpt-aero-verify
//!
//! Reusable [`proptest`] strategies for generating valid aerodynamic inputs,
//! plus the cross-crate invariant checks from `spec.txt` §2/§4:
//!
//! * symmetric geometry at zero AoA ⇒ zero rolling and yawing moments,
//! * reversing the freestream velocity vector reverses the lift vector,
//! * (Phase 2+) BEMT induction bounds and actuator-disk mass conservation.
//!
//! This crate is std-only (proptest needs the OS RNG); solver crates remain
//! `no_std`.

#![forbid(unsafe_code)]

use proptest::prelude::*;

/// A symmetric trapezoidal planform drawn from reasonable GA/UA V ranges.
#[derive(Debug, Clone, PartialEq)]
pub struct SymmetricWingParams {
    /// Span (m), 0.5 .. 30.
    pub span: f64,
    /// Root chord (m), 0.1 .. 5.
    pub root_chord: f64,
    /// Taper ratio tip/root, 0.3 .. 1.0.
    pub taper: f64,
    /// Tip leading-edge aft offset (m), 0 .. span/4.
    pub sweep: f64,
    /// Dihedral angle (rad), −0.35 .. 0.35.
    pub dihedral: f64,
}

/// Strategy for [`SymmetricWingParams`].
pub fn symmetric_wing() -> impl Strategy<Value = SymmetricWingParams> {
    (0.5_f64..30.0, 0.1_f64..5.0, 0.3_f64..1.0, 0.0_f64..1.0, -0.35_f64..0.35).prop_map(
        |(span, root_chord, taper, sweep_frac, dihedral)| {
            let max_sweep = span * 0.25;
            SymmetricWingParams {
                span,
                root_chord,
                taper,
                sweep: sweep_frac * max_sweep,
                dihedral,
            }
        },
    )
}

impl SymmetricWingParams {
    /// Materialize the planform.
    #[must_use]
    pub fn planform(&self) -> tpt_aero_core::WingPlanform {
        tpt_aero_core::WingPlanform::tapered(
            self.span,
            self.root_chord,
            self.root_chord * self.taper,
            self.sweep,
            self.dihedral,
        )
    }
}

/// Angle of attack in the valid VLM range, −15° .. 15° (radians).
pub fn angle_of_attack() -> impl Strategy<Value = f64> {
    (-15.0_f64..15.0).prop_map(|deg| deg.to_radians())
}

/// Freestream speed in m/s, 5 .. 200.
pub fn freestream_speed() -> impl Strategy<Value = f64> {
    5.0_f64..200.0
}

/// VLM discretization: 4..24 strips, 1..5 chordwise panels.
pub fn lattice_resolution() -> impl Strategy<Value = (usize, usize)> {
    (4_usize..24, 1_usize..5)
}
