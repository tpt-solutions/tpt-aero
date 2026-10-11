//! # tpt-aero-verify
//!
//! Reusable [`proptest`](mod@proptest) strategies for generating valid aerodynamic inputs,
//! plus the cross-crate invariant checks from `spec.txt` §2/§4:
//!
//! * symmetric geometry at zero AoA ⇒ zero rolling and yawing moments,
//! * reversing the freestream velocity vector reverses the lift vector,
//! * (Phase 2+) BEMT induction bounds and actuator-disk mass conservation.
//!
//! This crate is std-only (proptest needs the OS RNG); solver crates remain
//! `no_std`.
//!
//! # Example
//!
//! ```
//! use proptest::strategy::{Strategy, ValueTree};
//! use proptest::test_runner::TestRunner;
//!
//! let mut runner = TestRunner::deterministic();
//! let params = tpt_aero_verify::symmetric_wing()
//!     .new_tree(&mut runner)
//!     .unwrap()
//!     .current();
//! let wing = params.planform();
//! assert!(wing.span() > 0.0);
//! ```

#![forbid(unsafe_code)]

#[cfg(kani)]
mod kani_proofs;
#[cfg(not(kani))]
mod kani_proofs;

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
    (
        0.5_f64..30.0,
        0.1_f64..5.0,
        0.3_f64..1.0,
        0.0_f64..1.0,
        -0.35_f64..0.35,
    )
        .prop_map(|(span, root_chord, taper, sweep_frac, dihedral)| {
            let max_sweep = span * 0.25;
            SymmetricWingParams {
                span,
                root_chord,
                taper,
                sweep: sweep_frac * max_sweep,
                dihedral,
            }
        })
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

/// Flight condition for an aircraft trim or polar check: `(speed m/s, geopotential
/// altitude m)`. Speed 10 .. 80 m/s, altitude 0 .. 11 000 m (the troposphere).
pub fn flight_condition() -> impl Strategy<Value = (f64, f64)> {
    (10.0_f64..80.0, 0.0_f64..11_000.0)
}

/// Control-surface deflection in radians, ±20°. The trim solvers refuse larger
/// deflections than the linear VLM range supports.
pub fn control_deflection() -> impl Strategy<Value = f64> {
    (-20.0_f64..20.0).prop_map(|deg| deg.to_radians())
}

/// A valid section polar: a linear lift slope in `[4, 7]` per radian, a
/// parabolic drag, and a constant pitching moment, tabulated on −10°..15°.
///
/// Returns the slope and the built table so tests can compare the two.
pub fn section_polar() -> impl Strategy<Value = (f64, tpt_aero_core::AirfoilPolar)> {
    (4.0_f64..7.0, 0.005_f64..0.03).prop_map(|(slope, cd0)| {
        let alphas: Vec<f64> = (-10..=15).map(|d| f64::from(d).to_radians()).collect();
        let cl: Vec<f64> = alphas.iter().map(|&a| slope * a).collect();
        let cd: Vec<f64> = alphas.iter().map(|&a| cd0 + 0.05 * a * a).collect();
        let cm = vec![-0.05; alphas.len()];
        let polar = tpt_aero_core::AirfoilPolar::new(alphas, cl, cd, cm, 0.0, 1.0e6)
            .expect("generated polar is valid by construction");
        (slope, polar)
    })
}
