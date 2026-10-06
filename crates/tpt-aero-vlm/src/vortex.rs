//! Biot–Savart velocity induction for vortex elements.
//!
//! All functions are generic over [`Real`], so influence coefficients can be
//! assembled on `f64` or on dual numbers for exact first derivatives.

use tpt_aero_core::{Real, Vec3};

/// Four-π, kept as a helper so the formulas below read like the textbook.
#[inline]
fn four_pi<R: Real>() -> R {
    R::from_f64(4.0) * R::pi()
}

/// Velocity induced at `p` by a finite straight vortex segment `a → b` with
/// circulation `g` (Katz & Plotkin eq. 11.98).
///
/// Returns zero when `p` lies on the segment axis (the singular self-induction
/// is excluded by construction in the solver).
pub fn segment_velocity<R: Real>(a: &Vec3<R>, b: &Vec3<R>, p: &Vec3<R>, g: R) -> Vec3<R> {
    let r0 = *b - *a; // a → b
    let r1 = *p - *a;
    let r2 = *p - *b;
    let cross = r1.cross(&r2);
    let cross_sq = cross.dot(&cross);
    let d1 = r1.norm();
    let d2 = r2.norm();
    let eps = R::from_f64(1.0e-12);
    if cross_sq.le_f(eps * eps) || d1.le_f(eps) || d2.le_f(eps) {
        return Vec3::zero();
    }
    let coeff = (r0.dot(&(
        r1.scaled(d1.recip()) - r2.scaled(d2.recip())
    ))) / cross_sq * (g / four_pi::<R>());
    cross.scaled(coeff)
}

/// Velocity induced at `p` by a semi-infinite vortex starting at `a` and
/// extending to infinity in the unit direction `d`, circulation `g`
/// (Katz & Plotkin eq. 11.99 form).
pub fn semi_infinite_velocity<R: Real>(a: &Vec3<R>, d: &Vec3<R>, p: &Vec3<R>, g: R) -> Vec3<R> {
    let r = *p - *a;
    // Direction d × r (not r × d): for a +x vortex this yields +z velocity at
    // a +y offset, matching the right-hand rule (and the bound-segment
    // formula's sign convention above).
    let rx = d.cross(&r);
    let rx_sq = rx.dot(&rx);
    let rn = r.norm();
    let eps = R::from_f64(1.0e-12);
    if rx_sq.le_f(eps * eps) || rn.le_f(eps) {
        return Vec3::zero();
    }
    let cos_theta = r.dot(d) / rn;
    let coeff = (R::one() + cos_theta) / rx_sq * (g / four_pi::<R>());
    rx.scaled(coeff)
}

/// Velocity induced at `p` by a horseshoe vortex: bound segment `a → b` with
/// circulation `g`, trailing legs from `a` and `b` extending to infinity along
/// the unit freestream direction `d` (closed loop, opposite circulation on the
/// return leg).
pub fn horseshoe_velocity<R: Real>(
    a: &Vec3<R>,
    b: &Vec3<R>,
    d: &Vec3<R>,
    p: &Vec3<R>,
    g: R,
) -> Vec3<R> {
    let bound = segment_velocity(a, b, p, g);
    let leg_b = semi_infinite_velocity(b, d, p, g);
    let leg_a = semi_infinite_velocity(a, d, p, g.neg());
    bound + leg_b + leg_a
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(v: Vec3<f64>, x: f64, y: f64, z: f64, tol: f64) {
        assert!((v.x - x).abs() < tol, "x: {} vs {v:?}", v.x);
        assert!((v.y - y).abs() < tol, "y: {} vs {v:?}", v.y);
        assert!((v.z - z).abs() < tol, "z: {} vs {v:?}", v.z);
    }

    #[test]
    fn segment_velocity_direction() {
        // A straight vortex along +y through the origin; induced velocity at
        // a point on +x is −z (right-hand rule: Γ along +y, field circles it).
        let a = Vec3::new(0.0, -1.0, 0.0);
        let b = Vec3::new(0.0, 1.0, 0.0);
        let p = Vec3::new(1.0, 0.0, 0.0);
        let v = segment_velocity(&a, &b, &p, 1.0);
        // Exact finite-segment result: |v| = Γ/(4π d) (cosθ1 + cosθ2),
        // θ measured from the perpendicular; here cosθ1 = cosθ2 = 1/√2.
        let expected = 1.0 / (4.0 * core::f64::consts::PI) * core::f64::consts::SQRT_2;
        approx(v, 0.0, 0.0, -expected, 1e-12);
    }

    #[test]
    fn semi_infinite_matches_perpendicular_limit() {
        // Semi-infinite vortex along +x starting at origin, point on +y.
        let a = Vec3::zero();
        let d = Vec3::new(1.0, 0.0, 0.0);
        let p = Vec3::new(0.0, 1.0, 0.0);
        let v = semi_infinite_velocity(&a, &d, &p, 1.0);
        // |v| = Γ/(4π h)(1 + cosθ), θ = 90° here → Γ/(4π h); a +x circulation
        // gives +z velocity at a +y offset (right-hand rule), i.e. UPWASH
        // outboard of the right trailing leg.
        let expected = 1.0 / (4.0 * core::f64::consts::PI);
        approx(v, 0.0, 0.0, expected, 1e-12);
    }

    #[test]
    fn horseshoe_symmetry_on_midplane() {
        // On the mid-span plane the two trailing legs cancel in y and z.
        let a = Vec3::new(0.0, -1.0, 0.0);
        let b = Vec3::new(0.0, 1.0, 0.0);
        let d = Vec3::new(1.0, 0.0, 0.0);
        let p = Vec3::new(0.5, 0.0, 0.0);
        let v = horseshoe_velocity(&a, &b, &d, &p, 1.0);
        assert!(v.y.abs() < 1e-12);
        assert!(v.z < 0.0); // downwash between the legs
    }

    #[test]
    fn degenerate_points_return_zero() {
        let a = Vec3::zero();
        let b = Vec3::new(1.0, 0.0, 0.0);
        let d = Vec3::new(1.0, 0.0, 0.0);
        // p on the segment → cross ≈ 0 → zero.
        let p = Vec3::new(0.5, 0.0, 0.0);
        let v = segment_velocity(&a, &b, &p, 1.0);
        assert_eq!(v, Vec3::zero());
        // p on the leg line.
        let q = Vec3::new(-1.0, 0.0, 0.0);
        let v = semi_infinite_velocity(&a, &d, &q, 1.0);
        assert_eq!(v, Vec3::zero());
    }
}
