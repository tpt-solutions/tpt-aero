//! Property checks over the strategies in `tpt-aero-verify`: polar lookups stay
//! within the table's own range, and the standard atmosphere thins with altitude.

use proptest::prelude::*;

proptest! {
    #[test]
    fn polar_lookup_stays_within_table_extremes((_slope, polar) in tpt_aero_verify::section_polar(), deg in -10.0_f64..15.0) {
        let cl = polar.cl_at(deg.to_radians());
        let (lo, hi) = polar.alphas().iter().fold((f64::MAX, f64::MIN), |(lo, hi), &a| {
            let c = polar.cl_at(a);
            (lo.min(c), hi.max(c))
        });
        prop_assert!(cl >= lo - 1e-12 && cl <= hi + 1e-12, "cl {cl} outside [{lo}, {hi}]");
    }

    #[test]
    fn polar_lookup_matches_the_linear_slope_inside_the_grid((slope, polar) in tpt_aero_verify::section_polar(), deg in -9.9_f64..14.9) {
        let a = deg.to_radians();
        // Piecewise-linear interpolation of a linear table is exact to rounding.
        prop_assert!((polar.cl_at(a) - slope * a).abs() < 1e-9);
    }

    #[test]
    fn density_falls_with_altitude(a in 0.0_f64..10_000.0, step in 1.0_f64..1_000.0) {
        let lo = tpt_aero_atmo::at(a).unwrap().density;
        let hi = tpt_aero_atmo::at(a + step).unwrap().density;
        prop_assert!(hi < lo);
    }

    #[test]
    fn flight_conditions_are_inside_the_isa_range((speed, alt) in tpt_aero_verify::flight_condition()) {
        prop_assert!(speed > 0.0);
        prop_assert!(tpt_aero_atmo::at(alt).is_ok());
    }
}
