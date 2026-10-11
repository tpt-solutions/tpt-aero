//! Integration tests for the atmosphere crate through its public API: the ISA
//! table at reference altitudes, the offset and inversion round trips, the
//! discrete gust shape, and the turbulence field's realized intensity.

use tpt_aero_atmo::{
    altitude_from_pressure, at, at_isa_offset, dynamic_viscosity, kinematic_viscosity,
    turbulence::{simulate_axis, TurbulenceAxis, TurbulenceModel, WhiteNoise},
    DiscreteGust, GustShape, IsaError,
};

#[test]
fn reference_altitudes_match_the_standard_atmosphere() {
    // (altitude m, temperature K, pressure Pa): US Standard Atmosphere 1976.
    for &(h, t, p) in &[
        (0.0, 288.15, 101_325.0),
        (11_000.0, 216.65, 22_632.06),
        (20_000.0, 216.65, 5_474.889),
        (47_000.0, 270.65, 110.9063),
    ] {
        let s = at(h).expect("in range");
        assert!((s.temperature - t).abs() < 1e-6, "T at {h}");
        assert!(
            (s.pressure - p).abs() / p < 1e-4,
            "p at {h}: {}",
            s.pressure
        );
    }
}

#[test]
fn the_range_ends_at_the_geopotential_top() {
    assert!(at(84_852.0).is_ok());
    assert_eq!(at(86_000.0), Err(IsaError::AltitudeOutOfRange));
    assert_eq!(at(-1.0), Err(IsaError::AltitudeOutOfRange));
}

#[test]
fn offset_zero_equals_the_standard_day() {
    for h in [0.0, 5_000.0, 15_000.0, 60_000.0] {
        let a = at(h).expect("valid");
        let b = at_isa_offset(h, 0.0).expect("valid");
        // The offset model integrates the hydrostatic equation numerically, so it
        // agrees with the layer formula to the same 1e-4 the unit test uses.
        assert!((a.density / b.density - 1.0).abs() < 1e-4, "ρ at {h}");
    }
}

#[test]
fn pressure_altitude_round_trip_over_the_range() {
    for h in (0..=8).map(|k| k as f64 * 10_000.0 + 500.0) {
        let p = at(h).expect("valid").pressure;
        let back = altitude_from_pressure(p).expect("in range");
        assert!((back - h).abs() < 1.0, "{h} -> {back}");
    }
}

#[test]
fn sea_level_viscosity_and_kinematic_viscosity() {
    let s = at(0.0).expect("valid");
    let mu = dynamic_viscosity(s.temperature);
    assert!((mu - 1.789e-5).abs() < 2e-8, "μ = {mu}");
    let nu = kinematic_viscosity(&s);
    assert!((nu - 1.460e-5).abs() < 2e-8, "ν = {nu}");
}

#[test]
fn one_minus_cosine_gust_has_its_peak_at_half_the_gradient_distance() {
    let g = DiscreteGust::new(10.0, 30.0, GustShape::OneMinusCosine).expect("valid");
    assert_eq!(g.speed_at(0.0), 0.0);
    assert!((g.speed_at(30.0) - 10.0).abs() < 1e-9);
    assert_eq!(g.speed_at(60.0), 0.0);
    // The gradient peaks at x = h/2.
    let peak = g.gradient_at(15.0);
    assert!(peak > g.gradient_at(5.0) && peak > g.gradient_at(25.0));
}

#[test]
fn dryden_vertical_turbulence_realizes_its_design_intensity() {
    let mut axis = TurbulenceAxis::new(TurbulenceModel::Dryden, 1.2, 150.0, 40.0, true);
    let mut noise = WhiteNoise::new(2026);
    let series = simulate_axis(&mut axis, &mut noise, 0.05, 200_000);
    let tail = &series[20_000..];
    let mean = tail.iter().sum::<f64>() / tail.len() as f64;
    let var = tail.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / tail.len() as f64;
    assert!((var.sqrt() - 1.2).abs() < 0.12, "σ realized {}", var.sqrt());
}
