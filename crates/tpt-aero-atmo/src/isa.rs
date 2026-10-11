//! International Standard Atmosphere (ICAO/US Standard Atmosphere 1976
//! layer model), valid to 86 km geopotential.

/// ISA errors.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IsaError {
    /// Altitude outside the modeled range (0 .. 86 km geopotential).
    AltitudeOutOfRange,
}

impl core::fmt::Display for IsaError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            IsaError::AltitudeOutOfRange => write!(f, "altitude out of ISA range (0..86000 m)"),
        }
    }
}

impl core::error::Error for IsaError {}

/// One atmospheric state point.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct IsaState {
    /// Temperature (K).
    pub temperature: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Speed of sound (m/s).
    pub speed_of_sound: f64,
}

/// Layer base values (US Standard Atmosphere 1976):
/// (base geopotential altitude m, base temperature K, base pressure Pa,
///  lapse rate K/m).
const LAYERS: [(f64, f64, f64, f64); 7] = [
    (0.0, 288.15, 101_325.0, -6.5e-3),
    (11_000.0, 216.65, 22_632.06, 0.0),
    (20_000.0, 216.65, 5_474.889, 1.0e-3),
    (32_000.0, 228.65, 868.0187, 2.8e-3),
    (47_000.0, 270.65, 110.9063, 0.0),
    (51_000.0, 270.65, 66.93887, -2.8e-3),
    (71_000.0, 214.65, 3.956420, -2.0e-3),
];

#[allow(unused_imports)]
use tpt_aero_core::Real as _;

/// Specific gas constant for dry air (J/(kg·K)).
pub const R_AIR: f64 = 287.052_87;
/// Standard gravitational acceleration (m/s²).
pub const G0: f64 = 9.806_65;

/// The ISA state at a geopotential altitude (m).
///
/// # Errors
/// [`IsaError::AltitudeOutOfRange`] outside 0..86 km.
pub fn at(altitude_m: f64) -> Result<IsaState, IsaError> {
    if !(0.0..=86_000.0).contains(&altitude_m) {
        return Err(IsaError::AltitudeOutOfRange);
    }
    // Locate the layer (last base ≤ h).
    let mut i_layer = 0;
    for (i, layer) in LAYERS.iter().enumerate() {
        if altitude_m >= layer.0 {
            i_layer = i;
        }
    }
    let (h0, t0, p0, lapse) = LAYERS[i_layer];
    let dh = altitude_m - h0;
    let temperature = t0 + lapse * dh;
    let pressure = if lapse.abs() < 1.0e-12 {
        // Isothermal layer: exponential drop.
        p0 * (-G0 * dh / (R_AIR * t0)).exp()
    } else {
        // Gradient layer: power law.
        p0 * (temperature / t0).powf(-G0 / (lapse * R_AIR))
    };
    let density = pressure / (R_AIR * temperature);
    let speed_of_sound = (1.4 * R_AIR * temperature).sqrt();
    Ok(IsaState {
        temperature,
        pressure,
        density,
        speed_of_sound,
    })
}

/// The ISA state at a geopotential altitude (m) on a non-standard day whose
/// temperature is `delta_isa` kelvin above (positive) or below ISA.
///
/// The temperature profile is the ISA lapse-rate profile shifted by
/// `delta_isa` at every altitude. Pressure is integrated up from sea level
/// through the layers with that temperature, using the exact solution of the
/// hydrostatic equation `dp/dh = −ρg` for a linear temperature profile, so
/// the result stays self-consistent (`ρ = p / (R T)`). With `delta_isa = 0`
/// this agrees with [`at`].
///
/// # Errors
/// [`IsaError::AltitudeOutOfRange`] outside 0..86 km or for a non-finite
/// offset.
pub fn at_isa_offset(altitude_m: f64, delta_isa: f64) -> Result<IsaState, IsaError> {
    if !(0.0..=86_000.0).contains(&altitude_m) || !delta_isa.is_finite() {
        return Err(IsaError::AltitudeOutOfRange);
    }
    let mut pressure = LAYERS[0].2;
    let mut temperature = LAYERS[0].1 + delta_isa;
    for (i, &(h0, t0, _, lapse)) in LAYERS.iter().enumerate() {
        let h_top = LAYERS.get(i + 1).map_or(f64::INFINITY, |layer| layer.0);
        let h_end = altitude_m.min(h_top);
        let t_base = t0 + delta_isa;
        let dh = h_end - h0;
        temperature = t_base + lapse * dh;
        pressure = if lapse.abs() < 1.0e-12 {
            // Isothermal: exponential drop at the layer temperature.
            pressure * (-G0 * dh / (R_AIR * t_base)).exp()
        } else {
            // Gradient: p₂ = p₁ (T₂/T₁)^(−g/(L R)).
            pressure * (temperature / t_base).powf(-G0 / (lapse * R_AIR))
        };
        if altitude_m <= h_top {
            break;
        }
    }
    let density = pressure / (R_AIR * temperature);
    let speed_of_sound = (1.4 * R_AIR * temperature).sqrt();
    Ok(IsaState {
        temperature,
        pressure,
        density,
        speed_of_sound,
    })
}

/// Dynamic viscosity of air (Pa·s) at temperature `t` (K), from Sutherland's
/// law `μ = μ₀ (T/T₀)^{3/2} (T₀ + S) / (T + S)` with `μ₀ = 1.716e-5` Pa·s at
/// `T₀ = 273.15` K and `S = 110.4` K.
#[must_use]
pub fn dynamic_viscosity(t: f64) -> f64 {
    const MU0: f64 = 1.716e-5;
    const T0: f64 = 273.15;
    const S: f64 = 110.4;
    MU0 * (t / T0).powf(1.5) * (T0 + S) / (t + S)
}

/// Kinematic viscosity (m²/s) of an ISA state, `ν = μ / ρ`.
#[must_use]
pub fn kinematic_viscosity(state: &IsaState) -> f64 {
    dynamic_viscosity(state.temperature) / state.density
}

#[cfg(test)]
mod offset_tests {
    use super::*;

    #[test]
    fn zero_offset_matches_standard_atmosphere() {
        for &h in &[
            0.0_f64, 5_000.0, 11_000.0, 15_000.0, 25_000.0, 40_000.0, 60_000.0, 80_000.0,
        ] {
            let std = at(h).expect("valid");
            let off = at_isa_offset(h, 0.0).expect("valid");
            assert!((off.pressure / std.pressure - 1.0).abs() < 1e-4, "p at {h}");
            assert!((off.temperature - std.temperature).abs() < 1e-9, "T at {h}");
        }
    }

    #[test]
    fn hot_day_is_less_dense_and_keeps_hydrostatic_balance() {
        let cold = at_isa_offset(5_000.0, -15.0).expect("valid");
        let std = at(5_000.0).expect("valid");
        let hot = at_isa_offset(5_000.0, 15.0).expect("valid");
        assert!(hot.density < std.density && std.density < cold.density);
        // Hot air expands: pressure at a given altitude is higher than ISA.
        assert!(hot.pressure > std.pressure);
        // Ideal-gas consistency.
        assert!((hot.density * R_AIR * hot.temperature - hot.pressure).abs() < 1e-6);
        // Hydrostatic: a small step in altitude changes pressure by ρ g dh.
        let up = at_isa_offset(5_010.0, 15.0).expect("valid");
        let dp_dh = (up.pressure - hot.pressure) / 10.0;
        assert!(
            (dp_dh + G0 * hot.density).abs() / (G0 * hot.density) < 1e-3,
            "dp/dh = {dp_dh}"
        );
    }

    #[test]
    fn sutherland_viscosity_reference_values() {
        // Sea-level air: μ ≈ 1.789e-5 Pa·s at 288.15 K.
        let mu = dynamic_viscosity(288.15);
        assert!((mu - 1.789e-5).abs() < 2e-8, "μ = {mu}");
        let nu = kinematic_viscosity(&at(0.0).expect("valid"));
        assert!((nu - 1.46e-5).abs() < 1e-7, "ν = {nu}");
    }
}
/// Invert the ISA: geopotential altitude (m) from a pressure (Pa).
///
/// Useful for cockpit instrumentation; accurate within the modeled range.
///
/// # Errors
/// [`IsaError::AltitudeOutOfRange`] if the pressure is outside the range of
/// the layer model.
pub fn altitude_from_pressure(pressure_pa: f64) -> Result<f64, IsaError> {
    if !(3.0..=101_400.0).contains(&pressure_pa) {
        return Err(IsaError::AltitudeOutOfRange);
    }
    // Bisect on the monotone pressure profile — simple and layer-agnostic.
    let (mut lo, mut hi) = (0.0_f64, 86_000.0_f64);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if at(mid).map(|s| s.pressure).unwrap_or(0.0) > pressure_pa {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sea_level_reference() {
        let s = at(0.0).expect("valid");
        assert!((s.temperature - 288.15).abs() < 1e-9);
        assert!((s.pressure - 101_325.0).abs() < 1e-6);
        assert!((s.density - 1.225).abs() < 2e-4, "ρ = {}", s.density);
        assert!(
            (s.speed_of_sound - 340.29).abs() < 0.05,
            "a = {}",
            s.speed_of_sound
        );
    }

    #[test]
    fn tropopause_values() {
        let s = at(11_000.0).expect("valid");
        assert!((s.temperature - 216.65).abs() < 1e-9);
        assert!((s.pressure - 22_632.06).abs() < 0.01, "p = {}", s.pressure);
        // Density at the tropopause ≈ 0.3639 kg/m³.
        assert!((s.density - 0.3639).abs() < 1e-3, "ρ = {}", s.density);
    }

    #[test]
    fn stratopause_values() {
        let s = at(47_000.0).expect("valid");
        assert!((s.pressure - 110.9063).abs() < 0.01);
        // 20 km checkpoint: 5474.9 Pa at 216.65 K.
        let s = at(20_000.0).expect("valid");
        assert!((s.pressure - 5_474.889).abs() < 0.01);
        assert!((s.temperature - 216.65).abs() < 1e-9);
    }

    #[test]
    fn pressure_altitude_roundtrip() {
        for &h in &[0.0_f64, 1500.0, 5400.0, 11_000.0, 17_500.0, 32_000.0] {
            let p = at(h).expect("valid").pressure;
            let h_back = altitude_from_pressure(p).expect("in range");
            assert!((h_back - h).abs() < 0.1, "h {h} → {h_back}");
        }
    }

    #[test]
    fn rejects_out_of_range() {
        assert_eq!(at(-1.0), Err(IsaError::AltitudeOutOfRange));
        assert_eq!(at(90_000.0), Err(IsaError::AltitudeOutOfRange));
    }
}
