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

/// One atmospheric state point.
#[derive(Debug, Clone, Copy, PartialEq)]
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
        assert!((s.speed_of_sound - 340.29).abs() < 0.05, "a = {}", s.speed_of_sound);
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
