//! Per-station rotor arrays: a constant profile must match the linear case it
//! replaces, and a wrong-length array must be refused.

use tpt_aero_cli::run_rotor;

fn case(extra_rotor: &str) -> String {
    format!(
        "[rotor]\nradius = 1.0\nblades = 2\nchord = 0.1\nstations = 4\n\
         root_twist_deg = 10.0\ntip_twist_deg = 10.0\n{extra_rotor}\n\
         [operating]\nrpm = 3000.0\n"
    )
}

#[test]
fn constant_profiles_reproduce_the_linear_case() {
    let linear = run_rotor(&case("")).expect("linear case");
    let profiled = run_rotor(&case(
        "chord_profile = [0.1, 0.1, 0.1, 0.1]\ntwist_profile_deg = [10.0, 10.0, 10.0, 10.0]",
    ))
    .expect("profiled case");
    assert!(
        (linear.ct - profiled.ct).abs() < 1e-12,
        "{} vs {}",
        linear.ct,
        profiled.ct
    );
}

#[test]
fn a_profile_must_have_one_value_per_station() {
    let err = run_rotor(&case("chord_profile = [0.1, 0.1]")).unwrap_err();
    assert!(err.to_string().contains("one value per station"), "{err}");
}

#[test]
fn a_nonpositive_chord_in_a_profile_is_refused() {
    assert!(run_rotor(&case("chord_profile = [0.1, 0.0, 0.1, 0.1]")).is_err());
}
