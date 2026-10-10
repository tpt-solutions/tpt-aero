//! A complete light-aircraft configuration (wing + tailplane + fin): neutral
//! point and static margin, the full set of stability derivatives, the
//! short-period mode as the centre of gravity moves, level-flight trim across
//! speed, and the lateral-directional modes with and without dihedral.
//!
//! Wing-and-tail only (no fuselage or propeller effects), so the numbers are
//! indicative, not a substitute for a flight-test-validated model.
//!
//! Run: `cargo run --example aircraft_stability`

use tpt_aero::core::WingSection;
use tpt_aero::prelude::*;

fn fin_planform(height: f64, root: f64, tip: f64) -> WingPlanform {
    let section = |y: f64, chord: f64| WingSection {
        le: Vec3::new(0.0, y, 0.0),
        chord,
        twist: 0.0,
    };
    WingPlanform::from_sections(vec![section(0.0, root), section(height, tip)])
        .expect("two ordered sections")
}

fn main() {
    let (v, rho, alpha) = (50.0, 1.225, 3.0_f64.to_radians());
    let wing = WingPlanform::rectangular(10.0, 1.2);
    let build = |x_cg_over_c: f64| build_with(&wing, x_cg_over_c);
    fn build_with(wing: &WingPlanform, x_cg_over_c: f64) -> VlmModel {
        VlmModel::from_surfaces(&[
            Surface::horizontal(wing.clone(), 20, 3).with_control(
                "aileron",
                0.75,
                2.5,
                4.8,
                ControlSymmetry::Antisymmetric,
            ),
            Surface::horizontal(WingPlanform::rectangular(4.0, 0.6), 10, 3)
                .at(4.0, 0.0, 0.0)
                .with_control("elevator", 0.6, 0.0, 2.0, ControlSymmetry::Symmetric),
            Surface::vertical(fin_planform(1.5, 0.9, 0.5), 6, 3)
                .at(3.8, 0.0, 0.0)
                .with_control("rudder", 0.6, 0.0, 1.5, ControlSymmetry::Symmetric),
        ])
        .expect("valid configuration")
        .with_reference(ReferenceGeometry {
            area: wing.area(),
            chord: wing.mac(),
            span: wing.span(),
            moment_center: Vec3::new(x_cg_over_c * wing.mac(), 0.0, 0.0),
        })
    }

    let model = build(0.25);
    let d = longitudinal_fd(&model, alpha, v, rho, 1e-4);
    let c = wing.mac();
    let x_np = 0.25 * c - d.cm_alpha / d.cl_alpha * c;
    println!("=== Wing + tailplane + fin, V = {v} m/s, alpha = 3 deg ===\n");
    println!(
        "Neutral point: x = {x_np:.3} m = {:.1} % MAC",
        100.0 * x_np / c
    );

    let lat = lateral_fd(&model, alpha, v, rho).expect("solvable");
    println!("\nLongitudinal (per rad / per q c/2V):");
    println!(
        "  CL_alpha {:.3}   Cm_alpha {:.3}   CL_q {:.2}   Cm_q {:.2}",
        d.cl_alpha, d.cm_alpha, d.cl_q, d.cm_q
    );
    println!("Lateral-directional:");
    println!(
        "  Cy_beta {:.3}   Cl_beta {:.4}   Cn_beta {:.3}",
        lat.cy_beta, lat.cl_beta, lat.cn_beta
    );
    println!(
        "  Cl_p {:.3}   Cn_p {:.3}   Cl_r {:.3}   Cn_r {:.3}",
        lat.cl_p, lat.cn_p, lat.cl_r, lat.cn_r
    );

    println!("\nCG sweep (mass 600 kg, I_y 900 kg m^2):");
    println!(
        "{:>8} {:>14} {:>12} {:>10}",
        "CG [%MAC]", "static margin", "wn [rad/s]", "zeta"
    );
    for cg in [0.15, 0.25, 0.35, 0.45, 0.55] {
        let m = build(cg);
        let dd = longitudinal_fd(&m, alpha, v, rho, 1e-4);
        let margin = -dd.cm_alpha / dd.cl_alpha; // fraction of MAC ahead of the NP
        let params = PlantParameters {
            v_inf: v,
            chord: c,
            area: wing.area(),
            mass: 600.0,
            inertia_y: 900.0,
            density: rho,
        };
        let plant = LongitudinalPlant::from_derivatives(&dd, &params, 0.05);
        match plant.short_period() {
            Some((wn, zeta)) => println!(
                "{:>8.0} {:>13.1}% {wn:>12.2} {zeta:>10.2}",
                100.0 * cg,
                100.0 * margin
            ),
            None => println!(
                "{:>8.0} {:>13.1}% {:>12} {:>10}",
                100.0 * cg,
                100.0 * margin,
                "unstable",
                "-"
            ),
        }
    }

    println!(
        "
Level-flight trim (600 kg, CG at 25 % MAC, CD0 = 0.03):"
    );
    println!(
        "{:>8} {:>10} {:>12} {:>10} {:>11}",
        "V [m/s]", "alpha [deg]", "elevator [deg]", "CL", "thrust [N]"
    );
    let trim_model = build(0.25);
    for speed in [35.0, 45.0, 55.0, 70.0] {
        let req = TrimRequest {
            weight_n: 600.0 * 9.80665,
            v_inf: speed,
            density: rho,
            elevator: "elevator",
            cd0: 0.03,
        };
        match trim_level_flight(&trim_model, &req) {
            Ok(t) => println!(
                "{speed:>8.0} {:>10.2} {:>12.2} {:>10.3} {:>11.0}",
                t.alpha_rad.to_degrees(),
                t.elevator_rad.to_degrees(),
                t.cl,
                t.thrust_n
            ),
            Err(e) => println!("{speed:>8.0}  trim failed: {e}"),
        }
    }

    // Both longitudinal oscillations at the 50 m/s trim point.
    let req = TrimRequest {
        weight_n: 600.0 * 9.80665,
        v_inf: v,
        density: rho,
        elevator: "elevator",
        cd0: 0.03,
    };
    let t = trim_level_flight(&trim_model, &req).expect("trims");
    let dt = longitudinal_fd(&trim_model, t.alpha_rad, v, rho, 1e-4);
    let long_params = FullLongitudinalParameters {
        v_inf: v,
        chord: c,
        area: wing.area(),
        density: rho,
        mass: 600.0,
        inertia_y: 900.0,
        cl: t.cl,
        cd: t.cdi + 0.03,
        cd_alpha: drag_alpha_fd(&trim_model, t.alpha_rad, v, rho).expect("solvable"),
    };
    let modes = FullLongitudinalPlant::from_derivatives(&dt, &long_params)
        .and_then(|p| p.modes())
        .expect("oscillatory modes");
    println!(
        "
Longitudinal modes at 50 m/s trim (L/D {:.1}):",
        t.cl / long_params.cd
    );
    println!(
        "  short period: wn {:.2} rad/s, zeta {:.2}",
        modes.short_period_wn, modes.short_period_zeta
    );
    println!(
        "  phugoid:      wn {:.3} rad/s, zeta {:.3}  (Lanchester wn {:.3})",
        modes.phugoid_wn,
        modes.phugoid_zeta,
        2.0_f64.sqrt() * 9.80665 / v
    );

    println!(
        "
Control power per radian of command (positive = trailing edge down on +y):"
    );
    for name in ["elevator", "aileron", "rudder"] {
        let cd = control_derivatives(&trim_model, name, alpha, v, rho).expect("control exists");
        println!(
            "  {name:9} dCL {:+.3}  dCm {:+.3}  dCy {:+.3}  dCl {:+.3}  dCn {:+.3}",
            cd.cl, cd.cm, cd.cy, cd.c_roll, cd.c_yaw
        );
    }

    println!(
        "
Lateral-directional modes at 50 m/s (I_x 1000, I_z 1800 kg m^2):"
    );
    println!(
        "{:>12} {:>9} {:>9} {:>14} {:>16}",
        "dihedral", "DR wn", "DR zeta", "roll tau [s]", "spiral"
    );
    for dihedral in [0.0_f64, 3.0, 6.0] {
        let planform = WingPlanform::tapered(10.0, 1.2, 1.2, 0.0, dihedral.to_radians());
        let m = build_with(&planform, 0.25);
        let d = lateral_fd(&m, alpha, v, rho).expect("solvable");
        let params = LateralParameters {
            v_inf: v,
            span: m.reference().span,
            area: m.reference().area,
            density: rho,
            mass: 600.0,
            inertia_x: 1000.0,
            inertia_z: 1800.0,
            inertia_xz: 0.0,
        };
        let modes = LateralPlant::from_derivatives(&d, &params)
            .and_then(|p| p.modes())
            .expect("standard spectrum");
        let spiral = match modes.spiral_time_to_double() {
            Some(t) => format!("doubles in {t:.0} s"),
            None => "stable".to_string(),
        };
        println!(
            "{:>10.0} deg {:>8.2} {:>9.2} {:>14.3} {:>16}",
            dihedral,
            modes.dutch_roll_wn,
            modes.dutch_roll_zeta,
            modes.roll_time_constant().unwrap_or(f64::NAN),
            spiral
        );
    }
}
