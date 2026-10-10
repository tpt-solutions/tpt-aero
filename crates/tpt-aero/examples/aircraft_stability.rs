//! A complete light-aircraft configuration (wing + tailplane + fin): neutral
//! point and static margin, the full set of stability derivatives, and the
//! short-period mode as the centre of gravity moves.
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
    let build = |x_cg_over_c: f64| {
        VlmModel::from_surfaces(&[
            Surface::horizontal(wing.clone(), 20, 3),
            Surface::horizontal(WingPlanform::rectangular(4.0, 0.6), 10, 3).at(4.0, 0.0, 0.0),
            Surface::vertical(fin_planform(1.5, 0.9, 0.5), 6, 3).at(3.8, 0.0, 0.0),
        ])
        .expect("valid configuration")
        .with_reference(ReferenceGeometry {
            area: wing.area(),
            chord: wing.mac(),
            span: wing.span(),
            moment_center: Vec3::new(x_cg_over_c * wing.mac(), 0.0, 0.0),
        })
    };

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
}
