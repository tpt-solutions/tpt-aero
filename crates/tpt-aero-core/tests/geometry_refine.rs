//! `refine` keeps every original section as a station, so a kinked planform
//! keeps its area at any resolution.

use tpt_aero_core::{Vec3, WingPlanform, WingSection};

fn kinked() -> WingPlanform {
    // Inboard panel, then a kink at y = 2 with a different taper, then the tip.
    WingPlanform::from_sections(vec![
        WingSection {
            le: Vec3::new(0.0, -5.0, 0.0),
            chord: 1.0,
            twist: 0.0,
        },
        WingSection {
            le: Vec3::new(0.0, -2.0, 0.0),
            chord: 1.0,
            twist: 0.0,
        },
        WingSection {
            le: Vec3::new(0.5, 2.0, 0.0),
            chord: 0.4,
            twist: 0.0,
        },
        WingSection {
            le: Vec3::new(0.5, 5.0, 0.0),
            chord: 0.4,
            twist: 0.0,
        },
    ])
    .expect("valid")
}

#[test]
fn refine_preserves_area_and_original_stations() {
    let p = kinked();
    for n in [4, 5, 7, 11, 40] {
        let r = p.refine(n);
        assert_eq!(r.sections().len(), n, "n = {n}");
        assert!(
            (r.area() - p.area()).abs() < 1e-12,
            "n = {n}: {} vs {}",
            r.area(),
            p.area()
        );
        for kink in [-2.0, 2.0] {
            assert!(
                r.sections().iter().any(|s| (s.le.y - kink).abs() < 1e-12),
                "kink at {kink} lost for n = {n}"
            );
        }
    }
}
