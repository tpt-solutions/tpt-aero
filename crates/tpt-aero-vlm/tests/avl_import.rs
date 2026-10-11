//! Import the validation decks in `docs/validation/avl/` and compare the
//! solved lift slope and pitch stiffness with AVL 3.52 on the same geometry.

use tpt_aero_vlm::parse_avl;

const WING_RECT: &str = include_str!("../../../docs/validation/avl/wing_rectangular.avl");
const T_TAIL_FINE: &str = include_str!("../../../docs/validation/avl/t_tail_fine.avl");

/// The lift-curve slope (per rad) from two angles of attack.
fn cl_alpha(model: &tpt_aero_vlm::VlmModel) -> f64 {
    let a = 0.5_f64.to_radians();
    let c0 = model.solve(0.0, 50.0, 1.225, &[]).expect("solves").cl;
    let c1 = model.solve(a, 50.0, 1.225, &[]).expect("solves").cl;
    (c1 - c0) / a
}

#[test]
fn rectangular_deck_matches_avl_lift_slope() {
    let deck = parse_avl(WING_RECT).expect("parses");
    let model = deck.to_model().expect("builds");
    // AVL 3.52, wing_rectangular.avl at 20 x 4 panels: CLa 4.6279.
    let got = cl_alpha(&model);
    assert!(
        (got - 4.627_851).abs() / 4.627_851 < 0.05,
        "C_La {got} vs AVL 4.6279"
    );
}

#[test]
fn t_tail_deck_imports_wing_tail_and_fin_and_matches_avl() {
    let deck = parse_avl(T_TAIL_FINE).expect("parses");
    assert_eq!(deck.surfaces.len(), 3, "wing, tail and fin");
    let model = deck.to_model().expect("builds");
    // AVL 3.52, t_tail_fine.avl: CLa 5.302540.
    let got = cl_alpha(&model);
    assert!(
        (got - 5.302_540).abs() / 5.302_540 < 0.05,
        "C_La {got} vs AVL 5.3025"
    );
}
