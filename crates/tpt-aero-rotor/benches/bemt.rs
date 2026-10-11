use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use tpt_aero_rotor::{
    maps::performance_map, BemtSolver, HubLoss, OperatingCondition, RotorGeometry, SectionModel,
    TipLoss,
};

const RHO: f64 = 1.225;
const OMEGA: f64 = 300.0;

fn solver(stations: usize) -> BemtSolver {
    let geometry =
        RotorGeometry::rectangular(1.0, 0.1, 2, 0.08, stations).with_linear_twist(0.4, 0.1);
    let section = SectionModel::Linear {
        cl_slope: 2.0 * core::f64::consts::PI,
        cl_max: 1.4,
        cd0: 0.01,
        cd_cl2: 0.01,
    };
    BemtSolver::new(geometry, section, TipLoss::Prandtl, HubLoss::None).expect("valid rotor")
}

fn hover_solve(c: &mut Criterion) {
    let mut group = c.benchmark_group("bemt_hover_solve");
    for stations in [10, 20, 40] {
        let s = solver(stations);
        let op = OperatingCondition::hover(OMEGA, RHO);
        group.bench_with_input(BenchmarkId::from_parameter(stations), &s, |b, s| {
            b.iter(|| s.solve(&op).expect("converges"));
        });
    }
    group.finish();
}

fn advance_sweep(c: &mut Criterion) {
    let speeds: Vec<f64> = (0..10).map(|i| 5.0 + 5.0 * f64::from(i)).collect();
    let s = solver(20);
    c.bench_function("bemt_10_point_performance_map", |b| {
        b.iter(|| performance_map(&s, OMEGA, RHO, &speeds, 0.0).expect("converges"));
    });
}

criterion_group!(benches, hover_solve, advance_sweep);
criterion_main!(benches);
