use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use tpt_aero_core::WingPlanform;
use tpt_aero_vlm::{CamberModel, VlmModel};

const V_INF: f64 = 50.0;
const RHO: f64 = 1.225;

fn model(n_span: usize, n_chord: usize) -> VlmModel {
    let wing = WingPlanform::tapered(10.0, 1.4, 0.8, 0.5, 0.05);
    VlmModel::from_planform(
        &wing,
        n_span,
        n_chord,
        CamberModel::Naca4 { m: 0.02, p: 0.4 },
    )
    .expect("valid planform")
}

fn single_solve(c: &mut Criterion) {
    let mut group = c.benchmark_group("vlm_solve");
    for (ns, nc) in [(16, 4), (32, 4), (64, 4)] {
        let m = model(ns, nc);
        group.bench_with_input(BenchmarkId::from_parameter(ns * nc), &m, |b, m| {
            b.iter(|| m.solve(0.05, V_INF, RHO, &[]).expect("solvable"));
        });
    }
    group.finish();
}

fn alpha_sweep(c: &mut Criterion) {
    let alphas: Vec<f64> = (-4..=12).map(|d| f64::from(d).to_radians()).collect();
    let mut group = c.benchmark_group("vlm_17_angle_sweep");
    for (ns, nc) in [(16, 4), (32, 4)] {
        let m = model(ns, nc);
        group.bench_with_input(BenchmarkId::new("one_by_one", ns * nc), &m, |b, m| {
            b.iter(|| {
                for &a in &alphas {
                    m.solve(a, V_INF, RHO, &[]).expect("solvable");
                }
            });
        });
        group.bench_with_input(BenchmarkId::new("prepared_sweep", ns * nc), &m, |b, m| {
            b.iter(|| m.solve_sweep(&alphas, V_INF, RHO).expect("solvable"));
        });
    }
    group.finish();
}

criterion_group!(benches, single_solve, alpha_sweep);
criterion_main!(benches);
