use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use tpt_aero_panel::{solve_doublet_panels, solve_source_panels, source::Point2};

fn circle(n: usize) -> Vec<Point2> {
    (0..n)
        .map(|i| {
            let t = 2.0 * core::f64::consts::PI * i as f64 / n as f64;
            (t.cos(), t.sin())
        })
        .collect()
}

fn panel_solves(c: &mut Criterion) {
    let mut group = c.benchmark_group("panel_solve");
    for n in [64, 128, 256] {
        let pts = circle(n);
        group.bench_with_input(BenchmarkId::new("source", n), &pts, |b, pts| {
            b.iter(|| solve_source_panels(pts, 0.05, 50.0).expect("solvable"));
        });
        group.bench_with_input(BenchmarkId::new("doublet", n), &pts, |b, pts| {
            b.iter(|| solve_doublet_panels(pts, 0.05, 50.0, 1.0).expect("solvable"));
        });
    }
    group.finish();
}

criterion_group!(benches, panel_solves);
criterion_main!(benches);
