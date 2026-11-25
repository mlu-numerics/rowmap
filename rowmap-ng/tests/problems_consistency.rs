//! The Rust test problems agree with the Fortran ones: y0, f(t0,y0),
//! J(t0,y0)·v (v_k = sin k) and f_t(t0,y0), compared with the committed
//! `driver dump` output (tests/data/dump/).

mod problems;

use problems::{ALL, Problem, data_dir, read_table};
use rowmap::OdeSystem;

fn rel(a: &[f64], b: &[f64], scale: f64) -> f64 {
    let d = a
        .iter()
        .zip(b)
        .fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()));
    let s = b.iter().fold(scale, |m, y| m.max(y.abs()));
    d / s.max(1e-300)
}

#[test]
fn rust_problems_match_fortran_dumps() {
    for kind in ALL {
        let m = kind.coarse_m();
        let p = Problem::new(kind, m);
        let n = p.n();
        let data = read_table(&data_dir().join(format!("dump/{}_m{m}.txt", kind.name())));
        assert_eq!(data.len(), n, "{}", kind.name());
        let col = |c: usize| data.iter().map(|r| r[c]).collect::<Vec<_>>();
        let y0 = p.y0();
        let v: Vec<f64> = (1..=n).map(|k| (k as f64).sin()).collect();
        let mut f = vec![0.0; n];
        let mut z = vec![0.0; n];
        let mut ft = vec![0.0; n];
        p.rhs(0.0, &y0, &mut f);
        p.jac_vec(0.0, &y0, &v, &mut z);
        let _ = p.dfdt(0.0, &y0, &mut ft);
        // f is a difference of O(1/h^2) terms: scale by |J||y| ~ 8/h^2 max|y|
        let fscale =
            8.0 * (m as f64 + 1.0).powi(2) * y0.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
        let checks = [
            ("y0", rel(&y0, &col(0), 0.0)),
            ("f", rel(&f, &col(1), fscale)),
            ("jv", rel(&z, &col(2), 0.0)),
            ("ft", rel(&ft, &col(3), 0.0)),
        ];
        for (what, r) in checks {
            assert!(
                r < 1e-13,
                "{} m={m}: {what} differs, rel {r:e}",
                kind.name()
            );
        }
    }
}
