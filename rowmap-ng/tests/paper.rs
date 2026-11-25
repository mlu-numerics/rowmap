//! The paper problems on coarse grids: accuracy against the DOPRI5
//! reference, benchmarked against the Fortran runs (tests/data/fortran/*.toml;
//! for the non-autonomous problems the native form sft = 2). The Fortran code
//! is a benchmark, not an oracle (AGENTS.md §1): the bounds below catch
//! regressions, not small differences.

mod problems;

use problems::{ALL, Problem, data_dir, err_paper, reference};
use rowmap::{OdeSystem, Rowmap, Status};

fn fortran(kind: &str, m: usize, k: i32, sft: u32) -> (f64, usize) {
    let path = data_dir().join(format!("fortran/{kind}_m{m}_grk4t_sft{sft}_tol{k}.toml"));
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let t: toml::Table = text.parse().unwrap();
    let s = t["stats"].as_table().unwrap();
    (
        s["err"].as_float().unwrap(),
        s["nstep"].as_integer().unwrap() as usize,
    )
}

#[test]
fn coarse_paper_problems_vs_reference_and_fortran() {
    let mut report = Vec::new();
    for kind in ALL {
        let m = kind.coarse_m();
        let p = Problem::new(kind, m);
        let r = reference(kind, m);
        let sft = if p.autonomous() { 0 } else { 2 };
        for k in 2..=6 {
            let tol = 10f64.powi(-k);
            let mut y = p.y0();
            let sol = Rowmap::new(tol, tol)
                .max_krylov(66)
                .solve(&p, 0.0, &mut y, p.t_end())
                .unwrap_or_else(|e| panic!("{} tol {tol}: {e}", kind.name()));
            assert_eq!(sol.status, Status::Finished);
            assert_eq!(sol.t, p.t_end());
            let err = err_paper(&y, &r);
            let (ef, nf) = fortran(kind.name(), m, k, sft);
            if !p.autonomous() {
                assert_eq!(sol.stats.dfdt_evals, sol.stats.accepted, "{}", kind.name());
            }
            let ns = sol.stats.steps;
            report.push(format!(
                "{:7} tol 1e-{k}: err {err:.2e} (f77 {ef:.2e})  nstep {ns} (f77 {nf})",
                kind.name()
            ));
            assert!(
                err <= 10.0 * ef,
                "{} tol 1e-{k}: err {err:e} vs Fortran {ef:e}",
                kind.name()
            );
            // ±30%, at least ±3 steps (the initial step estimates differ:
            // Hairer–Wanner here, a heuristic in rowmap.F)
            let slack = (0.3 * nf as f64).max(3.0);
            assert!(
                (ns as f64 - nf as f64).abs() <= slack,
                "{} tol 1e-{k}: nstep {ns} vs Fortran {nf}",
                kind.name()
            );
        }
    }
    eprintln!("{}", report.join("\n"));
}

#[test]
fn error_decreases_with_tolerance() {
    // tolerance proportionality, checked over two decades (not monotone per step)
    for kind in ALL {
        let m = kind.coarse_m();
        let p = Problem::new(kind, m);
        let r = reference(kind, m);
        let err = |tol: f64| {
            let mut y = p.y0();
            Rowmap::new(tol, tol)
                .solve(&p, 0.0, &mut y, p.t_end())
                .unwrap();
            err_paper(&y, &r)
        };
        let (e3, e5) = (err(1e-3), err(1e-5));
        assert!(
            e5 < 0.2 * e3,
            "{}: err(1e-5) = {e5:e}, err(1e-3) = {e3:e}",
            kind.name()
        );
    }
}
