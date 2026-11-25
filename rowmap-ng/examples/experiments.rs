//! Paper experiments with the Rust integrator (AGENTS.md phases 4 and 6).
//!
//! ```sh
//! cargo run --release --example experiments -- [paper|coarse|all] [--reduced] [--only=<problem>]
//!     [--method=grk4t|ros34pw2] [--single-pass]
//! ```
//!
//! `--single-pass` uses one Gram–Schmidt pass (as rowmap.F) for the
//! compiler/memory comparison; files get the tag `_single`.
//!
//! For every autonomous paper problem and tol = 1e-2 … 1e-7, `rtol = atol =
//! tol`, max Krylov dimension 66 (as in the paper), the example writes
//! `results/<problem>_m<m>_grk4t_sft0_tol<k>.toml` (schema in
//! docs/data-format.md). With `--reduced`, it uses the schedule without the
//! stage-4 closure (docs/map.md §4) and writes `..._reduced_...` files.
//!
//! The CPU time is the mean over repeated runs (at least 0.5 s in total), as
//! in the Fortran driver.

#[path = "../tests/problems/mod.rs"]
mod problems;

use std::time::Instant;

use problems::{ALL, Problem, err_paper, reference};
use rowmap::{Method, OdeSystem, Rowmap};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let which = args
        .iter()
        .find(|a| !a.starts_with("--"))
        .map(String::as_str)
        .unwrap_or("all");
    let reduced = args.iter().any(|a| a == "--reduced");
    // --single-pass: one Gram–Schmidt pass, as rowmap.F (comparison only)
    let single = args.iter().any(|a| a == "--single-pass");
    // --only=<problem>: restrict to one problem (timing work)
    let only = args
        .iter()
        .find_map(|a| a.strip_prefix("--only=").map(str::to_owned));
    // --method=grk4t|ros34pw2 (default grk4t)
    let (method, mname) = match args.iter().find_map(|a| a.strip_prefix("--method=")) {
        Some("ros34pw2") => (Method::Ros34pw2, "ros34pw2"),
        Some("grk4t") | None => (Method::Grk4t, "grk4t"),
        Some(other) => panic!("unknown method {other}"),
    };
    assert!(
        !(reduced && method != Method::Grk4t),
        "--reduced applies to GRK4T only"
    );
    let out_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("results");
    std::fs::create_dir_all(&out_dir).unwrap();

    for kind in ALL {
        if only.as_deref().is_some_and(|o| o != kind.name()) {
            continue;
        }
        let autonomous = Problem::new(kind, kind.coarse_m()).autonomous();
        // the reduced schedule exists for the autonomous case only
        if reduced && !autonomous {
            continue;
        }
        // non-autonomous problems use the native form, paper strategy sft = 2
        let sft = if autonomous { 0 } else { 2 };
        let mut grids = Vec::new();
        if which == "paper" || which == "all" {
            grids.push(kind.paper_m());
        }
        if which == "coarse" || which == "all" {
            grids.push(kind.coarse_m());
        }
        for m in grids {
            let p = Problem::new(kind, m);
            let r = reference(kind, m);
            for k in 2..=7 {
                let tol = 10f64.powi(-k);
                let mut solver = Rowmap::new(tol, tol);
                solver.max_krylov(66).max_steps(1_000_000).method(method);
                if reduced {
                    solver.__schedule_without_stage4_closure();
                }
                if single {
                    solver.__reorthogonalize(false);
                }
                // repeat for a stable timing
                let mut reps = 0u32;
                let start = Instant::now();
                let (res, y) = loop {
                    let mut y = p.y0();
                    let res = solver.solve(&p, 0.0, &mut y, p.t_end());
                    reps += 1;
                    if start.elapsed().as_secs_f64() >= 0.5 || res.is_err() {
                        break (res, y);
                    }
                };
                let cpu = start.elapsed().as_secs_f64() / reps as f64;
                let (idid, err, nstep, nfeval, njv, kmean) = match &res {
                    Ok(s) => (
                        1,
                        err_paper(&y, &r),
                        s.stats.steps,
                        s.stats.f_evals,
                        s.stats.jv,
                        s.stats.krylov.map(|d| d.mean),
                    ),
                    Err(e) => {
                        eprintln!("{} m={m} tol=1e-{k}: {e}", kind.name());
                        (-1, f64::NAN, 0, 0, 0, [0.0; 4])
                    }
                };
                let sched = match (reduced, single) {
                    (true, false) => "no-stage4-closure",
                    (false, true) => "single-pass",
                    (true, true) => "no-stage4-closure-single-pass",
                    (false, false) => "paper",
                };
                let toml = format!(
                    "[run]\nproducer = \"rowmap-rs\"\nproblem = \"{}\"\nm = {m}\nn = {}\n\
                     method = \"{mname}\"\nschedule = \"{sched}\"\nrtol = {tol:e}\natol = {tol:e}\n\
                     max_krylov = 66\nsft = {sft}\nt_end = {:e}\nidid = {idid}\n\n[stats]\n\
                     err = {}\nnstep = {nstep}\nnfeval = {nfeval}\nnjv = {njv}\ncpu_time = {cpu:e}\n",
                    kind.name(),
                    p.n(),
                    p.t_end(),
                    if err.is_finite() {
                        format!("{err:e}")
                    } else {
                        "nan".into()
                    },
                );
                let mut tag = if reduced {
                    "grk4t_reduced".to_string()
                } else {
                    mname.to_string()
                };
                if single {
                    tag.push_str("_single");
                }
                let path = out_dir.join(format!("{}_m{m}_{tag}_sft{sft}_tol{k}.toml", kind.name()));
                std::fs::write(&path, toml).unwrap();
                println!(
                    "{:7} m={m:3} tol=1e-{k}  err={err:9.2e} nstep={nstep:6} nfeval={nfeval:7} \
                     cpu={cpu:8.4}s  kappa_mean=[{:.1} {:.1} {:.1} {:.1}]",
                    kind.name(),
                    kmean[0],
                    kmean[1],
                    kmean[2],
                    kmean[3]
                );
            }
        }
    }
}
