//! Where does the time go? Splits the solve time into problem evaluations
//! (f, J·v, f_t) and solver overhead (orthogonalization, small solves, ...).
//!
//! ```sh
//! cargo run --release --example profile -- [problem] [tol] [grk4t|ros34pw2]
//! ```

#[path = "../tests/problems/mod.rs"]
mod problems;

use std::cell::Cell;
use std::time::Instant;

use problems::{ALL, Problem};
use rowmap::{Method, OdeSystem, Rowmap, Unsupported};

struct Timed<'a> {
    p: &'a Problem,
    t_rhs: Cell<f64>,
    t_jv: Cell<f64>,
    n_rhs: Cell<usize>,
    n_jv: Cell<usize>,
}

impl OdeSystem for Timed<'_> {
    fn dim(&self) -> usize {
        self.p.dim()
    }
    fn autonomous(&self) -> bool {
        self.p.autonomous()
    }
    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) {
        let s = Instant::now();
        self.p.rhs(t, y, dy);
        self.t_rhs.set(self.t_rhs.get() + s.elapsed().as_secs_f64());
        self.n_rhs.set(self.n_rhs.get() + 1);
    }
    fn jac_vec(&self, t: f64, y: &[f64], v: &[f64], out: &mut [f64]) {
        let s = Instant::now();
        self.p.jac_vec(t, y, v, out);
        self.t_jv.set(self.t_jv.get() + s.elapsed().as_secs_f64());
        self.n_jv.set(self.n_jv.get() + 1);
    }
    fn dfdt(&self, t: f64, y: &[f64], ft: &mut [f64]) -> Result<(), Unsupported> {
        self.p.dfdt(t, y, ft)
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let tol: f64 = args.get(1).map(|s| s.parse().unwrap()).unwrap_or(1e-5);
    let method = match args.get(2).map(String::as_str) {
        Some("ros34pw2") => Method::Ros34pw2,
        _ => Method::Grk4t,
    };
    for kind in ALL {
        if args.first().is_some_and(|a| a != kind.name() && a != "all") {
            continue;
        }
        let p = Problem::new(kind, kind.paper_m());
        let sys = Timed {
            p: &p,
            t_rhs: 0.0.into(),
            t_jv: 0.0.into(),
            n_rhs: 0.into(),
            n_jv: 0.into(),
        };
        let mut y = p.y0();
        let start = Instant::now();
        let sol = Rowmap::new(tol, tol)
            .max_krylov(66)
            .max_steps(1_000_000)
            .method(method)
            .solve(&sys, 0.0, &mut y, p.t_end())
            .unwrap();
        let total = start.elapsed().as_secs_f64();
        let user = sys.t_rhs.get() + sys.t_jv.get();
        let k: Vec<String> = sol
            .stats
            .krylov
            .iter()
            .map(|d| format!("{:.1}", d.mean))
            .collect();
        println!(
            "{:7} n={:5} tol={tol:.0e} {method:?}: total {total:.4}s  f {:.4}s ({}x)  Jv {:.4}s ({}x)  \
             solver {:.4}s ({:.0}%)  steps {}  kappa [{}]",
            kind.name(),
            p.n(),
            sys.t_rhs.get(),
            sys.n_rhs.get(),
            sys.t_jv.get(),
            sys.n_jv.get(),
            total - user,
            100.0 * (total - user) / total,
            sol.stats.steps,
            k.join(" ")
        );
        // with `--features profile`: breakdown inside the solver
        #[cfg(feature = "profile")]
        {
            let mut t = rowmap::prof::take();
            t.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            for (name, s, calls) in t {
                println!(
                    "    {name:22} {s:8.4}s  {:5.1}%  {calls:8} calls",
                    100.0 * s / total
                );
            }
        }
    }
}
