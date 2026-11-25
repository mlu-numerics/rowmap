//! Order and breakdown tests for the MAP step (AGENTS.md §5.3, docs/map.md §4).
//!
//! The order tests measure the deviation from the exact ROW method (GRK4T
//! with exact linear solves, realized as a full-dimensional Krylov space):
//! Thm 2.1 states `u_K − u_ROW = O(h⁵)` at the minimal dimensions
//! 4/7/8/11. Omitting any single action must lower that order (docs/map.md,
//! table "Necessity").

use super::*;
use crate::grk4t::GRK4T;
use crate::krylov::tests::{DenseOp, Rng};
use crate::linalg::nrm2;

/// Nonlinear test system `y' = M y + φ(y)` with
/// `φ_i = 0.5 sin(y_i) y_{i+1} + 0.3 y_i²` (cyclic). The nonlinearity makes
/// `f''` point in generic directions, as the necessity argument requires.
pub(super) struct Nl {
    m: DenseOp,
}

impl Nl {
    pub(super) fn new(n: usize, seed: u64) -> Self {
        let mut rng = Rng::new(seed);
        let mut m = DenseOp::random(n, &mut rng);
        let s = 3.0 / (n as f64).sqrt();
        for (idx, v) in m.a.iter_mut().enumerate() {
            *v *= s;
            if idx % (n + 1) == 0 {
                *v -= 1.0;
            }
        }
        Nl { m }
    }
}

impl OdeSystem for Nl {
    fn dim(&self) -> usize {
        self.m.n
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        let n = self.m.n;
        self.m.apply(y, dy);
        for i in 0..n {
            dy[i] += 0.5 * y[i].sin() * y[(i + 1) % n] + 0.3 * y[i] * y[i];
        }
    }
    fn jac_vec(&self, _t: f64, y: &[f64], v: &[f64], out: &mut [f64]) {
        let n = self.m.n;
        self.m.apply(v, out);
        for i in 0..n {
            let ip = (i + 1) % n;
            out[i] += (0.5 * y[i].cos() * y[ip] + 0.6 * y[i]) * v[i] + 0.5 * y[i].sin() * v[ip];
        }
    }
    fn autonomous(&self) -> bool {
        true
    }
}

/// Linear test system `y' = M y`.
pub(super) struct Lin {
    m: DenseOp,
}

impl OdeSystem for Lin {
    fn dim(&self) -> usize {
        self.m.n
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        self.m.apply(y, dy);
    }
    fn jac_vec(&self, _t: f64, _y: &[f64], v: &[f64], out: &mut [f64]) {
        self.m.apply(v, out);
    }
    fn autonomous(&self) -> bool {
        true
    }
}

fn one_step(sys: &dyn OdeSystem, u: &[f64], h: f64, cfg: &StepConfig) -> StepOut {
    let n = u.len();
    let mut fu = vec![0.0; n];
    sys.rhs(0.0, u, &mut fu);
    let weights = vec![1e6; n]; // atol = rtol = 1e-6 at |u| ≈ 0 (only used by the residual test)
    let mut basis = Basis::new(n);
    step(
        sys,
        &GRK4T,
        &GRK4T_AUTONOMOUS,
        0.0,
        u,
        &fu,
        None,
        h,
        &weights,
        cfg,
        &mut basis,
    )
    .expect("step failed")
}

pub(super) fn forced(kappa1: usize, omit: &[(usize, usize)]) -> StepConfig {
    StepConfig {
        max_krylov: 1000,
        hook: TestHook {
            kappa1: Some(kappa1),
            omit: omit.to_vec(),
        },
        ..Default::default()
    }
}

/// Exact ROW: full-dimensional Krylov space (T = A).
fn full(n: usize) -> StepConfig {
    forced(n, &[])
}

const HS: [f64; 6] = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625];

/// Observed orders log2(d(h)/d(h/2)) of the deviation from exact ROW.
fn deviation_orders(sys: &dyn OdeSystem, u0: &[f64], cfg: &StepConfig) -> (Vec<f64>, Vec<f64>) {
    let n = u0.len();
    let d: Vec<f64> = HS
        .iter()
        .map(|&h| {
            let a = one_step(sys, u0, h, cfg).u_new;
            let b = one_step(sys, u0, h, &full(n)).u_new;
            let diff: Vec<f64> = a.iter().zip(&b).map(|(x, y)| x - y).collect();
            nrm2(&diff)
        })
        .collect();
    let p = d.windows(2).map(|w| (w[0] / w[1]).log2()).collect();
    (d, p)
}

fn setup() -> (Nl, Vec<f64>) {
    let sys = Nl::new(40, 11);
    let u0: Vec<f64> = Rng::new(12).vec(40);
    (sys, u0)
}

fn sci(v: &[f64]) -> String {
    v.iter()
        .map(|x| format!("{x:.1e}"))
        .collect::<Vec<_>>()
        .join(" ")
}

/// Mean of the last three observed orders (asymptotic regime).
fn tail_order(p: &[f64]) -> f64 {
    p[p.len() - 3..].iter().sum::<f64>() / 3.0
}

#[test]
fn minimal_dimensions_4_7_8_11() {
    let (sys, u0) = setup();
    let out = one_step(&sys, &u0, 0.1, &forced(4, &[]));
    assert_eq!(&out.dims, &[4, 7, 8, 11], "dims = {:?}", out.dims);
}

#[test]
fn order_5_deviation_at_minimal_dimensions() {
    let (sys, u0) = setup();
    let (d, p) = deviation_orders(&sys, &u0, &forced(4, &[]));
    eprintln!("minimal dims: d = {}\n  orders = {p:.2?}", sci(&d));
    assert!(tail_order(&p) > 4.7, "observed order {p:?}");
}

#[test]
fn omitting_w4_insertion_is_harmless() {
    // w4 ∈ K3 (f4 = f3): positive control for the omit hook
    let (sys, u0) = setup();
    let (_, p) = deviation_orders(&sys, &u0, &forced(4, &[(3, 0)]));
    assert!(tail_order(&p) > 4.7, "observed order {p:?}");
}

/// (description, kappa1, omitted actions, predicted order of u_K − u_ROW),
/// see docs/map.md §4 "Necessity". Stage indices are 0-based (stage 2 = 1).
type Omission = (&'static str, usize, &'static [(usize, usize)], f64);
const OMISSIONS: [Omission; 12] = [
    // needed for classical order 4
    ("kappa1 = 3", 3, &[], 4.0),
    ("stage 2: insert w2", 4, &[(1, 0)], 3.0),
    ("stage 2: A q_k1", 4, &[(1, 1)], 4.0),
    ("stage 2: A q_k1+1", 4, &[(1, 2)], 4.0),
    ("insert w3 and w4", 4, &[(2, 0), (3, 0)], 4.0),
    // not needed for classical order 4 (docs/map.md §4): order stays 5
    ("stage 3: insert w3", 4, &[(2, 0)], 5.0),
    ("stage 4: A q_k1+2", 4, &[(3, 1)], 5.0),
    ("stage 4: A q_k1+3", 4, &[(3, 2)], 5.0),
    ("stage 4: A q_k3", 4, &[(3, 3)], 5.0),
    ("stage 4: whole closure", 4, &[(3, 1), (3, 2), (3, 3)], 5.0),
    (
        "w3 + stage-4 closure",
        4,
        &[(2, 0), (3, 1), (3, 2), (3, 3)],
        5.0,
    ),
    ("stage 4: insert w4", 4, &[(3, 0)], 5.0),
];

#[test]
fn omissions_match_predicted_orders() {
    let (sys, u0) = setup();
    let mut failures = Vec::new();
    for (name, k1, omit, pred) in OMISSIONS {
        let (d, p) = deviation_orders(&sys, &u0, &forced(k1, omit));
        let q = tail_order(&p);
        eprintln!(
            "{name:22} predicted {pred}  observed {q:.2}   orders {p:.2?}  d = {}",
            sci(&d)
        );
        if (q - pred).abs() > 0.35 {
            failures.push(format!("{name}: predicted {pred}, observed {q:.2}"));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

// ------------------------------------------------------------------ breakdown

/// Diagonal operator with 3 distinct eigenvalues: every Krylov space has
/// dimension ≤ 3, so stage 1 breaks down exactly.
pub(super) fn three_eigenvalues(n: usize) -> Lin {
    let d: Vec<f64> = (0..n).map(|i| [-1.0, -10.0, -100.0][i % 3]).collect();
    Lin {
        m: DenseOp::diag(&d),
    }
}

/// In exact arithmetic K₁ = span{w, Aw, A²w} is invariant and the step is
/// the exact ROW step with dims [3,3,3,3]. In floating point the computed
/// q₃ deviates from the invariant subspace by forward rounding errors
/// (~1e-14), so A q₃ has a genuine remainder ~1e-12·‖A q₃‖ > τ: a few
/// extra, orthonormal and harmless vectors are kept (AGENTS.md §3.4).
/// Required: the result equals exact ROW to rounding, the stage-1 residual is
/// at rounding level (so stage 1 stops at max(3, KAPPA1_MIN)), and the
/// dimensions stay far below n.
#[test]
fn exact_breakdown_reproduces_row() {
    let n = 30;
    let sys = three_eigenvalues(n);
    let u0 = Rng::new(21).vec(n);
    for h in [0.01, 0.1, 1.0] {
        let k = one_step(&sys, &u0, h, &StepConfig::default());
        let r = one_step(&sys, &u0, h, &full(n));
        let diff: Vec<f64> = k.u_new.iter().zip(&r.u_new).map(|(a, b)| a - b).collect();
        assert!(
            nrm2(&diff) <= 1e-13 * nrm2(&r.u_new),
            "h = {h}: {:e}",
            nrm2(&diff)
        );
        assert_eq!(k.dims[0], KAPPA1_MIN, "h = {h}: dims {:?}", k.dims);
        assert!(k.dims[3] <= 8, "h = {h}: dims {:?}", k.dims);
        assert!(k.u_new.iter().all(|v| v.is_finite()));
    }
}

#[test]
#[ignore]
fn breakdown_remainder_report() {
    let n = 30;
    let sys = three_eigenvalues(n);
    let u0 = Rng::new(21).vec(n);
    let mut fu = vec![0.0; n];
    sys.rhs(0.0, &u0, &mut fu);
    let mut av = |v: &[f64], o: &mut [f64]| sys.m.apply(v, o);
    let mut b = Basis::new(n);
    b.insert(&fu, &mut av).unwrap();
    for j in 0..3 {
        let p = b.outside(
            &(0..b.dim())
                .map(|i| if i == j { 1.0 } else { 0.0 })
                .collect::<Vec<_>>(),
        );
        let mut aq = vec![0.0; n];
        sys.m.apply(b.q(j), &mut aq);
        eprintln!("j={j}: |p_j| / |A q_j| = {:.2e}", nrm2(&p) / nrm2(&aq));
        let o = b.promote(j, &mut av).unwrap();
        eprintln!("   promote -> {o:?}, dim {}", b.dim());
    }
}

// ------------------------------------------------------ ε-continuity sweep

/// One sweep point: perturbed operator `diag(3 eigenvalues) + ε G`.
struct SweepPoint {
    eps: f64,
    /// ‖u_K(ε) − u_ROW(ε)‖ / ‖u_ROW(ε)‖ — deviation from exact ROW.
    dev_row: f64,
    /// ‖u_K(ε) − u_K(0)‖ — change of the Krylov result.
    change_k: f64,
    /// ‖u_ROW(ε) − u_ROW(0)‖ — true change of the ROW result.
    change_row: f64,
    /// max |QᵀQ − I| after the step.
    orth: f64,
    dims: [usize; MAXS],
}

fn sweep(cfg: &StepConfig, h: f64) -> Vec<SweepPoint> {
    let n = 30;
    let base = three_eigenvalues(n);
    let g = DenseOp::random(n, &mut Rng::new(31));
    let u0 = Rng::new(32).vec(n);
    let run = |eps: f64, cfg: &StepConfig| -> (StepOut, f64) {
        let mut m = DenseOp {
            n,
            a: base.m.a.clone(),
        };
        for (a, b) in m.a.iter_mut().zip(&g.a) {
            *a += eps * b;
        }
        let sys = Lin { m };
        let mut fu = vec![0.0; n];
        sys.rhs(0.0, &u0, &mut fu);
        let mut basis = Basis::new(n);
        let w = vec![1e6; n];
        let out = step(
            &sys,
            &GRK4T,
            &GRK4T_AUTONOMOUS,
            0.0,
            &u0,
            &fu,
            None,
            h,
            &w,
            cfg,
            &mut basis,
        )
        .expect("step failed");
        let mut orth: f64 = 0.0;
        for i in 0..basis.dim() {
            for j in 0..basis.dim() {
                let e = crate::linalg::dot(basis.q(i), basis.q(j)) - if i == j { 1.0 } else { 0.0 };
                orth = orth.max(e.abs());
            }
        }
        (out, orth)
    };
    let diff =
        |a: &[f64], b: &[f64]| nrm2(&a.iter().zip(b).map(|(x, y)| x - y).collect::<Vec<_>>());
    let (k0, _) = run(0.0, cfg);
    let (r0, _) = run(0.0, &full(n));
    let mut pts = Vec::new();
    for e in 0..=28 {
        let eps = 10f64.powf(-16.0 + 0.5 * e as f64);
        let (k, orth) = run(eps, cfg);
        let (r, _) = run(eps, &full(n));
        pts.push(SweepPoint {
            eps,
            dev_row: diff(&k.u_new, &r.u_new) / nrm2(&r.u_new),
            change_k: diff(&k.u_new, &k0.u_new),
            change_row: diff(&r.u_new, &r0.u_new),
            orth,
            dims: k.dims,
        });
    }
    pts
}

/// Requirements of AGENTS.md §5.3 item 4, for the default configuration
/// (residual test) and for the forced minimal dimension κ₁ = 4.
#[test]
fn eps_continuity_sweep() {
    for (label, cfg) in [
        ("default", StepConfig::default()),
        ("kappa1=4", forced(4, &[])),
    ] {
        let pts = sweep(&cfg, 0.1);
        for p in &pts {
            assert!(
                p.orth < 1e-13,
                "{label} eps={:.1e}: orthogonality {:.1e}",
                p.eps,
                p.orth
            );
        }
        // no jump: the Krylov result changes like the exact ROW result, up to
        // the (small, h-dependent) Krylov deviation
        let floor = pts.iter().map(|p| p.dev_row).fold(0.0, f64::max);
        for p in &pts {
            assert!(
                p.change_k <= 2.0 * p.change_row + 4.0 * floor + 1e-13,
                "{label} eps={:.1e}: change {:.1e} vs ROW change {:.1e} (floor {:.1e})",
                p.eps,
                p.change_k,
                p.change_row,
                floor
            );
        }
        // dimensions grow monotonically with eps
        for w in pts.windows(2) {
            assert!(
                w[1].dims[3] + 1 >= w[0].dims[3],
                "{label}: dims not monotone {:?} -> {:?}",
                w[0].dims,
                w[1].dims
            );
        }
    }
}

/// Diagnostic: prints the sweep and writes `results/sweep_autonomous.toml`
/// (plot: `uv run plots/breakdown_sweep.py`). Run with
/// `cargo test sweep_report -- --ignored --nocapture`.
#[test]
#[ignore]
fn sweep_report() {
    use std::fmt::Write as _;
    let mut toml = String::from("# eps-continuity sweep, autonomous MAP (src/map/tests.rs)\n");
    for (label, cfg) in [
        ("default", StepConfig::default()),
        ("kappa1_4", forced(4, &[])),
    ] {
        let pts = sweep(&cfg, 0.1);
        eprintln!("--- {label}");
        let _ = writeln!(toml, "\n[{label}]");
        let col =
            |f: &dyn Fn(&SweepPoint) -> String| pts.iter().map(f).collect::<Vec<_>>().join(", ");
        let _ = writeln!(toml, "eps = [{}]", col(&|p| format!("{:e}", p.eps)));
        let _ = writeln!(toml, "dev_row = [{}]", col(&|p| format!("{:e}", p.dev_row)));
        let _ = writeln!(
            toml,
            "change_k = [{}]",
            col(&|p| format!("{:e}", p.change_k))
        );
        let _ = writeln!(
            toml,
            "change_row = [{}]",
            col(&|p| format!("{:e}", p.change_row))
        );
        let _ = writeln!(toml, "orth = [{}]", col(&|p| format!("{:e}", p.orth)));
        let _ = writeln!(toml, "dim4 = [{}]", col(&|p| p.dims[3].to_string()));
        for p in &pts {
            eprintln!(
                "eps {:8.1e}  dev_row {:8.1e}  change_k {:8.1e}  change_row {:8.1e}  orth {:8.1e}  dims {:?}",
                p.eps, p.dev_row, p.change_k, p.change_row, p.orth, p.dims
            );
        }
    }
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/results/sweep_autonomous.toml");
    std::fs::create_dir_all(concat!(env!("CARGO_MANIFEST_DIR"), "/results")).unwrap();
    std::fs::write(path, toml).unwrap();
}

/// Export the ε-sweep problem for the Fortran comparison
/// (`fortran/experiments/sweep.f90`, `plots/breakdown_sweep_f77.py`):
/// `results/sweep_input.txt` (n, h, ε list, D, G, u0) and
/// `results/sweep_rust.txt` (per ε: u_ROW, u_K with κ₁ = 4, u_K default).
/// Run with `cargo test sweep_export -- --ignored`.
#[test]
#[ignore]
fn sweep_export() {
    use std::fmt::Write as _;
    let n = 30;
    let h = 0.1;
    let base = three_eigenvalues(n);
    let g = DenseOp::random(n, &mut Rng::new(31));
    let u0 = Rng::new(32).vec(n);
    let eps: Vec<f64> = (0..=28)
        .map(|e| 10f64.powf(-16.0 + 0.5 * e as f64))
        .collect();
    let mut inp = format!(
        "# eps-sweep input (src/map/tests.rs::sweep_export)\n{n} {} {h:e}\n",
        eps.len()
    );
    for e in &eps {
        let _ = writeln!(inp, "{e:.17e}");
    }
    for i in 0..n {
        let _ = writeln!(inp, "{:.17e}", base.m.a[i * n + i]);
    }
    for v in &g.a {
        let _ = writeln!(inp, "{v:.17e}");
    }
    for v in &u0 {
        let _ = writeln!(inp, "{v:.17e}");
    }
    let mut out =
        String::from("# per eps: n lines u_ROW, n lines u_K(kappa1=4), n lines u_K(default)\n");
    for &e in &eps {
        let mut m = DenseOp {
            n,
            a: base.m.a.clone(),
        };
        for (a, b) in m.a.iter_mut().zip(&g.a) {
            *a += e * b;
        }
        let sys = Lin { m };
        for cfg in [full(n), forced(4, &[]), StepConfig::default()] {
            for v in one_step(&sys, &u0, h, &cfg).u_new {
                let _ = writeln!(out, "{v:.17e}");
            }
        }
    }
    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/results");
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(format!("{dir}/sweep_input.txt"), inp).unwrap();
    std::fs::write(format!("{dir}/sweep_rust.txt"), out).unwrap();
}

/// The Rust schedule equals the one the Lean theorems are about
/// (`lean/schedule.toml`, exported from `Rowmap.grk4tAutonomous`).
#[test]
fn schedule_matches_lean() {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/lean/schedule.toml");
    let text = std::fs::read_to_string(path).expect("lean/schedule.toml missing");
    let doc: toml::Table = text.parse().expect("invalid TOML");
    let table = |key: &str| -> Vec<StageKind> {
        doc[key]
            .as_array()
            .expect("array of stages")
            .iter()
            .map(|s| match s["kind"].as_str().unwrap() {
                "arnoldi" => StageKind::Arnoldi,
                "extend" => StageKind::Extend {
                    insert: s["insert"].as_bool().unwrap(),
                    rounds: s["rounds"].as_integer().unwrap() as usize,
                },
                k => panic!("unknown stage kind {k}"),
            })
            .collect()
    };
    assert_eq!(table("grk4t_autonomous"), GRK4T_AUTONOMOUS.to_vec());
    assert_eq!(table("grk4t_nonautonomous"), GRK4T_NONAUTONOMOUS.to_vec());
    assert_eq!(table("ros34pw2_autonomous"), ROS34PW2_AUTONOMOUS.to_vec());
    assert_eq!(
        table("ros34pw2_nonautonomous"),
        ROS34PW2_NONAUTONOMOUS.to_vec()
    );
}

#[test]
fn zero_right_hand_side() {
    // f(y0) = 0: w1 = 0, all stages vanish, the step returns y0 exactly
    let n = 20;
    let sys = three_eigenvalues(n);
    let u0 = vec![0.0; n];
    let out = one_step(&sys, &u0, 0.5, &StepConfig::default());
    assert_eq!(out.u_new, u0);
    assert_eq!(out.err, vec![0.0; n]);
    assert_eq!(&out.dims, &[0, 0, 0, 0]);
}

#[test]
fn default_config_reaches_tolerance_dims_and_cap() {
    // the residual test (2.4) drives kappa1; the cap reserves room for stages 2-4
    let (sys, u0) = setup();
    let out = one_step(&sys, &u0, 0.1, &StepConfig::default());
    assert!(out.dims[0] >= KAPPA1_MIN);
    let cfg = StepConfig {
        max_krylov: 15,
        ..Default::default()
    };
    let out = one_step(&sys, &u0, 0.1, &cfg);
    assert!(out.dims[3] <= 15, "dims {:?}", out.dims);
    assert_eq!(reserve(&GRK4T_AUTONOMOUS), 9);
}
