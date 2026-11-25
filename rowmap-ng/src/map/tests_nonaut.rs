//! Non-autonomous MAP: order tests for the candidate schedule
//! `GRK4T_NONAUTONOMOUS` (docs/map.md §7). Deviation from the exact
//! non-autonomous ROW method (full-dimensional Krylov space, exact f_t).
//!
//! The predicted orders in `OMISSIONS_NA` were written down from the
//! analysis in docs/map.md §7 *before* the first run.

use super::*;
use crate::grk4t::GRK4T;
use crate::krylov::tests::{DenseOp, Rng};
use crate::linalg::nrm2;

/// `y' = M y + φ(y) + c(t)∘y + g(t)` with `φ_i = 0.5 sin(y_i) y_{i+1} + 0.3 y_i²`,
/// `c_i(t) = 0.5 sin(t + i)`, `g_i(t) = cos(2t + 0.3 i)`. So f_t, f_ty and f_tt
/// are all nonzero.
pub(super) struct NlT {
    m: DenseOp,
}

impl NlT {
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
        NlT { m }
    }
}

impl OdeSystem for NlT {
    fn dim(&self) -> usize {
        self.m.n
    }
    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) {
        let n = self.m.n;
        self.m.apply(y, dy);
        for i in 0..n {
            let fi = i as f64;
            dy[i] += 0.5 * y[i].sin() * y[(i + 1) % n]
                + 0.3 * y[i] * y[i]
                + 0.5 * (t + fi).sin() * y[i]
                + (2.0 * t + 0.3 * fi).cos();
        }
    }
    fn jac_vec(&self, t: f64, y: &[f64], v: &[f64], out: &mut [f64]) {
        let n = self.m.n;
        self.m.apply(v, out);
        for i in 0..n {
            let ip = (i + 1) % n;
            out[i] += (0.5 * y[i].cos() * y[ip] + 0.6 * y[i] + 0.5 * (t + i as f64).sin()) * v[i]
                + 0.5 * y[i].sin() * v[ip];
        }
    }
    fn dfdt(&self, t: f64, y: &[f64], ft: &mut [f64]) -> Result<(), crate::system::Unsupported> {
        for i in 0..self.m.n {
            let fi = i as f64;
            ft[i] = 0.5 * (t + fi).cos() * y[i] - 2.0 * (2.0 * t + 0.3 * fi).sin();
        }
        Ok(())
    }
}

const T0: f64 = 0.3;

fn one_step(sys: &NlT, u: &[f64], h: f64, cfg: &StepConfig, sched: &[StageKind]) -> StepOut {
    let n = u.len();
    let mut fu = vec![0.0; n];
    let mut ft = vec![0.0; n];
    sys.rhs(T0, u, &mut fu);
    sys.dfdt(T0, u, &mut ft).unwrap();
    let w = vec![1e6; n];
    let mut basis = Basis::new(n);
    step(
        sys,
        &GRK4T,
        sched,
        T0,
        u,
        &fu,
        Some(&ft),
        h,
        &w,
        cfg,
        &mut basis,
    )
    .expect("step failed")
}

fn forced(kappa1: usize, omit: &[(usize, usize)]) -> StepConfig {
    StepConfig {
        max_krylov: 1000,
        hook: TestHook {
            kappa1: Some(kappa1),
            omit: omit.to_vec(),
        },
        ..Default::default()
    }
}

const HS: [f64; 6] = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625];

fn orders(sys: &NlT, u0: &[f64], cfg: &StepConfig, sched: &[StageKind]) -> (Vec<f64>, Vec<f64>) {
    let n = u0.len();
    let d: Vec<f64> = HS
        .iter()
        .map(|&h| {
            let a = one_step(sys, u0, h, cfg, sched).u_new;
            let b = one_step(sys, u0, h, &forced(n, &[]), sched).u_new;
            nrm2(&a.iter().zip(&b).map(|(x, y)| x - y).collect::<Vec<_>>())
        })
        .collect();
    let p = d.windows(2).map(|w| (w[0] / w[1]).log2()).collect();
    (d, p)
}

fn sci(v: &[f64]) -> String {
    v.iter()
        .map(|x| format!("{x:.1e}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn tail(p: &[f64]) -> f64 {
    p[p.len() - 3..].iter().sum::<f64>() / 3.0
}

fn setup() -> (NlT, Vec<f64>) {
    (NlT::new(40, 41), Rng::new(42).vec(40))
}

#[test]
fn nonautonomous_minimal_dimensions_4_9_10_15() {
    let (sys, u0) = setup();
    let out = one_step(&sys, &u0, 0.1, &forced(4, &[]), &GRK4T_NONAUTONOMOUS);
    assert_eq!(&out.dims, &[4, 9, 10, 15], "dims {:?}", out.dims);
    assert_eq!(reserve(&GRK4T_NONAUTONOMOUS), 11);
}

#[test]
fn nonautonomous_order_5_deviation() {
    let (sys, u0) = setup();
    let (d, p) = orders(&sys, &u0, &forced(4, &[]), &GRK4T_NONAUTONOMOUS);
    eprintln!(
        "non-autonomous minimal dims: d = {}\n  orders {p:.2?}",
        sci(&d)
    );
    assert!(tail(&p) > 4.7, "{p:?}");
}

#[test]
fn autonomous_schedule_is_not_enough_for_nonautonomous_problems() {
    // missing A²f_t at stage 2 ⇒ (I−P₂)A²w₂ = O(h) ⇒ S₂ = O(h³) ⇒ order 4
    let (sys, u0) = setup();
    let (_, p) = orders(&sys, &u0, &forced(4, &[]), &GRK4T_AUTONOMOUS);
    assert!((tail(&p) - 4.0).abs() < 0.35, "{p:?}");
}

/// (description, omitted actions, first prediction, observed order of
/// u_K − u_ROW). The first predictions (docs/map.md §7, written before the
/// first run) were wrong in 6 of 10 cases. The test asserts the *observed*
/// orders (regression); see docs/map.md §7 for the corrected explanation.
///
/// Stage 2 actions: 0 insert w₂; round 1: 1 = A q_κ₁, 2 = A q_{w₂};
/// round 2: 3 = A²q_κ₁ (A of the remainder of A q_κ₁), 4 = A² q_{w₂}.
/// Stage 4 actions: 0 insert w₄; round 1 (column order): 1, 2 = A of the two
/// round-2 vectors of stage 2, 3 = A q_{w₃}, 4 = A q_{w₄}.
type Omission = (&'static str, &'static [(usize, usize)], f64, f64);
const OMISSIONS_NA: [Omission; 10] = [
    ("stage 2: insert w2", &[(1, 0)], 2.0, 2.0),
    ("stage 2: A q_k1", &[(1, 1)], 4.0, 4.0),
    ("stage 2: A q_w2", &[(1, 2)], 3.0, 4.0),
    ("stage 2: A^2 q_k1 (round 2)", &[(1, 3)], 5.0, 4.0),
    ("stage 2: A^2 q_w2 (round 2)", &[(1, 4)], 4.0, 4.0),
    ("stage 3: insert w3", &[(2, 0)], 5.0, 4.0),
    ("stage 4: insert w4", &[(3, 0)], 3.0, 4.0),
    ("stage 4: A of round-2 vectors", &[(3, 1), (3, 2)], 5.0, 4.0),
    ("stage 4: A q_w3", &[(3, 3)], 4.0, 4.0),
    ("stage 4: A q_w4", &[(3, 4)], 4.0, 5.0),
];

#[test]
fn nonautonomous_omissions() {
    // two independent problem instances
    for (sys, u0) in [setup(), (NlT::new(40, 51), Rng::new(52).vec(40))] {
        omissions_for(&sys, &u0);
    }
}

fn omissions_for(sys: &NlT, u0: &[f64]) {
    let mut mismatches = Vec::new();
    for (name, omit, first, observed) in OMISSIONS_NA {
        let (d, p) = orders(sys, u0, &forced(4, omit), &GRK4T_NONAUTONOMOUS);
        let q = tail(&p);
        eprintln!(
            "{name:30} first prediction {first}  expected {observed}  measured {q:.2}  orders {p:.2?}  d = {}",
            sci(&d)
        );
        if (q - observed).abs() > 0.35 {
            mismatches.push(format!("{name}: expected {observed}, measured {q:.2}"));
        }
    }
    assert!(mismatches.is_empty(), "{mismatches:#?}");
}

// ------------------------------------------- breakdown with f_t = 0, f_t ∥ f

/// `y' = D y` (3 distinct eigenvalues), with a prescribed `f_t` passed
/// directly to `step` (the step only uses `f_t` through `w_i += h d_i f_t`).
struct Diag3 {
    m: DenseOp,
}

impl OdeSystem for Diag3 {
    fn dim(&self) -> usize {
        self.m.n
    }
    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        self.m.apply(y, dy);
    }
    fn jac_vec(&self, _t: f64, _y: &[f64], v: &[f64], out: &mut [f64]) {
        self.m.apply(v, out);
    }
}

fn step_with_ft(sys: &Diag3, u: &[f64], ft: &[f64], h: f64, cfg: &StepConfig) -> StepOut {
    let n = u.len();
    let mut fu = vec![0.0; n];
    sys.rhs(0.0, u, &mut fu);
    let mut basis = Basis::new(n);
    step(
        sys,
        &GRK4T,
        &GRK4T_NONAUTONOMOUS,
        0.0,
        u,
        &fu,
        Some(ft),
        h,
        &vec![1e6; n],
        cfg,
        &mut basis,
    )
    .expect("step failed")
}

fn diag3(n: usize) -> Diag3 {
    let d: Vec<f64> = (0..n).map(|i| [-1.0, -10.0, -100.0][i % 3]).collect();
    Diag3 {
        m: DenseOp::diag(&d),
    }
}

/// f_t = 0 and f_t = c·f(u): the f_t directions add nothing new (K₁ is
/// invariant and contains f). The step must equal exact ROW to rounding,
/// without errors or NaN, and with small dimensions.
#[test]
fn nonautonomous_exact_breakdown_ft_zero_and_parallel() {
    let n = 30;
    let sys = diag3(n);
    let u0 = Rng::new(61).vec(n);
    let mut fu = vec![0.0; n];
    sys.rhs(0.0, &u0, &mut fu);
    let zero = vec![0.0; n];
    let par: Vec<f64> = fu.iter().map(|v| -0.7 * v).collect();
    for (label, ft) in [("f_t = 0", &zero), ("f_t = -0.7 f", &par)] {
        for h in [0.01, 0.1] {
            let k = step_with_ft(&sys, &u0, ft, h, &StepConfig::default());
            let r = step_with_ft(&sys, &u0, ft, h, &forced(n, &[]));
            let d = nrm2(
                &k.u_new
                    .iter()
                    .zip(&r.u_new)
                    .map(|(a, b)| a - b)
                    .collect::<Vec<_>>(),
            );
            assert!(
                d <= 1e-13 * nrm2(&r.u_new),
                "{label}, h = {h}: deviation {d:e}"
            );
            assert!(k.dims[3] <= 10, "{label}, h = {h}: dims {:?}", k.dims);
        }
    }
}

/// ε-continuity with `f_t = f + ε g` (g generic, n = 30, 3 eigenvalues):
/// the Krylov result changes like the exact ROW result, no jumps.
#[test]
fn nonautonomous_eps_sweep_ft() {
    let n = 30;
    let sys = diag3(n);
    let u0 = Rng::new(71).vec(n);
    let g = Rng::new(72).vec(n);
    let mut fu = vec![0.0; n];
    sys.rhs(0.0, &u0, &mut fu);
    let diff =
        |a: &[f64], b: &[f64]| nrm2(&a.iter().zip(b).map(|(x, y)| x - y).collect::<Vec<_>>());
    let run = |eps: f64, cfg: &StepConfig| {
        let ft: Vec<f64> = fu.iter().zip(&g).map(|(f, gi)| f + eps * gi).collect();
        step_with_ft(&sys, &u0, &ft, 0.1, cfg).u_new
    };
    let (k0, r0) = (run(0.0, &StepConfig::default()), run(0.0, &forced(n, &[])));
    let mut floor: f64 = 0.0;
    let mut pts = Vec::new();
    for e in 0..=28 {
        let eps = 10f64.powf(-16.0 + 0.5 * e as f64);
        let (k, r) = (run(eps, &StepConfig::default()), run(eps, &forced(n, &[])));
        let dev = diff(&k, &r) / nrm2(&r);
        floor = floor.max(dev);
        pts.push((eps, diff(&k, &k0), diff(&r, &r0), dev));
    }
    for (eps, ck, cr, dev) in pts {
        assert!(dev < 1e-12, "eps {eps:.1e}: deviation from ROW {dev:e}");
        assert!(
            ck <= 2.0 * cr + 4.0 * floor + 1e-13,
            "eps {eps:.1e}: jump {ck:e} vs ROW {cr:e}"
        );
    }
}
