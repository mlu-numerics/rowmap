//! ROS34PW2 (order 3, public option) with the MAP: Krylov dimensions for
//! guaranteed order 3 (docs/map.md §8). Deviation from the exact ROS34PW2
//! method (full-dimensional Krylov space): order 4 (= p + 1) at the minimal
//! dimensions.
//!
//! The predicted orders were written down from the derivation in
//! docs/map.md §8 *before* the first run.

use super::tests::Nl;
use super::tests_nonaut::NlT;
use super::*;
use crate::grk4t::ROS34PW2;
use crate::krylov::tests::Rng;
use crate::linalg::nrm2;

const T0: f64 = 0.3;

fn one_step(
    sys: &dyn OdeSystem,
    u: &[f64],
    h: f64,
    cfg: &StepConfig,
    sched: &[StageKind],
) -> StepOut {
    let n = u.len();
    let mut fu = vec![0.0; n];
    sys.rhs(T0, u, &mut fu);
    let mut ft = vec![0.0; n];
    let ft_arg = if sys.autonomous() {
        None
    } else {
        sys.dfdt(T0, u, &mut ft).unwrap();
        Some(ft.as_slice())
    };
    let mut basis = Basis::new(n);
    step(
        sys,
        &ROS34PW2,
        sched,
        T0,
        u,
        &fu,
        ft_arg,
        h,
        &vec![1e6; n],
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

const HS: [f64; 8] = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125, 0.0015625];

/// Three more halvings, for the omission cases whose observed order is still
/// pre-asymptotic at h = 0.0016 (an h⁴ term dominates long). Not used for the
/// order-4 cases, whose deviation would reach rounding level.
const HS_LONG: [f64; 11] = [
    0.2,
    0.1,
    0.05,
    0.025,
    0.0125,
    0.00625,
    0.003125,
    0.0015625,
    7.8125e-4,
    3.90625e-4,
    1.953125e-4,
];

fn orders(sys: &dyn OdeSystem, u0: &[f64], cfg: &StepConfig, sched: &[StageKind]) -> Vec<f64> {
    orders_hs(sys, u0, cfg, sched, &HS)
}

fn orders_hs(
    sys: &dyn OdeSystem,
    u0: &[f64],
    cfg: &StepConfig,
    sched: &[StageKind],
    hs: &[f64],
) -> Vec<f64> {
    let n = u0.len();
    let d: Vec<f64> = hs
        .iter()
        .map(|&h| {
            let a = one_step(sys, u0, h, cfg, sched).u_new;
            let b = one_step(sys, u0, h, &forced(n, &[]), sched).u_new;
            nrm2(&a.iter().zip(&b).map(|(x, y)| x - y).collect::<Vec<_>>())
        })
        .collect();
    d.windows(2).map(|w| (w[0] / w[1]).log2()).collect()
}

fn tail(p: &[f64]) -> f64 {
    p[p.len() - 3..].iter().sum::<f64>() / 3.0
}

/// (description, schedule, kappa1, omitted actions, predicted order).
type Case = (
    &'static str,
    &'static [StageKind],
    usize,
    &'static [(usize, usize)],
    f64,
);

const AUTONOMOUS_CASES: [Case; 3] = [
    ("minimal 3/4/4/4", &ROS34PW2_AUTONOMOUS, 3, &[], 4.0),
    ("kappa1 = 2", &ROS34PW2_AUTONOMOUS, 2, &[], 3.0),
    (
        "stage 2: insert w2 omitted",
        &ROS34PW2_AUTONOMOUS,
        3,
        &[(1, 0)],
        3.0,
    ),
];

/// The first non-autonomous candidate inserted `w₄` too (dims 3/4/8/9). The
/// order tests showed it is not needed; kept as a tested variant.
const ROS34PW2_NA_WITH_W4: [StageKind; 4] = [
    StageKind::Arnoldi,
    StageKind::Extend {
        insert: true,
        rounds: 0,
    },
    StageKind::Extend {
        insert: true,
        rounds: 1,
    },
    StageKind::Extend {
        insert: true,
        rounds: 0,
    },
];

/// Stage 3 actions: 0 insert w₃; round 1 (column order): 1 = A q_κ₁,
/// 2 = A q_{w₂}, 3 = A q_{w₃}.
const NONAUTONOMOUS_CASES: [Case; 9] = [
    ("minimal 3/4/8/8", &ROS34PW2_NONAUTONOMOUS, 3, &[], 4.0),
    ("kappa1 = 2", &ROS34PW2_NONAUTONOMOUS, 2, &[], 3.0),
    ("autonomous schedule", &ROS34PW2_AUTONOMOUS, 3, &[], 2.0),
    (
        "stage 2: insert w2 omitted",
        &ROS34PW2_NONAUTONOMOUS,
        3,
        &[(1, 0)],
        3.0,
    ),
    (
        "stage 3: insert w3 omitted",
        &ROS34PW2_NONAUTONOMOUS,
        3,
        &[(2, 0)],
        2.0,
    ),
    (
        "stage 3: A q_k1 omitted",
        &ROS34PW2_NONAUTONOMOUS,
        3,
        &[(2, 1)],
        3.0,
    ),
    (
        "stage 3: A q_w2 omitted",
        &ROS34PW2_NONAUTONOMOUS,
        3,
        &[(2, 2)],
        3.0,
    ),
    (
        "stage 3: A q_w3 omitted",
        &ROS34PW2_NONAUTONOMOUS,
        3,
        &[(2, 3)],
        3.0,
    ),
    // first prediction: w₄ must be inserted (order 3 without). Measured: not
    // needed, so the minimal schedule omits it; the variant with w₄ also has 4.
    (
        "variant: w4 inserted (3/4/8/9)",
        &ROS34PW2_NA_WITH_W4,
        3,
        &[],
        4.0,
    ),
];

fn run_cases(sys: &dyn OdeSystem, u0: &[f64], cases: &[Case]) -> Vec<String> {
    let mut mismatches = Vec::new();
    for &(name, sched, k1, omit, pred) in cases {
        // the two stage-3 closure omissions converge slowly (see HS_LONG)
        let slow = matches!(omit, [(2, 1)] | [(2, 2)]);
        let p = if slow {
            orders_hs(sys, u0, &forced(k1, omit), sched, &HS_LONG)
        } else {
            orders(sys, u0, &forced(k1, omit), sched)
        };
        let q = tail(&p);
        eprintln!("{name:30} predicted {pred}  measured {q:.2}  orders {p:.2?}");
        if (q - pred).abs() > 0.35 {
            mismatches.push(format!("{name}: predicted {pred}, measured {q:.2}"));
        }
    }
    mismatches
}

#[test]
fn ros34pw2_minimal_dimensions() {
    let sys = Nl::new(40, 11);
    let u0 = Rng::new(12).vec(40);
    let out = one_step(&sys, &u0, 0.1, &forced(3, &[]), &ROS34PW2_AUTONOMOUS);
    assert_eq!(&out.dims, &[3, 4, 4, 4], "autonomous dims {:?}", out.dims);
    let sys = NlT::new(40, 41);
    let u0 = Rng::new(42).vec(40);
    let out = one_step(&sys, &u0, 0.1, &forced(3, &[]), &ROS34PW2_NONAUTONOMOUS);
    assert_eq!(
        &out.dims,
        &[3, 4, 8, 8],
        "non-autonomous dims {:?}",
        out.dims
    );
    assert_eq!(reserve(&ROS34PW2_AUTONOMOUS), 1);
    assert_eq!(reserve(&ROS34PW2_NONAUTONOMOUS), 5);
}

#[test]
fn ros34pw2_autonomous_orders() {
    let mut mismatches = Vec::new();
    for (sys, u0) in [
        (Nl::new(40, 11), Rng::new(12).vec(40)),
        (Nl::new(40, 21), Rng::new(22).vec(40)),
    ] {
        mismatches.extend(run_cases(&sys, &u0, &AUTONOMOUS_CASES));
    }
    assert!(mismatches.is_empty(), "{mismatches:#?}");
}

#[test]
fn ros34pw2_nonautonomous_orders() {
    let mut mismatches = Vec::new();
    for (sys, u0) in [
        (NlT::new(40, 41), Rng::new(42).vec(40)),
        (NlT::new(40, 51), Rng::new(52).vec(40)),
    ] {
        mismatches.extend(run_cases(&sys, &u0, &NONAUTONOMOUS_CASES));
    }
    assert!(mismatches.is_empty(), "{mismatches:#?}");
}
