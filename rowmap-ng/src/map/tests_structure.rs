//! Structure of H in the multiple Arnoldi process (docs/map.md §9).
//!
//! Checked on every schedule, with and without breakdown:
//! * (S1) zero pattern: `H[i, j] = 0` exactly for `i > c(j)` (closing row);
//! * (S2) a promoted column ends with `H[c(j), j] = ‖p_j‖ > 0`;
//! * (S3) monotone profile: `c(j) ≤ c(j+1)` (open columns count as `κ − 1`);
//! * (S4) row m has exactly `o(m)` nonzeros left of the diagonal, in the
//!   contiguous columns `m − o(m) .. m − 1` (generic data);
//! * (S5) `o` evolves as `o(m+1) = o(m) + 1 − [q_m promoted] − (remainders
//!   dropped in between)`, so the lower bandwidth `max o(m)` is 1 plus the
//!   accepted insertions after stage 1 that are followed by another vector.

use super::tests::{Nl, forced, three_eigenvalues};
use super::tests_nonaut::NlT;
use super::*;
use crate::grk4t::{GRK4T, ROS34PW2};

fn run(
    sys: &dyn OdeSystem,
    m: &RowMethod,
    sched: &[StageKind],
    nonaut: bool,
    kappa1: usize,
    omit: &[(usize, usize)],
) -> (Basis, StepOut) {
    let n = sys.dim();
    let t0 = 0.3;
    let u0: Vec<f64> = (0..n).map(|i| 0.3 * ((i as f64) * 0.7).sin()).collect();
    let mut fu = vec![0.0; n];
    sys.rhs(t0, &u0, &mut fu);
    let mut ft = vec![0.0; n];
    if nonaut {
        sys.dfdt(t0, &u0, &mut ft).unwrap();
    }
    let weights = vec![1e6; n];
    let mut basis = Basis::new(n);
    let out = step(
        sys,
        m,
        sched,
        t0,
        &u0,
        &fu,
        nonaut.then_some(&ft[..]),
        0.05,
        &weights,
        &forced(kappa1, omit),
        &mut basis,
    )
    .expect("step failed");
    (basis, out)
}

/// Checks S1–S4 and returns the lower bandwidth `max o(m)`.
fn check_structure(b: &Basis, generic: bool, label: &str) -> usize {
    let h = b.h_mat();
    let k = b.dim();
    let prof = b.profile();
    let open = b.open_at_entry();
    // S1
    for j in 0..k {
        for i in prof[j] + 1..k {
            assert_eq!(h.get(i, j), 0.0, "{label}: S1 violated at H[{i},{j}]");
        }
    }
    // S2
    for (j, c) in b.closed().iter().enumerate() {
        if let Some(c) = *c
            && c > j
            && generic
        {
            assert!(
                h.get(c, j) > 0.0,
                "{label}: S2 violated, H[{c},{j}] = {}",
                h.get(c, j)
            );
        }
    }
    // S3
    assert!(
        prof.windows(2).all(|w| w[0] <= w[1]),
        "{label}: S3 violated, profile {prof:?}"
    );
    // S4
    for m in 0..k {
        let nz: Vec<usize> = (0..m).filter(|&j| h.get(m, j) != 0.0).collect();
        let expected: Vec<usize> = (m - open[m]..m).collect();
        if generic {
            assert_eq!(nz, expected, "{label}: S4 violated in row {m}");
        } else {
            assert!(
                nz.iter().all(|j| expected.contains(j)),
                "{label}: S4 (pattern) violated in row {m}"
            );
        }
    }
    let sub: Vec<usize> = prof.iter().enumerate().map(|(j, &c)| c - j).collect();
    eprintln!("{label}: dim {k}, subdiagonals per column {sub:?}, o(m) {open:?}");
    open.iter().copied().max().unwrap_or(0)
}

#[test]
fn h_is_generalized_hessenberg() {
    let n = 40;
    let aut = Nl::new(n, 7);
    let non = NlT::new(n, 8);
    // Proved in lean/RowmapLean/HStructure.lean (`trace_*`, `bandwidths`):
    // o(0..κ) = 0, 1 (κ₁−1 times), then a schedule-specific tail, for every
    // κ₁ (no breakdown; GRK4T aut drops w₄ ∈ K₃). Bandwidth = max of the tail.
    let cases: [Case; 4] = [
        Case {
            label: "GRK4T aut",
            sys: &aut,
            m: &GRK4T,
            sched: &GRK4T_AUTONOMOUS,
            nonaut: false,
            tail: &[1, 2, 2, 2, 3, 3, 3],
            bw: 3,
        },
        Case {
            label: "GRK4T non-aut",
            sys: &non,
            m: &GRK4T,
            sched: &GRK4T_NONAUTONOMOUS,
            nonaut: true,
            tail: &[1, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4],
            bw: 4,
        },
        Case {
            label: "ROS34PW2 aut",
            sys: &aut,
            m: &ROS34PW2,
            sched: &ROS34PW2_AUTONOMOUS,
            nonaut: false,
            tail: &[1],
            bw: 1,
        },
        Case {
            label: "ROS34PW2 non-aut",
            sys: &non,
            m: &ROS34PW2,
            sched: &ROS34PW2_NONAUTONOMOUS,
            nonaut: true,
            tail: &[1, 2, 3, 3, 3],
            bw: 3,
        },
    ];
    for c in cases {
        for kappa1 in [c.m.order, 6, 9] {
            let (b, _) = run(c.sys, c.m, c.sched, c.nonaut, kappa1, &[]);
            let got = check_structure(&b, true, &format!("{} κ₁={kappa1}", c.label));
            assert_eq!(b.stats.dropped_promote, 0);
            let mut expected = vec![0];
            expected.extend(std::iter::repeat_n(1, kappa1 - 1));
            expected.extend_from_slice(c.tail);
            assert_eq!(
                b.open_at_entry(),
                &expected[..],
                "{} κ₁={kappa1}: o(m) ≠ Lean trace",
                c.label
            );
            assert_eq!(got, c.bw, "{} κ₁={kappa1}: bandwidth", c.label);
        }
    }
}

/// One schedule case of `h_is_generalized_hessenberg`.
struct Case<'a> {
    label: &'a str,
    sys: &'a dyn OdeSystem,
    m: &'a RowMethod,
    sched: &'a [StageKind],
    nonaut: bool,
    /// o(m) after stage 1, from the Lean theorems `trace_*`.
    tail: &'a [usize],
    bw: usize,
}

/// Breakdown: dropped remainders close a column without a positive entry and
/// lower the bandwidth by one; the pattern (S1, S3) is unaffected.
#[test]
fn structure_survives_breakdown() {
    let sys = three_eigenvalues(30);
    for kappa1 in [3, 4, 6] {
        let (b, _) = run(&sys, &GRK4T, &GRK4T_AUTONOMOUS, false, kappa1, &[]);
        check_structure(&b, false, &format!("breakdown κ₁={kappa1}"));
    }
}

/// The monotone profile is a property of the FIFO order of the closure
/// rounds: skipping a promotion (test hook only) leaves an old column open
/// while newer ones close, and S3 fails. S1 still holds.
#[test]
fn omitted_promotion_breaks_monotone_profile() {
    let sys = Nl::new(40, 7);
    // stage 2 (index 1): action 1 = promotion of the first pending column
    let (b, _) = run(&sys, &GRK4T, &GRK4T_AUTONOMOUS, false, 4, &[(1, 1)]);
    let prof = b.profile();
    assert!(
        prof.windows(2).any(|w| w[0] > w[1]),
        "expected a non-monotone profile: {prof:?}"
    );
}
