//! One step of the Krylov-W method with the multiple Arnoldi process
//! (AGENTS.md §3.5, docs/map.md).

#[cfg(test)]
use crate::dense::Mat;
use crate::dense::Singular;
use crate::givens::Givens;
use crate::grk4t::{MAXS, RowMethod};
use crate::krylov::{Basis, NonFinite, Outcome};
use crate::linalg::{all_finite, axpy};
use crate::system::OdeSystem;

/// How a stage builds its Krylov space.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum StageKind {
    /// Stage 1: Arnoldi from `w₁` until the residual test (2.4) holds.
    Arnoldi,
    /// Insert `w_i` (if `insert`), then `rounds` closure rounds `K ← K + A K`
    /// (each round promotes all remainders pending at its start).
    Extend { insert: bool, rounds: usize },
}

/// MAP schedule for GRK4T, autonomous case (docs/map.md §3, §5):
/// `K₂ = C(K₁ + ⟨w₂⟩)`, `K₃ = K₂ + ⟨w₃⟩`, `K₄ = C(K₃ + ⟨w₄⟩)`.
/// Mirrored in `lean/schedule.toml`.
pub(crate) const GRK4T_AUTONOMOUS: [StageKind; 4] = [
    StageKind::Arnoldi,
    StageKind::Extend {
        insert: true,
        rounds: 1,
    },
    StageKind::Extend {
        insert: true,
        rounds: 0,
    },
    StageKind::Extend {
        insert: true,
        rounds: 1,
    },
];

/// MAP schedule for GRK4T, non-autonomous case (docs/map.md §7):
/// `K₂ = C²(K₁ + ⟨w₂⟩)`, `K₃ = K₂ + ⟨w₃⟩`, `K₄ = C(K₃ + ⟨w₄⟩)`.
/// Increments ≤ +5, +1, +5 (paper §4: 4, 9, 10, 15).
///
/// Review gate passed 2026-10-31. Mirrored in `lean/schedule.toml`.
pub(crate) const GRK4T_NONAUTONOMOUS: [StageKind; 4] = [
    StageKind::Arnoldi,
    StageKind::Extend {
        insert: true,
        rounds: 2,
    },
    StageKind::Extend {
        insert: true,
        rounds: 0,
    },
    StageKind::Extend {
        insert: true,
        rounds: 1,
    },
];

/// MAP schedule for ROS34PW2, autonomous (docs/map.md §8): insert `w₂`
/// (its O(h²) direction must be in the space, `b₂ ≠ 0`); `w₃, w₄` then leave
/// `K₂` only at O(h³). Minimal dims 3, 4, 4, 4.
pub(crate) const ROS34PW2_AUTONOMOUS: [StageKind; 4] = [
    StageKind::Arnoldi,
    StageKind::Extend {
        insert: true,
        rounds: 0,
    },
    StageKind::Extend {
        insert: false,
        rounds: 0,
    },
    StageKind::Extend {
        insert: false,
        rounds: 0,
    },
];

/// MAP schedule for ROS34PW2, non-autonomous (docs/map.md §8): `w₃` carries
/// an O(h) `f_t` direction (`β₃' ≠ 0`): insert it and close (a full round:
/// `f_t` is generic). `w₄` then leaves `K₃` only at O(h³): the first
/// prediction (insert `w₄`) was not confirmed; the order tests show it is not
/// needed. Dims ≤ 3, 4, 8, 8.
pub(crate) const ROS34PW2_NONAUTONOMOUS: [StageKind; 4] = [
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
        insert: false,
        rounds: 0,
    },
];

/// Minimal stage-1 dimension for order 4 (paper Thm 2.1). In general the
/// minimal κ₁ equals the order of the method (`RowMethod::order`): the
/// stage-1 Arnoldi error is O(h^κ₁). The step uses `m.order`; this constant
/// documents the GRK4T value for the tests.
#[cfg(test)]
pub(crate) const KAPPA1_MIN: usize = 4;

/// θ of the stage-1 residual test (2.4), `h‖r‖ ≤ θ` (paper: 0.1).
pub(crate) const THETA: f64 = 0.1;

/// Worst-case number of vectors added after stage 1 (every insertion
/// accepted, no drops): the space to keep free when stage 1 stops.
/// Stage 1 ends with one pending remainder; an insertion adds a vector and a
/// pending remainder; a closure round adds one vector per pending remainder.
pub(crate) fn reserve(sched: &[StageKind]) -> usize {
    let (mut pending, mut total) = (1usize, 0usize);
    for st in &sched[1..] {
        if let StageKind::Extend { insert, rounds } = *st {
            if insert {
                pending += 1;
                total += 1;
            }
            total += rounds * pending;
        }
    }
    total
}

/// Internal hooks for the order tests (not public API).
#[derive(Debug, Clone, Default)]
pub(crate) struct TestHook {
    /// Use exactly this stage-1 dimension (no residual test, no cap).
    pub kappa1: Option<usize>,
    /// Skip actions `(stage, action)`. Action 0 is the insertion; action
    /// `1 + k` is the k-th promotion of the stage's closure rounds.
    pub omit: Vec<(usize, usize)>,
}

#[derive(Debug, Clone)]
pub(crate) struct StepConfig {
    pub max_krylov: usize,
    pub theta: f64,
    pub hook: TestHook,
}

impl Default for StepConfig {
    fn default() -> Self {
        StepConfig {
            max_krylov: 70,
            theta: THETA,
            hook: TestHook::default(),
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct StepOut {
    /// `u + h Σ b_i k_i`.
    pub u_new: Vec<f64>,
    /// `h Σ (b̂_i − b_i) k_i`.
    pub err: Vec<f64>,
    /// Krylov dimension after each stage.
    pub dims: [usize; MAXS],
    /// f evaluations done in the step (stages only).
    pub nfev: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum StepFail {
    /// `I − γhH` numerically singular.
    Singular,
    /// NaN/Inf in f or J·v.
    NonFinite,
}

impl From<Singular> for StepFail {
    fn from(_: Singular) -> Self {
        StepFail::Singular
    }
}

impl From<NonFinite> for StepFail {
    fn from(_: NonFinite) -> Self {
        StepFail::NonFinite
    }
}

/// Solve the projected stage system `(I − γh H) l = coords` by a dense LU
/// from scratch (coords padded with zeros). Superseded by the incremental
/// Givens QR (`crate::givens`); kept as the test oracle.
#[cfg(test)]
pub(crate) fn solve_projected(
    basis: &Basis,
    coords: &[f64],
    gh: f64,
) -> Result<Vec<f64>, Singular> {
    let k = basis.dim();
    let mut m = basis.h_mat();
    for v in m.v.iter_mut() {
        *v *= -gh;
    }
    for i in 0..k {
        m.set(i, i, m.get(i, i) + 1.0);
    }
    let mut rhs = coords.to_vec();
    rhs.resize(k, 0.0);
    Mat::solve(&m, &rhs)
}

/// One Krylov-W step from `(t, u)` with step size `h`.
///
/// `fu = f(t, u)`; `ft = ∂f/∂t (t, u)` for non-autonomous problems (`None`
/// for autonomous ones); `weights_i = 1/(atol_i + rtol_i |u_i|)`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn step<S: OdeSystem + ?Sized>(
    sys: &S,
    m: &RowMethod,
    sched: &[StageKind],
    t: f64,
    u: &[f64],
    fu: &[f64],
    ft: Option<&[f64]>,
    h: f64,
    weights: &[f64],
    cfg: &StepConfig,
    basis: &mut Basis,
) -> Result<StepOut, StepFail> {
    assert_eq!(sched.len(), m.s);
    let n = u.len();
    let gh = m.gamma * h;
    basis.clear();
    // one incremental QR of I − γhH for all stages (docs/map.md §9.4)
    let mut giv = Givens::new(gh);
    let mut av = |v: &[f64], out: &mut [f64]| sys.jac_vec(t, u, v, out);
    let mut k: Vec<Vec<f64>> = Vec::with_capacity(m.s);
    let mut fvals: Vec<Vec<f64>> = Vec::with_capacity(m.s);
    let mut dims = [0usize; MAXS];
    let mut nfev = 0;
    let cap = cfg.max_krylov.saturating_sub(reserve(sched)).max(1);

    for i in 0..m.s {
        // f_i (stage reuse detected from the coefficients, e.g. f₄ = f₃)
        let fi = if i == 0 {
            fu.to_vec()
        } else if let Some(r) = m.reuses_f(i) {
            fvals[r].clone()
        } else {
            let mut ui = u.to_vec();
            for j in 0..i {
                if m.a[i][j] != 0.0 {
                    axpy(h * m.a[i][j], &k[j], &mut ui);
                }
            }
            let mut out = vec![0.0; n];
            {
                let _g = crate::prof::section("f (user)");
                sys.rhs(t + m.alpha(i) * h, &ui, &mut out);
            }
            nfev += 1;
            if !all_finite(&out) {
                return Err(StepFail::NonFinite);
            }
            out
        };
        let mut w = fi.clone();
        fvals.push(fi);
        for j in 0..i {
            axpy(m.c(i, j), &k[j], &mut w);
        }
        if let Some(ft) = ft {
            axpy(h * m.d(i), ft, &mut w);
        }

        // Remainders are needed only if a later closure round can promote them
        // (or for the stage-1 residual); otherwise only H columns are computed
        // and pending remainders are discarded (phase 7, saves sweeps).
        let later_rounds = sched[i + 1..]
            .iter()
            .any(|s| matches!(s, StageKind::Extend { rounds, .. } if *rounds > 0));
        let l = match sched[i] {
            StageKind::Arnoldi => {
                basis.set_keep_remainders(true);
                let (_, coords) = basis.insert(&w, &mut av)?;
                loop {
                    giv.update(basis);
                    let d = basis.dim();
                    let pend = basis.pending_cols();
                    // no pending remainder: invariant subspace (or w = 0), exact
                    let Some(&newest) = pend.last() else {
                        break giv.solve(basis, &coords)?;
                    };
                    let stop = match cfg.hook.kappa1 {
                        Some(k1) => d >= k1,
                        None => {
                            // κ₁ ≥ order of the method; the residual needs only
                            // the open (trailing) components of l, no full back
                            // substitution (docs/map.md §9.4)
                            d >= cap
                                || (d >= m.order && {
                                    let (first, lt) = giv.solve_tail(basis, &coords)?;
                                    h * gh * basis.residual_wrms(first, &lt, weights) <= cfg.theta
                                })
                        }
                    };
                    if stop {
                        break giv.solve(basis, &coords)?;
                    }
                    basis.promote(newest, &mut av)?;
                }
            }
            StageKind::Extend { insert, rounds } => {
                let omit = |a: usize| cfg.hook.omit.contains(&(i, a));
                basis.set_keep_remainders(rounds > 0 || later_rounds);
                let coords = if insert && !omit(0) {
                    basis.insert(&w, &mut av)?.1
                } else {
                    basis.coords(&w)
                };
                let mut action = 1;
                for r in 0..rounds {
                    basis.set_keep_remainders(r + 1 < rounds || later_rounds);
                    for j in basis.pending_cols() {
                        if !omit(action) {
                            let _: Outcome = basis.promote(j, &mut av)?;
                        }
                        action += 1;
                    }
                }
                giv.update(basis);
                giv.solve(basis, &coords)?
            }
        };

        if !later_rounds {
            basis.freeze_pending();
        }

        let mut ki = basis.combine(&l);
        for j in 0..i {
            axpy(-m.c(i, j), &k[j], &mut ki);
        }
        k.push(ki);
        dims[i] = basis.dim();
    }

    let mut u_new = u.to_vec();
    let mut err = vec![0.0; n];
    for i in 0..m.s {
        if m.b[i] != 0.0 {
            axpy(h * m.b[i], &k[i], &mut u_new);
        }
        if m.e(i) != 0.0 {
            axpy(h * m.e(i), &k[i], &mut err);
        }
    }
    Ok(StepOut {
        u_new,
        err,
        dims,
        nfev,
    })
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod tests_nonaut;
#[cfg(test)]
mod tests_ros34pw2;
#[cfg(test)]
mod tests_structure;
