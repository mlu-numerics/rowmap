//! Krylov basis for the multiple Arnoldi process (AGENTS.md §3.3, §3.4).
//!
//! The basis replaces the `kl`/`asv`/`reortn` bookkeeping of rowmap.F by an
//! explicit invariant:
//!
//! * `q_1..q_κ` are orthonormal (two Gram–Schmidt passes everywhere);
//! * for every basis vector `q_j` the product `A q_j` was computed when `q_j`
//!   was added, so `H = QᵀAQ` is always complete;
//! * the remainder `p_j = (I − QQᵀ) A q_j` is kept as *pending* until it is
//!   promoted to a basis vector or dropped as numerically dependent, so
//!   `A q_j = Q H[:,j] + p_j` (pending), `= Q H[:,j]` (promoted) or
//!   `= Q H[:,j] + O(τ‖A q_j‖)` (dropped).
//!
//! The same dependency test (`DEP_TOL`) decides for inserted right-hand
//! sides and for promoted remainders.

#[cfg(test)]
use crate::dense::Mat;
use crate::linalg::{all_finite, axpy, dot, nrm2, q_times, qt_v, scale, v_minus_qc, wrms};

/// Relative dependency threshold τ (AGENTS.md §3.4): a candidate `v₀` whose
/// component outside the current space, after two orthogonalization passes,
/// satisfies `‖v₂‖ ≤ τ‖v₀‖` lies in the space to working precision and is
/// dropped. After two passes an accepted vector is orthonormal to working
/// precision even if it is mostly rounding noise, and such a direction is
/// harmless for the Galerkin approximation, so τ is chosen just above the
/// rounding level of the projection (≈ ε√κ), not at √ε. Subject to tuning
/// with the ε-continuity sweep (phase 3).
pub(crate) const DEP_TOL: f64 = 1e3 * f64::EPSILON;

/// `f` or `J·v` returned NaN/Inf.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct NonFinite;

/// Result of adding a candidate to the basis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Outcome {
    /// A new basis vector was created.
    Added,
    /// The candidate lies in the current space (numerically); nothing added.
    Dependent,
}

#[derive(Debug, Clone)]
struct Pending {
    /// `(I − QQᵀ) A q_j`, kept orthogonal to the current basis.
    p: Vec<f64>,
    /// `‖A q_j‖`, the reference for the dependency test at promotion.
    norm_aq: f64,
    /// `‖p‖` right after its last full orthogonalization. If `p` later shrinks
    /// below `REORTH_ALPHA · norm_ref` (cancellation against new vectors),
    /// its relative orthogonality has degraded by that factor and it is
    /// re-orthogonalized before promotion.
    norm_ref: f64,
    /// `false`: frozen, never promoted (phase 7). `p` is then either the raw
    /// `A q_j` or a remainder that was orthogonal to the basis when frozen. It
    /// is no longer updated, only used for the H entries of later basis
    /// vectors: `H[new, j] = q_newᵀ A q_j = q_newᵀ p` (valid since `q_new` is
    /// orthogonal to the basis at freeze time). Without these entries H would
    /// not be `QᵀAQ`.
    promotable: bool,
}

/// Kahan–Parlett-type criterion for selective re-orthogonalization: another
/// Gram–Schmidt pass is needed when a pass (or later projections) removed
/// more than `1 − REORTH_ALPHA` of the norm. Without it, orthogonality was
/// lost near breakdown (‖QᵀQ − I‖ up to 1 in the ε-sweep, phase 7).
const REORTH_ALPHA: f64 = 0.5;

/// DGKS criterion (Daniel–Gragg–Kaufman–Stewart) for the *second* pass: one
/// classical Gram–Schmidt pass leaves an orthogonality error of about
/// `ε‖v₀‖/‖v₁‖`; if the norm did not drop below `‖v₀‖/√2`, that is at
/// rounding level and the second pass is skipped (phase 7). Otherwise more
/// passes follow (REORTH_ALPHA). Near breakdown the drop is large, so the
/// extra passes always happen there.
const DGKS_ETA: f64 = std::f64::consts::FRAC_1_SQRT_2;

/// Maximal number of Gram–Schmidt passes in one orthogonalization.
const MAX_PASSES: usize = 4;

/// Counters for statistics and tests.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct BasisStats {
    pub jv: usize,
    pub dropped_insert: usize,
    pub dropped_promote: usize,
}

#[derive(Debug, Clone)]
pub(crate) struct Basis {
    n: usize,
    /// Number of basis vectors.
    k: usize,
    /// `q_1..q_k`, contiguous column-major `n × k` (blocked kernels, phase 7).
    q: Vec<f64>,
    /// `h[j]` = column j of H (length = current dimension).
    h: Vec<Vec<f64>>,
    pending: Vec<Option<Pending>>,
    /// Closing row of column j (docs/map.md §9): `Some(c)` once the remainder
    /// of column j was promoted (into `q_c`) or dropped (when the dimension
    /// was `c + 1`); then `A q_j ∈ span(q_0..q_c)` and `H[i, j] = 0` for
    /// `i > c`. `None`: the column is open (pending or frozen).
    closed: Vec<Option<usize>>,
    /// `open_at_entry[m]` = number of open columns `j < m` when `q_m` entered,
    /// = number of nonzeros of row m of H left of the diagonal (§9).
    open_at_entry: Vec<usize>,
    /// See [`Basis::set_keep_remainders`].
    keep: bool,
    /// See [`Basis::set_reorth`].
    reorth: bool,
    pub stats: BasisStats,
}

/// Re-orthogonalization (selective second Gram–Schmidt pass) on by default
/// (decision Helmut, 2026-10-31, after the single-pass experiment). A single
/// classical pass, like the single MGS pass of rowmap.F, can be selected with
/// [`Basis::set_reorth`] (`Rowmap::__reorthogonalize(false)`) for
/// comparisons only. With a single pass, 13 tests fail: at (near) breakdown
/// the normalized rounding remainder is not orthogonal to Q (‖QᵀQ − I‖ = 1,
/// garbage step results), `w₄ ∈ K₃` is no longer detected, and on the paper
/// problems VDPL loses a factor 30 in ERR and HEAT needs 2.5× the steps at
/// tol 1e-7 (docs/report).
pub(crate) const REORTH_DEFAULT: bool = true;

impl Basis {
    pub fn new(n: usize) -> Self {
        Basis {
            n,
            k: 0,
            q: Vec::new(),
            h: Vec::new(),
            pending: Vec::new(),
            closed: Vec::new(),
            open_at_entry: Vec::new(),
            keep: true,
            reorth: REORTH_DEFAULT,
            stats: BasisStats::default(),
        }
    }

    /// Remove all vectors (start of a new step); statistics are kept, and the
    /// storage is reused.
    pub fn clear(&mut self) {
        self.k = 0;
        self.q.clear();
        self.h.clear();
        self.pending.clear();
        self.closed.clear();
        self.open_at_entry.clear();
        self.keep = true;
    }

    pub fn dim(&self) -> usize {
        self.k
    }

    /// Lower profile of H (docs/map.md §9): for each column the last row that
    /// can be nonzero, `closed[j]` or the last row if the column is open.
    #[cfg(test)]
    pub fn profile(&self) -> Vec<usize> {
        let last = self.k.saturating_sub(1);
        self.closed.iter().map(|c| c.unwrap_or(last)).collect()
    }

    /// See the field `open_at_entry`.
    #[cfg(test)]
    pub fn open_at_entry(&self) -> &[usize] {
        &self.open_at_entry
    }

    /// Closing rows (`None` = open), see the field `closed`.
    #[cfg(test)]
    pub fn closed(&self) -> &[Option<usize>] {
        &self.closed
    }

    /// Closing row of column j (`None` = open).
    pub fn closed_row(&self, j: usize) -> Option<usize> {
        self.closed[j]
    }

    /// Column j of H as stored (rows `0..len`, zero below).
    pub fn h_col(&self, j: usize) -> &[f64] {
        &self.h[j]
    }

    /// `‖Σ_{j ≥ first} l_{j−first} p_j‖_wrms` over the pending remainders: the
    /// Galerkin residual (without the factor γh) from the trailing solution
    /// components only (docs/map.md §9.4). With one open column (stage 1) it
    /// is `|l| ‖p‖_wrms`, one pass over the vector.
    pub fn residual_wrms(&self, first: usize, l: &[f64], weights: &[f64]) -> f64 {
        let _g = crate::prof::section("residual");
        let terms: Vec<(f64, &[f64])> = (first..self.k)
            .zip(l)
            .filter_map(|(j, &lj)| {
                self.pending[j].as_ref().map(|p| {
                    debug_assert!(p.promotable, "residual with frozen remainders");
                    (lj, &p.p[..])
                })
            })
            .collect();
        match terms.as_slice() {
            [] => 0.0,
            [(lj, p)] => lj.abs() * wrms(p, weights),
            _ => {
                let mut out = vec![0.0; self.n];
                for (lj, p) in &terms {
                    axpy(*lj, p, &mut out);
                }
                wrms(&out, weights)
            }
        }
    }

    /// Switch the re-orthogonalization (second Gram–Schmidt pass, DGKS /
    /// Kahan–Parlett, selective re-orthogonalization at promotion) on or off.
    /// Default [`REORTH_DEFAULT`].
    pub fn set_reorth(&mut self, on: bool) {
        self.reorth = on;
    }

    #[cfg(test)]
    pub fn q(&self, j: usize) -> &[f64] {
        &self.q[j * self.n..(j + 1) * self.n]
    }

    /// `true` if column `j` has a promotable remainder.
    pub fn is_pending(&self, j: usize) -> bool {
        self.pending
            .get(j)
            .is_some_and(|p| p.as_ref().is_some_and(|p| p.promotable))
    }

    /// Indices of all columns whose remainder is still pending.
    pub fn pending_cols(&self) -> Vec<usize> {
        (0..self.dim()).filter(|&j| self.is_pending(j)).collect()
    }

    /// `H = QᵀAQ` as a dense matrix (tests).
    #[cfg(test)]
    pub fn h_mat(&self) -> Mat {
        let k = self.dim();
        let mut m = Mat::zeros(k);
        for (j, col) in self.h.iter().enumerate() {
            for (i, &x) in col.iter().enumerate() {
                m.set(i, j, x);
            }
        }
        m
    }

    /// `Q l`.
    pub fn combine(&self, l: &[f64]) -> Vec<f64> {
        let _g = crate::prof::section("combine Ql");
        assert_eq!(l.len(), self.dim());
        let mut out = vec![0.0; self.n];
        q_times(&self.q, self.n, self.k, l, &mut out);
        out
    }

    /// `Qᵀw` (projection without insertion, for stages that do not insert).
    pub fn coords(&self, w: &[f64]) -> Vec<f64> {
        let mut c = vec![0.0; self.k];
        qt_v(&self.q, self.n, self.k, w, &mut c);
        c
    }

    /// `Σ_{pending j} l_j p_j` — the part of `A Q l` outside the space.
    /// The Galerkin residual of `(I − γhA) k̃ = w` with `k̃ = Q l` and
    /// `w ∈ K` is `−γh` times this vector (AGENTS.md §3.3). Tests; the
    /// integrator uses [`Basis::residual_wrms`].
    #[cfg(test)]
    pub fn outside(&self, l: &[f64]) -> Vec<f64> {
        let _g = crate::prof::section("residual");
        assert_eq!(l.len(), self.dim());
        let mut out = vec![0.0; self.n];
        for (pj, &lj) in self.pending.iter().zip(l) {
            if let Some(pj) = pj {
                // only up-to-date remainders (the residual is used in stage 1,
                // where all remainders are kept)
                debug_assert!(pj.promotable, "residual with frozen remainders");
                axpy(lj, &pj.p, &mut out);
            }
        }
        out
    }

    /// Classical Gram–Schmidt with selective re-orthogonalization, in blocked
    /// sweeps. Leaves `v ← (I − QQᵀ)v₀` and returns the accumulated
    /// coefficients `Qᵀv₀`, `‖v₀‖` and the final `‖v‖`.
    ///
    /// Pass 1 always; pass 2 only if pass 1 reduced the norm below
    /// `DGKS_ETA·‖v₀‖`; further passes while a pass removes more than
    /// `1 − REORTH_ALPHA` of the norm (near dependence), at most MAX_PASSES.
    /// Each update sweep also computes the next pass's coefficients and the
    /// new norm, so a pass costs about one sweep over Q (phase 7).
    ///
    /// A non-finite `‖v₀‖` means non-finite entries (or overflow of the
    /// unscaled sum of squares, which the caller disambiguates).
    fn orth(&self, v: &mut [f64]) -> (Vec<f64>, f64, f64) {
        let _g = crate::prof::section("orth (CGS, DGKS)");
        let (n, k) = (self.n, self.k);
        let mut c = vec![0.0; k];
        let mut ci = vec![0.0; k];
        let n0sq = qt_v(&self.q, n, k, v, &mut ci);
        let n0 = if n0sq.is_finite() && n0sq > 1e-280 {
            n0sq.sqrt()
        } else {
            nrm2(v)
        };
        if k == 0 || !n0.is_finite() || n0 == 0.0 {
            return (c, n0, n0);
        }
        let max_passes = if self.reorth { MAX_PASSES } else { 1 };
        let mut next = vec![0.0; k];
        let mut prev = n0;
        let mut now = n0;
        for pass in 0..max_passes {
            let _p = crate::prof::section(if pass == 0 {
                "  orth pass 1"
            } else {
                "  orth pass 2+"
            });
            for (a, b) in c.iter_mut().zip(&ci) {
                *a += b;
            }
            // The next pass's coefficients are fused into the update sweep only
            // after pass 1, where a further pass is likely; after later passes a
            // further one is rare, and its coefficients are computed on demand.
            let fuse = pass == 0 && max_passes > 1;
            let nsq = v_minus_qc(
                &self.q,
                n,
                k,
                &ci,
                v,
                if fuse { Some(&mut next) } else { None },
            );
            now = if nsq.is_finite() && nsq > 1e-280 {
                nsq.sqrt()
            } else {
                nrm2(v)
            };
            let eta = if pass == 0 { DGKS_ETA } else { REORTH_ALPHA };
            if now >= eta * prev || pass + 1 == max_passes {
                break;
            }
            if fuse {
                std::mem::swap(&mut ci, &mut next);
            } else {
                qt_v(&self.q, n, k, v, &mut ci);
            }
            prev = now;
        }
        (c, n0, now)
    }

    /// Append the unit vector `qn` (orthogonal to the basis): update the
    /// pending remainders and H, then compute `A qn` and its remainder.
    fn push_unit(
        &mut self,
        qn: Vec<f64>,
        av: &mut dyn FnMut(&[f64], &mut [f64]),
    ) -> Result<(), NonFinite> {
        // new row of H for the existing columns: nonzero exactly in the open
        // columns (closed ones have A q_j in the old span, ⊥ qn)
        self.open_at_entry
            .push(self.closed.iter().filter(|c| c.is_none()).count());
        {
            let _g = crate::prof::section("pending update");
            let passes = if self.reorth { 2 } else { 1 };
            for (col, pj) in self.h.iter_mut().zip(self.pending.iter_mut()) {
                let mut hij = 0.0;
                if let Some(pj) = pj {
                    if pj.promotable {
                        // keep the remainder orthogonal (two passes with reorth)
                        for _pass in 0..passes {
                            let c = dot(&qn, &pj.p);
                            axpy(-c, &qn, &mut pj.p);
                            hij += c;
                        }
                    } else {
                        // frozen: only the H entry q_newᵀ A q_j
                        hij = dot(&qn, &pj.p);
                    }
                }
                col.push(hij);
            }
        }
        // new column
        let mut z = vec![0.0; self.n];
        {
            let _g = crate::prof::section("Jv (user)");
            av(&qn, &mut z);
        }
        self.stats.jv += 1;
        self.q.extend_from_slice(&qn);
        self.k += 1;
        self.closed.push(None);
        if !self.keep {
            // This remainder will never be promoted and no residual is needed:
            // the H column `Qᵀ(A q)` in one sweep (its fused `‖A q‖²` gives the
            // finiteness check), and the raw `A q` is kept frozen for the H
            // entries of later vectors. No update or reorthogonalization sweeps
            // (phase 7).
            let mut col = vec![0.0; self.k];
            let s = qt_v(&self.q, self.n, self.k, &z, &mut col);
            if !s.is_finite() && !all_finite(&z) {
                return Err(NonFinite);
            }
            self.h.push(col);
            self.pending.push(Some(Pending {
                p: z,
                norm_aq: s.sqrt(),
                norm_ref: 0.0,
                promotable: false,
            }));
            return Ok(());
        }
        // ‖A q‖ and the finiteness check come out of the first sweep of `orth`
        let (col, norm_aq, norm_ref) = self.orth(&mut z);
        if !norm_aq.is_finite() && !all_finite(&z) {
            return Err(NonFinite);
        }
        self.h.push(col);
        self.pending.push(Some(Pending {
            p: z,
            norm_aq,
            norm_ref,
            promotable: true,
        }));
        Ok(())
    }

    /// Whether new basis vectors keep a promotable remainder `(I − QQᵀ)Aq`
    /// (needed if it may be promoted later or for the stage-1 residual).
    /// Default `true`. The MAP step switches it off once no promotion can
    /// follow (phase 7).
    pub fn set_keep_remainders(&mut self, keep: bool) {
        self.keep = keep;
    }

    /// No promotion follows: freeze all remainders. They are no longer
    /// updated, but still give the H entries of later basis vectors.
    pub fn freeze_pending(&mut self) {
        for p in self.pending.iter_mut().flatten() {
            p.promotable = false;
        }
    }

    /// Insert a right-hand side `w` (stage vector). Returns the outcome and
    /// the coordinates `Qᵀw` in the basis *after* the operation (exact
    /// representation `w = Q·coords` up to `τ‖w‖` if dependent).
    pub fn insert(
        &mut self,
        w: &[f64],
        av: &mut dyn FnMut(&[f64], &mut [f64]),
    ) -> Result<(Outcome, Vec<f64>), NonFinite> {
        assert_eq!(w.len(), self.n);
        let mut v = w.to_vec();
        let (mut coords, nw, nv) = self.orth(&mut v);
        if !nw.is_finite() && !all_finite(w) {
            return Err(NonFinite);
        }
        if nw == 0.0 {
            self.stats.dropped_insert += 1;
            return Ok((Outcome::Dependent, vec![0.0; self.dim()]));
        }
        if nv <= DEP_TOL * nw {
            self.stats.dropped_insert += 1;
            return Ok((Outcome::Dependent, coords));
        }
        scale(1.0 / nv, &mut v);
        self.push_unit(v, av)?;
        coords.push(nv);
        Ok((Outcome::Added, coords))
    }

    /// Promote the pending remainder of column `j` to a basis vector, or drop
    /// it if `A q_j` already lies in the space. Panics if `j` is not pending
    /// (a scheduling error: the caller skips roles whose vector was dropped).
    pub fn promote(
        &mut self,
        j: usize,
        av: &mut dyn FnMut(&[f64], &mut [f64]),
    ) -> Result<Outcome, NonFinite> {
        assert!(self.is_pending(j), "promote: column is not pending");
        let Pending {
            mut p,
            norm_aq,
            norm_ref,
            ..
        } = self.pending[j].take().expect("checked above");
        // Selective re-orthogonalization: `p` is orthogonal to working
        // precision *relative to `norm_ref`*. If later projections in
        // `push_unit` shrank it below REORTH_ALPHA · norm_ref, its relative
        // orthogonality degraded by that factor: re-orthogonalize. (An
        // unconditional pass here doubled the cost; dropping it entirely lost
        // orthogonality near breakdown.)
        let mut np = nrm2(&p);
        if self.reorth && np < REORTH_ALPHA * norm_ref {
            let (c, _, now) = self.orth(&mut p);
            for (hij, ci) in self.h[j].iter_mut().zip(c) {
                *hij += ci;
            }
            np = now;
        }
        if np <= DEP_TOL * norm_aq {
            self.stats.dropped_promote += 1;
            self.closed[j] = Some(self.k - 1);
            return Ok(Outcome::Dependent);
        }
        scale(1.0 / np, &mut p);
        // column j is still open while q_new enters: H[new, j] = ‖p_j‖ is its
        // last nonzero; it closes right after
        self.push_unit(p, av)?;
        let new = self.dim() - 1;
        self.h[j][new] = np;
        self.closed[j] = Some(new);
        Ok(Outcome::Added)
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// Small deterministic generator (no dependencies).
    pub(crate) struct Rng(u64);
    impl Rng {
        pub fn new(seed: u64) -> Self {
            Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1)
        }
        pub fn next(&mut self) -> f64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            (self.0 >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
        }
        pub fn vec(&mut self, n: usize) -> Vec<f64> {
            (0..n).map(|_| self.next()).collect()
        }
    }

    /// Dense test operator (row-major).
    pub(crate) struct DenseOp {
        pub n: usize,
        pub a: Vec<f64>,
    }
    impl DenseOp {
        pub fn random(n: usize, rng: &mut Rng) -> Self {
            DenseOp {
                n,
                a: rng.vec(n * n),
            }
        }
        pub fn diag(d: &[f64]) -> Self {
            let n = d.len();
            let mut a = vec![0.0; n * n];
            for i in 0..n {
                a[i * n + i] = d[i];
            }
            DenseOp { n, a }
        }
        pub fn apply(&self, x: &[f64], y: &mut [f64]) {
            for (i, yi) in y.iter_mut().enumerate() {
                *yi = dot(&self.a[i * self.n..(i + 1) * self.n], x);
            }
        }
        pub fn norm(&self) -> f64 {
            nrm2(&self.a)
        }
    }

    /// Check all basis invariants; returns the max deviation found.
    pub(crate) fn check_invariants(b: &Basis, op: &DenseOp, dropped_ok: f64) -> f64 {
        let k = b.dim();
        let an = op.norm();
        let mut worst: f64 = 0.0;
        // orthonormality
        for i in 0..k {
            for j in 0..k {
                let e = dot(b.q(i), b.q(j)) - if i == j { 1.0 } else { 0.0 };
                worst = worst.max(e.abs());
            }
        }
        assert!(worst < 1e-13, "‖QᵀQ − I‖ = {worst:e}");
        // H = QᵀAQ and A q_j = Q H[:,j] + p_j
        let hm = b.h_mat();
        let mut z = vec![0.0; op.n];
        for j in 0..k {
            op.apply(b.q(j), &mut z);
            for i in 0..k {
                let e = (dot(b.q(i), &z) - hm.get(i, j)).abs() / an;
                assert!(e < 1e-13, "H[{i},{j}] deviates by {e:e}");
                worst = worst.max(e);
            }
            // frozen columns: H = QᵀAQ is checked above; there is no remainder
            // to reconstruct A q_j from
            if matches!(b.pending.get(j), Some(Some(pj)) if !pj.promotable) {
                continue;
            }
            let col: Vec<f64> = (0..k).map(|i| hm.get(i, j)).collect();
            let mut r = b.combine(&col);
            if let Some(Some(pj)) = b.pending.get(j) {
                axpy(1.0, &pj.p, &mut r);
                for i in 0..k {
                    let e = dot(b.q(i), &pj.p).abs() / an;
                    assert!(e < 1e-13, "pending p_{j} not orthogonal to q_{i}: {e:e}");
                }
            }
            axpy(-1.0, &z, &mut r);
            let e = nrm2(&r) / an;
            let tol = if b.is_pending(j) {
                1e-13
            } else {
                1e-13 + dropped_ok
            };
            assert!(e < tol, "A q_{j} reconstruction error {e:e}");
            worst = worst.max(e);
        }
        worst
    }

    #[test]
    fn arnoldi_invariants() {
        let mut rng = Rng::new(1);
        let op = DenseOp::random(40, &mut rng);
        let mut av = |x: &[f64], y: &mut [f64]| op.apply(x, y);
        let mut b = Basis::new(40);
        let w = rng.vec(40);
        let (o, coords) = b.insert(&w, &mut av).unwrap();
        assert_eq!(o, Outcome::Added);
        assert!((coords[0] - nrm2(&w)).abs() < 1e-14);
        for _ in 0..9 {
            let last = b.dim() - 1;
            assert_eq!(b.promote(last, &mut av).unwrap(), Outcome::Added);
        }
        assert_eq!(b.dim(), 10);
        assert_eq!(b.pending_cols(), vec![9]);
        assert_eq!(b.stats.jv, 10);
        check_invariants(&b, &op, 0.0);
        // Arnoldi ⇒ H upper Hessenberg
        let h = b.h_mat();
        for j in 0..10 {
            for i in j + 2..10 {
                assert!(h.get(i, j).abs() < 1e-13 * op.norm());
            }
        }
    }

    #[test]
    fn map_like_sequence_invariants() {
        // stage-1 Arnoldi, then insert a new rhs, promote out of order
        let mut rng = Rng::new(2);
        let op = DenseOp::random(50, &mut rng);
        let mut av = |x: &[f64], y: &mut [f64]| op.apply(x, y);
        let mut b = Basis::new(50);
        b.insert(&rng.vec(50), &mut av).unwrap();
        for _ in 0..3 {
            let last = b.dim() - 1;
            b.promote(last, &mut av).unwrap();
        }
        let k1 = b.dim(); // 4
        let w2 = rng.vec(50);
        let (o, coords) = b.insert(&w2, &mut av).unwrap();
        assert_eq!(o, Outcome::Added);
        // coordinates reproduce w2
        let mut r = b.combine(&coords);
        axpy(-1.0, &w2, &mut r);
        assert!(nrm2(&r) < 1e-13 * nrm2(&w2));
        b.promote(k1 - 1, &mut av).unwrap(); // A q_{κ1}
        b.promote(k1, &mut av).unwrap(); // A q_{κ1+1}
        assert_eq!(b.dim(), k1 + 3);
        check_invariants(&b, &op, 0.0);
    }

    #[test]
    fn exact_breakdown_invariant_subspace() {
        // A has 3 distinct eigenvalues ⇒ Krylov space of any w has dim ≤ 3
        let d: Vec<f64> = (0..30).map(|i| [1.0, 2.0, 5.0][i % 3]).collect();
        let op = DenseOp::diag(&d);
        let mut av = |x: &[f64], y: &mut [f64]| op.apply(x, y);
        let mut b = Basis::new(30);
        let mut rng = Rng::new(3);
        b.insert(&rng.vec(30), &mut av).unwrap();
        assert_eq!(b.promote(0, &mut av).unwrap(), Outcome::Added);
        assert_eq!(b.promote(1, &mut av).unwrap(), Outcome::Added);
        assert_eq!(b.promote(2, &mut av).unwrap(), Outcome::Dependent);
        assert_eq!(b.dim(), 3);
        assert!(b.pending_cols().is_empty());
        assert_eq!(b.stats.dropped_promote, 1);
        check_invariants(&b, &op, 1e-12);
        // any vector built from A^k w is in the span
        let mut z = vec![0.0; 30];
        op.apply(b.q(2), &mut z);
        let (o, _) = b.insert(&z, &mut av).unwrap();
        assert_eq!(o, Outcome::Dependent);
        assert_eq!(b.stats.dropped_insert, 1);
    }

    #[test]
    fn zero_and_dependent_insert() {
        let mut rng = Rng::new(4);
        let op = DenseOp::random(20, &mut rng);
        let mut av = |x: &[f64], y: &mut [f64]| op.apply(x, y);
        let mut b = Basis::new(20);
        let (o, c) = b.insert(&[0.0; 20], &mut av).unwrap();
        assert_eq!((o, c.len()), (Outcome::Dependent, 0));
        let w = rng.vec(20);
        b.insert(&w, &mut av).unwrap();
        let mut w3 = w.clone();
        scale(-3.0, &mut w3);
        let (o, c) = b.insert(&w3, &mut av).unwrap();
        assert_eq!(o, Outcome::Dependent);
        assert!((c[0] + 3.0 * nrm2(&w)).abs() < 1e-13 * nrm2(&w));
        assert_eq!(b.dim(), 1);
    }

    #[test]
    fn near_dependent_insert_keeps_orthonormality() {
        // w2 = w1 + ε g: accepted down to ~τ, always orthonormal
        let mut rng = Rng::new(5);
        let op = DenseOp::random(30, &mut rng);
        let w1 = rng.vec(30);
        let g = rng.vec(30);
        for e in [1e-2, 1e-6, 1e-10, 1e-12, 1e-13, 1e-15, 1e-17] {
            let mut av = |x: &[f64], y: &mut [f64]| op.apply(x, y);
            let mut b = Basis::new(30);
            b.insert(&w1, &mut av).unwrap();
            let mut w2 = w1.clone();
            axpy(e, &g, &mut w2);
            let (o, _) = b.insert(&w2, &mut av).unwrap();
            if e >= 1e-12 {
                assert_eq!(o, Outcome::Added, "ε = {e:e}");
            }
            if e <= 1e-15 {
                assert_eq!(o, Outcome::Dependent, "ε = {e:e}");
            }
            check_invariants(&b, &op, 1e-12);
        }
    }

    /// Diagnostic (run with `--ignored --nocapture`): where the dependency
    /// test flips, and orthonormality across the threshold zone.
    #[test]
    #[ignore]
    fn threshold_zone_report() {
        let mut rng = Rng::new(6);
        let op = DenseOp::random(200, &mut rng);
        let w1 = rng.vec(200);
        let g = rng.vec(200);
        for k in 0..=40 {
            let e = 10f64.powf(-8.0 - 0.25 * k as f64);
            let mut av = |x: &[f64], y: &mut [f64]| op.apply(x, y);
            let mut b = Basis::new(200);
            b.insert(&w1, &mut av).unwrap();
            for _ in 0..5 {
                let l = b.dim() - 1;
                b.promote(l, &mut av).unwrap();
            }
            let mut w2 = w1.clone();
            axpy(e, &g, &mut w2);
            let (o, _) = b.insert(&w2, &mut av).unwrap();
            let worst = check_invariants(&b, &op, 1e-12);
            println!(
                "eps = {e:9.2e}  {o:?}  dim = {}  worst invariant dev = {worst:.1e}",
                b.dim()
            );
        }
    }

    #[test]
    fn non_finite_detected() {
        let mut av = |_x: &[f64], y: &mut [f64]| y.fill(f64::NAN);
        let mut b = Basis::new(3);
        assert_eq!(b.insert(&[1.0, 0.0, 0.0], &mut av), Err(NonFinite));
        let mut b = Basis::new(3);
        let mut ok = |x: &[f64], y: &mut [f64]| y.copy_from_slice(x);
        assert_eq!(
            b.insert(&[f64::INFINITY, 0.0, 0.0], &mut ok),
            Err(NonFinite)
        );
    }
}
