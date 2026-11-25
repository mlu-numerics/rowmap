# rowmap-ng — ROWMAP rewritten in Rust

A clean, from-scratch Rust implementation of ROWMAP, a Krylov-W method of
order 4(3) for large stiff ODEs `y' = f(t, y)`. It is based on the ROW method
GRK4T and solves the stage equations in a special *Multiple Arnoldi Process*
(MAP), which preserves order 4 for small Krylov dimensions independent of n.

Sources in this repository:

- `references/rpha9639.ps` — Weiner, Schmitt, Podhaisky: *ROWMAP – a ROW-code
  with Krylov techniques for large stiff ODEs* (Report 39, Halle 1996).
  **The paper is the specification.**
- `fortran/rowmap.F` — the original 1997 Fortran code. It is a *hint*, not a
  reference: it is known to be not 100 % correct (see "Known Fortran issues").
  Do not modify `rowmap.F`; it is the baseline for the comparison runs.

Status (2026-10-31): phases 0–7 done (5b, the order-3 method, is still
optional and undecided).
- Phase 1: `docs/fortran-paper-runs.md`.
- Phase 2: core numerics in `src/linalg.rs`, `src/grk4t.rs`, `src/dense.rs`,
  `src/krylov.rs`.
- Phase 3: MAP step (`src/map.rs`) with order, breakdown and ε-sweep tests;
  `docs/map.md` (autonomous); Lean proof in `lean/`.
- Phase 4: public integrator `Rowmap` (`src/solver.rs`). Rust test problems
  in `tests/problems/`; results in `docs/rust-runs.md`.
- Phase 5: non-autonomous MAP (`docs/map.md` §7, review gate passed). Schedule
  `GRK4T_NONAUTONOMOUS` (4/9/10/15) is used by the integrator; DIFKON and
  DIFFU2 are in `docs/rust-runs.md` §4.
- Phase 6: `docs/comparison.md` (summary of all paper experiments).
- Phase 7: crate docs, README, performance (2.1–2.9× Fortran CPU).

Possible next steps: 5b (order-3 method); FOOD3 (if the true coefficients
appear); the open questions below.

Phase 5b (2026-10-31): the order-3 method is **ROS34PW2**, public as
`Method::Ros34pw2` (decision Helmut), with GRK4T still the default. See
`docs/map.md` §8.
- **Minimal Krylov dimensions for classical order 3:** 3/4/4/4 (autonomous)
  and 3/4/8/8 (non-autonomous). Confirmed by order tests (one prediction
  wrong: inserting `w₄` is not needed) and in Lean, including tightness.
- **ROS3P was examined first and dropped:** its embedded estimate vanishes
  on linear problems. Its code was removed; one paragraph in `docs/map.md`
  §8 records why.
- The step-size controller and the initial step use the method's orders:
  exponent 1/(q+1) and 1/(p+1).
- Work-precision report: `docs/report/report.tex` (LuaLaTeX, `latexmk
  -lualatex`). Figures come from `uv run plots/work_precision.py --native-only
  --out docs/report/fig`, the equal-accuracy table from `uv run
  tools/wp_ratio.py --tex docs/report/ratios.tex`.
- Performance: the kernels pick AVX2/FMA at run time. Blocked CGS uses
  selective DGKS re-orthogonalization, and remainders are frozen when no later
  closure round needs them. Since then, an incremental Givens QR. At equal ERR,
  Rust GRK4T takes 0.92–1.56× the Fortran CPU time. The rest of the gap is
  algorithmic: the second GS pass, and the full schedule with no early stage
  stop. Profile with
  `cargo run --release --features profile --example profile -- <problem> <tol>`.
- Single GS pass (decision Helmut: two-pass stays the default). The single
  pass is kept only for comparison (`REORTH_DEFAULT`, experiment flag
  `--single-pass`). With it, Rust is 0.75–0.96× Fortran CPU *per J·v* on 6 of
  7 problems (NILIDI 1.09; `uv run tools/jv_cost.py`, report §5). Single pass
  fails 13 tests: breakdown gives garbage, and `w₄ ∈ K₃` is no longer
  detected. On the paper problems, VDPL is 30× worse in ERR and HEAT needs
  2.5× the steps at tol 1e-7.

**Open for Helmut:**
1. `docs/map.md` was reviewed as "plausible". Still open: a rigorous
   derivation of the omission orders (§4, §7.4).
2. The schedule question: the stage-4 closure is not needed for classical
   order 4. The data in `docs/rust-runs.md` §2 shows equal step counts,
   ERR within a factor of about 2, and 10–50% less CPU.
3. Exact-breakdown dimensions (§3.4 "Finding").

Working problem set: BRUSS, HEAT, NILIDI, FOOD, VDPL (autonomous), DIFKON,
DIFFU2 (non-autonomous). **FOOD3 is parked**: the printed coefficients are
suspicious and source [3] is unavailable. The extra problem ROTADV (GRK4T vs
LSTAB) was dropped after two tries, because LSTAB showed no advantage.

Quick reference:
- `make -C fortran/experiments`: build the drivers.
- `fortran/run_experiments.sh [paper|coarse|all]`: run all Fortran experiments
  and write `tests/data/fortran/*.toml`.
- `uv run tools/crosscheck.py`: check the Python and Fortran discretizations
  against each other.
- `uv run tools/reference.py <problem> [m]`: compute a reference solution
  (one-off).
- `uv run plots/work_precision.py`: plot to `plots/out/`.
- `uv run tools/fortran_tables.py`: generate markdown tables.
- `cargo run --release --example experiments -- [paper|coarse|all] [--reduced]`:
  Rust paper experiments, written to `results/*.toml`.
- `uv run tools/compare_runs.py [--reduced-vs-paper] [--markdown FILE]`:
  compare Rust with Fortran, or the two schedules.
- ε-sweep vs Fortran: run `cargo test sweep_export -- --ignored`, then
  `fortran/experiments/build/sweep results/sweep_input.txt results/sweep_f77.txt`,
  then `uv run plots/breakdown_sweep_f77.py`.
- `cargo test`: all Rust tests. Diagnostics run with
  `cargo test <name> -- --ignored --nocapture` (`threshold_zone_report`,
  `breakdown_remainder_report`, `sweep_report`; the last one writes
  `results/sweep_autonomous.toml`, plotted by `uv run plots/breakdown_sweep.py`).
- `cd lean && lake build`: Lean proofs. `lake env lean --run DumpSchedule.lean
  schedule.toml` re-exports the schedule. `rg -n "sorry|admit|^axiom"
  lean/RowmapLean` must be empty. See `lean/README.md`.

---

## 1. Decisions (from the interview)

| Topic | Decision |
|---|---|
| Language / numbers | Rust, `f64` only, no generics over the scalar type. |
| Dependencies | **None** in the library. Hand-written BLAS-1/2 kernels on `&[f64]`. Dev-deps for tests are OK only if really needed. |
| Problem description | A **trait** (`OdeSystem`). The user always supplies the Jacobian-vector product. There is **no** finite-difference J·v. |
| ∂f/∂t | Non-autonomous problems **must** supply `dfdt`. No finite differences. |
| Autonomy | Part of the problem (`autonomous()`), default `false`. It sets the Krylov increments: autonomous `+3,+1,+3`; non-autonomous `+5,+1,+5`. |
| Order-3 method | **ROS34PW2**, public as `Method::Ros34pw2` (decision Helmut, 2026-10-31). Its own derived and Lean-checked minimal Krylov dimensions: 3/4/4/4 (autonomous), 3/4/8/8 (non-autonomous), `docs/map.md` §8. ROS3P was examined and dropped (estimate vanishes on linear problems). |
| Base method | `enum Method { Grk4t (default), Ros34pw2 }`. LSTAB was added as an option and then dropped (no advantage, §3.1). The coefficients live in one place (`src/grk4t.rs`). |
| Output | `solve` returns `Result`. An optional per-step callback gets a `Step` view with order-3 Hermite dense output and can stop the integration. |
| Options | Primary: `rtol`, `atol` (each a scalar or a vector). Secondary, through builder methods with defaults: `max_krylov`, `max_steps`, `initial_step`, `max_step`. Fixed constants: step factors 0.25 / 2 / 0.8, θ = 0.1. |
| Errors | `Result<_, RowmapError>` (enum). The library never prints. |
| Near-breakdown | Two orthogonalization passes (CGS2 / "twice is enough"), plus one documented relative dependency criterion for **all** candidate vectors (see §3.4). A dropped vector counts as "done" and the stage dimension shrinks. Nothing gets replaced. |
| Robustness vs speed | **Robust and correct first.** Speed is secondary for now. Prefer simple dense small-matrix solves over clever incremental updates. |
| Fortran role | **Benchmark, not oracle.** Store the Fortran endpoint values and statistics. Rust must reach comparable accuracy against an independent high-accuracy reference, with similar step and evaluation counts. Deviations get investigated and documented; they do not auto-fail tests. |
| Problem sizes | Coarse grids in `cargo test`; the paper's sizes in `examples/`. |
| Order-test hooks | Fixed step size and forced minimal Krylov dimensions are **internal only** (`pub(crate)` / `#[doc(hidden)]`). |
| Breakdown tests | Exact-breakdown cases **and** an ε-continuity sweep (§5.3, items 3–4). |
| Formal proof | **Lean 4 + Mathlib** in `lean/`. It proves the *structural* part of the order argument: the MAP schedule, including the skip rule, yields the required subspace identities. The analytic step "identities ⇒ order 4" is proved on paper in `docs/map.md` and enters Lean as **one named, documented hypothesis**. Autonomous proof alongside phase 3. The non-autonomous proof is part of the phase-5 review gate. See §3.8. |
| Schedule single source | The MAP schedule is written in both Rust (`const`) and Lean. A test compares the Lean-exported table with the Rust one. |
| Crate | `rowmap` at the repository root. MIT license (as the original). Git repository. |

---

## 2. Public API (target sketch)

```rust
pub trait OdeSystem {
    /// Dimension n.
    fn dim(&self) -> usize;
    /// f(t, y) -> dy
    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]);
    /// out = J(t, y) · v with J = ∂f/∂y (exact; the order relies on it).
    fn jac_vec(&self, t: f64, y: &[f64], v: &[f64], out: &mut [f64]);
    /// true ⇒ f does not depend on t. Enables smaller Krylov increments.
    /// Declaring true for a t-dependent f gives wrong results.
    fn autonomous(&self) -> bool { false }
    /// ft = ∂f/∂t (t, y). Required unless `autonomous()` is true.
    /// The default returns `Err(Unsupported)`, and the solver then fails
    /// up front with `RowmapError::MissingDfdt`. It never fails mid-run.
    fn dfdt(&self, t: f64, y: &[f64], ft: &mut [f64]) -> Result<(), Unsupported> { Err(Unsupported) }
}

pub enum Tol { Scalar(f64), Vector(Vec<f64>) }   // From<f64>, From<Vec<f64>>

let mut solver = Rowmap::new(rtol, atol);        // impl Into<Tol>
solver.max_krylov(70).max_steps(100_000);        // optional
let out = solver.solve(&sys, t0, &mut y, tend)?; // y updated in place
let out = solver.solve_with(&sys, t0, &mut y, tend, |s: &Step| {
    s.interpolate(t_out, &mut buf);              // Hermite, order 3
    ControlFlow::Continue(())
})?;
// out: Solution { t, h_next, status: Finished | Interrupted, stats: Stats }
```

`Step` exposes `t_old, t_new, y_old, y_new, f_old, f_new` and `interpolate(s, out)`.

`Stats`: steps (accepted, rejected), f evaluations, J·v products, ∂f/∂t
evaluations, and Krylov dimensions per stage (min / max / mean). It also counts
dropped (dependent) vectors per kind: insertion vs. Arnoldi, and stage-4
insertions in particular.

`RowmapError`: `InvalidInput(..)`, `MissingDfdt`, `StepSizeTooSmall { t, h }`,
`MaxStepsExceeded { t }`, `SingularStageMatrix { t }` (singular after 3
consecutive step halvings), `NonFinite { t }` (NaN/Inf in f or J·v). On error
the caller still has the last accepted `y`; the error carries `t`.

---

## 3. Algorithm

### 3.1 ROW scheme (GRK4T, Hairer–Wanner ROS4 coefficients)

Stage `i = 1..4`, with `c_ij = γ_ij/γ` and `k̃_i = k_i + Σ_{j<i} c_ij k_j`:

```
(I − γhA) k̃_i = w_i,   w_i = f_i + Σ_{j<i} c_ij k_j  [+ h d_i f_t  if non-autonomous]
f_i = f(t + α_i h, u + h Σ_{j<i} a_ij k_j),   A = f_y(t_m, u_m)
u_{m+1} = u + h Σ b_i k_i,   err = h Σ (ê_i − b_i) k_i
```

Facts the code relies on: `b₃ = 0` and `f₄ = f₃` (a₄ⱼ = a₃ⱼ,
so only 3 f evaluations per step plus one at the new point). Check them in a
unit test. In the same test, verify the order conditions (order 4, embedded
order 3) from the coefficients themselves; do not trust the Fortran table
blindly.

Stability, computed from the `rowcoe` coefficients (2026-10-31; reproduce
this in a unit test):

| Method | γ | \|R(∞)\| | max \|R(iy)\| | A(α) | embedded \|R̂(∞)\| / max \|R̂(iy)\| |
|---|---|---|---|---|---|
| GRK4T | 0.231 | 0.45 | 1.028 (y≈4.4) | ≈89.35° | 2.60 / 2.60 |

GRK4T is not L-stable and is slightly unstable near the imaginary axis.
**Known limitation, documented, accepted.** An LSTAB option was considered
and dropped (2026-10-31). Even on a purely skew-symmetric advection problem
with discontinuous data, LSTAB had no work-precision advantage
(`docs/fortran-paper-runs.md`). In the Krylov-W method, stability also depends
on the stage-1 residual test (2.4) and on κ₁. The embedded method's
|R̂(∞)| ≈ 2.6 makes the error estimate pessimistic on stiff components. That
costs efficiency, not correctness.

### 3.2 Krylov-W approximation

`k̃_i = Q_i (I − γh H_i)⁻¹ Q_iᵀ w_i` with `H_i = Q_iᵀ A Q_i`. This is valid only
if `w_i ∈ K_i = range Q_i` (Lemma 2.1(1)). The design keeps that as an
**explicit invariant**.

### 3.3 Krylov basis data structure (replaces Fortran `kl` / `asv` / `reortn` bookkeeping)

```
Basis {
  q:       Vec<Vec<f64>>,      // orthonormal basis q_1..q_κ
  h:       dense κ×κ,          // H = QᵀAQ, every column known
  pending: map j -> p_j,       // p_j = (I − QQᵀ) A q_j, not yet promoted
}
```

- Every basis vector gets `A q_j` computed right when it is added (one J·v), so
  `H` is always the full `QᵀAQ`. The remainder `p_j` stays *pending* until it is
  promoted.
- **Adding** a vector `q_new` (an inserted w or a promoted p): for every pending
  `p_j`, set `H[new, j] = q_newᵀ p_j` and `p_j ← p_j − H[new, j] q_new`, with a
  re-orthogonalization pass. Then compute `A q_new` and orthogonalize it, which
  gives column `new` of `H` and a new pending `p_new`.
- **Promoting** `p_j`: normalize it into a new basis vector. `H[new, j] = ‖p_j‖`
  has already been accounted for.
- **Galerkin residual, exact and cheap:** `r = (I − γhA) Q l − w = −γh Σ_{pending j} l_j p_j`.
  Compute its weighted RMS norm from the actual vector, not from the Fortran
  triangle-inequality bound.
- **Small system** `(I − γh H) l = Qᵀw`: build it densely and solve with LU with
  partial pivoting (or Householder QR), from scratch for each stage and each
  stage-1 iteration. κ ≤ ~70, so the O(κ³) cost is negligible next to n·κ. A
  pivot below `uround·‖I − γhH‖` counts as singular: reject the step and halve
  h.
  **Correction (2026-10-31, Helmut):** this ignores the point of the algorithm.
  H is *generalized Hessenberg*: each vector entering the basis has a fixed
  structural consequence for H, which is proved in Lean
  (`HStructure.lean`, `docs/map.md` §9). The Fortran incremental Givens QR
  (`kryarn` / `ztocs`) exploits exactly this. It gives one factorization per
  step, and a stage-1 residual without back substitution. **Done:**
  `src/givens.rs` (incremental Givens QR, driven by the closing rows that the
  basis records). The dense LU stays as the test oracle (`givens::tests`).

### 3.4 Dependency ("lucky breakdown") criterion — the critical part

The same routine handles every candidate: inserted `w_i`, and promoted `p_j`
(themselves products of orthogonalization).

```
v0 = candidate
v1 = v0 − Q(Qᵀv0)          (pass 1, CGS or MGS)
v2 = v1 − Q(Qᵀv1)          (pass 2)
dependent  ⇔  ‖v2‖ ≤ τ · ‖v0‖
otherwise accept q = v2/‖v2‖ (both passes' coefficients go into H or l)
```

- One documented τ for everything. Start near `τ = 1e3·ε`. The sweep near
  "√ε range" from the interview is a tuning *upper* bound to explore, not a
  default. Rationale: after two passes, an accepted vector is orthonormal to
  working precision even when it is mostly noise, and an orthonormal noise
  direction does **not** break the Galerkin formula. What breaks things is a
  vector that is *not* normalized or *not* orthogonal (the Fortran
  `h(m,j) > 1e-14` skip leaves an unnormalized column in Q). So drop only
  vectors that are numerically in the span. If the ‖v2‖/‖v1‖ ratio shows
  further severe cancellation, run a third pass or drop. Decide this when
  testing, and document it here.
- **Finding (phase 3): exact breakdown is fuzzy in floating point.** For
  `diag(−1, −10, −100)` the remainder at the theoretical breakdown is about
  `5e-12·‖A q₃‖`, above τ. It is not first-pass rounding but *forward error*:
  the computed q₃ is off the invariant subspace by ~1e-14, because the earlier
  Krylov ratios (~0.08) amplify rounding. So `A q₃` really does leave the
  space by that much.
  - Kept: τ stays at `1e3·ε`. The extra vectors are orthonormal and harmless.
    The result equals exact ROW to 1e-13. The stage-1 residual is at rounding
    level, so stage 1 stops at `KAPPA1_MIN`, and dims are `[4,5,5,6]` instead
    of `[3,3,3,3]`.
  - The ε-sweep shows full continuity: no jump at any dimension switch,
    orthogonality ≤ 1e-15.
  - Possible refinement (not done; decision Helmut): track a per-vector
    forward-error estimate and drop remainders below their own noise
    estimate. That would give clean dimensions at exact breakdown, with the
    same accuracy.
- For a `p_j`, "dependent" means `A q_j ∈ K` already. Its H column is complete,
  and it is removed from `pending`. **Invariant to maintain and prove:** after
  stage i, for every basis vector `q_l` with `l ≤ κ_{i−1}`, `A q_l` is in `K_i`
  (exactly, or dropped as dependent). This is Lemma 2.1(2) in the version that
  tolerates breakdown.
- When `w_i` is dependent, `Q_iᵀw_i` already represents it. Record it in
  `Stats`. In the autonomous case, `w₄` is **expected** to be dependent
  (`f₄ = f₃`); with the generic insertion routine this happens automatically.
- `w₁ = 0` (e.g. at an equilibrium) gives `k₁ = 0`. Handle it without dividing
  by zero.
- Breakdown at stage 1 before κ₁ = 4 means an invariant subspace, so the stage
  solution is exact. Stop the iteration. This is not an error.
- **Proof obligation (roadmap phase 3 / 5):** a written argument, in
  `docs/map.md`, that the skip rule preserves the conditions used in Thm 2.1:
  `w_i ∈ K_i` and `(T_i^l − A^l) w_k = 0` for the required `l`. Do it for the
  autonomous MAP and for the non-autonomous MAP.

### 3.5 MAP schedule (autonomous, paper §2)

| Stage | Actions (each one may be dropped as dependent) | κ |
|---|---|---|
| 1 | Arnoldi from `w₁` until residual test (2.4) holds, with `κ₁ ≥ 4` and `κ₁ ≤ max_krylov − 7` | κ₁ |
| 2 | insert `w₂`; promote `p(q_κ₁)`; promote `p(q_κ₁₊₁)` | ≤ κ₁+3 |
| 3 | insert `w₃` | ≤ κ₂+1 |
| 4 | promote `p(q_κ₁₊₂)`, `p(q_κ₁₊₃)`, `p(q_κ₃)`; `w₄ ∈ K₃` expected | ≤ κ₃+3 |

The actions are referenced by *role*, not by index. When a vector does not
exist because it was dropped, its follow-up promotion is skipped. The schedule
is data: a `const` table per case. That keeps the dimension rules reviewable,
and the order tests can use it directly. Write the machinery method-generic:
the number of stages, the stage roles and the coefficients come from the
method definition, not from hard-coded GRK4T assumptions. The `f₄ = f₃`
shortcut, for example, is a property flag of the method. That keeps phase 5b
(order 3) a matter of adding data, a derivation and a proof.

Stage-1 stopping (paper 2.4): `h · ‖r‖_wrms ≤ θ = 0.1`, with
`scal_i = atol_i + rtol_i |u_i(t_m)|`. Check this against the Fortran
`dfkt < ktol/h` once the scalings line up, and document any difference. If the
cap is reached first, continue. Order 4(3) still holds; the step-size control
takes over (paper p. 7).

### 3.6 Non-autonomous MAP (to be derived; phase 5)

The paper states only the increments `+5, +1, +5` (minimal dimensions
4, 9, 10, 15) and refers to [21] for the construction. The Fortran code inserts
`w_i` whenever it is not in the span, including at stage 4. There
`w₄ − w₃ ≡ h(d₄−d₃) f_t (mod K₃)`, and for GRK4T `d₃ ≈ 0.551`,
`d₄ ≈ −0.055`, so `w₄ ∉ K₃` in general. The task:

1. Derive the vector sequence with the Taylor argument of Thm 2.1, including
   `f_t` terms. A likely ingredient is including `f_t`, or its component
   orthogonal to `w₁`, explicitly in the basis.
2. Write it up in `docs/map.md`, with exactly which `A`-powers are needed for
   which vectors and why the counts come out as 5/1/5.
3. **Review gate:** Helmut reviews the derivation before implementation.
4. Verify it with fixed-step order tests at minimal dimensions. As a negative
   test, show that one dimension fewer lowers the order.

### 3.8 Formal verification in Lean (structural part)

**Split of the proof of Thm 2.1:**

- **(A) Analytic, on paper (`docs/map.md`).** A precise list `L` of
  conditions such that, if they hold for the stage vectors, the Krylov-W
  solution has local error O(h⁵), and the embedded solution O(h⁴). The
  conditions are memberships like `Aˡ w_k ∈ K_i`, equivalently
  `(P_i A)ˡ w_k = Aˡ w_k` for any projection `P_i` onto `K_i`. They also
  include `w_i ∈ K_i` and `A(K_{i−1}) ⊆ K_i`, i.e. Lemma 2.1(2). Derive `L`
  from the Taylor argument of Thm 2.1 for autonomous and non-autonomous
  problems separately. In the Lean code `L` appears as the hypothesis
  `OrderFromIdentities` with a doc comment pointing to the paper proof. **`L`
  must be minimal:** each item is justified by a term of the expansion.
- **(B) Structural, in Lean.** For any field `𝕜`, vector space `V`, linear map
  `A : V →ₗ[𝕜] V`, and arbitrary vectors `w₁..w₄` (and `f_t`): running the MAP
  schedule produces subspaces `K₁ ⊆ … ⊆ K₄` with
  `dim K_i ≤ κ₁ + Δ_i` (4/7/8/11 resp. 4/9/10/15) that satisfy every item of
  `L`.
  - Subspace-level model: an action adds `span{v}` to the current space.
    "Dropped" ⇔ `v ∈ K`, which leaves the space unchanged. Model exactly the
    Rust semantics of which follow-up actions are skipped when a vector is
    dropped (roles, §3.5).
  - The model should be independent of orthogonalization and of the choice of
    representatives modulo `K`. If the proof needs a specific representative,
    i.e. it really needs the inner product, that is a finding: document it,
    and switch to `EuclideanSpace ℝ (Fin n)`.
  - Stage 1 is modelled as `K₁ = span{w₁, A w₁, …, A^{κ₁−1} w₁}` with κ₁ ≥ 4,
    including early invariance (breakdown), where `A K₁ ⊆ K₁`.
- **(C) Structural necessity, in Lean.** For each stage increment, one
  concrete counterexample over ℚ (explicit matrix and vectors, discharged by
  `decide` / `norm_num`) where one action fewer violates an item of `L`. Its
  effect on the order is shown by the negative numerical order tests (§5.3).
- **Not covered:** floating point, the τ threshold (covered by the ε sweep),
  and the stage-1 residual test (which only affects stability).

**Setup:** a `lean/` Lake project pinned to a fixed Mathlib release, with
`lean-toolchain` committed and a ~5 GB Mathlib cache, via `lake exe cache get`
(elan is in `~/.elan`). `lake build` must succeed with **no `sorry`**. Check
this with a script: `rg sorry lean/` must be empty.

*Implemented (phase 3):* `OrderFromIdentities` is **not** stated as a Lean
`Prop`. Without formalizing the ODE setting, any placeholder would be vacuous.
The boundary is documented in `Map.lean` instead. The main theorems depend
only on `propext`, `Classical.choice` and `Quot.sound` (checked with
`#print axioms`).

**Rust↔Lean consistency:** `lake env lean --run DumpSchedule.lean
schedule.toml` writes `lean/schedule.toml`, which is committed. It runs in the
interpreter, because a `lake exe` would compile all of Mathlib natively. The
Rust test `schedule_matches_lean` asserts that the Rust `const` schedule
tables match it. Regenerate it whenever the schedule changes.

### 3.7 Step-size control and the rest

- `err = ‖h Σ(ê_i − b_i) k_i‖_wrms` with `scal` from `u_m`;
  `h_new = h · min(2, max(0.25, 0.8 · err^{-1/4}))`. Accept if `err < 1`.
  After a rejection, do not allow growth.
- Initial step: Hairer–Wanner style `hinit` (2 f evaluations), not the
  dimensionally odd Fortran heuristic. Then `min(h, max_step)`.
- The last step lands exactly on `tend`: stretch when `t + 1.005h ≥ tend`. Fail
  with `StepSizeTooSmall` when `0.1|h| ≤ ε|t|`.
- `f(u_{m+1})` is reused as `f₁` of the next step and for the Hermite output.
  f_t is evaluated once per step at `(t_m, u_m)`.
- Check that all f and J·v outputs are finite, and fail with `NonFinite` if not.
  Optionally, a rejection plus halving first.

---

## 4. Known Fortran issues (to verify; the Rust code must not copy them)

1. Non-autonomous stage 4: ad-hoc insertion of `w₄` (`nrmv/nrmrhs ≥ 1e-8`),
   with no theory behind it.
2. Inconsistent thresholds: `1e-14` (reortn scaling, which leaves the column
   unnormalized), `1e-8` (stage-4 insertion), `uround` (singularity).
   Orthogonalization is a single MGS pass only.
3. Stage-1 cap reserves 12 / 16 vectors instead of the 7 / 11 the paper needs.
4. "Repeatedly singular" returns `idid = −5`; the documentation says `−4`.
5. The `solout` interrupt message goes to unit `*` instead of `lun`.
6. `driver.F` "average krylov dim." is `njacv / nsteps`, not the mean dimension.
7. The residual norm in `krdfkt` is a triangle-inequality bound, not the true
   residual.
8. The initial step formula `0.1 · rms(1/scal)^{-1/4}` is not dimensionally
   consistent.
9. `driver.F` (NILIDI demo): `fex` uses `pi = 3.141592654`, while `jacv` and
   `init` use the full π. So f and J·v belong to slightly different problems.
   The experiment drivers use exact π everywhere.
10. The default step limit of 10000 (`iwork(1)`) is too small for VDPL at tight
    tolerances. The experiments set 10⁶.

Add to this list whenever the comparison runs reveal more.

---

## 5. Testing strategy

### 5.1 Fortran baseline (`fortran/experiments/`)

- Build with `gfortran -O2 … -lopenblas`. `rowmap.F` needs BLAS only. Use one
  driver per problem, with `ijacv = 1` (exact J·v, to match Rust) and `ifdt = 1`
  where non-autonomous.
- Each run writes `tests/data/fortran/<problem>_<grid>_tol<k>.toml`: the
  settings, the endpoint `y` (or `y` at a few points, for paper grids), the
  stats, and a checksum of `f(t0, y0)`. The checksum catches discretization
  mismatches between the Fortran, Rust and Python versions.
- A script `fortran/run_experiments.sh` regenerates everything. The generated
  data is committed.

### 5.2 Reference solutions

- **Computed once with DOPRI5** (Dormand–Prince 5(4), explicit; the problems
  are only mildly stiff), with `atol = rtol = 1e-10`. Use SciPy `RK45` (which
  is DOPRI5) via `uv run tools/reference.py`. Run it once, and **commit** the
  result to the repo. It is regenerated only when a problem definition
  changes.
- Format: `tests/data/reference/<problem>_<grid>.txt`, plain text: a `#`
  header (problem, grid, n, t_end, solver, tol, date), then one value
  `y_end[i]` per line in `es24.16e3` format. Fortran list-directed `read`,
  Rust and numpy all read it trivially. This is *input* data, so it is
  deliberately not TOML.
- If DOPRI5 becomes too expensive on a problem (watch VDPL over t ∈ [0, 1000]
  and HEAT on the paper grid), note it here and decide case by case.
- Error measure as in the paper: `ERR = sqrt(1/n Σ ((u_i − u_ex,i)/(1 + |u_ex,i|))²)`.
  The **producer computes it** (the Fortran driver, the Rust example) by
  reading the reference file.

### 5.2a Data format and plots

- **All result files are TOML.** That covers Rust examples and experiment
  runs, the Fortran drivers and the Lean schedule export, so one reader
  handles everything. The exception is the reference solutions, which are
  plain-text input data (§5.2).
- The library itself stays dependency-free. Examples and tests may use the
  dev-dependencies `serde` + `toml`.
- The Fortran drivers write TOML by hand with `FORMAT`; Lean writes it by hand
  too. Floats are written with full precision: `es24.16e3` in Fortran,
  `{:e}` / round-trip in Rust.
- One schema for every producer (`producer = "rowmap-rs" | "rowmap-f77"`):

  ```toml
  [run]
  producer = "rowmap-rs"   # version/commit, date
  problem = "nilidi"
  grid = [69, 69]
  n = 4761
  method = "grk4t"
  rtol = 1e-4
  atol = 1e-4
  max_krylov = 66

  t_end = 1.0

  [stats]                  # exactly these, from every producer
  err = 3.2e-5             # endpoint error ERR (paper formula) vs reference
  nstep = 7                # number of steps (accepted + rejected, as Fortran iwork(5))
  nfeval = 28              # number of f evaluations
  cpu_time = 0.84          # seconds, integration only
  ```

  No solution vector goes in the output files. More detailed statistics (J·v, Krylov dimensions,
  dropped vectors) remain available in Rust `Stats`. They are printed by
  tests and examples when needed, but are not part of the output files.

  The schema is fixed in `docs/data-format.md`. Its first version goes in
  with phase 1.
- **Plots: Python, in this repository**, under `plots/`. Use uv scripts with
  PEP 723 inline dependencies (matplotlib, numpy), so `uv run plots/<name>.py`
  needs no setup. They read the TOML with `tomllib`.
  - `plots/work_precision.py`: Figs. 3.1–3.6 and 4.7–4.8 style, lg(ERR) vs
    time, overlaying Fortran, Rust and the method variants.
  - `plots/order.py`: log-log local-error slopes from the order tests,
    including the negative cases.
  - `plots/breakdown_sweep.py`: the ε sweep (‖u(ε) − u(0)‖,
    orthogonality, dimensions vs ε).
  - The order and ε-sweep plots read their own small TOML files written by
    the respective tests: `results/order_*.toml` with `h` / `local_err`
    arrays, and `results/sweep_*.toml`.
  - Output goes to `plots/out/*.pdf` + `.png`. The generated images are not
    committed. The `.toml` data is.

### 5.3 Rust test layers

1. **Unit tests:** BLAS kernels; the GRK4T order conditions and `b₃ = 0`,
   `f₄ = f₃`; the small dense solver; orthogonality `‖QᵀQ − I‖ ≈ ε` after every
   basis operation (debug assertions plus tests).
2. **Order tests** (internal fixed-step / fixed-dimension hook): local error
   `∝ h⁵` at the minimal dimensions, 4/7/8/11 autonomous and 4/9/10/15
   non-autonomous, on small nonlinear problems with n ≫ κ and exact J·v.
   **Negative tests:** one dimension less in a stage that matters gives lower
   order. Use the full dimension (κ = n) as a ROW cross-check: the result
   equals plain GRK4T with direct solves.
3. **Exact breakdown:** linear `y' = My + g` with few distinct eigenvalues,
   where `f` lies in a low-dimensional invariant subspace. Also
   `f(y0) = 0`, and the non-autonomous path with `f_t = 0` and `f_t ∥ f`. The
   result must be exact (to rounding) or equal the full ROW result, with
   reduced dimensions recorded in `Stats`, and no error or NaN.
4. **ε-continuity sweep (near breakdown):** perturb each exact-breakdown case
   by `ε·g`, where g is a fixed generic direction, for
   `ε = 1e-16, 1e-15, …, 1e-2`. Run it on one fixed-step step and on a short
   adaptive run.
   - Requirements: the step result `u_{m+1}(ε)` is continuous in ε,
     `‖u(ε) − u(0)‖ ≲ C·ε`, with no jump when a vector switches between
     dropped and kept at the τ threshold.
   - `‖QᵀQ − I‖ = O(ε_mach)` and `‖QᵀAQ − H‖` small, for every ε.
   - The Krylov dimensions change monotonically and only near τ.
   - The local error stays order 4 across the whole ε range.
   - Run the same sweep in both the autonomous and the non-autonomous path
     (`f_t = f + ε·g`). Run the same cases through `rowmap.F` and document its
     behaviour (expected: trouble from the 1e-14 / 1e-8 thresholds).
5. **Integration tests:** the paper problems on coarse grids. Error against the
   reference ≲ tolerance-consistent bound; stats within a documented factor of
   Fortran.
6. **Examples** (`examples/`): paper sizes and work-precision tables for
   tol = 1e-2 … 1e-7. Not run in CI by default.

### 5.4 Paper problems

Autonomous: BRUSSELATOR (N=5000), HEAT (2500), NILIDI (4761, driver exists),
FOOD (5000), FOOD3 (2000), VDPL (800). Non-autonomous: DIFKON (4761),
DIFFU2 (4761). FOOD3 is parked. An extra imaginary-spectrum problem (ROTADV)
was tried and dropped, because LSTAB showed no advantage
(docs/fortran-paper-runs.md). GRK4T is now the only method.

The paper text is partly garbled (PS→text). **Open:** confirm the exact
definitions of DIFKON (exact solution, exponent), FOOD3 (`b₁ = b₂`?,
coefficients), and the grid layouts (inner points; Neumann discretization) with
Helmut or with [16] / [3] / [12]. Record the final definitions in
`docs/problems.md`, because Fortran, Rust and Python must agree exactly.

---

## 6. Repository layout (target)

```
Cargo.toml            # crate `rowmap`, no deps
src/lib.rs            # public API, re-exports
src/system.rs         # OdeSystem, Tol
src/solver.rs         # Rowmap, step loop, step-size control, Stats, errors
src/grk4t.rs          # GRK4T coefficients (single source)
src/krylov.rs         # Basis, orthogonalization, dependency test, MAP schedule
src/dense.rs          # small dense LU / QR
src/linalg.rs         # dot, axpy, nrm2, wrms
tests/                # order, breakdown, problems, fortran comparison
tests/problems/       # shared problem definitions (Rust)
tests/data/           # fortran/*.toml (generated), reference/*.txt (DOPRI5 1e-10, computed once), committed
results/              # Rust example/experiment outputs (*.toml, committed)
plots/                # uv Python scripts (matplotlib) reading TOML; out/ ignored
docs/data-format.md   # TOML schema
examples/             # paper experiments, full size
docs/map.md           # MAP derivations + breakdown proofs
docs/problems.md      # exact test problem definitions
lean/                 # Lake project: MAP structural proofs, schedule.toml export
fortran/              # original code (unchanged) + experiments/
tools/reference.py    # one-off DOPRI5 (SciPy RK45) reference solutions (uv)
```

---

## 7. Roadmap

Each phase ends with green `cargo test`, `cargo clippy -- -D warnings` and
`cargo fmt --check`, plus an update of this file's status line.

- [x] **0. Setup:** `cargo init --lib`, license, `.gitignore`, CI-like `just`
  or shell script (optional).
- [x] **1. Fortran baseline: run *all* paper experiments with `rowmap.F`.**
  Done 2026-10-31: 168 runs (8 problems, ROTADV dropped), all `idid = 1`. sft=4/5 use the patched copy,
  not the `solout` trick.
  Build; drivers for all 8 paper problems, on the paper grids and on coarse
  grids, with tol = 1e-2…1e-7, `rtol = atol = tol`, max Krylov dim 66 as in
  the paper, and GRK4T.
  - Reproduce every work-precision series of Figs. 3.1–3.6 and 4.7–4.8:
    ERR (paper formula) vs CPU time, plus nstep and nfeval. Write them as
    tables in
    `docs/fortran-paper-runs.md`, and note where the shape deviates from the
    paper's figures.
  - Figs. 4.7/4.8 strategies: sft=2 is the code itself. sft=3 is an augmented
    autonomous system in the driver. sft=4 resets `unew(n+1) = tnew` in
    `solout` (rowmap.F copies `uu → u` after `solout`). Caveat: `fm = f(uu)` was already
    evaluated with the uncorrected `uu(n+1)`, so this is only an approximation
    of sft=4. A patched copy gives the clean version. sft=5 needs internal
    changes, so it goes in a *patched copy* under `fortran/experiments/`,
    never in `rowmap.F`, or is skipped with a note.
  - VODPK, DOP853 and ROS4 curves are not reproduced (codes not available).
  - Additionally: LSTAB on an imaginary-spectrum problem (ROTADV). Tried and
    dropped.
  - TOML output (schema in `docs/data-format.md`); first plots
    (`plots/work_precision.py`); confirm the problem
  definitions (§5.4 open points); `docs/problems.md`; DOPRI5 references (1e-10, computed once, committed).
  Record the Fortran anomalies found in §4.
- [x] **2. Core numerics:** linalg kernels, GRK4T coefficients with
  order-condition, `f₄ = f₃` and stability-function tests (§3.1 table), dense LU, `Basis` with the CGS2 dependency test and orthogonality tests.
  Done 2026-10-31.
  - `RowMethod` is method-generic. Stage reuse (`f₄ = f₃`) is *detected*
    from the coefficient rows, not hard-coded.
  - `τ = DEP_TOL = 1e3·ε`. The threshold-zone diagnostic
    (`cargo test threshold_zone_report -- --ignored --nocapture`) shows the
    flip at a relative new component ≈ τ, with ‖QᵀQ − I‖ and ‖QᵀAQ − H‖ at
    about 1e-15 on both sides.
  - `promote` does one extra re-orthogonalization pass before its
    dependency test.
- [x] **3. Autonomous MAP + single step:** schedule table, stage solves, exact
  residual; fixed-step hook; order tests (positive and negative); exact
  breakdown tests and the ε-continuity sweep (used to tune τ); `docs/map.md`
  autonomous part (list `L` + paper proof). **Lean:** set up `lean/`, then the
  autonomous structural theorem (B), the necessity counterexamples (C),
  `schedule.toml` and the Rust consistency test.
  Done 2026-10-31.
  - **Representation-free MAP:** `K₂ = C(K₁+⟨w₂⟩)`, `K₄ = C(K₃+⟨w₄⟩)` with
    `C(S) = S + AS`. "Promote all pending remainders" is exactly C. The
    schedule is generic: `Extend { insert, rounds }`.
  - **Order tests:** order 5 of `u_K − u_ROW` at 4/7/8/11 (observed 4.99).
  - **Omission tests (12 cases):** each needed action lowers the order as
    predicted. **Finding:** the stage-4 closure and the insertion of w₃ are
    *not* needed for classical order 4. Dims 4/7/7/8 still give 4.95.
    Reason: the only new O(h²) direction is f''(f,f), already in K₂ via w₂.
    `docs/map.md` §4 is corrected; the question goes to Helmut.
  - **Exact breakdown and ε-sweep:** see §3.4 "Finding".
  - **Lean:** F1–F5 plus the tightness example (dims 4/7/8/11 attained over
    ℚ). Instead of one counterexample per omitted action, (C) is shown as
    tightness plus the minimality of the closure.
  - **Not done:** the ε-sweep through `rowmap.F` (it needs a Fortran driver
    for the test operator). Moved to phase 4.
- [x] **4. Integrator:** step-size control, hinit, last step, callback +
  Hermite, Stats, errors; integration tests on the autonomous paper problems
  vs reference and Fortran. Also, carried over from phase 3:
  - the ε-sweep through `rowmap.F`;
  - a comparison of the paper schedule (4/7/8/11) with the reduced schedule
    (without the stage-4 closure) on the *stiff* paper problems, where the
    classical expansion is not valid. That gives Helmut data for the
    schedule question.

  Done 2026-10-31 (`docs/rust-runs.md`):
  - The step counts match Fortran within a few percent.
  - ERR is equal or better in most cases (BRUSS/FOOD up to 25× better at
    tight tolerances).
  - CPU is 5–10× slower (phase 7).
  - The reduced schedule leaves step counts equal, changes ERR by a factor of
    about 0.4–3.6, and saves 10–50% CPU.
  - `rowmap.F` is also continuous in the ε-sweep (null result).
  - Non-autonomous problems returned `InvalidInput` until phase 5. They are
    supported since phase 5.
  - Tests: `tests/api.rs`, `tests/paper.rs` (coarse grids vs reference and
    Fortran), and `tests/problems_consistency.rs` (Rust problems vs Fortran
    dumps in `tests/data/dump/`).
- [x] **5. Non-autonomous MAP:** derivation in `docs/map.md` (list `L` with
  f_t terms) → Lean structural proof + counterexamples for 5/1/5 → **review
  gate (Helmut): paper proof of `L` + green `lake build`** → implementation → order tests at 4/9/10/15 (positive and
  negative) → breakdown tests (`f_t = 0`, `f_t ∥ f`) → DIFKON / DIFFU2.
  Status 2026-10-31, up to the gate:
  - **Derivation:** f_t enters w₂ at O(h). That requires a second closure
    round at stage 2 (A²f_t), and w₄ ∉ K₃ requires inserting w₄. Together
    this reproduces the paper's +5/+1/+5.
  - **Order tests:** order 5 at 4/9/10/15. The autonomous schedule gives only
    order 4 on this problem.
  - **Omissions:** all actions are needed except A q_{w₄} (6 of 10 first
    predictions were wrong; a rigorous derivation of the omission orders is
    open).
  - **Breakdown tests** (`f_t = 0`, `f_t ∥ f`, ε-sweep) pass, done ahead of
    the gate via the step hook.
  - **Lean:** structural facts, +5/+1/+5, tightness 4/9/10/15, both
    schedules exported. The planned per-action counterexamples were replaced
    by the tightness argument, as in phase 3.
  - **Review passed** ("looks plausible", 2026-10-31). After the gate:
    integrator wiring (f_t once per accepted point); DIFKON/DIFFU2 match or
    beat Fortran sft = 2 in accuracy with slightly fewer steps
    (`docs/rust-runs.md` §4); tests in `tests/paper.rs` (all 7 problems) and
    `tests/api.rs::nonautonomous_forced_problem`.
- [ ] **5b. Order-3 method (optional, decided later):**
  1. Select a candidate (ROS3P / ROS3Pw / ROS34PW2) from its coefficients:
     order conditions, stability function, embedded method.
  2. Derive the condition list L and the minimal Krylov schedule, autonomous
     and non-autonomous, in `docs/map.md`.
  3. Lean structural proof plus necessity counterexamples.
  4. **Review gate (Helmut).**
  5. Implement it as a second method table, with order-3 tests at the minimal
     dimensions (positive and negative).
  6. Work-precision comparison with GRK4T on the paper problems.
  7. Only then decide whether it becomes a public option.

  Status 2026-10-31: steps 1–3 done for ROS3P (Helmut asked for ROS3P).
  The method table and the order-3 tests at minimal dimensions (positive and
  negative) from step 5 are done as well, via the step hook. Step 3 used a
  tightness example instead of counterexamples, as before. **At the review
  gate (step 4).** Still to do: integrator selection, then steps 6 and 7.
- [x] **6. Paper experiments:** examples at full size, work-precision tables,
  comparison write-up against Fortran (`docs/comparison.md`).
  Done 2026-10-31: `docs/comparison.md` plus `docs/tables/summary.md`
  (`tools/compare_runs.py --summary`). Rust has the same steps as Fortran,
  equal or better ERR (geometric-mean ratio < 1 in 10 of 14 groups), and was
  3.4–8.9× slower at that point (2.1–2.9× after phase 7).
- [x] **7. Polish:** rustdoc with a usage example, API review, then (later)
  performance work: MGS vs CGS2 cost, allocation reuse, optional parallel
  kernels.
  Done 2026-10-31:
  - Crate documentation with a doctest (1D heat equation, checked against
    the exact solution); `README.md`.
  - API review: no changes; "forward integration only" is now documented.
  - **Performance: 2.5–4× faster**, CPU ratio to Fortran now 2.1–2.9×.
    1. The unconditional safety pass in `promote` was replaced by
       **selective re-orthogonalization**: re-orthogonalize only if the
       remainder shrank below ½ of its norm at creation, and in `orth2` add a
       3rd/4th pass if a pass removed more than ½. **Dropping the safety pass
       entirely broke orthogonality near breakdown** (‖QᵀQ − I‖ up to 1 in
       the ε-sweep; the test caught it).
    2. Reductions with 8 independent accumulators (SIMD without fast-math).
    3. `nrm2` uses a plain sum of squares, with the scaled version only on
       overflow/underflow.
  - Not done (minor or not needed): allocation reuse, BLAS, parallel
    kernels, MGS instead of CGS2 (rejected for robustness).

## 8. Working rules for agents

- Edit files with the Edit / Write tools only (no `sed -i`, no interpreter
  heredocs). Use `rg` / `fd`, and `uv` for Python.
- Never change `fortran/rowmap.F`. Put experiments in `fortran/experiments/`.
- Correctness before speed. Every threshold is a named constant with a
  comment explaining why.
- If the paper and the Fortran code disagree, the paper wins. Record the
  disagreement in §4.
- Commit only when asked.
