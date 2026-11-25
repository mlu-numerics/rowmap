# The multiple Arnoldi process (MAP): construction and order proof

Status (2026-10-31): autonomous part §1–6 and non-autonomous part §7 written.
Helmut reviewed them ("looks plausible"), which passed the phase-5 gate. The
non-autonomous schedule is now used by the integrator. Still open: a rigorous
derivation of the omission orders (§4, §7.4).

This document has two parts:
- the *analytic* argument, written on paper: the conditions L imply order 4;
- the *structural* facts that the Lean proof (`lean/`) checks: the MAP
  delivers L.

The Rust implementation realizes the same construction with an orthonormal
basis (§5).

## 1. Setting

Autonomous problem `y' = f(y)`. At `u = u_m`, `A = f'(u)`. GRK4T has 4 stages,
`b₃ = 0`, `a₄ⱼ = a₃ⱼ` (so `f₄ = f₃`), and `a₄₃ = 0`. Stage `i` of the ROW
method (paper (1.2)) in the `k̃` form:

```
(I − γhA) k̃_i = w_i,   k̃_i = k_i + Σ_{j<i} c_ij k_j,
w_i = f(u + h Σ_{j<i} a_ij k_j) + Σ_{j<i} c_ij k_j .
```

The Krylov-W method replaces `A` in stage `i` by `T_i = P_i A`, where `P_i` is
the orthogonal projector onto a subspace `K_i`, with `K₁ ⊆ K₂ ⊆ K₃ ⊆ K₄`.

**Lemma 1 (projected solve).** If `w ∈ K` and `P` is the orthogonal projector
onto `K`, with `T = PA`, then `x = Q(I − γhH)⁻¹Qᵀw` (`H = QᵀAQ`) solves
`(I − γhT)x = w`, and `x ∈ K`. For `γh‖A‖ < 1`,
`x = Σ_{l≥0} (γh)^l T^l w`.

**Lemma 2 (telescoping).** For every `w` and `l ≥ 1`,

```
A^l w − T^l w = Σ_{j=1}^{l} T^{l−j} (I − P) A^j w .
```

Proof: `A^l − T^l = Σ_j T^{l−j}(A − T)A^{j−1}` and `A − T = (I − P)A`. ∎

So `(T^l − A^l) w` is small as soon as the parts of `A^j w` outside `K`,
`(I − P)A^j w` for `j ≤ l`, are small. All conditions below have this form.

## 2. Error propagation

Let `k_i^R` be the stages of the exact ROW method from the same `u`. As in the
paper (2.2),

```
k_i − k_i^R = S_i + O(h) Σ_{j<i} (k_j − k_j^R),   S_i := Σ_{l≥1} (γh)^l (T_i^l − A^l) w_i .
```

The `c_ij k_j` terms appear on both sides of the stage equation. Their errors
therefore enter only through `[(I − γhT_i)⁻¹ − I]`, i.e. with a factor `h`,
and the stage arguments add a further factor `h`.

**Sufficient for order 4** (local error `u_{m+1} − u^R_{m+1} = O(h⁵)`):
`S₁, S₂, S₄ = O(h⁴)` and `S₃ = O(h³)`.

- Because `b₃ = 0`, `k₃` enters `u_{m+1}` only through `k₄`, with a factor
  `h`.
- Because `a₄₃ = 0`, `k₃` does not enter the argument of `f₄`.
- The embedded solution then has order 3 as well.

By Lemma 2 and with `‖T_i‖ ≤ ‖A‖`, it suffices that

```
(I − P_i) A^j w_i = O(h^{4−j})   for i ∈ {1,2,4}, j = 1,2,3,
(I − P_i) A^j w_i = O(h^{3−j})   for i = 3,     j = 1,2,
```

together with `w_i ∈ K_i` (Lemma 1).

## 3. The construction (representation-free)

With the *closure* `C(S) := S + A·S`:

```
K₁ = span{w₁, A w₁, …, A^{κ₁−1} w₁},  κ₁ ≥ 4   (or A K₁ ⊆ K₁: breakdown)
K₂ = C(K₁ + ⟨w₂⟩)
K₃ = K₂ + ⟨w₃⟩
K₄ = C(K₃ + ⟨w₄⟩)
```

This is exactly the paper's MAP:
- `A K₁ ⊆ K₁ + ⟨A q_κ₁⟩`, so stage 2 adds `w₂`, `A q_κ₁` and `A q_κ₁₊₁`.
- Stage 4 adds `A q_κ₁₊₂`, `A q_κ₁₊₃` and `A q_κ₃`, because the images of all
  older vectors are already in `K₂`.
- `w₄ ∈ K₃` (see F1), so it adds nothing.

**Structural facts (checked in Lean, §6):**

| | Statement |
|---|---|
| F1 | `w_i ∈ K_i`; in particular `w₄ ∈ K₃` (from `f₄ = f₃` and `k_j ∈ K_j`) |
| F2 | `A^j w₁ ∈ K₁` for `j ≤ κ₁ − 1` (`j ≤ 3`) |
| F3 | `A(K₁ + ⟨w₂⟩) ⊆ K₂`, hence `A w₂ ∈ K₂` |
| F4 | `A(K₃ + ⟨w₄⟩) ⊆ K₄`, hence `A K₃ ⊆ K₄` and `A² w₂ ∈ K₄` |
| F5 | `dim K₂ ≤ dim K₁ + 3`, `dim K₃ ≤ dim K₂ + 1`, `dim K₄ ≤ dim K₃ + 3` |

The minimal dimensions are therefore 4, 7, 8, 11 (paper (2.3)).

## 4. Analytic part: the facts imply the conditions of §2

Notation: `v ≡ O(h^p) (mod K)` means `‖(I − P)v‖ = O(h^p)`.

**Taylor expansions.** These are general: the `O(h²)` remainders contain
`f''(·,·)` terms in *arbitrary* directions.

```
(E1) w₁ = f
(E2) w₂ = (1+c₂₁) w₁ + h z₂ + O(h²),                 z₂ ∈ ⟨A w₁⟩
(E3) w₃ = λ₁ w₁ + λ₂ w₂ + h z₃ + O(h²),              z₃ ∈ ⟨A w₁, A w₂⟩
(E4) w₄ = μ₁ w₁ + μ₂ w₂ + μ₃ w₃ + h z₄ + O(h²),      z₄ ∈ ⟨A w₁, A w₂, A w₃⟩
```

For (E2): `k₁ = w₁ + γh A w₁ + O(h²)` (exact for κ₁ ≥ 2), and
`f(u + h a₂₁ k₁) = f + h a₂₁ A f + O(h²)`. Similarly for (E3) and (E4), with
`k̃_i = w_i + γh T_i w_i + O(h²)`.

**i = 1.** `(I−P₁)A^j w₁ = 0` for `j ≤ 3` (F2, κ₁ ≥ 4).

**i = 2.**
- `j=1`: exactly 0 (F3).
- `j=2`: `A²w₂ = (1+c₂₁)A²w₁ + h A z₂ + O(h²)`. Here `A²w₁` and
  `A z₂ ∈ ⟨A²w₁⟩` lie in `K₁`, so the result is `O(h²)`.
- `j=3`: `A³w₂ = (1+c₂₁)A³w₁ + O(h)` with `A³w₁ ∈ K₁` (κ₁ ≥ 4), so the result
  is `O(h)`.

**i = 3** (needs `O(h^{3−j})`).
- `j=1`: `A w₃ = λ₁Aw₁ + λ₂Aw₂ + h A z₃ + O(h²)`. Here `Aw₁, Aw₂ ∈ K₂`;
  `A z₃ ∈ ⟨A²w₁, A²w₂⟩`; `A²w₁ ∈ K₁`; and `A²w₂ ≡ O(h²)` (i=2, j=2). So the
  result is `O(h²)`.
- `j=2`: `A²w₃ = λ₁A²w₁ + λ₂A²w₂ + O(h) ≡ O(h)`.

**i = 4.**
- `j=1`: exactly 0 (F4).
- `j=2`: `A²w₄ = Σ μ_k A²w_k + h A² z₄ + O(h²)`.
  - `A²w₁ ∈ K₁`, and `A²w₂ ∈ K₄` (F4).
  - `A²w₃ = λ₁A²w₁ + λ₂A²w₂ + h A²z₃ + O(h²)`, with
    `A²z₃ ∈ ⟨A³w₁, A³w₂⟩ ≡ O(h)`, so `A²w₃ ≡ O(h²)`.
  - `A²z₄ ∈ ⟨A³w₁, A³w₂, A³w₃⟩ ≡ O(h)` (using i=2, j=3 and (E3)).
  - So `A²w₄ ≡ O(h²)`.
- `j=3`: `A³w₄ = Σ μ_k A³w_k + O(h) ≡ O(h)`.

All conditions of §2 hold. **Hence ROWMAP has order 4** (and the embedded
solution order 3) for `κ₁ ≥ 4`, independently of n. ∎

### Which actions are necessary? (corrected 2026-10-31, confirmed numerically)

**Key observation.** At `O(h²)`, the stage vectors `w₂, w₃, w₄` leave the
Krylov space `K₁` in only *one* new direction, `g = f''(f,f)`:

- `w₂ = … + (h²a₂₁²/2) g + O(h³)`
- `w₃ = … + (h²α₃²/2) g + O(h³)`
- likewise for `w₄`

All other `h²`-terms are `A²f`-type and lie in `K₁`. Because `w₂` is
inserted at stage 2, `g ∈ K₂ + O(h)`, and `A g ∈ K₂ + O(h)` by the stage-2
closure. Hence

```
(I − P₂) w₃ = O(h³),   (I − P₂) w₄ = O(h³),   (I − P) A w₄ = O(h³) already in K₂-terms.
```

| Omitted | Consequence | Predicted order of `u_K − u_ROW` | Observed |
|---|---|---|---|
| κ₁ = 3 | S₁ = O(h³) | h⁴ | 3.97 |
| insert w₂ | (I−P₂)w₂ = O(h²) (the direction g is missing) | h³ | 2.97 |
| A q_κ₁ (stage 2) | component of w₂ along q_κ₁ is O(h²) (g generic) ⇒ (I−P₂)Aw₂ = O(h²) | h⁴ | 4.02 |
| A q_κ₁₊₁ (stage 2) | A g missing ⇒ (I−P₂)Aw₂ = O(h²) | h⁴ | 3.97 |
| insert w₃ **and** w₄ | (I−P)w₄ = O(h³) ⇒ k₄ error O(h³) | h⁴ | 3.95 |
| insert w₃ only (w₄ inserted at stage 4) | k₃ error O(h³) ⇒ k₄ error O(h⁴), b₃ = 0 | h⁵ | 4.95 |
| any or all stage-4 closure vectors | (I−P₄)Aw₄ = O(h³) ⇒ S₄ = O(h⁴) | h⁵ | 4.95–4.99 |
| w₃ + whole stage-4 closure (dims 4/7/7/8) | as above | h⁵ | 4.95 |

So for **classical order 4** (‖hA‖ → 0) the paper's construction is
sufficient but not minimal. These suffice:

```
K₁ (κ₁ ≥ 4),  K₂ = C(K₁ + ⟨w₂⟩),  w₄ (or w₃) inserted   →  dims 4, 7, 7–8, 8
```

The stage-4 closure (`A K₃ ⊆ K₄`, Lemma 2.1(2) in the paper) and the
insertion of w₃ are not needed for the *classical* order. They may matter for
stiff problems (`‖hA‖` large, where the expansion in `hA` is not valid; cf.
paper Remark 2.3) or for B-consistency. **Open question for Helmut.** The
implementation keeps the paper's schedule for now; the comparison on stiff
problems is roadmap phase 4.

The proof in §4 remains a valid *sufficiency* proof for the paper's schedule.
My first table (all actions necessary) was wrong for stages 3 and 4.

**Linear f.** Here `g = 0`, and more actions become superfluous: the
components are only `O(h^{κ₁−1})`. That is why the test problem must be
nonlinear.

Tests: `src/map/tests.rs` (`omissions_match_predicted_orders`, 12 cases on a
nonlinear system with n = 40, h = 0.2 … 0.00625).

## 5. Breakdown and the Rust realization

**Breakdown.** All definitions in §3 are subspace sums. A vector that already
lies in the space leaves the sum unchanged. F1–F5 and §4 therefore hold
*verbatim* when vectors are dropped (exact breakdown).

**Numerical breakdown.** With the relative threshold `τ = 1e3·ε`, a dropped
remainder perturbs `T_i` by `O(τ‖A‖)`. That is an h-independent error floor of
relative size `τ‖γhA‖`, far below any tolerance. It is checked by the
ε-continuity sweep.

**Basis realization** (`src/krylov.rs`). Invariant: for every non-pending
basis vector `q_j`, `A q_j ∈ K`. So "promote all pending remainders" computes
exactly `C(K)`, and the closure *rounds* of a stage are iterated closures.
The schedule (`src/map.rs`, mirrored in `lean/schedule.toml`) is:

| Stage | Action |
|---|---|
| 1 | Arnoldi from `w₁` until the residual test (2.4) holds, with `κ₁ ≥ 4` and `κ₁ ≤ max_krylov − reserve` |
| 2 | insert `w₂`, then 1 closure round |
| 3 | insert `w₃`, no closure |
| 4 | insert `w₄` (dependent in exact arithmetic), then 1 closure round |

The generic "insert, then r closure rounds" form is meant to cover the
non-autonomous case and an order-3 method too (phases 5 and 5b).

## 6. What Lean checks

`lean/RowmapLean/Map.lean`, over any field `𝕜`, a finite-dimensional `V`, a
linear `A`, and arbitrary vectors `w₁..w₄` with `w₄ ∈ K₃`:

- **F1–F4** for the subspaces defined in §3;
- **F5** as `finrank` inequalities;
- **Tightness** (`Tightness.lean`): the shift operator on ℚ¹¹ with
  `w₁ = e₀, w₂ = e₆, w₃ = e₉` attains dims 4, 7, 8, 11. Since `C(S)` is the
  smallest subspace containing `S` and `AS`, no construction with fewer
  vectors satisfies F1–F4 in general. This replaces the per-action
  counterexamples planned originally. It concerns the *structural* facts.
  Whether the *order* needs every vector is the analytic question of §4: for
  the stage-4 actions it does not, classically.

The analytic step (§4) is **not formalized** and is deliberately not stated as
a Lean `Prop`, because without the ODE setting it would be vacuous. All main
theorems depend only on `propext`, `Classical.choice` and `Quot.sound`.

## 7. Non-autonomous problems (phase 5)

### 7.1 Setting

The stage vectors are `w_i = f_i + Σ_{j<i} c_ij k_j + h d_i f_t`, with
`f_i = f(t + α_i h, u + h Σ a_ij k_j)`, `d_i = γ + Σ_j γ_ij`, `A = f_y(t,u)`,
and `φ := f_t(t,u)`. The error-propagation argument of §2 is unchanged.
Sufficient conditions: `w_i ∈ K_i`, plus `(I−P_i)A^j w_i = O(h^{4−j})` for
`i ∈ {1,2,4}` and `O(h^{3−j})` for `i = 3`.

### 7.2 What is new: f_t enters at O(h)

```
w₁ = f + hγφ
w₂ = (1+c₂₁) w₁ + h β₂' φ + h(·) A w₁ + O(h²),     β₂' = α₂ + γ₂₁   (GRK4T: 0.191)
```

So the new direction of `w₂` relative to `K₁` is `h β₂' (I−P₁)φ + O(h²)`.
In the autonomous case it was `O(h²)`. Consequences:

- Stage 2, `j = 2`: `(I−P₂)A²w₂ = h β₂' (I−P₂)A²φ + O(h²)` must be `O(h²)`,
  so `A²φ ∈ K₂ + O(h)` is needed. One closure round gives `Aφ`; **a second
  round gives `A²φ`.**
- `K₂` contains `φ` to leading order (from `w₁`, `w₂`). Together with `w₃`,
  observed numerically: `(I−P₃)φ = O(h²)`.
- Stage 4: since `f₄ = f₃`,
  `w₄ = w₃ + Σ(c₄ⱼ − c₃ⱼ)k_j + c₄₃k₃ + h(d₄ − d₃)φ` with `d₃ ≈ 0.551`,
  `d₄ ≈ −0.055`. So `w₄ ∉ K₃`, and the out-of-space part is
  `h(d₄−d₃)(I−P₃)φ = O(h³)`. **`w₄` must be inserted**, followed by one closure
  round.

### 7.3 Construction

```
K₁ = Krylov(A, w₁, κ₁),  K₂ = C²(K₁ + ⟨w₂⟩),  K₃ = K₂ + ⟨w₃⟩,  K₄ = C(K₃ + ⟨w₄⟩)
schedule: [Arnoldi, insert + 2 rounds, insert, insert + 1 round]
```

The worst-case increments are +5, +1, +5, giving minimal dims **4, 9, 10, 15**.
This is exactly the paper's §4 statement, derived here, not assumed. The
autonomous schedule is the same with 1 round at stage 2 and `w₄ ∈ K₃`.

### 7.4 Evidence (`src/map/tests_nonaut.rs`)

Test problem: `y' = M y + φ(y) + c(t)∘y + g(t)`, n = 40, nonlinear, with
`f_t`, `f_ty`, `f_tt` ≠ 0. Two independent instances.

- Dims with κ₁ = 4: **[4, 9, 10, 15]**.
- Deviation from exact non-autonomous ROW: **order 5** (observed 4.95 → 5.00).
- The *autonomous* schedule on this problem: order 4 only. The second
  stage-2 round is needed.
- Omissions (both instances agree to ±0.1). The first predictions were
  written before the first run; **6 of 10 were wrong**:

| Omitted | First prediction | Observed |
|---|---|---|
| stage 2: insert w₂ | 2 | 2 |
| stage 2: A q_κ₁ | 4 | 4 |
| stage 2: A q_{w₂} | 3 | 4 |
| stage 2, round 2: A²q_κ₁ | 5 | **4** |
| stage 2, round 2: A²q_{w₂} | 4 | 4 |
| stage 3: insert w₃ | 5 | **4** |
| stage 4: insert w₄ | 3 | 4 |
| stage 4: A of the two round-2 vectors | 5 | **4** |
| stage 4: A q_{w₃} | 4 | 4 |
| stage 4: A q_{w₄} | 4 | **5** |

So in the non-autonomous case, **the paper's construction is minimal for
classical order 4 except for one vector** (`A q_{w₄}`; dims 4/9/10/14 would
suffice). The autonomous case was different: the whole stage-4 closure was
superfluous there. My heuristic size estimates were too coarse. The
consistent reading of the observations is that `φ` is captured in `K₃` up to
`O(h²)`, through both `w₂` and `w₃`. **A rigorous derivation of the
omission orders is open.** Sufficiency, the property the implementation
relies on, is supported by the order-5 test and by the structural facts.

- Breakdown: `f_t = 0` and `f_t = −0.7 f` on a 3-eigenvalue operator give
  exact ROW to 1e-13. The ε-sweep with `f_t = f + εg` is continuous with
  deviation < 1e-12.

### 7.5 Lean (`Map.lean` section `GRK4TNonautonomous`, `TightnessNonautonomous.lean`)

- `mapSpaces_grk4t_nonautonomous`: the generic semantics of the exported
  schedule gives exactly `K₁ … K₄` above.
- Memberships `w_i ∈ K_i`. Closure facts `A(K₁+⟨w₂⟩) ⊆ K₂`,
  **`A²(K₁+⟨w₂⟩) ⊆ K₂`**, `A(K₃+⟨w₄⟩) ⊆ K₄` (so `A K₃ ⊆ K₄`), and
  `A²w₂ ∈ K₂`.
- `finrank K₂ ≤ finrank K₁ + 5`, `K₃ ≤ K₂ + 1`, `K₄ ≤ K₃ + 5`, with no
  assumption on `w₄` (κ ≥ 1).
- Tightness: the shift on ℚ¹⁵ (`w₁ = e₀, w₂ = e₇, w₃ = e₁₁, w₄ = e₁₃ ∉ K₃`)
  attains 4, 9, 10, 15.
- `schedule.toml` contains both tables; `schedule_matches_lean` checks both.

### 7.6 After the review (done 2026-10-31)

1. `GRK4T_NONAUTONOMOUS` and `dfdt` are wired into the integrator. `f_t` is
   evaluated once per accepted point `(t_m, u_m)` and reused after rejections.
2. DIFKON and DIFFU2 on coarse and paper grids: `docs/rust-runs.md` §4. The
   Rust code matches or beats Fortran sft = 2 in accuracy, with slightly
   fewer steps.
3. Tests: `tests/paper.rs` (now all 7 problems; f_t count = accepted steps)
   and `tests/api.rs::nonautonomous_forced_problem` (exact solution, stiff).

## 8. ROS34PW2, order 3 (phase 5b; public option `Method::Ros34pw2`)

### 8.1 The method

ROS34PW2 (Rang & Angermann, BIT 45 (2005) 761–787). Coefficients as in
OrdinaryDiffEq.jl v6.66.0 (`ROS34PW2Tableau`), given there in the original
α, Γ, b, b̂ form.

Checked in `src/grk4t.rs`:
- 4 stages; order 3 (residual < 1e-14), not order 4; embedded order 2;
- stiffly accurate (b = last row of α + Γ);
- L-stable (`R(∞) = 0`), A-stable, `|R̂(∞)| = 0.48`;
- no zero weight and no stage reuse;
- `R̂ ≠ R`, so the embedded estimate also works on linear problems.

**Why not ROS3P** (examined first and dropped, decision Helmut). ROS3P has
`β₂₁ = α₂₁ + γ₂₁ = 0`, so `k₂ = k₁` on linear problems. Its estimate
`(b̂ − b)·k` then vanishes identically. This is unavoidable for ROS3P: the
order-2 conditions force `R̂ = R` on linear problems for every embedding. Its
step-size control is blind on linear and nearly linear problems (observed:
the error stays at about 1e-2 for tol 1e-3 … 1e-7).

### 8.2 Conditions for classical order 3

The local error must be `O(h⁴)`. All weights are nonzero, so every stage
needs `S_i = O(h³)`. With Lemma 2.1: `(I−P_i)A^j w_i = O(h^{3−j})` and
`w_i ∈ K_i` up to `O(h³)`.

### 8.3 Derivation

**Stage 1:** `κ₁ ≥ 3` gives `S₁ = O(h³)`. In general, **κ₁ ≥ order**: the
step uses `RowMethod::order` (4 for GRK4T).

**Autonomous:**
- `w₂` leaves `K₁` at `O(h²)`, in the direction `f''(f,f)`. Since `b₂ ≠ 0`
  that is not allowed, so **insert `w₂`**. No closure is needed:
  `(I−P)Aw₂ = O(h²)`.
- `w₃` and `w₄` have the same `O(h²)` direction `f''(f,f)`, which is already
  in `K₂`. So they leave `K₂` only at `O(h³)`: **nothing at stages 3 and
  4**.
- **Dims 3, 4, 4, 4.**

**Non-autonomous:**
- `β₂' = α₂₁ + γ₂₁ = 0`, so stage 2 is as in the autonomous case.
- `β₃' = α₃ + γ₃₁ + γ₃₂ ≈ −0.118 ≠ 0`, so `w₃` has an `O(h)` direction
  `f_t`: **insert `w₃`**.
- `(I−P)Aw₃ = O(h²)` needs `A f_t ∈ K₃` up to `O(h)`. Since `f_t` is
  generic, all of its components matter, so a **full closure round** is
  needed (+3).
- `w₄`: first predicted to need insertion (`b₄ ≠ 0`). The order tests show
  it does not: `f_t` is captured in `K₃` up to `O(h²)` (through `w₁`, `w₃`
  and the closure), so `w₄` leaves `K₃` only at `O(h³)`.
- **Dims 3, 4, 8, 8.**

Schedules (`src/map.rs`, mirrored in `lean/schedule.toml`):

```
ROS34PW2_AUTONOMOUS    = [Arnoldi, (insert, 0), (no insert, 0), (no insert, 0)]   dims ≤ κ₁, κ₁+1, κ₁+1, κ₁+1
ROS34PW2_NONAUTONOMOUS = [Arnoldi, (insert, 0), (insert, 1),    (no insert, 0)]   dims ≤ κ₁, κ₁+1, κ₁+5, κ₁+5
```

### 8.4 Evidence (`src/map/tests_ros34pw2.rs`)

Measured: the order of the deviation from exact ROS34PW2 (full-dimensional
space), with h = 0.2 … 0.0016, and down to 2·10⁻⁴ for two slowly converging
cases. Two instances each of `Nl` and `NlT`. The expected order at the
minimal dims is 4. The predictions were written before the first run.

| Case | Predicted | Measured (two instances) |
|---|---|---|
| autonomous, 3/4/4/4 | 4 | 3.97, 3.97 |
| autonomous, κ₁ = 2 | 3 | 2.99, 2.99 |
| autonomous, w₂ not inserted | 3 | 3.00, 3.00 |
| non-autonomous, 3/4/8/8 | 4 | 3.99, 3.99 |
| non-autonomous, κ₁ = 2 | 3 | 2.99, 2.99 |
| non-autonomous with the autonomous schedule | 2 | 2.01, 2.01 |
| non-autonomous, w₂ not inserted | 3 | 2.99, 2.99 |
| non-autonomous, w₃ not inserted | 2 | 1.99, 1.99 |
| non-autonomous, A q_κ₁ omitted | 3 | 3.04, 3.01 (slow convergence) |
| non-autonomous, A q_{w₂} omitted | 3 | 3.02, 3.02 (slow convergence) |
| non-autonomous, A q_{w₃} omitted | 3 | 2.99, 3.00 |
| non-autonomous, w₄ also inserted (3/4/8/9) | first prediction: needed | 3.99, 3.99 (not needed) |

The only miss, the `w₄` insertion, is corrected: the minimal schedule omits
it.

API test (`tests/api.rs::ros34pw2_method`): on the very stiff spectrum
λ = −1 … −981 (autonomous and non-autonomous), the error is proportional to
the tolerance (tol 1e-3 → 1e-7: error 9.7e-5 → 1.2e-7).

### 8.5 Lean (`lean/RowmapLean/Ros34pw2.lean`)

- The generic semantics of the exported schedules gives
  `[K₁, K₁+⟨w₂⟩, K₁+⟨w₂⟩, K₁+⟨w₂⟩]` and `[K₁, K₂, C(K₂+⟨w₃⟩), C(K₂+⟨w₃⟩)]`.
- `w₂ ∈ K₂`, `w₃ ∈ K₃`; the closure `A(K₂+⟨w₃⟩) ⊆ K₃`.
- `finrank K₂ ≤ finrank K₁ + 1` and `finrank K₃ ≤ finrank K₂ + 4`
  (non-autonomous).
- Tightness: the shift on ℚ¹¹ (`w₁ = e₀, w₂ = e₄, w₃ = e₇`, κ = 3) attains
  3, 4, 8.
- Standard axioms only. `schedule.toml` has four tables (GRK4T and
  ROS34PW2, each autonomous and non-autonomous), checked by
  `schedule_matches_lean`.

## 9. The structure of H (generalized Hessenberg)

Every vector that enters the basis has a fixed structural consequence for
`H = QᵀAQ`. The MAP is a generalized Arnoldi process: plain Arnoldi gives an
upper Hessenberg H, and each further insertion adds one subdiagonal, for as
long as later vectors follow. This section states the structure exactly. It
is proved in Lean (`lean/RowmapLean/HStructure.lean`) and checked against the
implementation (`src/map/tests_structure.rs`).

### 9.1 The process

The basis `q₀, q₁, …` is built by a sequence of operations. Each operation is
*accepted* (a new vector enters) or *dropped* (its candidate lies in the
current space, §3.4):

- **ins**: insert a stage vector `w`; the candidate is `w`.
- **prom**: promote the remainder of the *oldest open column* `f`; the
  candidate is `A q_f`.

Column `j` is *open* from the moment `q_j` enters until its remainder is
promoted or dropped. Then it is *closed* at its **closing row** `c(j)`:

- accepted promotion: `c(j)` is the index of the new vector (`q_{c(j)}` is
  the normalized remainder of `A q_j`);
- dropped promotion: `c(j) = k − 1`, the last row at that moment.

**FIFO.** The schedule promotes in the order in which columns opened. A
closure round promotes the pending columns at its start in index order, and
the columns that open during a round wait for the next round. Stage 1 has a
single open column. So the open columns always form a queue.

Let `k` be the dimension, `f` the number of closed columns, and `o(m)` the
number of open columns at the moment `q_m` enters (for a promotion, this
includes the column being closed).

### 9.2 Theorem (`H_structure`)

Let the vectors `q₀ … q_{k−1}` be mutually orthogonal and satisfy

- (flag) if `q_m` is the promotion of column `j`, then `A q_j ∈ span{q₀ … q_m}`;
- (drop) if the promotion of `j` was dropped, then `A q_j ∈ span{q₀ … q_{c(j)}}`.

Then, for every sequence of operations:

- **(S1) Zero pattern.** For every closed column j, `H[i, j] = ⟨q_i, A q_j⟩ = 0`
  for all `i > c(j)`.
- **(S2) Last entry.** An accepted promotion gives `H[c(j), j] = ‖p_j‖ > 0`,
  the norm of the remainder. This holds by construction; it is not in Lean.
- **(S3) Staircase.** The closed columns are exactly `0 … f−1`, the open ones
  `f … k−1`, and `c` is nondecreasing: `i ≤ j < f ⇒ c(i) ≤ c(j)`.
- **(S4) Rows.** Row m vanishes left of column `m − o(m)`. It has at most
  `o(m)` nonzeros below the diagonal, in the contiguous columns
  `m − o(m) … m − 1`. For generic data, exactly these are nonzero (Rust test).
- **(S5) Counts.** The sequence `o(0), o(1), …` is a function of the
  operation sequence alone. Write L for the number of open columns:

  | operation | entry `o(m)` | L afterwards |
  |---|---|---|
  | ins, accepted | L | L + 1 |
  | ins, dropped | — | L |
  | prom, accepted | L | L |
  | prom, dropped | — | L − 1 |

The inserted vectors themselves never enter the statement. The structure
depends only on the order of the operations and on which ones were accepted.
Frozen remainders (phase 7) are open columns that are never promoted, so
they change nothing.

*Proof sketch.* S1: by (flag) or (drop), `A q_j ∈ span{q₀ … q_{c(j)}}`, and
each `q_i` with `i > c(j)` is orthogonal to that span. S3–S5: an invariant
of the operation semantics, by induction over the operations. For S4, the
invariant is `j < m − o(m) ⇔ (j closed ∧ c(j) < m)`. FIFO is essential: a
skipped promotion (test hook only) leaves an old column open while newer ones
close, and the staircase breaks (`omitted_promotion_breaks_monotone_profile`).

### 9.3 The four schedules (`trace_*`, `bandwidths`)

Without breakdown (all promotions accepted; for GRK4T autonomous, w₄ ∈ K₃ is
dropped), Lean proves for **every κ₁**:

| schedule | `o(0 … κ−1)` | lower bandwidth |
|---|---|---|
| GRK4T autonomous | 0, 1×(κ₁−1), 1, 2, 2, 2, 3, 3, 3 | 3 |
| GRK4T non-autonomous | 0, 1×(κ₁−1), 1, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4 | 4 |
| ROS34PW2 autonomous | 0, 1×(κ₁−1), 1 | 1 (upper Hessenberg) |
| ROS34PW2 non-autonomous | 0, 1×(κ₁−1), 1, 2, 3, 3, 3 | 3 |

The bandwidth is 1 plus the number of accepted insertions after stage 1 that
are followed by another vector. With breakdown, each dropped remainder lowers
L by one from then on (S5), so the bandwidth can only shrink. For ROS34PW2
autonomous, `w₂` is the last vector, so no row sees the extra open column.
`h_is_generalized_hessenberg` checks these sequences against the
implementation for κ₁ ∈ {order, 6, 9}, together with S1–S4.

### 9.4 Consequences for the stage systems

All stages of a step solve with the same matrix `M = I − γhH`. Only H grows,
and the right-hand side changes. M has the profile of H. Hence:

1. **Givens QR without fill-in below the profile.** Eliminating column j
   needs `c(j) − j` rotations on rows `j … c(j)`. The rotations of earlier
   columns act on rows `≤ c(j−1) ≤ c(j)` (S3), so they create no entries
   below the staircase.
2. **Closed columns are final.** Once column j is closed, all its entries are
   known (S1), and so are its rotations. Later rows are zero in column j, so
   the rotations never have to be redone. The factorization grows column by
   column in closing order, across the stage-1 iterations and across stages
   2–4. Only the open columns (the last ones, S3) have to be treated at each
   solve.
3. **Residual without back substitution.** The Galerkin residual is
   `r = −γh Σ_{open j} l_j p_j`, and the open columns are the *last*
   components of `l`, exactly the first ones a back substitution produces.
   In stage 1 there is exactly one open column (o = 1):
   `r = −γh l_κ p_κ`, with `l_κ = g_κ / r_κκ` read directly off the rotated
   right-hand side, so `‖r‖ = γh |l_κ| ‖p_κ‖`. One Arnoldi iteration then
   costs O(κ) for the new column's rotations instead of a fresh O(κ³) LU
   plus an n-vector residual. The full back substitution is needed only once,
   when the test (2.4) passes. (For the same reason, Fortran issue §4.7 is
   moot in stage 1: with one open column the triangle-inequality bound is
   exact.)
4. **Right-hand sides.** `Qᵀw₁ = ‖w₁‖e₀`. An inserted `w_i` has coordinates
   only up to its own row (later vectors are orthogonal to it). Each stage
   needs the stored rotations applied to the new right-hand side, O(κ·o),
   plus one back substitution.

The current implementation (dense LU from scratch, AGENTS.md §3.3) does not
use any of this yet.
