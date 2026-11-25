/-
Copyright (c) 2026 Helmut Podhaisky. All rights reserved.
Released under MIT license as described in the file LICENSE.
Authors: Helmut Podhaisky
-/
import RowmapLean.Map
import RowmapLean.Ros34pw2

/-!
# The generalized Hessenberg structure of `H = QᵀAQ` in the MAP

See `docs/map.md` §9. The multiple Arnoldi process builds the basis by a
sequence of operations; each one may be accepted (a new basis vector enters)
or dropped (the candidate lies in the current space):

* `ins` — insert a stage vector `w` (the candidate is `w`);
* `prom` — promote the remainder of the *oldest open column* `f` (the
  candidate is `A q_f`). Promotion is FIFO: the closure rounds of
  `src/map.rs` promote the pending columns in index order, and stage 1 has a
  single open column.

Column `j` is *open* until its remainder is promoted or dropped; then it is
*closed* at row `c j` (the index of the vector created from it, or the last
row if dropped). The process state (`St`) records `c`, the origin of every
basis vector, and `o m`, the number of open columns when vector `m` entered.

Main results (`H_structure`), for any vectors `q` that are mutually orthogonal
and satisfy the Gram–Schmidt flag property (a vector created by promoting
`j` spans, together with its predecessors, `A q_j`):

* (S1) zero pattern: `⟪q i, A q j⟫ = 0` for every closed column `j` and every
  row `i > c j`;
* (S3) the closed columns are exactly `0, …, f-1`, the open ones `f, …, k-1`,
  and `c` is monotone: the profile of `H` is a staircase;
* (S4) row `m` vanishes left of column `m - o m`;
* (S5) the sequence `o 0, …, o (k-1)` is `trace 0 ops`, a function of the
  operation sequence alone.

For the four schedules (no breakdown; for GRK4T autonomous `w₄ ∈ K₃` is
dropped), `trace` gives the lower bandwidth exactly:
GRK4T 3 (autonomous) / 4 (non-autonomous), ROS34PW2 1 / 3, for every `κ₁`.

The inserted vectors themselves never enter: the structure depends only on
the order of operations and on which ones were accepted.
-/

namespace Rowmap.HStructure

open Submodule

/-- A basis operation of the MAP. -/
inductive Op where
  /-- Insert a stage vector. -/
  | ins
  /-- Promote the remainder of the oldest open column. -/
  | prom
  deriving DecidableEq, Repr

/-- State of the process. Columns `< f` are closed, columns `f … k-1` open. -/
structure St where
  /-- Dimension: basis vectors `q 0 … q (k-1)`. -/
  k : ℕ
  /-- Number of closed columns (the open ones form the FIFO queue `f … k-1`). -/
  f : ℕ
  /-- Closing row of column `j` (meaningful for `j < f`). -/
  c : ℕ → ℕ
  /-- Column `j` was closed by a dropped promotion. -/
  dropped : ℕ → Bool
  /-- Origin of vector `m`: `some j` = promotion of column `j`, `none` = insertion. -/
  origin : ℕ → Option ℕ
  /-- Number of open columns when vector `m` entered. -/
  o : ℕ → ℕ

/-- The empty basis. -/
def St.init : St := ⟨0, 0, fun _ => 0, fun _ => false, fun _ => none, fun _ => 0⟩

/-- One operation; the flag says whether it was accepted. A promotion with
no open column is skipped (its role does not exist). -/
def step (s : St) : Op × Bool → St
  | (.ins, false) => s
  | (.ins, true) => insAcc s
  | (.prom, acc) => if s.f < s.k then (if acc then promAcc s else promDrop s) else s
where
  /-- Accepted insertion: vector `k` enters, `o k = k - f` columns are open. -/
  insAcc (s : St) : St :=
    ⟨s.k + 1, s.f, s.c, s.dropped, Function.update s.origin s.k none,
     Function.update s.o s.k (s.k - s.f)⟩
  /-- Accepted promotion of column `f`: vector `k` enters and closes column `f`
  at row `k`. -/
  promAcc (s : St) : St :=
    ⟨s.k + 1, s.f + 1, Function.update s.c s.f s.k, s.dropped,
     Function.update s.origin s.k (some s.f), Function.update s.o s.k (s.k - s.f)⟩
  /-- Dropped promotion: column `f` closes at the last row `k - 1`. -/
  promDrop (s : St) : St :=
    ⟨s.k, s.f + 1, Function.update s.c s.f (s.k - 1),
     Function.update s.dropped s.f true, s.origin, s.o⟩

/-- Run a sequence of operations. -/
def run (s : St) (ops : List (Op × Bool)) : St := ops.foldl step s

/-- The invariant of the process. -/
structure Inv (s : St) : Prop where
  f_le : s.f ≤ s.k
  c_ge : ∀ j < s.f, j ≤ s.c j
  c_lt : ∀ j < s.f, s.c j < s.k
  mono : ∀ i j, i ≤ j → j < s.f → s.c i ≤ s.c j
  origin_c : ∀ j < s.f, s.dropped j = false → s.origin (s.c j) = some j
  row : ∀ m < s.k, ∀ j < m, (j < m - s.o m ↔ (j < s.f ∧ s.c j < m))

theorem inv_init : Inv St.init where
  f_le := le_rfl
  c_ge := fun _ h => absurd h (Nat.not_lt_zero _)
  c_lt := fun _ h => absurd h (Nat.not_lt_zero _)
  mono := fun _ _ _ h => absurd h (Nat.not_lt_zero _)
  origin_c := fun _ h => absurd h (Nat.not_lt_zero _)
  row := fun _ h => absurd h (Nat.not_lt_zero _)

open Function in
theorem inv_insAcc (s : St) (hs : Inv s) : Inv (step.insAcc s) := by
  obtain ⟨hf, hge, hlt, hmono, hor, hrow⟩ := hs
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩ <;> simp only [step.insAcc]
  · omega
  · exact hge
  · intro j hj; have := hlt j hj; omega
  · exact hmono
  · intro j hj hd
    rw [update_of_ne (by have := hlt j hj; omega)]
    exact hor j hj hd
  · intro m hm j hjm
    by_cases hmk : m = s.k
    · subst hmk
      rw [update_self]
      by_cases hjf : j < s.f
      · have := hlt j hjf; omega
      · omega
    · rw [update_of_ne hmk]
      exact hrow m (by omega) j hjm

open Function in
theorem inv_promAcc (s : St) (hs : Inv s) (hfk : s.f < s.k) : Inv (step.promAcc s) := by
  obtain ⟨hf, hge, hlt, hmono, hor, hrow⟩ := hs
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩ <;> simp only [step.promAcc]
  · omega
  · intro j hj
    by_cases h : j = s.f
    · subst h; rw [update_self]; omega
    · rw [update_of_ne h]; exact hge j (by omega)
  · intro j hj
    by_cases h : j = s.f
    · subst h; rw [update_self]; omega
    · rw [update_of_ne h]; have := hlt j (by omega); omega
  · intro i j hij hj
    by_cases hj' : j = s.f
    · subst hj'
      rw [update_self]
      by_cases hi : i = s.f
      · subst hi; rw [update_self]
      · rw [update_of_ne hi]; have := hlt i (by omega); omega
    · have hi : i ≠ s.f := by omega
      rw [update_of_ne hj', update_of_ne hi]
      exact hmono i j hij (by omega)
  · intro j hj hd
    by_cases h : j = s.f
    · subst h; rw [update_self, update_self]
    · rw [update_of_ne h]
      have := hlt j (by omega)
      rw [update_of_ne (by omega)]
      exact hor j (by omega) hd
  · intro m hm j hjm
    by_cases hmk : m = s.k
    · subst hmk
      rw [update_self]
      by_cases h : j = s.f
      · subst h; rw [update_self]; omega
      · rw [update_of_ne h]
        by_cases hjf : j < s.f
        · have := hlt j hjf; omega
        · omega
    · rw [update_of_ne hmk]
      have hr := hrow m (by omega) j hjm
      by_cases h : j = s.f
      · subst h; rw [update_self]; omega
      · rw [update_of_ne h]; omega

open Function in
theorem inv_promDrop (s : St) (hs : Inv s) (hfk : s.f < s.k) : Inv (step.promDrop s) := by
  obtain ⟨hf, hge, hlt, hmono, hor, hrow⟩ := hs
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩ <;> simp only [step.promDrop]
  · omega
  · intro j hj
    by_cases h : j = s.f
    · subst h; rw [update_self]; omega
    · rw [update_of_ne h]; exact hge j (by omega)
  · intro j hj
    by_cases h : j = s.f
    · subst h; rw [update_self]; omega
    · rw [update_of_ne h]; exact hlt j (by omega)
  · intro i j hij hj
    by_cases hj' : j = s.f
    · subst hj'
      rw [update_self]
      by_cases hi : i = s.f
      · subst hi; rw [update_self]
      · rw [update_of_ne hi]; have := hlt i (by omega); omega
    · have hi : i ≠ s.f := by omega
      rw [update_of_ne hj', update_of_ne hi]
      exact hmono i j hij (by omega)
  · intro j hj hd
    by_cases h : j = s.f
    · subst h; rw [update_self] at hd; exact absurd hd (by simp)
    · rw [update_of_ne h] at hd; rw [update_of_ne h]
      exact hor j (by omega) hd
  · intro m hm j hjm
    have hr := hrow m hm j hjm
    by_cases h : j = s.f
    · subst h; rw [update_self]; omega
    · rw [update_of_ne h]; omega

theorem inv_step (s : St) (op : Op × Bool) (hs : Inv s) : Inv (step s op) := by
  rcases op with ⟨_ | _, _ | _⟩
  · exact hs
  · exact inv_insAcc s hs
  · by_cases hfk : s.f < s.k
    · simp only [step, hfk, ite_true, Bool.false_eq_true, ite_false]
      exact inv_promDrop s hs hfk
    · simp only [step, hfk, ite_false]; exact hs
  · by_cases hfk : s.f < s.k
    · simp only [step, hfk, ite_true]
      exact inv_promAcc s hs hfk
    · simp only [step, hfk, ite_false]; exact hs

theorem inv_run (ops : List (Op × Bool)) : ∀ s, Inv s → Inv (run s ops) := by
  induction ops with
  | nil => intro s hs; exact hs
  | cons op ops ih => intro s hs; exact ih _ (inv_step s op hs)

/-! ## The open-column counts as a function of the operations -/

/-- The values `o m` of the vectors entering during `ops`, starting with `L`
open columns. -/
def trace : ℕ → List (Op × Bool) → List ℕ
  | _, [] => []
  | L, (.ins, true) :: ops => L :: trace (L + 1) ops
  | L, (.ins, false) :: ops => trace L ops
  | 0, (.prom, _) :: ops => trace 0 ops
  | L + 1, (.prom, true) :: ops => (L + 1) :: trace (L + 1) ops
  | L + 1, (.prom, false) :: ops => trace L ops

theorem run_cons (s : St) (op : Op × Bool) (ops : List (Op × Bool)) :
    run s (op :: ops) = run (step s op) ops := rfl

theorem map_update_range (o : ℕ → ℕ) (k v : ℕ) :
    (List.range (k + 1)).map (Function.update o k v) = (List.range k).map o ++ [v] := by
  rw [List.range_succ, List.map_append]
  congr 1
  · apply List.map_congr_left
    intro m hm
    rw [Function.update_of_ne (by simp at hm; omega)]
  · simp

theorem o_trace (ops : List (Op × Bool)) : ∀ s : St, s.f ≤ s.k →
    (List.range (run s ops).k).map (run s ops).o
      = (List.range s.k).map s.o ++ trace (s.k - s.f) ops := by
  induction ops with
  | nil => intro s _; simp [run, trace]
  | cons op ops ih =>
    intro s hf
    rcases op with ⟨_ | _, _ | _⟩
    · rw [run_cons, show step s (.ins, false) = s from rfl, ih s hf, trace]
    · rw [run_cons, show step s (.ins, true) = step.insAcc s from rfl,
        ih _ (by simp only [step.insAcc]; omega)]
      simp only [step.insAcc, map_update_range, trace, List.append_assoc,
        List.singleton_append]
      rw [show s.k + 1 - s.f = s.k - s.f + 1 by omega]
    · by_cases hfk : s.f < s.k
      · obtain ⟨L, hL⟩ : ∃ L, s.k - s.f = L + 1 := ⟨s.k - s.f - 1, by omega⟩
        have hstep : step s (.prom, false) = step.promDrop s := by simp [step, hfk]
        rw [run_cons, hstep, ih _ (by simp only [step.promDrop]; omega)]
        simp only [step.promDrop]
        rw [hL, trace, show s.k - (s.f + 1) = L by omega]
      · have hstep : step s (.prom, false) = s := by simp [step, hfk]
        rw [run_cons, hstep, ih s hf, show s.k - s.f = 0 by omega, trace]
    · by_cases hfk : s.f < s.k
      · obtain ⟨L, hL⟩ : ∃ L, s.k - s.f = L + 1 := ⟨s.k - s.f - 1, by omega⟩
        have hstep : step s (.prom, true) = step.promAcc s := by simp [step, hfk]
        rw [run_cons, hstep, ih _ (by simp only [step.promAcc]; omega)]
        simp only [step.promAcc, map_update_range, List.append_assoc, List.singleton_append]
        rw [hL, trace, show s.k + 1 - (s.f + 1) = L + 1 by omega]
      · have hstep : step s (.prom, true) = s := by simp [step, hfk]
        rw [run_cons, hstep, ih s hf, show s.k - s.f = 0 by omega, trace]

/-! ## The structure theorem -/

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- A vector in the span of `q 0 … q c` is orthogonal to every `q i` that is
orthogonal to these. -/
theorem inner_eq_zero_of_mem_span (q : ℕ → E) (i c : ℕ)
    (horth : ∀ l, l ≤ c → inner ℝ (q i) (q l) = 0) {x : E}
    (hx : x ∈ span ℝ (q '' Set.Iic c)) : inner ℝ (q i) x = 0 := by
  have hle : span ℝ (q '' Set.Iic c) ≤ LinearMap.ker (innerₛₗ ℝ (q i)) := by
    rw [span_le]
    rintro _ ⟨l, hl, rfl⟩
    simpa [LinearMap.mem_ker, innerₛₗ_apply_apply] using horth l hl
  simpa [LinearMap.mem_ker, innerₛₗ_apply_apply] using hle hx

/-- **Structure of `H = QᵀAQ`** for any state reachable by the process.
Hypotheses on the vectors: mutual orthogonality, the Gram–Schmidt flag
property of promoted vectors, and the meaning of a dropped promotion. -/
theorem H_structure_of_inv (s : St) (hs : Inv s) (A : E →ₗ[ℝ] E) (q : ℕ → E)
    (horth : ∀ i l, l < i → i < s.k → inner ℝ (q i) (q l) = 0)
    (hprom : ∀ m < s.k, ∀ j, s.origin m = some j → A (q j) ∈ span ℝ (q '' Set.Iic m))
    (hdrop : ∀ j < s.f, s.dropped j = true → A (q j) ∈ span ℝ (q '' Set.Iic (s.c j))) :
    -- (S1) zero pattern below the closing row
    (∀ j < s.f, ∀ i < s.k, s.c j < i → inner ℝ (q i) (A (q j)) = 0) ∧
    -- (S3) closed columns = prefix, open columns = suffix, monotone closing rows
    s.f ≤ s.k ∧ (∀ i j, i ≤ j → j < s.f → s.c i ≤ s.c j) ∧
    -- (S4) row `m` vanishes left of column `m - o m`
    (∀ m < s.k, ∀ j < m - s.o m, inner ℝ (q m) (A (q j)) = 0) := by
  have S1 : ∀ j < s.f, ∀ i < s.k, s.c j < i → inner ℝ (q i) (A (q j)) = 0 := by
    intro j hj i hi hci
    apply inner_eq_zero_of_mem_span q i (s.c j)
      (fun l hl => horth i l (by omega) hi)
    cases hd : s.dropped j
    · exact hprom (s.c j) (hs.c_lt j hj) j (hs.origin_c j hj hd)
    · exact hdrop j hj hd
  refine ⟨S1, hs.f_le, hs.mono, ?_⟩
  intro m hm j hj
  have h := (hs.row m hm j (by omega)).mp hj
  exact S1 j h.1 m hm h.2

/-- The structure theorem for a run of the process from the empty basis,
with the open-column counts (S5). -/
theorem H_structure (ops : List (Op × Bool)) (A : E →ₗ[ℝ] E) (q : ℕ → E)
    (horth : ∀ i l, l < i → i < (run St.init ops).k → inner ℝ (q i) (q l) = 0)
    (hprom : ∀ m < (run St.init ops).k, ∀ j, (run St.init ops).origin m = some j →
      A (q j) ∈ span ℝ (q '' Set.Iic m))
    (hdrop : ∀ j < (run St.init ops).f, (run St.init ops).dropped j = true →
      A (q j) ∈ span ℝ (q '' Set.Iic ((run St.init ops).c j))) :
    (∀ j < (run St.init ops).f, ∀ i < (run St.init ops).k, (run St.init ops).c j < i →
        inner ℝ (q i) (A (q j)) = 0) ∧
    (run St.init ops).f ≤ (run St.init ops).k ∧
    (∀ i j, i ≤ j → j < (run St.init ops).f → (run St.init ops).c i ≤ (run St.init ops).c j) ∧
    (∀ m < (run St.init ops).k, ∀ j < m - (run St.init ops).o m,
        inner ℝ (q m) (A (q j)) = 0) ∧
    (List.range (run St.init ops).k).map (run St.init ops).o = trace 0 ops := by
  have hs := inv_run ops _ inv_init
  obtain ⟨h1, h2, h3, h4⟩ := H_structure_of_inv _ hs A q horth hprom hdrop
  refine ⟨h1, h2, h3, h4, ?_⟩
  simpa [St.init] using o_trace ops St.init le_rfl

/-! ## The schedules: operation sequences and bandwidths -/

/-- Operations of one stage when every promotion is accepted, starting with
`L` open columns; the insertion (if any) is accepted iff `acc`. A closure
round promotes all `L'` open columns (FIFO), and an accepted promotion keeps
the number of open columns. Returns the operations and the new count. -/
def stageOps (κ L : ℕ) : StageKind → Bool → List (Op × Bool) × ℕ
  | .arnoldi, _ => ((.ins, true) :: List.replicate (κ - 1) (.prom, true), 1)
  | .extend ins rounds, acc =>
      let L' := if ins && acc then L + 1 else L
      ((if ins then [(.ins, acc)] else []) ++ List.replicate (rounds * L') (.prom, true), L')

/-- Operations of a schedule; `acc` gives the acceptance of each stage's insertion. -/
def schedOps (κ : ℕ) : List StageKind → List Bool → ℕ → List (Op × Bool)
  | st :: sts, a :: as, L => (stageOps κ L st a).1 ++ schedOps κ sts as (stageOps κ L st a).2
  | _, _, _ => []

theorem trace_prom_replicate (n : ℕ) (t : List (Op × Bool)) :
    trace 1 (List.replicate n (.prom, true) ++ t) = List.replicate n 1 ++ trace 1 t := by
  induction n with
  | zero => simp
  | succ n ih => simp [List.replicate_succ, trace, ih]

theorem trace_arnoldi (κ : ℕ) (t : List (Op × Bool)) :
    trace 0 (((.ins, true) :: List.replicate (κ - 1) (.prom, true)) ++ t)
      = 0 :: List.replicate (κ - 1) 1 ++ trace 1 t := by
  simp [trace, trace_prom_replicate]

/-- GRK4T, autonomous (no breakdown; `w₄ ∈ K₃` is dropped). -/
theorem trace_grk4t_aut (κ : ℕ) :
    trace 0 (schedOps κ grk4tAutonomous [true, true, true, false] 0)
      = 0 :: List.replicate (κ - 1) 1 ++ [1, 2, 2, 2, 3, 3, 3] := by
  simp only [grk4tAutonomous, schedOps, stageOps]
  rw [List.append_assoc, trace_arnoldi]
  rfl

/-- GRK4T, non-autonomous (no breakdown). -/
theorem trace_grk4t_nonaut (κ : ℕ) :
    trace 0 (schedOps κ grk4tNonautonomous [true, true, true, true] 0)
      = 0 :: List.replicate (κ - 1) 1 ++ [1, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4] := by
  simp only [grk4tNonautonomous, schedOps, stageOps]
  rw [List.append_assoc, trace_arnoldi]
  rfl

/-- ROS34PW2, autonomous: `H` is upper Hessenberg. -/
theorem trace_ros34pw2_aut (κ : ℕ) :
    trace 0 (schedOps κ ros34pw2Autonomous [true, true, true, true] 0)
      = 0 :: List.replicate (κ - 1) 1 ++ [1] := by
  simp only [ros34pw2Autonomous, schedOps, stageOps]
  rw [List.append_assoc, trace_arnoldi]
  rfl

/-- ROS34PW2, non-autonomous. -/
theorem trace_ros34pw2_nonaut (κ : ℕ) :
    trace 0 (schedOps κ ros34pw2Nonautonomous [true, true, true, true] 0)
      = 0 :: List.replicate (κ - 1) 1 ++ [1, 2, 3, 3, 3] := by
  simp only [ros34pw2Nonautonomous, schedOps, stageOps]
  rw [List.append_assoc, trace_arnoldi]
  rfl

/-- The lower bandwidth `max o m` of `H` for the four schedules, every `κ₁`. -/
theorem bandwidths (κ : ℕ) :
    (∀ x ∈ trace 0 (schedOps κ grk4tAutonomous [true, true, true, false] 0), x ≤ 3) ∧
    3 ∈ trace 0 (schedOps κ grk4tAutonomous [true, true, true, false] 0) ∧
    (∀ x ∈ trace 0 (schedOps κ grk4tNonautonomous [true, true, true, true] 0), x ≤ 4) ∧
    4 ∈ trace 0 (schedOps κ grk4tNonautonomous [true, true, true, true] 0) ∧
    (∀ x ∈ trace 0 (schedOps κ ros34pw2Autonomous [true, true, true, true] 0), x ≤ 1) ∧
    1 ∈ trace 0 (schedOps κ ros34pw2Autonomous [true, true, true, true] 0) ∧
    (∀ x ∈ trace 0 (schedOps κ ros34pw2Nonautonomous [true, true, true, true] 0), x ≤ 3) ∧
    3 ∈ trace 0 (schedOps κ ros34pw2Nonautonomous [true, true, true, true] 0) := by
  rw [trace_grk4t_aut, trace_grk4t_nonaut, trace_ros34pw2_aut, trace_ros34pw2_nonaut]
  refine ⟨?_, by simp, ?_, by simp, ?_, by simp, ?_, by simp⟩ <;>
  · intro x hx
    simp only [List.mem_cons, List.mem_append, List.mem_replicate, List.not_mem_nil,
      or_false] at hx
    omega

end Rowmap.HStructure
