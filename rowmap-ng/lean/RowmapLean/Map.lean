/-
Copyright (c) 2026 Helmut Podhaisky. All rights reserved.
Released under MIT license as described in the file LICENSE.
Authors: Helmut Podhaisky
-/
import Mathlib

/-!
# Structural facts of the ROWMAP multiple Arnoldi process (MAP)

See `docs/map.md`. The MAP is described *representation-free* by subspaces:
with the closure `C(S) = S + A S`,

  K₁ = span{w₁, A w₁, …, A^{κ-1} w₁},
  K₂ = C(K₁ + ⟨w₂⟩),  K₃ = K₂ + ⟨w₃⟩,  K₄ = C(K₃ + ⟨w₄⟩).

The schedule is data (`StageKind`, mirrored in `src/map.rs` and exported to
`schedule.toml`); `mapSpaces` gives it a generic semantics, and the theorems
below are about `mapSpaces A κ w grk4tAutonomous`.

Everything holds over any field and for arbitrary vectors `w i`; dropped
(dependent) vectors need no special treatment, because a vector already in
the space does not change a subspace sum.

The analytic step "these facts imply classical order 4" is *not* formalized;
it is the hypothesis `OrderFromIdentities` (docs/map.md §4).
-/

namespace Rowmap

open Submodule Module

/-- How a stage builds its space (mirrors `StageKind` in `src/map.rs`). -/
inductive StageKind where
  /-- Stage 1: Arnoldi (Krylov space of dimension κ). -/
  | arnoldi
  /-- Insert `w_i` (if `insert`), then `rounds` closure rounds. -/
  | extend (insert : Bool) (rounds : ℕ)
  deriving DecidableEq, Repr

/-- The MAP schedule for GRK4T, autonomous case. -/
def grk4tAutonomous : List StageKind :=
  [.arnoldi, .extend true 1, .extend true 0, .extend true 1]

variable {𝕜 V : Type*} [Field 𝕜] [AddCommGroup V] [Module 𝕜 V]
variable (A : V →ₗ[𝕜] V)

/-- Krylov space `span{w, A w, …, A^{κ-1} w}`. -/
def krylov (w : V) (κ : ℕ) : Submodule 𝕜 V :=
  span 𝕜 (Set.range fun j : Fin κ => (⇑A)^[j] w)

/-- Closure `C(S) = S + A S` ("promote all pending remainders"). -/
def close (S : Submodule 𝕜 V) : Submodule 𝕜 V := S ⊔ S.map A

/-- `r` closure rounds. -/
def closeIter : ℕ → Submodule 𝕜 V → Submodule 𝕜 V
  | 0, S => S
  | r + 1, S => closeIter r (close A S)

/-- The space after an `extend` stage, starting from `S`, with stage vector `w`. -/
def stageSpace (S : Submodule 𝕜 V) (w : V) : StageKind → Submodule 𝕜 V
  | .arnoldi => S
  | .extend ins r => closeIter A r (if ins then S ⊔ span 𝕜 {w} else S)

/-- Run the stages `i, i+1, …` from the space `S`. -/
def run (w : ℕ → V) : Submodule 𝕜 V → ℕ → List StageKind → List (Submodule 𝕜 V)
  | _, _, [] => []
  | S, i, st :: rest => stageSpace A S (w i) st :: run w (stageSpace A S (w i) st) (i + 1) rest

/-- Spaces `K₁, K₂, …` of a schedule starting with an Arnoldi stage of
dimension `κ` (stage vectors `w 1, w 2, …`). -/
def mapSpaces (κ : ℕ) (w : ℕ → V) : List StageKind → List (Submodule 𝕜 V)
  | .arnoldi :: rest => krylov A (w 1) κ :: run A w (krylov A (w 1) κ) 2 rest
  | _ => []

section GRK4T

variable (κ : ℕ) (w : ℕ → V)

/-- `K₁`. -/
abbrev K1 : Submodule 𝕜 V := krylov A (w 1) κ
/-- `K₁ + ⟨w₂⟩`. -/
abbrev S2 : Submodule 𝕜 V := K1 A κ w ⊔ span 𝕜 {w 2}
/-- `K₂ = C(K₁ + ⟨w₂⟩)`. -/
abbrev K2 : Submodule 𝕜 V := close A (S2 A κ w)
/-- `K₃ = K₂ + ⟨w₃⟩`. -/
abbrev K3 : Submodule 𝕜 V := K2 A κ w ⊔ span 𝕜 {w 3}
/-- `K₄ = C(K₃ + ⟨w₄⟩)`. -/
abbrev K4 : Submodule 𝕜 V := close A (K3 A κ w ⊔ span 𝕜 {w 4})

/-- The generic schedule semantics yields exactly `K₁, …, K₄`. -/
theorem mapSpaces_grk4t :
    mapSpaces A κ w grk4tAutonomous = [K1 A κ w, K2 A κ w, K3 A κ w, K4 A κ w] := by
  simp [mapSpaces, run, stageSpace, closeIter, grk4tAutonomous]

/-! ### F1: `w_i ∈ K_i` -/

theorem w1_mem (hκ : 1 ≤ κ) : w 1 ∈ K1 A κ w :=
  subset_span ⟨⟨0, hκ⟩, rfl⟩

theorem w2_mem : w 2 ∈ K2 A κ w :=
  mem_sup_left (mem_sup_right (mem_span_singleton_self _))

theorem w3_mem : w 3 ∈ K3 A κ w :=
  mem_sup_right (mem_span_singleton_self _)

theorem w4_mem : w 4 ∈ K4 A κ w :=
  mem_sup_left (mem_sup_right (mem_span_singleton_self _))

theorem K1_le_K2 : K1 A κ w ≤ K2 A κ w := le_sup_left.trans le_sup_left

theorem K2_le_K3 : K2 A κ w ≤ K3 A κ w := le_sup_left

/-- `w₄ ∈ K₃` for GRK4T: `f₄ = f₃` gives `w₄ = w₃ + Σ_j (c₄ⱼ − c₃ⱼ) k_j + c₄₃ k₃`,
and the projected stage solutions satisfy `k_j ∈ K_j` (docs/map.md, Lemma 1). -/
theorem w4_mem_K3 (k1 k2 k3 : V) (hk1 : k1 ∈ K1 A κ w) (hk2 : k2 ∈ K2 A κ w)
    (hk3 : k3 ∈ K3 A κ w) (c1 c2 c3 : 𝕜) (h : w 4 = w 3 + c1 • k1 + c2 • k2 + c3 • k3) :
    w 4 ∈ K3 A κ w := by
  rw [h]
  have h1 : k1 ∈ K3 A κ w := K2_le_K3 A κ w (K1_le_K2 A κ w hk1)
  have h2 : k2 ∈ K3 A κ w := K2_le_K3 A κ w hk2
  exact add_mem (add_mem (add_mem (w3_mem A κ w) (smul_mem _ _ h1)) (smul_mem _ _ h2))
    (smul_mem _ _ hk3)

/-! ### F2: Krylov powers -/

theorem pow_w1_mem (j : ℕ) (hj : j < κ) : (⇑A)^[j] (w 1) ∈ K1 A κ w :=
  subset_span ⟨⟨j, hj⟩, rfl⟩

/-! ### F3, F4: closure properties -/

theorem F3 : (S2 A κ w).map A ≤ K2 A κ w := le_sup_right

theorem A_w2_mem : A (w 2) ∈ K2 A κ w :=
  F3 A κ w (mem_map_of_mem (mem_sup_right (mem_span_singleton_self _)))

theorem F4 : (K3 A κ w ⊔ span 𝕜 {w 4}).map A ≤ K4 A κ w := le_sup_right

theorem A_K3_le_K4 : (K3 A κ w).map A ≤ K4 A κ w :=
  (map_mono le_sup_left).trans (F4 A κ w)

theorem A_A_w2_mem : A (A (w 2)) ∈ K4 A κ w :=
  A_K3_le_K4 A κ w (mem_map_of_mem (K2_le_K3 A κ w (A_w2_mem A κ w)))

/-! ### F5: dimension increments (κ₂ ≤ κ₁+3, κ₃ ≤ κ₂+1, κ₄ ≤ κ₃+3) -/

/-- Krylov shift: `A K₁ ⊆ K₁ + ⟨A^κ w⟩`. -/
theorem map_krylov_le (v : V) :
    (krylov A v κ).map A ≤ krylov A v κ ⊔ span 𝕜 {(⇑A)^[κ] v} := by
  rw [krylov, map_span, span_le]
  rintro _ ⟨_, ⟨j, rfl⟩, rfl⟩
  simp only [SetLike.mem_coe]
  have hs : A ((⇑A)^[(j : ℕ)] v) = (⇑A)^[(j : ℕ) + 1] v :=
    (Function.iterate_succ_apply' _ _ _).symm
  rw [hs]
  by_cases h : (j : ℕ) + 1 < κ
  · exact mem_sup_left (subset_span ⟨⟨(j : ℕ) + 1, h⟩, rfl⟩)
  · have : (j : ℕ) + 1 = κ := by omega
    rw [this]
    exact mem_sup_right (mem_span_singleton_self _)

theorem map_span_singleton (v : V) : (span 𝕜 {v}).map A = span 𝕜 {A v} := by
  rw [map_span, Set.image_singleton]

variable [FiniteDimensional 𝕜 V]

theorem finrank_sup_span_le (S : Submodule 𝕜 V) (s : Finset V) :
    finrank 𝕜 ↥(S ⊔ span 𝕜 (s : Set V)) ≤ finrank 𝕜 S + s.card :=
  (finrank_add_le_finrank_add_finrank _ _).trans
    (Nat.add_le_add_left (finrank_span_finset_le_card s) _)

theorem finrank_le_add_three (S T : Submodule 𝕜 V) (a b c : V)
    (h : T ≤ S ⊔ span 𝕜 {a, b, c}) : finrank 𝕜 T ≤ finrank 𝕜 S + 3 := by
  classical
  have hT : T ≤ S ⊔ span 𝕜 (({a, b, c} : Finset V) : Set V) := by simpa using h
  exact (finrank_mono hT).trans
    ((finrank_sup_span_le S _).trans (Nat.add_le_add_left Finset.card_le_three _))

omit [FiniteDimensional 𝕜 V] in
theorem K2_le : K2 A κ w ≤ K1 A κ w ⊔ span 𝕜 {w 2, (⇑A)^[κ] (w 1), A (w 2)} := by
  have hK1 := map_krylov_le A κ (w 1)
  apply sup_le
  · exact sup_le le_sup_left (le_sup_of_le_right (span_mono (by simp)))
  · rw [Submodule.map_sup, map_span_singleton]
    apply sup_le
    · exact hK1.trans (sup_le le_sup_left (le_sup_of_le_right (span_mono (by simp))))
    · exact le_sup_of_le_right (span_mono (by simp))

/-- `κ₂ ≤ κ₁ + 3`. -/
theorem finrank_K2 : finrank 𝕜 (K2 A κ w) ≤ finrank 𝕜 (K1 A κ w) + 3 :=
  finrank_le_add_three _ _ _ _ _ (K2_le A κ w)

/-- `κ₃ ≤ κ₂ + 1`. -/
theorem finrank_K3 : finrank 𝕜 (K3 A κ w) ≤ finrank 𝕜 (K2 A κ w) + 1 := by
  classical
  have h := finrank_sup_span_le (K2 A κ w) ({w 3} : Finset V)
  rw [Finset.coe_singleton, Finset.card_singleton] at h
  exact h

omit [FiniteDimensional 𝕜 V] in
theorem A_K3_le : (K3 A κ w).map A ≤
    K3 A κ w ⊔ span 𝕜 {(⇑A)^[κ + 1] (w 1), A (A (w 2)), A (w 3)} := by
  have hshift : A ((⇑A)^[κ] (w 1)) = (⇑A)^[κ + 1] (w 1) :=
    (Function.iterate_succ_apply' _ _ _).symm
  set R := K3 A κ w ⊔ span 𝕜 {(⇑A)^[κ + 1] (w 1), A (A (w 2)), A (w 3)} with hR
  have hK3R : K3 A κ w ≤ R := le_sup_left
  have hspan : ∀ v ∈ ({(⇑A)^[κ + 1] (w 1), A (A (w 2)), A (w 3)} : Set V),
      span 𝕜 {v} ≤ R := fun v hv =>
    le_sup_of_le_right (span_mono (Set.singleton_subset_iff.mpr hv))
  -- A S₂ ⊆ K₂ ⊆ K₃
  have hS2 : (S2 A κ w).map A ≤ R := (F3 A κ w).trans ((K2_le_K3 A κ w).trans hK3R)
  have hK1 : (K1 A κ w).map A ≤ R := (map_mono le_sup_left).trans hS2
  -- A (A K₁) ⊆ A K₁ + ⟨A^{κ+1} w₁⟩
  have hAAK1 : ((K1 A κ w).map A).map A ≤ R := by
    refine (map_mono (map_krylov_le A κ (w 1))).trans ?_
    rw [Submodule.map_sup, map_span_singleton, hshift]
    exact sup_le hK1 (hspan _ (by simp))
  have hAS2 : ((S2 A κ w).map A).map A ≤ R := by
    rw [S2, Submodule.map_sup, map_span_singleton, Submodule.map_sup, map_span_singleton]
    exact sup_le hAAK1 (hspan _ (by simp))
  have hK2 : (K2 A κ w).map A ≤ R := by
    rw [K2, close, Submodule.map_sup]
    exact sup_le hS2 hAS2
  rw [K3, Submodule.map_sup, map_span_singleton]
  exact sup_le hK2 (hspan _ (by simp))

/-- `κ₄ ≤ κ₃ + 3`, given `w₄ ∈ K₃` (`w4_mem_K3`). -/
theorem finrank_K4 (hw4 : w 4 ∈ K3 A κ w) :
    finrank 𝕜 (K4 A κ w) ≤ finrank 𝕜 (K3 A κ w) + 3 := by
  have hS : K3 A κ w ⊔ span 𝕜 {w 4} = K3 A κ w :=
    sup_eq_left.mpr (span_le.mpr (by simpa using hw4))
  refine finrank_le_add_three _ _ ((⇑A)^[κ + 1] (w 1)) (A (A (w 2))) (A (w 3)) ?_
  rw [K4, close, hS]
  exact sup_le le_sup_left (A_K3_le A κ w)

end GRK4T

/-! ## Non-autonomous schedule (docs/map.md §7)

`K₂ = C²(K₁ + ⟨w₂⟩)`, `K₃ = K₂ + ⟨w₃⟩`, `K₄ = C(K₃ + ⟨w₄⟩)`. Here `w₄ ∉ K₃` in
general (`w₄ − w₃ ∋ h(d₄ − d₃) f_t`), so no hypothesis on `w₄`. -/

/-- The MAP schedule for GRK4T, non-autonomous case. -/
def grk4tNonautonomous : List StageKind :=
  [.arnoldi, .extend true 2, .extend true 0, .extend true 1]

section GRK4TNonautonomous

variable (κ : ℕ) (w : ℕ → V)

/-- `K₂ = C²(K₁ + ⟨w₂⟩)`. -/
abbrev K2n : Submodule 𝕜 V := close A (close A (S2 A κ w))
/-- `K₃ = K₂ + ⟨w₃⟩`. -/
abbrev K3n : Submodule 𝕜 V := K2n A κ w ⊔ span 𝕜 {w 3}
/-- `K₄ = C(K₃ + ⟨w₄⟩)`. -/
abbrev K4n : Submodule 𝕜 V := close A (K3n A κ w ⊔ span 𝕜 {w 4})

theorem mapSpaces_grk4t_nonautonomous :
    mapSpaces A κ w grk4tNonautonomous = [K1 A κ w, K2n A κ w, K3n A κ w, K4n A κ w] := by
  simp [mapSpaces, run, stageSpace, closeIter, grk4tNonautonomous]

theorem S2_le_K2n : S2 A κ w ≤ K2n A κ w := le_sup_left.trans le_sup_left

theorem AS2_le_K2n : (S2 A κ w).map A ≤ K2n A κ w := le_sup_right.trans le_sup_left

/-- Second closure round: `A²(K₁ + ⟨w₂⟩) ⊆ K₂` (needed for `A² f_t`). -/
theorem AAS2_le_K2n : ((S2 A κ w).map A).map A ≤ K2n A κ w :=
  (map_mono le_sup_right).trans le_sup_right

theorem n_w2_mem : w 2 ∈ K2n A κ w :=
  S2_le_K2n A κ w (mem_sup_right (mem_span_singleton_self _))

theorem n_w3_mem : w 3 ∈ K3n A κ w := mem_sup_right (mem_span_singleton_self _)

theorem n_w4_mem : w 4 ∈ K4n A κ w := mem_sup_left (mem_sup_right (mem_span_singleton_self _))

theorem n_F4 : (K3n A κ w ⊔ span 𝕜 {w 4}).map A ≤ K4n A κ w := le_sup_right

theorem n_A_K3_le_K4 : (K3n A κ w).map A ≤ K4n A κ w := (map_mono le_sup_left).trans (n_F4 A κ w)

theorem n_AA_w2_mem : A (A (w 2)) ∈ K2n A κ w :=
  AAS2_le_K2n A κ w (mem_map_of_mem (mem_map_of_mem (mem_sup_right (mem_span_singleton_self _))))

/-- `C²(S) = S ⊔ AS ⊔ (AS ⊔ A²S)`. -/
theorem close_close_eq (S : Submodule 𝕜 V) :
    close A (close A S) = S ⊔ S.map A ⊔ (S.map A ⊔ (S.map A).map A) := by
  simp only [close, Submodule.map_sup]

/-- `A C²(S) = AS ⊔ A²S ⊔ (A²S ⊔ A³S)`. -/
theorem map_close_close_eq (S : Submodule 𝕜 V) :
    (close A (close A S)).map A =
      S.map A ⊔ (S.map A).map A ⊔ ((S.map A).map A ⊔ ((S.map A).map A).map A) := by
  rw [close_close_eq]
  simp only [Submodule.map_sup]

/-- `S ≤ T ⊔ ⟨s⟩` implies `A S ≤ A T ⊔ ⟨A s⟩`. -/
theorem map_le_sup_span {S T : Submodule 𝕜 V} {s : Set V} (h : S ≤ T ⊔ span 𝕜 s) :
    S.map A ≤ T.map A ⊔ span 𝕜 (A '' s) := by
  refine (map_mono h).trans ?_
  rw [Submodule.map_sup, map_span]

/-- Dimension bound with a list of extra vectors. -/
theorem finrank_le_add_list [FiniteDimensional 𝕜 V] (S T : Submodule 𝕜 V) (l : List V)
    (h : T ≤ S ⊔ span 𝕜 {x | x ∈ l}) : finrank 𝕜 T ≤ finrank 𝕜 S + l.length := by
  classical
  have hT : T ≤ S ⊔ span 𝕜 (l.toFinset : Set V) := by simpa using h
  exact (finrank_mono hT).trans
    ((finrank_sup_span_le S _).trans (Nat.add_le_add_left (List.toFinset_card_le l) _))

theorem K2n_le : K2n A κ w ≤ K1 A κ w ⊔ span 𝕜
    {x | x ∈ [w 2, (⇑A)^[κ] (w 1), A (w 2), (⇑A)^[κ + 1] (w 1), A (A (w 2))]} := by
  have hshift : A ((⇑A)^[κ] (w 1)) = (⇑A)^[κ + 1] (w 1) :=
    (Function.iterate_succ_apply' _ _ _).symm
  set R := K1 A κ w ⊔ span 𝕜
    {x | x ∈ [w 2, (⇑A)^[κ] (w 1), A (w 2), (⇑A)^[κ + 1] (w 1), A (A (w 2))]} with hR
  have hsp : ∀ v, v ∈ [w 2, (⇑A)^[κ] (w 1), A (w 2), (⇑A)^[κ + 1] (w 1), A (A (w 2))] →
      span 𝕜 {v} ≤ R := fun v hv =>
    le_sup_of_le_right (span_mono (Set.singleton_subset_iff.mpr hv))
  have hS : S2 A κ w ≤ R := sup_le le_sup_left (hsp _ (by simp))
  have hK1 : (K1 A κ w).map A ≤ R :=
    (map_krylov_le A κ (w 1)).trans (sup_le le_sup_left (hsp _ (by simp)))
  have hAS : (S2 A κ w).map A ≤ R := by
    rw [S2, Submodule.map_sup, map_span_singleton]
    exact sup_le hK1 (hsp _ (by simp))
  have hAAK1 : ((K1 A κ w).map A).map A ≤ R := by
    refine (map_mono (map_krylov_le A κ (w 1))).trans ?_
    rw [Submodule.map_sup, map_span_singleton, hshift]
    exact sup_le hK1 (hsp _ (by simp))
  have hAAS : ((S2 A κ w).map A).map A ≤ R := by
    rw [S2, Submodule.map_sup, map_span_singleton, Submodule.map_sup, map_span_singleton]
    exact sup_le hAAK1 (hsp _ (by simp))
  -- C²(S) = S ⊔ AS ⊔ (AS ⊔ AAS)
  rw [K2n, close_close_eq]
  exact sup_le (sup_le hS hAS) (sup_le hAS hAAS)

/-- `κ₂ ≤ κ₁ + 5`. -/
theorem finrank_K2n [FiniteDimensional 𝕜 V] :
    finrank 𝕜 (K2n A κ w) ≤ finrank 𝕜 (K1 A κ w) + 5 :=
  finrank_le_add_list _ _ _ (K2n_le A κ w)

/-- `κ₃ ≤ κ₂ + 1`. -/
theorem finrank_K3n [FiniteDimensional 𝕜 V] :
    finrank 𝕜 (K3n A κ w) ≤ finrank 𝕜 (K2n A κ w) + 1 := by
  classical
  have h := finrank_sup_span_le (K2n A κ w) ({w 3} : Finset V)
  rw [Finset.coe_singleton, Finset.card_singleton] at h
  exact h

/-- `A^κ w₁ ∈ A K₁` (for κ ≥ 1). -/
theorem pow_kappa_mem (hκ : 1 ≤ κ) : (⇑A)^[κ] (w 1) ∈ (K1 A κ w).map A := by
  obtain ⟨k, rfl⟩ : ∃ k, κ = k + 1 := ⟨κ - 1, by omega⟩
  rw [Function.iterate_succ_apply']
  exact mem_map_of_mem (pow_w1_mem A (k + 1) w k (by omega))

theorem A_K3n_le (hκ : 1 ≤ κ) : (K3n A κ w).map A ≤
    K3n A κ w ⊔ span 𝕜 {x | x ∈ [(⇑A)^[κ + 2] (w 1), A (A (A (w 2))), A (w 3)]} := by
  have h1 : A ((⇑A)^[κ] (w 1)) = (⇑A)^[κ + 1] (w 1) :=
    (Function.iterate_succ_apply' _ _ _).symm
  have h2 : A ((⇑A)^[κ + 1] (w 1)) = (⇑A)^[κ + 2] (w 1) :=
    (Function.iterate_succ_apply' _ _ _).symm
  set R := K3n A κ w ⊔ span 𝕜 {x | x ∈ [(⇑A)^[κ + 2] (w 1), A (A (A (w 2))), A (w 3)]}
    with hR
  have hK3R : K3n A κ w ≤ R := le_sup_left
  have hK2R : K2n A κ w ≤ R := le_sup_left.trans hK3R
  have hsp : ∀ v, v ∈ [(⇑A)^[κ + 2] (w 1), A (A (A (w 2))), A (w 3)] → span 𝕜 {v} ≤ R :=
    fun v hv => le_sup_of_le_right (span_mono (Set.singleton_subset_iff.mpr hv))
  have hAS : (S2 A κ w).map A ≤ R := (AS2_le_K2n A κ w).trans hK2R
  have hAAS : ((S2 A κ w).map A).map A ≤ R := (AAS2_le_K2n A κ w).trans hK2R
  -- A³ S₂ = A³ K₁ ⊔ ⟨A³ w₂⟩, and A³ K₁ ≤ A² K₁ ⊔ ⟨A^{κ+2} w₁⟩ with A² K₁ ≤ A² S₂
  have hAAK1 : ((K1 A κ w).map A).map A ≤ R :=
    (map_mono (map_mono le_sup_left)).trans hAAS
  have hx2 : (⇑A)^[κ + 1] (w 1) ∈ R := by
    rw [← h1]
    exact hAAS (mem_map_of_mem (map_mono le_sup_left (pow_kappa_mem A κ w hκ)))
  have hAAAK1 : (((K1 A κ w).map A).map A).map A ≤ R := by
    have hA2 : ((K1 A κ w).map A).map A ≤ (K1 A κ w).map A ⊔ span 𝕜 {(⇑A)^[κ + 1] (w 1)} := by
      have := map_le_sup_span A (map_krylov_le A κ (w 1))
      rwa [Set.image_singleton, h1] at this
    have := map_le_sup_span A hA2
    rw [Set.image_singleton, h2] at this
    exact this.trans (sup_le hAAK1 (hsp _ (by simp)))
  have hAAAS : (((S2 A κ w).map A).map A).map A ≤ R := by
    rw [S2, Submodule.map_sup, map_span_singleton, Submodule.map_sup, map_span_singleton,
      Submodule.map_sup, map_span_singleton]
    exact sup_le hAAAK1 (hsp _ (by simp))
  have hAK2 : (K2n A κ w).map A ≤ R := by
    rw [K2n, map_close_close_eq]
    exact sup_le (sup_le hAS hAAS) (sup_le hAAS hAAAS)
  rw [K3n, Submodule.map_sup, map_span_singleton]
  exact sup_le hAK2 (hsp _ (by simp))

/-- `κ₄ ≤ κ₃ + 5` (no assumption on `w₄`). -/
theorem finrank_K4n [FiniteDimensional 𝕜 V] (hκ : 1 ≤ κ) :
    finrank 𝕜 (K4n A κ w) ≤ finrank 𝕜 (K3n A κ w) + 5 := by
  refine finrank_le_add_list _ _
    [w 4, A (w 4), (⇑A)^[κ + 2] (w 1), A (A (A (w 2))), A (w 3)] ?_
  have hA := A_K3n_le A κ w hκ
  set R := K3n A κ w ⊔ span 𝕜
    {x | x ∈ [w 4, A (w 4), (⇑A)^[κ + 2] (w 1), A (A (A (w 2))), A (w 3)]} with hR
  have hsp : ∀ v, v ∈ [w 4, A (w 4), (⇑A)^[κ + 2] (w 1), A (A (A (w 2))), A (w 3)] →
      span 𝕜 {v} ≤ R := fun v hv =>
    le_sup_of_le_right (span_mono (Set.singleton_subset_iff.mpr hv))
  have hK3 : (K3n A κ w).map A ≤ R := hA.trans (sup_le le_sup_left
    (span_le.mpr (by
      intro v hv
      exact (hsp v (by simp only [List.mem_cons] at hv ⊢; tauto)) (mem_span_singleton_self _))))
  rw [K4n, close, Submodule.map_sup, map_span_singleton]
  exact sup_le (sup_le le_sup_left (hsp _ (by simp))) (sup_le hK3 (hsp _ (by simp)))

end GRK4TNonautonomous

/-!
### Boundary of the formalization

The analytic step is not formalized: "a Krylov-W method based on GRK4T whose
subspaces satisfy F1–F4 with κ ≥ 4 has classical order 4" (docs/map.md §4,
`OrderFromIdentities` in AGENTS.md). Stating it meaningfully would need the
ODE/Taylor setting. A placeholder `Prop` would be vacuous, so none is given.
The step is covered by the paper proof and by the numerical order tests in
`src/map/tests.rs`.
-/

end Rowmap
