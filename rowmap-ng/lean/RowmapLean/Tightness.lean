/-
Copyright (c) 2026 Helmut Podhaisky. All rights reserved.
Released under MIT license as described in the file LICENSE.
Authors: Helmut Podhaisky
-/
import RowmapLean.Map

/-!
# Tightness of the dimension bounds (structural necessity)

A concrete example over `ℚ` where the bounds of F5 are attained:
`dim K₁ = 4, dim K₂ = 7, dim K₃ = 8, dim K₄ = 11`.

Since `K₂ = C(K₁ + ⟨w₂⟩)` is the *smallest* subspace that contains
`K₁ + ⟨w₂⟩` and its image under `A` (F1, F3), and likewise for `K₄` (F4),
no construction with fewer vectors can satisfy F1–F4 in general.
This is the structural necessity. Whether the *order* needs every vector is
an analytic question (docs/map.md §4: classically, the stage-4 closure is not
needed).

Example: `V = Fin 11 → ℚ`, `A` the shift `e_i ↦ e_{i+1}`, `w₁ = e₀`,
`w₂ = e₆`, `w₃ = e₉`, `w₄ = 0`, `κ = 4`.
-/

namespace Rowmap.Tightness

open Submodule Module

/-- Standard basis vector `e_i` of `ℚ¹¹`. -/
abbrev e (i : ℕ) (h : i < 11 := by omega) : Fin 11 → ℚ := Pi.single (⟨i, h⟩ : Fin 11) 1

/-- Shift matrix: `M i j = 1` iff `i = j + 1`. -/
def shiftMat : Matrix (Fin 11) (Fin 11) ℚ := fun i j => if (i : ℕ) = (j : ℕ) + 1 then 1 else 0

/-- The shift `e_j ↦ e_{j+1}` (and `e₁₀ ↦ 0`). -/
def shift : (Fin 11 → ℚ) →ₗ[ℚ] (Fin 11 → ℚ) := Matrix.mulVecLin shiftMat

theorem shift_e (i : ℕ) (hi : i + 1 < 11) : shift (e i) = e (i + 1) hi := by
  funext j
  simp [shift, shiftMat, Pi.single_apply, Fin.ext_iff]

theorem iterate_shift_e0 (j : ℕ) (hj : j < 11) : (⇑shift)^[j] (e 0) = e j hj := by
  induction j with
  | zero => rfl
  | succ k ih => rw [Function.iterate_succ_apply', ih (by omega), shift_e k hj]

/-- Stage vectors of the example. -/
def wv : ℕ → Fin 11 → ℚ
  | 1 => e 0
  | 2 => e 6
  | 3 => e 9
  | _ => 0

/-- `m` distinct standard basis vectors in `S` give `m ≤ dim S`. -/
theorem card_le_finrank {m : ℕ} (σ : Fin m → Fin 11) (hσ : Function.Injective σ)
    (S : Submodule ℚ (Fin 11 → ℚ)) (hS : ∀ i, Pi.single (σ i) (1 : ℚ) ∈ S) :
    m ≤ finrank ℚ S := by
  let b : Fin m → Fin 11 → ℚ := fun i => Pi.single (σ i) 1
  have hli : LinearIndependent ℚ b := by
    have h := (Pi.basisFun ℚ (Fin 11)).linearIndependent.comp σ hσ
    convert h using 1
    funext i
    simp [b, Pi.basisFun_apply]
  have hle : span ℚ (Set.range b) ≤ S := span_le.mpr (by rintro _ ⟨i, rfl⟩; exact hS i)
  have h1 := finrank_span_eq_card hli
  have h2 := finrank_mono hle
  rw [Fintype.card_fin] at h1
  omega

local notation "𝐊₁" => K1 shift 4 wv
local notation "𝐊₂" => K2 shift 4 wv
local notation "𝐊₃" => K3 shift 4 wv
local notation "𝐊₄" => K4 shift 4 wv

theorem e_mem_K1 (i : ℕ) (hi : i < 4) : e i ∈ 𝐊₁ := by
  have := pow_w1_mem shift 4 wv i hi
  rwa [show wv 1 = e 0 from rfl, iterate_shift_e0 i (by omega)] at this

theorem A_mem_K2 {v : Fin 11 → ℚ} (hv : v ∈ S2 shift 4 wv) : shift v ∈ 𝐊₂ :=
  F3 shift 4 wv (mem_map_of_mem hv)

theorem A_mem_K4 {v : Fin 11 → ℚ} (hv : v ∈ 𝐊₃) : shift v ∈ 𝐊₄ :=
  A_K3_le_K4 shift 4 wv (mem_map_of_mem hv)

theorem e6_mem_S2 : e 6 ∈ S2 shift 4 wv := mem_sup_right (mem_span_singleton_self _)

theorem e_mem_K2 (i : ℕ) (hi : i < 4) : e i ∈ 𝐊₂ := K1_le_K2 shift 4 wv (e_mem_K1 i hi)

theorem e4_mem_K2 : e 4 ∈ 𝐊₂ := by
  rw [← shift_e 3 (by omega)]
  exact A_mem_K2 (mem_sup_left (e_mem_K1 3 (by omega)))

theorem e6_mem_K2 : e 6 ∈ 𝐊₂ := w2_mem shift 4 wv

theorem e7_mem_K2 : e 7 ∈ 𝐊₂ := by
  rw [← shift_e 6 (by omega)]
  exact A_mem_K2 e6_mem_S2

theorem e9_mem_K3 : e 9 ∈ 𝐊₃ := w3_mem shift 4 wv

theorem K2_K3 {v : Fin 11 → ℚ} (h : v ∈ 𝐊₂) : v ∈ 𝐊₃ := K2_le_K3 shift 4 wv h

theorem K3_K4 {v : Fin 11 → ℚ} (h : v ∈ 𝐊₃) : v ∈ 𝐊₄ :=
  (le_sup_left.trans le_sup_left : 𝐊₃ ≤ 𝐊₄) h

/-- Membership of all 11 basis vectors in `K₄`. -/
theorem all_mem_K4 (i : Fin 11) : Pi.single i (1 : ℚ) ∈ 𝐊₄ := by
  have h5 : e 5 ∈ 𝐊₄ := by rw [← shift_e 4 (by omega)]; exact A_mem_K4 (K2_K3 e4_mem_K2)
  have h8 : e 8 ∈ 𝐊₄ := by rw [← shift_e 7 (by omega)]; exact A_mem_K4 (K2_K3 e7_mem_K2)
  have h10 : e 10 ∈ 𝐊₄ := by rw [← shift_e 9 (by omega)]; exact A_mem_K4 e9_mem_K3
  fin_cases i
  · exact K3_K4 (K2_K3 (e_mem_K2 0 (by omega)))
  · exact K3_K4 (K2_K3 (e_mem_K2 1 (by omega)))
  · exact K3_K4 (K2_K3 (e_mem_K2 2 (by omega)))
  · exact K3_K4 (K2_K3 (e_mem_K2 3 (by omega)))
  · exact K3_K4 (K2_K3 e4_mem_K2)
  · exact h5
  · exact K3_K4 (K2_K3 e6_mem_K2)
  · exact K3_K4 (K2_K3 e7_mem_K2)
  · exact h8
  · exact K3_K4 e9_mem_K3
  · exact h10

theorem finrank_K1_ge : 4 ≤ finrank ℚ 𝐊₁ :=
  card_le_finrank (fun i : Fin 4 => ⟨i, by omega⟩) (fun a b h => by simpa [Fin.ext_iff] using h)
    _ (fun i => e_mem_K1 i i.isLt)

theorem finrank_K1_le : finrank ℚ 𝐊₁ ≤ 4 := by
  have h := finrank_range_le_card (R := ℚ) (fun j : Fin 4 => (⇑shift)^[j] (wv 1))
  rw [Fintype.card_fin] at h
  exact h

/-- The example attains every bound of F5. -/
theorem tight :
    finrank ℚ 𝐊₁ = 4 ∧ finrank ℚ 𝐊₂ = 7 ∧ finrank ℚ 𝐊₃ = 8 ∧ finrank ℚ 𝐊₄ = 11 := by
  have h1 : finrank ℚ 𝐊₁ = 4 := le_antisymm finrank_K1_le finrank_K1_ge
  -- lower bounds from explicit independent basis vectors
  have l2 : 7 ≤ finrank ℚ 𝐊₂ := by
    refine card_le_finrank ![0, 1, 2, 3, 4, 6, 7] (by decide) _ ?_
    intro i
    fin_cases i
    · exact e_mem_K2 0 (by omega)
    · exact e_mem_K2 1 (by omega)
    · exact e_mem_K2 2 (by omega)
    · exact e_mem_K2 3 (by omega)
    · exact e4_mem_K2
    · exact e6_mem_K2
    · exact e7_mem_K2
  have l3 : 8 ≤ finrank ℚ 𝐊₃ := by
    refine card_le_finrank ![0, 1, 2, 3, 4, 6, 7, 9] (by decide) _ ?_
    intro i
    fin_cases i
    · exact K2_K3 (e_mem_K2 0 (by omega))
    · exact K2_K3 (e_mem_K2 1 (by omega))
    · exact K2_K3 (e_mem_K2 2 (by omega))
    · exact K2_K3 (e_mem_K2 3 (by omega))
    · exact K2_K3 e4_mem_K2
    · exact K2_K3 e6_mem_K2
    · exact K2_K3 e7_mem_K2
    · exact e9_mem_K3
  have l4 : 11 ≤ finrank ℚ 𝐊₄ := card_le_finrank id Function.injective_id _ all_mem_K4
  -- upper bounds from F5 (w₄ = 0 ∈ K₃)
  have u2 := finrank_K2 shift 4 wv
  have u3 := finrank_K3 shift 4 wv
  have u4 := finrank_K4 shift 4 wv (zero_mem _)
  omega

end Rowmap.Tightness
