/-
Copyright (c) 2026 Helmut Podhaisky. All rights reserved.
Released under MIT license as described in the file LICENSE.
Authors: Helmut Podhaisky
-/
import RowmapLean.Map

/-!
# Tightness of the non-autonomous dimension bounds

The shift `e_i ↦ e_{i+1}` on `ℚ¹⁵` with `w₁ = e₀`, `w₂ = e₇`, `w₃ = e₁₁`,
`w₄ = e₁₃` (so `w₄ ∉ K₃`) and `κ = 4` attains
`dim K = 4, 9, 10, 15` (paper §4: κ₂ = κ₁+5, κ₃ = κ₂+1, κ₄ = κ₃+5).
-/

namespace Rowmap.TightnessNA

open Submodule Module

/-- Standard basis vector `e_i` of `ℚ¹⁵`. -/
abbrev e (i : ℕ) (h : i < 15 := by omega) : Fin 15 → ℚ := Pi.single (⟨i, h⟩ : Fin 15) 1

def shiftMat : Matrix (Fin 15) (Fin 15) ℚ := fun i j => if (i : ℕ) = (j : ℕ) + 1 then 1 else 0

def shift : (Fin 15 → ℚ) →ₗ[ℚ] (Fin 15 → ℚ) := Matrix.mulVecLin shiftMat

theorem shift_e (i : ℕ) (hi : i + 1 < 15) : shift (e i) = e (i + 1) hi := by
  funext j
  simp [shift, shiftMat, Pi.single_apply, Fin.ext_iff]

theorem iterate_shift_e0 (j : ℕ) (hj : j < 15) : (⇑shift)^[j] (e 0) = e j hj := by
  induction j with
  | zero => rfl
  | succ k ih => rw [Function.iterate_succ_apply', ih (by omega), shift_e k hj]

def wv : ℕ → Fin 15 → ℚ
  | 1 => e 0
  | 2 => e 7
  | 3 => e 11
  | 4 => e 13
  | _ => 0

theorem card_le_finrank {m : ℕ} (σ : Fin m → Fin 15) (hσ : Function.Injective σ)
    (S : Submodule ℚ (Fin 15 → ℚ)) (hS : ∀ i, Pi.single (σ i) (1 : ℚ) ∈ S) :
    m ≤ finrank ℚ S := by
  let b : Fin m → Fin 15 → ℚ := fun i => Pi.single (σ i) 1
  have hli : LinearIndependent ℚ b := by
    have h := (Pi.basisFun ℚ (Fin 15)).linearIndependent.comp σ hσ
    convert h using 1
    funext i
    simp [b, Pi.basisFun_apply]
  have hle : span ℚ (Set.range b) ≤ S := span_le.mpr (by rintro _ ⟨i, rfl⟩; exact hS i)
  have h1 := finrank_span_eq_card hli
  have h2 := finrank_mono hle
  rw [Fintype.card_fin] at h1
  omega

local notation "𝐊₁" => K1 shift 4 wv
local notation "𝐒₂" => S2 shift 4 wv
local notation "𝐊₂" => K2n shift 4 wv
local notation "𝐊₃" => K3n shift 4 wv
local notation "𝐊₄" => K4n shift 4 wv

theorem e_mem_K1 (i : ℕ) (hi : i < 4) : e i ∈ 𝐊₁ := by
  have := pow_w1_mem shift 4 wv i hi
  rwa [show wv 1 = e 0 from rfl, iterate_shift_e0 i (by omega)] at this

theorem K1_S2 {v : Fin 15 → ℚ} (h : v ∈ 𝐊₁) : v ∈ 𝐒₂ := mem_sup_left h
theorem e7_S2 : e 7 ∈ 𝐒₂ := mem_sup_right (mem_span_singleton_self _)

theorem AS2 {v : Fin 15 → ℚ} (h : v ∈ 𝐒₂) : shift v ∈ 𝐊₂ := AS2_le_K2n shift 4 wv (mem_map_of_mem h)
theorem AAS2 {v : Fin 15 → ℚ} (h : v ∈ 𝐒₂) : shift (shift v) ∈ 𝐊₂ :=
  AAS2_le_K2n shift 4 wv (mem_map_of_mem (mem_map_of_mem h))
theorem S2_K2 {v : Fin 15 → ℚ} (h : v ∈ 𝐒₂) : v ∈ 𝐊₂ := S2_le_K2n shift 4 wv h
theorem K2_K3 {v : Fin 15 → ℚ} (h : v ∈ 𝐊₂) : v ∈ 𝐊₃ := mem_sup_left h
theorem K3_K4 {v : Fin 15 → ℚ} (h : v ∈ 𝐊₃) : v ∈ 𝐊₄ :=
  (le_sup_left.trans le_sup_left : 𝐊₃ ≤ 𝐊₄) h
theorem AK3 {v : Fin 15 → ℚ} (h : v ∈ 𝐊₃) : shift v ∈ 𝐊₄ :=
  n_A_K3_le_K4 shift 4 wv (mem_map_of_mem h)
theorem e13_K4 : e 13 ∈ 𝐊₄ := n_w4_mem shift 4 wv
theorem e14_K4 : e 14 ∈ 𝐊₄ := by
  rw [← shift_e 13 (by omega)]
  exact n_F4 shift 4 wv (mem_map_of_mem (mem_sup_right (mem_span_singleton_self _)))

/-- Members of `K₂`: `e₀..e₅, e₇, e₈, e₉`. -/
theorem K2_members : ∀ i ∈ ([0, 1, 2, 3, 4, 5, 7, 8, 9] : List ℕ), ∀ h : i < 15, e i h ∈ 𝐊₂ := by
  intro i hi h
  simp only [List.mem_cons, List.not_mem_nil, or_false] at hi
  rcases hi with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · exact S2_K2 (K1_S2 (e_mem_K1 0 (by omega)))
  · exact S2_K2 (K1_S2 (e_mem_K1 1 (by omega)))
  · exact S2_K2 (K1_S2 (e_mem_K1 2 (by omega)))
  · exact S2_K2 (K1_S2 (e_mem_K1 3 (by omega)))
  · rw [← shift_e 3 (by omega)]; exact AS2 (K1_S2 (e_mem_K1 3 (by omega)))
  · rw [← shift_e 4 (by omega), ← shift_e 3 (by omega)]
    exact AAS2 (K1_S2 (e_mem_K1 3 (by omega)))
  · exact S2_K2 e7_S2
  · rw [← shift_e 7 (by omega)]; exact AS2 e7_S2
  · rw [← shift_e 8 (by omega), ← shift_e 7 (by omega)]; exact AAS2 e7_S2

theorem e11_K3 : e 11 ∈ 𝐊₃ := n_w3_mem shift 4 wv

theorem all_mem_K4 (i : Fin 15) : Pi.single i (1 : ℚ) ∈ 𝐊₄ := by
  have k2 := K2_members
  have h6 : e 6 ∈ 𝐊₄ := by
    rw [← shift_e 5 (by omega)]; exact AK3 (K2_K3 (k2 5 (by simp) (by omega)))
  have h10 : e 10 ∈ 𝐊₄ := by
    rw [← shift_e 9 (by omega)]; exact AK3 (K2_K3 (k2 9 (by simp) (by omega)))
  have h12 : e 12 ∈ 𝐊₄ := by rw [← shift_e 11 (by omega)]; exact AK3 e11_K3
  fin_cases i
  · exact K3_K4 (K2_K3 (k2 0 (by simp) (by omega)))
  · exact K3_K4 (K2_K3 (k2 1 (by simp) (by omega)))
  · exact K3_K4 (K2_K3 (k2 2 (by simp) (by omega)))
  · exact K3_K4 (K2_K3 (k2 3 (by simp) (by omega)))
  · exact K3_K4 (K2_K3 (k2 4 (by simp) (by omega)))
  · exact K3_K4 (K2_K3 (k2 5 (by simp) (by omega)))
  · exact h6
  · exact K3_K4 (K2_K3 (k2 7 (by simp) (by omega)))
  · exact K3_K4 (K2_K3 (k2 8 (by simp) (by omega)))
  · exact K3_K4 (K2_K3 (k2 9 (by simp) (by omega)))
  · exact h10
  · exact K3_K4 e11_K3
  · exact h12
  · exact e13_K4
  · exact e14_K4

/-- The example attains every bound: dims 4, 9, 10, 15. -/
theorem tight :
    finrank ℚ 𝐊₁ = 4 ∧ finrank ℚ 𝐊₂ = 9 ∧ finrank ℚ 𝐊₃ = 10 ∧ finrank ℚ 𝐊₄ = 15 := by
  have k2 := K2_members
  have l1 : 4 ≤ finrank ℚ 𝐊₁ :=
    card_le_finrank (fun i : Fin 4 => ⟨i, by omega⟩) (fun a b h => by simpa [Fin.ext_iff] using h)
      _ (fun i => e_mem_K1 i i.isLt)
  have u1 : finrank ℚ 𝐊₁ ≤ 4 := by
    have h := finrank_range_le_card (R := ℚ) (fun j : Fin 4 => (⇑shift)^[j] (wv 1))
    rw [Fintype.card_fin] at h
    exact h
  have l2 : 9 ≤ finrank ℚ 𝐊₂ := by
    refine card_le_finrank ![0, 1, 2, 3, 4, 5, 7, 8, 9] (by decide) _ ?_
    intro i
    fin_cases i <;> exact k2 _ (by simp) (by omega)
  have l3 : 10 ≤ finrank ℚ 𝐊₃ := by
    refine card_le_finrank ![0, 1, 2, 3, 4, 5, 7, 8, 9, 11] (by decide) _ ?_
    intro i
    fin_cases i
    all_goals first
      | exact e11_K3
      | exact K2_K3 (k2 _ (by simp) (by omega))
  have l4 : 15 ≤ finrank ℚ 𝐊₄ := card_le_finrank id Function.injective_id _ all_mem_K4
  have u2 := finrank_K2n shift 4 wv
  have u3 := finrank_K3n shift 4 wv
  have u4 := finrank_K4n shift 4 wv (by omega)
  omega

end Rowmap.TightnessNA
