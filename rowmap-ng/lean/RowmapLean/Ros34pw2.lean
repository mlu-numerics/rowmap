/-
Copyright (c) 2026 Helmut Podhaisky. All rights reserved.
Released under MIT license as described in the file LICENSE.
Authors: Helmut Podhaisky
-/
import RowmapLean.Map
import RowmapLean.Tightness

/-!
# MAP schedules for ROS34PW2 (order 3, phase 5b; docs/map.md §8)

- Autonomous: `K₂ = K₃ = K₄ = K₁ + ⟨w₂⟩`, increments `+1, +0, +0`
  (dims 3, 4, 4, 4).
- Non-autonomous: `K₂ = K₁ + ⟨w₂⟩`, `K₃ = K₄ = C(K₂ + ⟨w₃⟩)`, increments
  `+1, +4, +0` (dims 3, 4, 8, 8).

Plus a tightness example for the non-autonomous bound (shift on ℚ¹¹).
-/

namespace Rowmap

open Submodule Module

/-- ROS34PW2, autonomous. -/
def ros34pw2Autonomous : List StageKind :=
  [.arnoldi, .extend true 0, .extend false 0, .extend false 0]

/-- ROS34PW2, non-autonomous. -/
def ros34pw2Nonautonomous : List StageKind :=
  [.arnoldi, .extend true 0, .extend true 1, .extend false 0]

variable {𝕜 V : Type*} [Field 𝕜] [AddCommGroup V] [Module 𝕜 V]
variable (A : V →ₗ[𝕜] V) (κ : ℕ) (w : ℕ → V)

/-- ROS34PW2 `K₃ = K₄`, non-autonomous: `C(K₁ + ⟨w₂⟩ + ⟨w₃⟩)`. -/
abbrev P3n : Submodule 𝕜 V := close A (S2 A κ w ⊔ span 𝕜 {w 3})

theorem mapSpaces_ros34pw2 :
    mapSpaces A κ w ros34pw2Autonomous = [K1 A κ w, S2 A κ w, S2 A κ w, S2 A κ w] := by
  simp [mapSpaces, run, stageSpace, closeIter, ros34pw2Autonomous]

theorem mapSpaces_ros34pw2_nonautonomous :
    mapSpaces A κ w ros34pw2Nonautonomous =
      [K1 A κ w, S2 A κ w, P3n A κ w, P3n A κ w] := by
  simp [mapSpaces, run, stageSpace, closeIter, ros34pw2Nonautonomous]

theorem ros34pw2_w2_mem : w 2 ∈ S2 A κ w := mem_sup_right (mem_span_singleton_self _)

theorem ros34pw2_n_w3_mem : w 3 ∈ P3n A κ w :=
  mem_sup_left (mem_sup_right (mem_span_singleton_self _))

/-- Closure: `A(K₂ + ⟨w₃⟩) ⊆ K₃` (non-autonomous; needed for `A f_t`). -/
theorem ros34pw2_n_closure : (S2 A κ w ⊔ span 𝕜 {w 3}).map A ≤ P3n A κ w := le_sup_right

theorem P3n_le : P3n A κ w ≤
    S2 A κ w ⊔ span 𝕜 {x | x ∈ [w 3, (⇑A)^[κ] (w 1), A (w 2), A (w 3)]} := by
  set R := S2 A κ w ⊔ span 𝕜 {x | x ∈ [w 3, (⇑A)^[κ] (w 1), A (w 2), A (w 3)]} with hR
  have hsp : ∀ v, v ∈ [w 3, (⇑A)^[κ] (w 1), A (w 2), A (w 3)] → span 𝕜 {v} ≤ R :=
    fun v hv => le_sup_of_le_right (span_mono (Set.singleton_subset_iff.mpr hv))
  have hS : S2 A κ w ≤ R := le_sup_left
  have hK1 : (K1 A κ w).map A ≤ R :=
    (map_krylov_le A κ (w 1)).trans (sup_le (le_sup_left.trans hS) (hsp _ (by simp)))
  apply sup_le
  · exact sup_le hS (hsp _ (by simp))
  · rw [Submodule.map_sup, map_span_singleton, S2, Submodule.map_sup, map_span_singleton]
    exact sup_le (sup_le hK1 (hsp _ (by simp))) (hsp _ (by simp))

variable [FiniteDimensional 𝕜 V]

/-- `κ₂ ≤ κ₁ + 1` (both cases; autonomous `κ₃ = κ₄ = κ₂`). -/
theorem finrank_P2 : finrank 𝕜 (S2 A κ w) ≤ finrank 𝕜 (K1 A κ w) + 1 := by
  classical
  have h := finrank_sup_span_le (K1 A κ w) ({w 2} : Finset V)
  rw [Finset.coe_singleton, Finset.card_singleton] at h
  exact h

/-- Non-autonomous: `κ₃ ≤ κ₂ + 4` (and `κ₄ = κ₃`). -/
theorem finrank_P3n : finrank 𝕜 (P3n A κ w) ≤ finrank 𝕜 (S2 A κ w) + 4 :=
  finrank_le_add_list _ _ _ (P3n_le A κ w)

end Rowmap

/-! ### Tightness of the non-autonomous ROS34PW2 bounds: dims 3, 4, 8 -/

namespace Rowmap.TightnessRos34pw2

open Submodule Module Rowmap.Tightness

/-- `w₁ = e₀`, `w₂ = e₄`, `w₃ = e₇` on ℚ¹¹ with the shift (κ = 3). -/
def wq : ℕ → Fin 11 → ℚ
  | 1 => e 0
  | 2 => e 4
  | 3 => e 7
  | _ => 0

theorem e_mem_K1 (i : ℕ) (hi : i < 3) : e i ∈ K1 shift 3 wq := by
  have := pow_w1_mem shift 3 wq i hi
  rwa [show wq 1 = e 0 from rfl, iterate_shift_e0 i (by omega)] at this

theorem tight :
    finrank ℚ (K1 shift 3 wq) = 3 ∧ finrank ℚ (S2 shift 3 wq) = 4 ∧
      finrank ℚ (P3n shift 3 wq) = 8 := by
  -- members of K₂ + ⟨w₃⟩: e₀ e₁ e₂ (K₁), e₄ (w₂), e₇ (w₃)
  have inS : ∀ {v}, v ∈ S2 shift 3 wq → v ∈ S2 shift 3 wq ⊔ span ℚ {wq 3} :=
    fun h => mem_sup_left h
  have k1 : ∀ i (hi : i < 3), e i ∈ S2 shift 3 wq := fun i hi => mem_sup_left (e_mem_K1 i hi)
  have e4 : e 4 ∈ S2 shift 3 wq := mem_sup_right (mem_span_singleton_self _)
  have e7 : e 7 ∈ S2 shift 3 wq ⊔ span ℚ {wq 3} := mem_sup_right (mem_span_singleton_self _)
  have hS : ∀ {v}, v ∈ S2 shift 3 wq ⊔ span ℚ {wq 3} → v ∈ P3n shift 3 wq :=
    fun h => mem_sup_left h
  have hA : ∀ {v}, v ∈ S2 shift 3 wq ⊔ span ℚ {wq 3} → shift v ∈ P3n shift 3 wq :=
    fun h => ros34pw2_n_closure shift 3 wq (mem_map_of_mem h)
  have l1 : 3 ≤ finrank ℚ (K1 shift 3 wq) :=
    card_le_finrank (fun i : Fin 3 => ⟨i, by omega⟩) (fun a b h => by simpa [Fin.ext_iff] using h)
      _ (fun i => e_mem_K1 i i.isLt)
  have u1 : finrank ℚ (K1 shift 3 wq) ≤ 3 := by
    have h := finrank_range_le_card (R := ℚ) (fun j : Fin 3 => (⇑shift)^[j] (wq 1))
    rw [Fintype.card_fin] at h
    exact h
  have l2 : 4 ≤ finrank ℚ (S2 shift 3 wq) := by
    refine card_le_finrank ![0, 1, 2, 4] (by decide) _ ?_
    intro i
    fin_cases i
    · exact k1 0 (by omega)
    · exact k1 1 (by omega)
    · exact k1 2 (by omega)
    · exact e4
  -- K₃ ∋ e₀ e₁ e₂ e₄ e₇ and A of e₂ e₄ e₇: e₃ e₅ e₈
  have l3 : 8 ≤ finrank ℚ (P3n shift 3 wq) := by
    refine card_le_finrank ![0, 1, 2, 3, 4, 5, 7, 8] (by decide) _ ?_
    intro i
    fin_cases i
    · exact hS (inS (k1 0 (by omega)))
    · exact hS (inS (k1 1 (by omega)))
    · exact hS (inS (k1 2 (by omega)))
    · change e 3 ∈ _
      rw [← shift_e 2 (by omega)]; exact hA (inS (k1 2 (by omega)))
    · exact hS (inS e4)
    · change e 5 ∈ _
      rw [← shift_e 4 (by omega)]; exact hA (inS e4)
    · exact hS e7
    · change e 8 ∈ _
      rw [← shift_e 7 (by omega)]; exact hA e7
  have u2 := finrank_P2 shift 3 wq
  have u3 := finrank_P3n shift 3 wq
  omega

end Rowmap.TightnessRos34pw2
