# RowmapLean — structural proofs for the ROWMAP multiple Arnoldi process

Lean 4 + Mathlib formalization of the *structural* part of the order proof
(docs/map.md, AGENTS.md §3.8).

- `RowmapLean/Map.lean`:
  - the MAP schedule as data (`StageKind`, `grk4tAutonomous`) with a generic
    subspace semantics (`mapSpaces`);
  - the facts F1–F5 for GRK4T, autonomous: `w_i ∈ K_i`, the Krylov powers,
    `A(K₁+⟨w₂⟩) ⊆ K₂`, `A(K₃+⟨w₄⟩) ⊆ K₄`, `w₄ ∈ K₃`, and the dimension
    increments `+3, +1, +3`.
- `RowmapLean/Tightness.lean`: an example over ℚ (shift operator on ℚ¹¹)
  that attains `dim K = 4, 7, 8, 11`.
- `RowmapLean/HStructure.lean`: the generalized Hessenberg structure of
  `H = QᵀAQ` (docs/map.md §9). `H_structure` covers the MAP as a FIFO process
  of insert/promote operations (each one accepted or dropped). For mutually
  orthogonal vectors with the Gram–Schmidt flag property, it proves the zero
  pattern below the closing rows, the staircase (monotone closing rows,
  closed columns form a prefix), the row supports, and the open-column counts
  `o(m) = trace 0 ops`. `trace_*` and `bandwidths` give the exact sequences
  and the lower bandwidths 3/4 (GRK4T) and 1/3 (ROS34PW2) for every κ₁.
- `DumpSchedule.lean`: exports the schedule to `schedule.toml`. The Rust test
  `schedule_matches_lean` checks it against `src/map.rs`.

Not formalized: the analytic step "F1–F4 with κ₁ ≥ 4 ⇒ classical order 4"
(paper proof in docs/map.md §4, numerical tests in `src/map/tests.rs`).

```sh
lake exe cache get                               # once: Mathlib cache (~7 GB)
lake build                                       # must succeed without `sorry`
lake env lean --run DumpSchedule.lean schedule.toml
rg -n "sorry|admit|^axiom" RowmapLean            # must be empty
```
