# ROWMAP in Rust vs the original Fortran code: paper experiments

Roadmap phase 6 (2026-10-31). This document summarizes the results. The
details are in:

- `docs/fortran-paper-runs.md`: the Fortran baseline, and how it compares with
  the paper's figures;
- `docs/rust-runs.md`: Rust runs, the schedule comparison, the ε-sweep through
  rowmap.F, and the non-autonomous problems;
- `docs/tables/`: generated tables (`summary.md`, `rust-vs-f77.md`,
  `reduced-vs-paper.md`);
- `plots/out/wp_<problem>_m<m>.png`: work-precision diagrams in the style of
  Figs. 3.1–3.6 and 4.7–4.8, overlaying Fortran (and its sft variants) with
  Rust.

## Setup

Both codes use the same settings as in the paper:
- 7 problems (FOOD3 parked, see `docs/problems.md`), on paper grids and
  coarse grids;
- tol = 1e-2 … 1e-7 with `rtol = atol = tol`, GRK4T, max Krylov dimension 66;
- exact J·v and exact f_t;
- error: the paper's ERR against one-off DOPRI5 references at 1e-10;
- CPU: mean over repeated runs (≥ 0.5 s in total).

The Rust code uses the paper's MAP schedules: autonomous 4/7/8/11 and the
derived non-autonomous 4/9/10/15. For non-autonomous problems the Fortran
comparison is sft = 2, the native form.

## Summary (`docs/tables/summary.md`)

Ratios Rust / Fortran, geometric mean over the six tolerances. State after
the phase-7 performance work; before it, the CPU ratio was 3.4–8.9:

| problem | m | sft | ERR ratio (geo. mean) | ERR ratio range | nstep ratio | CPU ratio |
|---|---|---|---|---|---|---|
| bruss | 50 | 0 | 0.26 | 0.05–0.98 | 0.99 | 2.2 |
| bruss | 10 | 0 | 0.21 | 0.08–0.38 | 0.96 | 2.5 |
| heat | 50 | 0 | 0.85 | 0.70–0.96 | 1.00 | 2.2 |
| heat | 10 | 0 | 1.08 | 1.00–1.21 | 1.03 | 2.7 |
| nilidi | 69 | 0 | 0.72 | 0.19–1.75 | 0.75 | 2.2 |
| nilidi | 10 | 0 | 1.44 | 0.70–6.44 | 0.87 | 2.6 |
| food | 50 | 0 | 0.62 | 0.33–0.81 | 1.00 | 2.5 |
| food | 10 | 0 | 0.24 | 0.04–1.59 | 0.99 | 2.7 |
| vdpl | 20 | 0 | 0.73 | 0.10–1.81 | 1.01 | 2.9 |
| vdpl | 8 | 0 | 0.85 | 0.29–1.48 | 0.98 | 2.1 |
| difkon | 69 | 2 | 0.81 | 0.45–1.20 | 0.84 | 2.4 |
| difkon | 10 | 2 | 0.53 | 0.37–0.73 | 0.84 | 2.2 |
| diffu2 | 69 | 2 | 0.59 | 0.29–1.37 | 0.84 | 2.6 |
| diffu2 | 10 | 2 | 1.30 | 0.65–6.16 | 0.85 | 2.6 |

## Conclusions

1. **Same method, same step-size control.** Step counts agree within a few
   percent. The Rust code takes a few fewer steps on NILIDI, DIFKON and
   DIFFU2, all of which need only a handful of steps, because of the
   Hairer–Wanner initial step.
2. **Accuracy equal or better.** The geometric-mean ERR ratio is < 1 in 11 of
   14 groups. BRUSS and FOOD are much better at tight tolerances (up to 25×).
   The few worse single runs are within a factor of about 6, at crude
   tolerances or on very few steps.
   The most likely reason for the better accuracy is the more careful Krylov
   basis: two Gram–Schmidt passes everywhere, the exact stage-1 residual, and
   no ad-hoc thresholds. That is plausible but not separately verified.
3. **Non-autonomous problems.** The derived schedule (`docs/map.md` §7)
   matches or beats Fortran sft = 2. The paper's finding (sft = 2 is better
   than the augmented system, sft = 3/4/5) is reproduced in the Fortran runs,
   and the Rust curve lies on the sft = 2 shape.
4. **CPU: Rust is 2.1–2.9× slower** (before phase 7: 3.4–8.9×). Phase 7 made
   two changes: a selective instead of unconditional re-orthogonalization
   before promotion, and reductions with independent accumulators so the
   compiler can use SIMD, plus a fast path for `nrm2`. The remaining gap
   comes mostly from the deliberate robustness choice of at least two
   Gram–Schmidt passes, against Fortran's single MGS pass. Allocations per
   stage and missing BLAS are minor.
5. **Robustness.** All runs of both codes succeed. The ε-sweep near an exact
   breakdown is continuous for both codes (`docs/rust-runs.md` §3).

## Not reproduced

- The VODPK, DOP853 and ROS4 curves of the paper: those codes are not
  available.
- FOOD3: the printed coefficients are suspicious and source [3] is
  unavailable (parked).
- Absolute times: the 1996 RS6000 timings are about 1000× larger.

## Reproduce

```sh
fortran/run_experiments.sh all                                # Fortran (~2 min)
cargo run --release --example experiments -- all              # Rust (~1.5 min)
cargo run --release --example experiments -- all --reduced    # schedule comparison
uv run tools/compare_runs.py --summary --markdown docs/tables/summary.md
uv run tools/compare_runs.py --markdown docs/tables/rust-vs-f77.md
uv run plots/work_precision.py                                # plots/out/
```
