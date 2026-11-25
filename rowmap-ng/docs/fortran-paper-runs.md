# Paper experiments with the original Fortran ROWMAP

Phase 1 of the roadmap (AGENTS.md). All runs used `fortran/rowmap.F`
(unchanged) or, for sft = 4/5, the patched copy
`fortran/experiments/rowmap_sft.F`.

Settings, as in the paper: GRK4T (`iwork(2)=2`), `rtol = atol = tol`,
tol = 1e-2 … 1e-7, max Krylov dim 66. Exact J·v (`ijacv=1`) and exact f_t
(`ifdt=1`). The paper used finite-difference J·v, so f-evaluation counts are
not directly comparable. Step limit 10⁶.

- **Error:** the paper's ERR at t_end against the DOPRI5 reference (1e-10; FOOD3:
  Radau), see docs/problems.md.
- **CPU time:** mean over repeated runs (at least 0.5 s accumulated).
- **Reproduce:** `fortran/run_experiments.sh`, then `uv run plots/work_precision.py`.
- **Tables:** every number is in `docs/fortran-paper-runs.tables.md`
  (generated); the raw data is in `tests/data/fortran/*.toml`.

All 168 runs (84 paper grid, 84 coarse grid; 8 problems) returned `idid = 1`.
The dropped ROTADV runs also all succeeded.

## Reference accuracy (checked 2026-10-31)

Coarse-grid references were compared with DOP853 at 1e-13: ERR is 1e-11 to
7e-11 for all problems except FOOD (2.6e-8). FOOD3 Radau
against Radau/BDF at 1e-12 gives about 3e-13. All of this is far below the
errors measured below.

## Comparison with the paper's figures

Only the *shapes* and ERR ranges can be compared: absolute times are about
1000× smaller than on the 1996 RS6000, and no VODPK, DOP853 or ROS4 curves
exist here.

| Problem | ERR range here (tol 1e-2 → 1e-7) | Paper | Remark |
|---|---|---|---|
| BRUSS (Fig. 3.1) | 10^-3.2 → 10^-7.0 | ≈ same range | smooth, regular |
| HEAT (Fig. 3.2) | 10^-5.3 → 10^-10.3 | down to ≈ 10^-11 | error far below tol (linear, decaying) |
| NILIDI (Fig. 3.3) | 10^-3.8 → 10^-8.1 | ≈ 10^-2 … 10^-8 | regular |
| FOOD (Fig. 3.4) | 10^-0.7 → 10^-4.9 | ≈ 10^0 … 10^-6 | ERR ≈ 20–100 × tol: global error growth (oscillatory dynamics, t_end = 3) |
| FOOD3 (Fig. 3.5) | 10^-5.3 → 10^-8.8 | ≈ 10^-2 … 10^-12 | **range differs from the paper** (see below) |
| VDPL (Fig. 3.6) | 10^-0.0 → 10^-7.2 | — | tol = 1e-2 useless (ERR ≈ 0.9), 1e-3 poor. Matches the paper's "reliable only for tol ≤ 1e-3" |
| DIFKON (Fig. 4.7) | sft=2 best; sft=3 worst by up to 2 orders; sft=4/5 in between | same ranking | **paper result reproduced** |
| DIFFU2 (Fig. 4.8) | same ranking; sft=3 worse by up to 3 orders | same ranking | **paper result reproduced** |

### Points to look at

1. **FOOD3.** With the coefficients exactly as printed, the solution
   collapses to the local reaction equilibrium within ~1e-9 time units
   (u¹ → 1 + 0.2xyz, u² → 1e-6). The step counts are *identical* on the
   5³, 7³ and 10³ grids, and the ERR range is much narrower than in Fig. 3.5.
   That suggests the paper's FOOD3 coefficients differ from what is printed
   (decided: literal; open question for Helmut / Brown–Hindmarsh [3]).
   **Parked on 2026-10-31.** [3] could not be obtained online (SIAM paywall,
   OSTI unreachable). The 1989 follow-up paper has a different food-web model.
2. **FOOD.** The global error stays 20–100× above the tolerance at all tolerances.
   This is a property of the problem (error propagation), not a failure.
3. **ROTADV** (extra problem, GRK4T vs LSTAB) — **dropped**.
   - First try (smooth hump, ν = 1e-3, one revolution): LSTAB was worse by
     ~2× in ERR and ~1.7× in steps.
   - Hard variant (ν = 0, discontinuous disk, four revolutions, m = 50): at
     tol 1e-2 GRK4T needed **325 steps against 67** for LSTAB, and 409 against
     343 at 1e-3. That is GRK4T's stability limit at work: the step control
     holds its steps small.
   - Even so, GRK4T was still faster in CPU time (0.105 s vs 0.128 s) at the
     same ERR (≈ 0.1), because LSTAB's larger γ costs more Krylov work per
     step. At tol ≤ 1e-4, accuracy limits the step size and GRK4T is better by
     ~2× in ERR and needs ~1.8× fewer steps.
   - Conclusion: no work-precision advantage for LSTAB on this problem. As
     agreed (one try), the problem was removed, and the LSTAB option was
     dropped from the design: GRK4T only.

### Fortran findings (added to AGENTS.md §4)

- `driver.F` uses an inconsistent π (f vs J·v).
- The default step limit of 10000 is too small for VDPL at tight tolerances.
- No failures and no messages from ROWMAP (e.g. "Matrix is singular") in any
  of the runs (driver stdout is kept as `*.log` only when non-empty; none
  exist). Rejected steps are not recorded in the output files. The ad-hoc thresholds of AGENTS.md §4 were not visibly triggered by
  these problems, so the breakdown tests (phase 3) need targeted problems.
