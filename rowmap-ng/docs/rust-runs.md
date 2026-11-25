# Rust integrator: paper experiments, schedule comparison, ε-sweep vs Fortran

Phase 4 (2026-10-31). Settings as for the Fortran runs: GRK4T,
`rtol = atol = tol`, tol = 1e-2 … 1e-7, max Krylov dim 66, exact J·v. Error =
paper ERR against the DOPRI5 reference. Sections 1–3 cover the autonomous
problems (phase 4). Section 4 covers the non-autonomous ones (phase 5).

Reproduce:

```sh
cargo run --release --example experiments -- all             # results/*_grk4t_*.toml
cargo run --release --example experiments -- all --reduced   # results/*_grk4t_reduced_*.toml
uv run tools/compare_runs.py --markdown docs/tables/rust-vs-f77.md
uv run tools/compare_runs.py --reduced-vs-paper --markdown docs/tables/reduced-vs-paper.md
uv run plots/work_precision.py
```

## 1. Rust vs Fortran (`docs/tables/rust-vs-f77.md`)

- **Step-size control:** `nstep` and `nfeval` agree within a few percent on
  all 60 runs. The only differences come from the initial step (Hairer–Wanner
  here, a heuristic in rowmap.F), visible on NILIDI with its 3–7 steps.
- **Accuracy:** ERR is equal or better in most runs.
  - BRUSS and FOOD at tight tolerances: Rust is 3–25× more accurate (e.g.
    BRUSS m=50, tol 1e-7: 5.0e-9 vs 9.5e-8). The likely cause is the more
    careful Krylov basis (two Gram–Schmidt passes, exact residual), but that
    is not verified.
  - Exceptions: NILIDI m=69 at 1e-4/1e-5 (1.6–1.8× worse), and VDPL m=20 at
    1e-5 (9× worse). In the VDPL case the Fortran value is the outlier: it
    drops from 1.9e-4 to 7.9e-7 between neighbouring tolerances.
- **CPU:** Rust was 5–10× slower in phase 4. After phase 7 it is 2.1–2.9×
  slower (`docs/comparison.md`). The phase-4 causes were:
  - about 4 orthogonalization passes per new vector (CGS2 in
    `push_unit`, plus the safety pass in `promote`), against Fortran's single
    MGS pass;
  - naive dot products without BLAS/SIMD.
- All 60 runs succeed. No rejections are caused by singular stage matrices.

## 2. Paper schedule vs reduced schedule (`docs/tables/reduced-vs-paper.md`)

The reduced schedule omits the stage-4 closure (dims ≤ κ₁+3, +1, +0;
docs/map.md §4). It is not needed for the *classical* order.

- **Step counts:** identical (±2% on VDPL). The error estimate does not react.
- **ERR:** mostly within a factor of 0.8–1.6 either way.
  - Worse with the reduced schedule: FOOD m=50 at tol 1e-2 (3.6×), BRUSS m=50
    at 1e-4 (2.1×), HEAT m=50 at 1e-4 (2.2×).
  - Better with the reduced schedule: FOOD m=50 at 1e-4…1e-6 (0.35–0.44×),
    VDPL m=20 at 1e-2 (0.2×).
- **CPU:** 10–50% lower (BRUSS about −45%).
- No failures and no instability.

Interpretation: the stage-4 closure buys a little accuracy in the
pre-asymptotic regime (crude tolerances, stiff components), but it is not
essential on these problems. **Decision open (Helmut).** The paper schedule
stays the default.

## 3. ε-sweep through rowmap.F (`plots/out/sweep_f77.png`)

Problem `y' = (diag(−1,−10,−100) + εG) y` with n = 30, one step h = 0.1, and
ε = 1e-16 … 1e-2. The Fortran run takes one step with `rtol = atol = 1`; its
later stages stop early once their residual is small (which the Rust code
never does), so it uses 6 J·v products.

- **rowmap.F is continuous too**: `‖u_F(ε) − u_F(0)‖` follows the exact ROW
  change, and the deviation from exact ROW grows like ε. It is the same size
  as with Rust (slightly smaller around ε ≈ 1e-11).
- The ad-hoc thresholds (AGENTS.md §4, items 1–2) cause **no visible damage
  on this test**. That is a null result for the design's expectation. They may
  still matter for problems where a remainder falls *between* 1e-14 and the
  rounding level relative to a large `‖A‖`. Such a case was not constructed.

## 4. Non-autonomous problems: DIFKON, DIFFU2 (phase 5)

The Rust code uses the native non-autonomous form, comparable to Fortran
sft = 2, with the derived schedule `GRK4T_NONAUTONOMOUS` (4/9/10/15,
docs/map.md §7). f_t is evaluated once per accepted step. Tables are in
`docs/tables/rust-vs-f77.md` (rows `difkon`, `diffu2`). Plots:
`plots/out/wp_difkon_m69.png`, `wp_diffu2_m69.png`, with Fortran sft = 2–5
overlaid.

- **Steps:** slightly fewer than Fortran sft = 2 at crude tolerances (4 vs 6,
  from the initial step), converging at tight tolerances (e.g. DIFFU2 m=69,
  tol 1e-7: 94 vs 99).
- **ERR:** equal or better in 21 of 24 runs (DIFKON m=69 at 1e-7: 7.1e-9 vs
  1.6e-8; DIFFU2 m=69 at 1e-4/1e-5: about 3× better). Exceptions:
  - DIFFU2 m=10 at tol 1e-2: 1.1e-3 vs 1.8e-4. Still at tolerance level; the
    run has only 4 steps.
  - DIFFU2 m=10 at 1e-5/1e-6: 1.3–1.4× worse.
- **Work precision:** the same shape as Fortran sft = 2. It is clearly better
  than the augmented-system strategies sft = 3/4/5 in accuracy, and shifted
  only by the general CPU factor (7–10× in phase 5, 2–3× after phase 7).
- All 24 runs succeed.
