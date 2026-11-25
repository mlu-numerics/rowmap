# Data formats

## Run results (TOML) — every producer (Fortran driver, Rust examples)

One file per run. Fortran: `tests/data/fortran/<problem>_m<m>_<method>_sft<s>_tol<k>.toml`;
Rust: `results/<same name>.toml`.

```toml
[run]
producer = "rowmap-f77"      # or "rowmap-rs"
problem = "nilidi"           # see docs/problems.md
m = 69                       # grid parameter
n = 4761                     # ODE dimension (without an added time component)
method = "grk4t"             # "grk4t" or "ros34pw2" (Rust only)
schedule = "paper"           # Rust only: "paper", "no-stage4-closure" or "single-pass" (comparison runs)
rtol = 1.0e-4
atol = 1.0e-4
max_krylov = 66
sft = 0                      # 0 autonomous; 2..5 non-autonomous strategies (paper §4)
t_end = 1.0
idid = 1                     # Fortran idid; Rust: 1 = ok, <0 = error

[stats]
err = 3.2e-5                 # endpoint ERR (paper formula) vs reference; nan if failed / no reference
nstep = 7                    # steps, accepted + rejected (Fortran iwork(5))
nfeval = 22                  # f evaluations (Fortran iwork(7))
njv = 120                    # diagnostic: J·v products (Fortran iwork(8)); optional
cpu_time = 0.045             # seconds, integration only
```

Rules:
- Floats use full precision (`es24.16e3` in Fortran, round-trip in Rust).
  TOML allows leading zeros in the exponent.
- No solution vectors in run files.
- A failed run has `idid < 0` and `err = nan`. Plots omit it, like the missing
  marks in the paper.

## Reference solutions (plain text, input data)

`tests/data/reference/<problem>_m<m>.txt`: `#` header lines (problem, m, n,
t0, t_end, solver, tolerance, nfev, cpu, date), then one `y_end[i]` per line.
Fortran (list-directed read), Rust and numpy can all read it. Computed once by
`tools/reference.py` (DOPRI5 = SciPy RK45, 1e-10; FOOD3: Radau) and committed.

## Discretization dump (`driver dump`, for cross-checks only)

`#` header, then per line: `y0  f(t0,y0)  J(t0,y0)·v  f_t(t0,y0)` with
`v_k = sin(k)`. Used by `tools/crosscheck.py`, not committed.
