# Test problems — exact definitions

These definitions are binding for all three implementations: Fortran
(`fortran/experiments/problems.f90`), Python (`tools/problems.py`) and Rust
(`tests/problems/`). They follow the paper (§3, §4). Deviations and repairs are
marked **[fix]**. Agreed with Helmut on 2026-10-31.

## Common conventions

- All problems are semi-discretized with second-order central differences.
  `Δ` denotes the Laplacian; `D_x`, `D_y` denote central first differences
  `(u_{i+1} − u_{i−1}) / (2h)`.
- **Dirichlet grids:** `m × m` inner points, `h = L/(m+1)`, `x_i = x_0 + i·h`,
  `i = 1..m`. For all Dirichlet problems here, the boundary values of the
  exact or initial solution are **zero**, so the boundary contributions vanish.
- **Neumann grids** (vertex-centered + ghost mirror): `m` points per direction
  including the boundary, `h = 1/(m−1)`, `x_i = (i−1)·h`, `i = 1..m`.
  `∂u/∂n = 0` by mirroring: `u_0 := u_2`, `u_{m+1} := u_{m−1}`.
- **Ordering:** x index fastest. 2D: `k = i + (j−1)·m`. 3D:
  `k = i + (j−1)·m + (l−1)·m²`. Systems with two components are stored as
  blocks: `y = [u¹(all points), u²(all points)]`.
- `x` belongs to index `i`, `y` to `j`, `z` to `l`.
- ERR (paper): `sqrt(1/n Σ ((u_i − u_ref,i)/(1 + |u_ref,i|))²)`, taken over the
  n ODE components only, i.e. without an added time component (sft=3..5).

## Autonomous problems (paper §3)

| Name | Ω, grid | paper m → n | t_end |
|---|---|---|---|
| BRUSS  | [0,1]², Neumann, 2 comp.  | 50 → 5000 | 1 |
| HEAT   | [0,1]², Dirichlet          | 50 → 2500 | 10.24 |
| NILIDI | [0,π/3]², Dirichlet        | 69 → 4761 | 1 |
| FOOD   | [0,1]², Neumann, 2 comp.  | 50 → 5000 | 3 |
| FOOD3  | [0,1]³, Neumann, 2 comp.  | 10 → 2000 | 10 |
| VDPL   | [0,1]², Neumann, 2 comp.  | 20 → 800  | 1000 |

**BRUSS** — `α = 0.02`
```
u_t = 1 + u²v − 4u + αΔu
v_t = 3u − u²v + αΔv              [fix] paper prints α(u_xx+u_yy) in the v-equation
u(0) = 0.5 + y,  v(0) = 1 + 5x
```

**HEAT**
```
u_t = Δu,  u(0) = 16x(1−x)y(1−y),  u = 0 on ∂Ω      [fix] boundary value missing in the paper
```

**NILIDI** (= `fortran/driver.F`, but with exact π everywhere)
```
u_t = e^u Δu + u(18e^u − 1),  u(0) = sin(3x)sin(3y),  exact PDE solution e^{−t} sin(3x) sin(3y)
J·v = e^u (Δv + (Δu) v) + (18 e^u (1+u) − 1) v
```

**FOOD** — `d₁ = 0.05, d₂ = 1, b₁ = 1, b₂ = −1000, a₁ = 0.1, a₂ = −100`
```
u¹_t = d₁Δu¹ + u¹(b₁ − a₁u²)
u²_t = d₂Δu² + u²(b₂ − a₂u¹)
u¹(0) = 10 − 5cos(πx)cos(10πy),  u²(0) = 17 + 5cos(10πx)cos(πy)
```

**FOOD3** — **PARKED (2026-10-31).** It stays in the harness and data, but
it is not used in Rust tests or comparisons until the true coefficients are
known. The printed coefficients look wrong (see `docs/fortran-paper-runs.md`).
The source [3] (Brown/Hindmarsh, SIAM J. Numer. Anal. 23, 1986) could not be
obtained. The 1989 follow-up (Appl. Math. Comp. 31) contains a different
multi-species food web: a_ii = −1, a_ij = −½·10⁻⁶ (prey–predator),
10⁴ (predator–prey), b = ±(1 + αxy), 2-D.

Exactly as printed: `d₁ = 0.05, d₂ = 1, b₁ = b₂ = b(x,y,z) = (1 + 0.2xyz)(10⁶ − 1 + 10⁻⁶)`,
`a₁₁ = 10⁶, a₁₂ = 1, a₂₁ = 10⁶ − 1, a₂₂ = 10⁶`
```
u¹_t = d₁Δu¹ + u¹(b − a₁₁u¹ − a₁₂u²)
u²_t = d₂Δu² + u²(b − a₂₁u¹ − a₂₂u²)
u¹(0) = 500 + 250cos(πx)cos(3πy)cos(10πz),  u²(0) = 200 + 150cos(10πx)cos(πy)cos(3πz)
```
Very stiff initial transient (`f¹ ≈ −2.5·10¹¹` at t=0; equilibrium near
`u¹ ≈ 1, u² ≈ 10⁻⁶`). The reference therefore uses **Radau**, not DOPRI5.

**VDPL** — `d₁ = 0.05, d₂ = 1`
```
u¹_t = d₁Δu¹ + u²
u²_t = d₂Δu² + 100(1 − (u¹)²)u² − u¹
u¹(0) = 1 − cos(πx)cos(2πy),  u²(0) = 1 + 2cos(2πx)cos(πy)
```

## Non-autonomous problems (paper §4)

Both are on [0,1]², Dirichlet, paper `m = 69 → n = 4761`, `t_end = 1`. Let
`s = sin(πx)sin(πy)`, `cx = cos(πx)sin(πy)`, `cy = sin(πx)cos(πy)`.
The source term is taken from the **continuous** exact solution at the grid
points. The semi-discrete solution therefore differs from it by O(h²), and the
reference is DOPRI5 as for all other problems.

**DIFKON** — `ε = 0.1`
```
u_t = εΔu − y D_x u + x D_y u + g(t,x,y)
u_exact = s (1 + 4xy E),  E = exp(−2επ² t)
g = g₀ + E g₁
g₀ = 2επ² s + π (y cx − x cy)
g₁ = −8επ (y cx + x cy) + 4π x y (y cx − x cy) + 4 s (y² − x²)
∂f/∂t = g_t = −2επ² E g₁
```

**DIFFU2**
```
u_t = Δu + q(t,x,y)
u_exact = s (1 + 4xy sin t)
q   = 4xy s cos t + 2π² s (1 + 4xy sin t) − 8π sin t (y cx + x cy)
q_t = −4xy s sin t + 8π² xy s cos t − 8π cos t (y cx + x cy)
```
(`g` and `q` are checked symbolically by `tools/check_sources.py`.)

## Dropped: ROTADV (extra problem, 2026-10-31)

Rotating advection `u_t = −y D_x u + x D_y u (+ ν Δu)` on [−1,1]² was meant to
show GRK4T's imaginary-axis weakness against LSTAB. It was tried twice: a
smooth hump with ν = 1e-3, then the hard variant (ν = 0, a discontinuous disk,
four revolutions). LSTAB never won in work-precision; see
`docs/fortran-paper-runs.md`. It was removed from the harness, as agreed.

## Reference solutions

`tools/reference.py` computes each one once and stores it in
`tests/data/reference/<problem>_m<m>.txt`. Solvers: DOPRI5 (SciPy `RK45`),
`rtol = atol = 1e-10`; FOOD3 uses Radau with an analytic sparse Jacobian. The
Python discretization is cross-checked against the Fortran one
(`f(t0,y0)`, `J·v`, `f_t`) to ~1e-13 relative before the references are
trusted.
