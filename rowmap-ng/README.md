# rowmap — ROWMAP in Rust

[![crates.io](https://img.shields.io/crates/v/rowmap.svg)](https://crates.io/crates/rowmap)
[![docs.rs](https://docs.rs/rowmap/badge.svg)](https://docs.rs/rowmap)

A from-scratch Rust implementation of **ROWMAP**, a Krylov-W method of order
4(3) for large stiff ODEs `y' = f(t, y)` (Weiner, Schmitt, Podhaisky, 1996).
It is the successor of the original Fortran 77 code in the
[parent directory](https://github.com/mlu-numerics/rowmap).

- **GRK4T** stages, with linear systems solved in Krylov subspaces built by
  the **multiple Arnoldi process** (MAP). Order 4 holds for n-independent
  Krylov dimensions 4/7/8/11 (autonomous) and 4/9/10/15 (non-autonomous).
- No Jacobian matrix: only an exact Jacobian-vector product `J·v` is needed,
  plus `∂f/∂t` for non-autonomous problems.
- Optional order-3 method **ROS34PW2** (L-stable, stiffly accurate), often
  cheaper at crude tolerances.
- `f64`, no dependencies. Robustness first: two Gram–Schmidt passes, one
  relative dependency test for (near) Krylov breakdown, and the exact stage
  residual.

## Installation

```sh
cargo add rowmap
```

or in `Cargo.toml`:

```toml
[dependencies]
rowmap = "0.1"
```

Rust 1.85 or newer (edition 2024).

## Usage

You describe the problem by implementing the trait `OdeSystem`, then call
`Rowmap::solve`:

| Method | Meaning |
|---|---|
| `dim()` | dimension n |
| `rhs(t, y, dy)` | `dy = f(t, y)` |
| `jac_vec(t, y, v, out)` | `out = ∂f/∂y (t, y) · v`. It must be **exact** (no difference quotient), because the order 4 depends on it. |
| `autonomous()` | `true` if `f` does not depend on `t` (default `false`). This allows smaller Krylov spaces. |
| `dfdt(t, y, ft)` | `ft = ∂f/∂t (t, y)`. Required for non-autonomous problems. |

### Example: a non-autonomous stiff problem

The Prothero–Robinson problem `y_i' = -λ_i (y_i - sin t) + cos t` has the
exact solution `y = sin t`.

```rust
use rowmap::{OdeSystem, Rowmap, RowmapError, Unsupported};

struct ProtheroRobinson {
    lambda: Vec<f64>,
}

impl OdeSystem for ProtheroRobinson {
    fn dim(&self) -> usize {
        self.lambda.len()
    }
    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) {
        for i in 0..y.len() {
            dy[i] = -self.lambda[i] * (y[i] - t.sin()) + t.cos();
        }
    }
    fn jac_vec(&self, _t: f64, _y: &[f64], v: &[f64], out: &mut [f64]) {
        for i in 0..v.len() {
            out[i] = -self.lambda[i] * v[i];
        }
    }
    // f depends on t, so autonomous() keeps its default `false` and
    // dfdt is required
    fn dfdt(&self, t: f64, _y: &[f64], ft: &mut [f64]) -> Result<(), Unsupported> {
        for i in 0..ft.len() {
            ft[i] = self.lambda[i] * t.cos() - t.sin();
        }
        Ok(())
    }
}

fn main() -> Result<(), RowmapError> {
    let sys = ProtheroRobinson { lambda: (1..=1000).map(|k| k as f64).collect() };
    let mut y = vec![0.0; 1000]; // initial value, overwritten with y(t_end)

    let sol = Rowmap::new(1e-6, 1e-6) // rtol, atol
        .solve(&sys, 0.0, &mut y, 1.0)?;

    println!("y(1) = {}, exact {}", y[999], 1f64.sin());
    println!("{} steps, {} f-evals, {} J·v", sol.stats.steps, sol.stats.f_evals, sol.stats.jv);
    Ok(())
}
```

The crate documentation contains a second, autonomous example (the 1D heat
equation).

### Options

`Rowmap::new(rtol, atol)` creates the integrator. Tolerances are scalars
(`f64`) or one per component (`Vec<f64>`, `&[f64]`). The error is measured in
the weighted RMS norm with weights `atol_i + rtol_i |y_i|`. The other
settings have defaults that fit most problems:

| Setter | Default | Meaning |
|---|---|---|
| `.method(Method::Grk4t)` | `Grk4t` | `Grk4t` (order 4) or `Ros34pw2` (order 3) |
| `.max_krylov(m)` | 70 | maximal Krylov dimension (memory n × m) |
| `.max_steps(n)` | 100 000 | maximal number of steps (accepted + rejected) |
| `.initial_step(h)` | estimated | initial step size |
| `.max_step(h)` | unlimited | maximal step size |

Only forward integration is supported (`t_end >= t0`).

### Dense output and early stop

`solve_with` calls a closure after every accepted step. `Step::interpolate`
evaluates a cubic Hermite interpolant on the step `[t_old, t_new]`. If the
closure returns `ControlFlow::Break(())`, the integration stops
(`Solution::status == Status::Interrupted`).

```rust
use std::ops::ControlFlow;

let mut grid = (1..=10).map(|k| k as f64 / 10.0).peekable(); // output times
let mut out = vec![0.0; sys.dim()];
Rowmap::new(1e-6, 1e-6).solve_with(&sys, 0.0, &mut y, 1.0, |step| {
    while let Some(&s) = grid.peek() {
        if s > step.t_new {
            break;
        }
        step.interpolate(s, &mut out);
        println!("t = {s:.1}: y[0] = {}", out[0]);
        grid.next();
    }
    ControlFlow::Continue(())
})?;
```

### Result and errors

`Solution` holds the final time `t`, the predicted next step size `h_next`,
the `status` and `stats`. The statistics include steps (accepted/rejected),
f evaluations, J·v products and the Krylov dimensions per stage.

On failure, `RowmapError` says why (`InvalidInput`, `MissingDfdt`,
`StepSizeTooSmall`, `MaxStepsExceeded`, `SingularStageMatrix`, `NonFinite`).
`y` then holds the last accepted solution.

## How it works

A ROW (Rosenbrock-type) method needs, in every stage, the solution of
`(I - hγJ) k = r`. ROWMAP never forms `J`. It solves these systems in Krylov
subspaces and builds the subspaces of the later stages from those of the
earlier stages (*multiple Arnoldi process*). With this construction the
Krylov approximation does not lower the order of the method, already for
small Krylov dimensions that do not depend on n. The theory behind it is in
the original report:

> R. Weiner, B.A. Schmitt, H. Podhaisky: *ROWMAP – a ROW-code with Krylov
> techniques for large stiff ODEs*. Report 39, FB Mathematik und Informatik,
> Universität Halle, 1996
> ([PostScript](https://github.com/mlu-numerics/rowmap/blob/main/rowmap-ng/references/rpha9639.ps)).

More detail in this repository:

- [`docs/map.md`](https://github.com/mlu-numerics/rowmap/blob/main/rowmap-ng/docs/map.md):
  construction of the MAP, minimal Krylov dimensions, ROS34PW2
- [`docs/comparison.md`](https://github.com/mlu-numerics/rowmap/blob/main/rowmap-ng/docs/comparison.md):
  Rust compared with the Fortran code on all paper problems
- [`docs/problems.md`](https://github.com/mlu-numerics/rowmap/blob/main/rowmap-ng/docs/problems.md):
  definitions of the test problems
- [`AGENTS.md`](https://github.com/mlu-numerics/rowmap/blob/main/rowmap-ng/AGENTS.md):
  design decisions and development log

## Verification

- **Order:** the deviation from the exact ROW method behaves like h⁵ at the
  minimal Krylov dimensions. Tests omit single MAP actions and check that
  the order drops as predicted (`src/map/tests*.rs`, `docs/map.md`).
- **Lean 4 / Mathlib** (`lean/`): the structural facts of the MAP (subspace
  memberships, closure properties, dimension increments, and examples that
  attain them).
- **Paper experiments:** all paper problems were run with the original
  Fortran code and with this code (`docs/comparison.md`). Step counts are the
  same, and accuracy is equal or better. At equal accuracy, Rust GRK4T takes
  about 0.9–1.6× the Fortran CPU time.

Run the tests and the paper experiments from the `rowmap-ng/` directory:

```sh
cargo test --release
# rewrites results/*.toml; options: [paper|coarse|all] [--only=<problem>]
#   [--method=grk4t|ros34pw2] [--reduced] [--single-pass]
cargo run --release --example experiments -- paper
```

## Layout

| Path | Contents |
|---|---|
| `src/` | library (`solver.rs` API, `map.rs` MAP step, `krylov.rs` basis) |
| `tests/`, `examples/experiments.rs` | integration tests, paper experiments |
| `docs/` | MAP construction and proofs, problem definitions, comparisons, report |
| `lean/` | formal proofs |
| `fortran/` | original rowmap.F (unchanged) and experiment drivers |
| `results/`, `tools/`, `plots/` | run results, reference solutions, comparison tables, plots (uv scripts) |
| `AGENTS.md` | design decisions and roadmap |

The crate on crates.io contains only `src/`, this README and the license.

## License

MIT, see [LICENSE](LICENSE).
