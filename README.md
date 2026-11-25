# Rowmap

A ROW-code with Krylov techniques for large stiff ODEs

This repository contains two implementations:

| Directory | Contents |
|---|---|
| `.` (this directory) | the original Fortran 77 code `rowmap.F` (1997) with a demo driver |
| [`rowmap-ng/`](rowmap-ng/) | **rowmap-ng**, a from-scratch Rust implementation, available as the crate [`rowmap`](https://crates.io/crates/rowmap) |

## Fortran code

For compilation, see makefile. 

Running the demo yields:

```sh 
make
./demo

❯ ./demo

Statistics: 
number of steps     =    7     number of rejected steps  =    0
number of f evals.  =  343     number of jac-vector op.  =  321
average krylov dim. =  45.86
minimum for lwork   = 390831   minimum for liwork        =   90
```

## rowmap-ng (Rust)

The Rust version implements the same method (GRK4T with the multiple Arnoldi
process), plus an optional order-3 method ROS34PW2. It needs only a
Jacobian-vector product, has no dependencies, and is verified against the
Fortran code on all problems of the original paper.

```sh
cargo add rowmap
```

```rust
use rowmap::{OdeSystem, Rowmap};

let sol = Rowmap::new(1e-6, 1e-6)       // rtol, atol (scalars or vectors)
    .solve(&sys, t0, &mut y, t_end)?;    // y: initial value -> solution
```

See [`rowmap-ng/README.md`](rowmap-ng/README.md) for the full documentation
and [docs.rs/rowmap](https://docs.rs/rowmap) for the API reference.

## Reference

R. Weiner, B.A. Schmitt, H. Podhaisky: *ROWMAP – a ROW-code with Krylov
techniques for large stiff ODEs*. Report 39, FB Mathematik und Informatik,
Universität Halle, 1996.
