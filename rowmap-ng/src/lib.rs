//! ROWMAP: a Krylov-W method of order 4(3) for large stiff ODEs
//! `y' = f(t, y)`.
//!
//! The method is GRK4T (Kaps–Rentrop). Its linear stage systems are solved
//! in Krylov subspaces built by a *multiple Arnoldi process* (MAP). That
//! process keeps the order 4 already for small, n-independent Krylov
//! dimensions: 4/7/8/11 for autonomous problems, 4/9/10/15 for
//! non-autonomous ones. The Jacobian is only needed as a product `J·v`, which
//! the user supplies exactly.
//!
//! Based on R. Weiner, B.A. Schmitt, H. Podhaisky: *ROWMAP – a ROW-code with
//! Krylov techniques for large stiff ODEs* (1996). Design, proofs and
//! experiments are in the [repository]: `AGENTS.md`, `docs/map.md`,
//! `docs/comparison.md`.
//!
//! [repository]: https://github.com/mlu-numerics/rowmap/tree/main/rowmap-ng
//!
//! # Example
//!
//! The 1D heat equation `u_t = u_xx` on (0, 1), zero boundary values,
//! semi-discretized with central differences:
//!
//! ```
//! use rowmap::{OdeSystem, Rowmap};
//!
//! struct Heat {
//!     n: usize,
//! }
//!
//! impl Heat {
//!     fn lap(&self, u: &[f64], out: &mut [f64]) {
//!         let h2 = 1.0 / ((self.n + 1) as f64).powi(2);
//!         for i in 0..self.n {
//!             let l = if i > 0 { u[i - 1] } else { 0.0 };
//!             let r = if i + 1 < self.n { u[i + 1] } else { 0.0 };
//!             out[i] = (l - 2.0 * u[i] + r) / h2;
//!         }
//!     }
//! }
//!
//! impl OdeSystem for Heat {
//!     fn dim(&self) -> usize {
//!         self.n
//!     }
//!     fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
//!         self.lap(y, dy);
//!     }
//!     // exact Jacobian-vector product (the problem is linear: J = Laplacian)
//!     fn jac_vec(&self, _t: f64, _y: &[f64], v: &[f64], out: &mut [f64]) {
//!         self.lap(v, out);
//!     }
//!     fn autonomous(&self) -> bool {
//!         true
//!     }
//! }
//!
//! let sys = Heat { n: 200 };
//! let pi = std::f64::consts::PI;
//! let x = |i: usize| (i + 1) as f64 / 201.0;
//! let mut y: Vec<f64> = (0..200).map(|i| (pi * x(i)).sin()).collect();
//!
//! let sol = Rowmap::new(1e-6, 1e-6).solve(&sys, 0.0, &mut y, 0.1)?;
//!
//! // exact solution of the semi-discrete problem: sin(πx) e^{λ t}
//! let lam = -4.0 * 201.0f64.powi(2) * (pi / (2.0 * 201.0)).sin().powi(2);
//! let err = (0..200).map(|i| (y[i] - (pi * x(i)).sin() * (lam * 0.1).exp()).abs()).fold(0.0, f64::max);
//! assert!(err < 1e-5);
//! println!("{} steps, {} f-evaluations, {} J·v", sol.stats.steps, sol.stats.f_evals, sol.stats.jv);
//! # Ok::<(), rowmap::RowmapError>(())
//! ```
//!
//! Non-autonomous problems (the default of [`OdeSystem::autonomous`]) must
//! also implement [`OdeSystem::dfdt`]. Dense output and early termination
//! are available through [`Rowmap::solve_with`].

// Index loops that mirror the formulas (Σ_j a_ij x_j) are clearer in numerical
// kernels than iterator chains.
#![allow(clippy::needless_range_loop)]

mod dense;
mod givens;
mod grk4t;
mod krylov;
mod linalg;
mod map;
#[doc(hidden)]
pub mod prof;
mod solver;
mod system;

pub use solver::{DimStats, Method, Rowmap, RowmapError, Solution, Stats, Status, Step, Tol};
pub use system::{OdeSystem, Unsupported};
