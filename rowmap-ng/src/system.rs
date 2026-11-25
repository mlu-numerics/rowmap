//! The user-facing problem description (AGENTS.md §2).

/// Returned by the default [`OdeSystem::dfdt`]: the system does not provide
/// `∂f/∂t`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Unsupported;

/// An initial value problem `y' = f(t, y)` with an exact Jacobian-vector
/// product.
pub trait OdeSystem {
    /// Dimension n.
    fn dim(&self) -> usize;

    /// `dy = f(t, y)`.
    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]);

    /// `out = J(t, y) · v` with `J = ∂f/∂y`. Must be exact (up to rounding):
    /// the order 4 of the method relies on it.
    fn jac_vec(&self, t: f64, y: &[f64], v: &[f64], out: &mut [f64]);

    /// `true` if `f` does not depend on `t`. This enables the smaller Krylov
    /// increments. Declaring `true` for a t-dependent `f` gives wrong results.
    fn autonomous(&self) -> bool {
        false
    }

    /// `ft = ∂f/∂t (t, y)`. Required unless [`autonomous`](Self::autonomous)
    /// returns `true`. The default returns `Err(Unsupported)`, and the solver
    /// then fails up front.
    fn dfdt(&self, _t: f64, _y: &[f64], _ft: &mut [f64]) -> Result<(), Unsupported> {
        Err(Unsupported)
    }
}
