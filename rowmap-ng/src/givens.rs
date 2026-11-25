//! Incremental Givens QR of the projected stage matrix `M = I − γhH`
//! (docs/map.md §9.4).
//!
//! All stages of a step share `M`; only the basis grows. By the structure
//! theorem (§9.2, Lean `H_structure`) a *closed* column of `H` is final and
//! zero below its closing row, and the closed columns form a prefix. So each
//! closed column is rotated to triangular form exactly once (`update`), and
//! its rotations are kept. A solve applies the stored rotations to the
//! right-hand side and to the few *open* trailing columns, triangularizes
//! that small trailing block, and back-substitutes.
//!
//! The code reads the closing rows that the basis records at run time; it
//! makes no assumption about the schedule (autonomous or non-autonomous,
//! GRK4T or ROS34PW2) and stays correct if the staircase is violated (test
//! hook): a closed column behind an open one is simply treated as open.

use crate::dense::{SING_TOL, Singular};
use crate::krylov::Basis;

/// Plane rotation acting on rows `(i, i+1)`: `[c s; −s c]`.
#[derive(Debug, Clone, Copy)]
struct Rot {
    i: usize,
    c: f64,
    s: f64,
}

impl Rot {
    #[inline]
    fn apply(&self, x: &mut [f64]) {
        let (a, b) = (x[self.i], x[self.i + 1]);
        x[self.i] = self.c * a + self.s * b;
        x[self.i + 1] = -self.s * a + self.c * b;
    }

    /// Rotation on rows `(i, i+1)` that zeroes `x[i+1]`.
    #[inline]
    fn zeroing(x: &mut [f64], i: usize) -> Rot {
        let (a, b) = (x[i], x[i + 1]);
        let r = a.hypot(b);
        let rot = if r == 0.0 {
            Rot { i, c: 1.0, s: 0.0 }
        } else {
            Rot {
                i,
                c: a / r,
                s: b / r,
            }
        };
        x[i] = r;
        x[i + 1] = 0.0;
        rot
    }
}

/// Incremental QR state for one step.
#[derive(Debug, Clone)]
pub(crate) struct Givens {
    gh: f64,
    /// Finalized columns `0 .. r.len()`: column j of R (rows `0..=j`).
    r: Vec<Vec<f64>>,
    /// Rotations of the finalized columns, in the order of application.
    rots: Vec<Rot>,
    /// Largest |M_ij| seen: reference for the singularity test (as the LU).
    mmax: f64,
}

/// The open trailing part after the stored rotations, triangularized.
struct Tail {
    /// Number of finalized columns.
    nf: usize,
    /// Open columns `nf..k` of R (length k each, upper triangular overall).
    cols: Vec<Vec<f64>>,
    /// Rotated right-hand side (length k).
    b: Vec<f64>,
}

impl Givens {
    pub fn new(gh: f64) -> Self {
        Givens {
            gh,
            r: Vec::new(),
            rots: Vec::new(),
            mmax: 0.0,
        }
    }

    /// Column j of `M = I − γhH`, rows `0..k` (H entries beyond the stored
    /// column length are zero).
    fn column(&mut self, basis: &Basis, j: usize, k: usize) -> Vec<f64> {
        let h = basis.h_col(j);
        let mut x = vec![0.0; k];
        for (xi, &hi) in x.iter_mut().zip(h) {
            *xi = -self.gh * hi;
        }
        x[j] += 1.0;
        for &v in &x {
            self.mmax = self.mmax.max(v.abs());
        }
        x
    }

    /// Finalize the closed columns that extend the closed prefix.
    pub fn update(&mut self, basis: &Basis) {
        let _g = crate::prof::section("givens update");
        let k = basis.dim();
        while self.r.len() < k {
            let j = self.r.len();
            if basis.closed_row(j).is_none() {
                break;
            }
            let mut x = self.column(basis, j, k);
            for rot in &self.rots {
                rot.apply(&mut x);
            }
            if let Some(last) = (j + 1..k).rev().find(|&i| x[i] != 0.0) {
                for i in (j..last).rev() {
                    let rot = Rot::zeroing(&mut x, i);
                    self.rots.push(rot);
                }
            }
            x.truncate(j + 1);
            self.r.push(x);
        }
    }

    /// Apply the stored rotations to `rhs` and to the open columns, and
    /// triangularize the trailing block.
    fn tail(&mut self, basis: &Basis, rhs: &[f64]) -> Result<Tail, Singular> {
        let k = basis.dim();
        let nf = self.r.len();
        let mut b = rhs.to_vec();
        b.resize(k, 0.0);
        for rot in &self.rots {
            rot.apply(&mut b);
        }
        let mut cols: Vec<Vec<f64>> = (nf..k)
            .map(|j| {
                let mut x = self.column(basis, j, k);
                for rot in &self.rots {
                    rot.apply(&mut x);
                }
                x
            })
            .collect();
        for jj in 0..cols.len() {
            let j = nf + jj;
            for i in (j..k - 1).rev() {
                if cols[jj][i + 1] == 0.0 {
                    continue;
                }
                let rot = Rot::zeroing(&mut cols[jj], i);
                for c in cols[jj + 1..].iter_mut() {
                    rot.apply(c);
                }
                rot.apply(&mut b);
            }
        }
        // singularity test on all pivots, relative to max |M_ij|
        let tol = SING_TOL * self.mmax;
        let pivots_ok = self.r.iter().enumerate().all(|(j, c)| c[j].abs() > tol)
            && cols
                .iter()
                .enumerate()
                .all(|(jj, c)| c[nf + jj].abs() > tol);
        if !pivots_ok {
            return Err(Singular);
        }
        Ok(Tail { nf, cols, b })
    }

    /// Solve `M l = rhs` (rhs padded with zeros to the current dimension).
    pub fn solve(&mut self, basis: &Basis, rhs: &[f64]) -> Result<Vec<f64>, Singular> {
        let _g = crate::prof::section("givens solve");
        let Tail { nf, cols, mut b } = self.tail(basis, rhs)?;
        let k = b.len();
        for j in (0..k).rev() {
            let col: &[f64] = if j < nf { &self.r[j] } else { &cols[j - nf] };
            b[j] /= col[j];
            let bj = b[j];
            for i in 0..j {
                b[i] -= col[i] * bj;
            }
        }
        Ok(b)
    }

    /// Only the trailing components `l[nf..k]` of the solution (the open
    /// columns), without the full back substitution: what the Galerkin
    /// residual `−γh Σ_open l_j p_j` needs. Returns `(nf, l[nf..k])`.
    pub fn solve_tail(
        &mut self,
        basis: &Basis,
        rhs: &[f64],
    ) -> Result<(usize, Vec<f64>), Singular> {
        let _g = crate::prof::section("givens tail");
        let Tail { nf, cols, b } = self.tail(basis, rhs)?;
        let k = b.len();
        let mut l = b[nf..].to_vec();
        for j in (nf..k).rev() {
            let col = &cols[j - nf];
            l[j - nf] /= col[j];
            let lj = l[j - nf];
            for i in nf..j {
                l[i - nf] -= col[i] * lj;
            }
        }
        Ok((nf, l))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::krylov::tests::{DenseOp, Rng};

    /// The former dense-LU stage solve (test oracle).
    fn lu_solve(basis: &Basis, rhs: &[f64], gh: f64) -> Vec<f64> {
        crate::map::solve_projected(basis, rhs, gh).unwrap()
    }

    /// MAP-like sequence (Arnoldi, insertions, closure rounds): the Givens
    /// solve equals the dense LU after every operation, also for the tail.
    #[test]
    fn matches_dense_lu() {
        let n = 30;
        let mut rng = Rng::new(5);
        let op = DenseOp::random(n, &mut rng);
        let mut av = |v: &[f64], o: &mut [f64]| op.apply(v, o);
        let gh = 0.37;
        let mut basis = Basis::new(n);
        let mut giv = Givens::new(gh);
        let w1 = rng.vec(n);
        let (_, c1) = basis.insert(&w1, &mut av).unwrap();
        let check = |basis: &Basis, giv: &mut Givens, c: &[f64]| {
            giv.update(basis);
            let l = giv.solve(basis, c).unwrap();
            let lref = lu_solve(basis, c, gh);
            let scale = lref.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            for (a, b) in l.iter().zip(&lref) {
                assert!((a - b).abs() <= 1e-12 * scale, "solve {a} vs {b}");
            }
            let (nf, lt) = giv.solve_tail(basis, c).unwrap();
            for (a, b) in lt.iter().zip(&lref[nf..]) {
                assert!((a - b).abs() <= 1e-12 * scale, "tail {a} vs {b}");
            }
        };
        for _ in 0..4 {
            check(&basis, &mut giv, &c1);
            let j = *basis.pending_cols().last().unwrap();
            basis.promote(j, &mut av).unwrap();
        }
        // two insertions, then two closure rounds (FIFO)
        for _ in 0..2 {
            let (_, c) = basis.insert(&rng.vec(n), &mut av).unwrap();
            check(&basis, &mut giv, &c);
        }
        for _ in 0..2 {
            for j in basis.pending_cols() {
                basis.promote(j, &mut av).unwrap();
                check(&basis, &mut giv, &c1);
            }
        }
    }
}
