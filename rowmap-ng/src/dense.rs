//! Small dense linear algebra for the projected stage systems
//! `(I − γh H) l = Qᵀw` (dimension κ ≤ max_krylov, typically < 100).
//!
//! The integrator uses the incremental Givens QR (`crate::givens`, which
//! exploits the generalized Hessenberg structure, docs/map.md §9). The dense
//! LU here, rebuilt from scratch for every solve, is kept as the test oracle.

/// Square matrix, column-major.
#[cfg(test)]
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Mat {
    pub k: usize,
    pub v: Vec<f64>,
}

/// The matrix is numerically singular (pivot below the relative threshold).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Singular;

/// Relative pivot threshold: a pivot `|u_jj| ≤ SING_TOL · ‖M‖_max` counts as
/// singular. Chosen a few orders above machine epsilon, so that the step is
/// rejected (and h halved) before the projected solve amplifies rounding errors
/// by more than ~1e12.
pub(crate) const SING_TOL: f64 = 1e3 * f64::EPSILON;

#[cfg(test)]
impl Mat {
    pub fn zeros(k: usize) -> Self {
        Mat {
            k,
            v: vec![0.0; k * k],
        }
    }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.v[i + j * self.k]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, x: f64) {
        self.v[i + j * self.k] = x;
    }

    /// `y = M x`.
    #[cfg(test)]
    pub fn mul_vec(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; self.k];
        for j in 0..self.k {
            for i in 0..self.k {
                y[i] += self.get(i, j) * x[j];
            }
        }
        y
    }

    /// Solve `M x = b` by LU with partial pivoting (M is copied).
    pub fn solve(&self, b: &[f64]) -> Result<Vec<f64>, Singular> {
        let k = self.k;
        assert_eq!(b.len(), k);
        let mut a = self.v.clone();
        let mut x = b.to_vec();
        let at = |i: usize, j: usize| i + j * k;
        let mnorm = a.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        if k == 0 {
            return Ok(x);
        }
        if mnorm == 0.0 || !mnorm.is_finite() {
            return Err(Singular);
        }
        for j in 0..k {
            // pivot search in column j
            let (p, pv) = (j..k)
                .map(|i| (i, a[at(i, j)].abs()))
                .fold((j, -1.0), |best, c| if c.1 > best.1 { c } else { best });
            if pv <= SING_TOL * mnorm {
                return Err(Singular);
            }
            if p != j {
                for c in 0..k {
                    a.swap(at(j, c), at(p, c));
                }
                x.swap(j, p);
            }
            let d = a[at(j, j)];
            for i in j + 1..k {
                let f = a[at(i, j)] / d;
                if f != 0.0 {
                    a[at(i, j)] = f;
                    for c in j + 1..k {
                        a[at(i, c)] -= f * a[at(j, c)];
                    }
                    x[i] -= f * x[j];
                }
            }
        }
        for j in (0..k).rev() {
            let mut s = x[j];
            for c in j + 1..k {
                s -= a[at(j, c)] * x[c];
            }
            x[j] = s / a[at(j, j)];
        }
        Ok(x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solve_with_pivoting() {
        // needs a row swap (zero leading entry)
        let mut m = Mat::zeros(3);
        let rows = [[0.0, 2.0, 1.0], [1.0, 1.0, 0.0], [3.0, 0.0, 1.0]];
        for (i, r) in rows.iter().enumerate() {
            for (j, &x) in r.iter().enumerate() {
                m.set(i, j, x);
            }
        }
        let xt = [1.0, -2.0, 0.5];
        let b = m.mul_vec(&xt);
        let x = m.solve(&b).unwrap();
        for (a, e) in x.iter().zip(xt) {
            assert!((a - e).abs() < 1e-14);
        }
    }

    #[test]
    fn singular_detected() {
        let mut m = Mat::zeros(2);
        m.set(0, 0, 1.0);
        m.set(0, 1, 2.0);
        m.set(1, 0, 2.0);
        m.set(1, 1, 4.0 + 1e-15);
        assert_eq!(m.solve(&[1.0, 1.0]), Err(Singular));
        assert_eq!(Mat::zeros(2).solve(&[1.0, 1.0]), Err(Singular));
        assert_eq!(Mat::zeros(0).solve(&[]), Ok(vec![]));
    }
}
