//! Rust versions of the test problems (docs/problems.md), shared by the
//! integration tests and the examples. Must agree with
//! `fortran/experiments/problems.f90` (checked against the committed dumps in
//! `tests/data/dump/` by `tests/problems_consistency.rs`).

#![allow(dead_code)]

use rowmap::{OdeSystem, Unsupported};
use std::f64::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Heat,
    Nilidi,
    Bruss,
    Food,
    Vdpl,
    Difkon,
    Diffu2,
}

pub const AUTONOMOUS: [Kind; 5] = [
    Kind::Heat,
    Kind::Nilidi,
    Kind::Bruss,
    Kind::Food,
    Kind::Vdpl,
];
pub const ALL: [Kind; 7] = [
    Kind::Heat,
    Kind::Nilidi,
    Kind::Bruss,
    Kind::Food,
    Kind::Vdpl,
    Kind::Difkon,
    Kind::Diffu2,
];

impl Kind {
    pub fn name(self) -> &'static str {
        match self {
            Kind::Heat => "heat",
            Kind::Nilidi => "nilidi",
            Kind::Bruss => "bruss",
            Kind::Food => "food",
            Kind::Vdpl => "vdpl",
            Kind::Difkon => "difkon",
            Kind::Diffu2 => "diffu2",
        }
    }
    pub fn paper_m(self) -> usize {
        match self {
            Kind::Heat | Kind::Bruss | Kind::Food => 50,
            Kind::Nilidi | Kind::Difkon | Kind::Diffu2 => 69,
            Kind::Vdpl => 20,
        }
    }
    pub fn coarse_m(self) -> usize {
        match self {
            Kind::Vdpl => 8,
            _ => 10,
        }
    }
    fn dirichlet(self) -> bool {
        matches!(
            self,
            Kind::Heat | Kind::Nilidi | Kind::Difkon | Kind::Diffu2
        )
    }
}

/// A semi-discretized test problem.
pub struct Problem {
    pub kind: Kind,
    pub m: usize,
    /// points per component
    np: usize,
    h: f64,
    xs: Vec<f64>,
}

impl Problem {
    pub fn new(kind: Kind, m: usize) -> Self {
        let (h, xs): (f64, Vec<f64>) = if kind.dirichlet() {
            let len = if kind == Kind::Nilidi { PI / 3.0 } else { 1.0 };
            let h = len / (m + 1) as f64;
            (h, (1..=m).map(|i| i as f64 * h).collect())
        } else {
            let h = 1.0 / (m - 1) as f64;
            (h, (0..m).map(|i| i as f64 * h).collect())
        };
        Problem {
            kind,
            m,
            np: m * m,
            h,
            xs,
        }
    }

    pub fn n(&self) -> usize {
        if self.kind.dirichlet() {
            self.np
        } else {
            2 * self.np
        }
    }

    pub fn t_end(&self) -> f64 {
        match self.kind {
            Kind::Heat => 10.24,
            Kind::Food => 3.0,
            Kind::Vdpl => 1000.0,
            _ => 1.0,
        }
    }

    pub fn y0(&self) -> Vec<f64> {
        let m = self.m;
        let mut y = vec![0.0; self.n()];
        let np = self.np;
        for j in 0..m {
            for i in 0..m {
                let (x, yy) = (self.xs[i], self.xs[j]);
                let k = i + j * m;
                match self.kind {
                    Kind::Heat => y[k] = 16.0 * x * (1.0 - x) * yy * (1.0 - yy),
                    Kind::Nilidi => y[k] = (3.0 * x).sin() * (3.0 * yy).sin(),
                    Kind::Bruss => {
                        y[k] = 0.5 + yy;
                        y[np + k] = 1.0 + 5.0 * x;
                    }
                    Kind::Food => {
                        y[k] = 10.0 - 5.0 * (PI * x).cos() * (10.0 * PI * yy).cos();
                        y[np + k] = 17.0 + 5.0 * (10.0 * PI * x).cos() * (PI * yy).cos();
                    }
                    Kind::Vdpl => {
                        y[k] = 1.0 - (PI * x).cos() * (2.0 * PI * yy).cos();
                        y[np + k] = 1.0 + 2.0 * (2.0 * PI * x).cos() * (PI * yy).cos();
                    }
                    Kind::Difkon => y[k] = (PI * x).sin() * (PI * yy).sin() * (1.0 + 4.0 * x * yy),
                    Kind::Diffu2 => y[k] = (PI * x).sin() * (PI * yy).sin(),
                }
            }
        }
        y
    }

    /// Dirichlet Laplacian (zero boundary values), m×m inner points.
    fn lap_dir(&self, u: &[f64], out: &mut [f64]) {
        let m = self.m;
        let h2 = self.h * self.h;
        for j in 0..m {
            for i in 0..m {
                let k = i + j * m;
                let mut s = -4.0 * u[k];
                if i > 0 {
                    s += u[k - 1];
                }
                if i + 1 < m {
                    s += u[k + 1];
                }
                if j > 0 {
                    s += u[k - m];
                }
                if j + 1 < m {
                    s += u[k + m];
                }
                out[k] = s / h2;
            }
        }
    }

    /// Neumann Laplacian, vertex grid with ghost mirror.
    fn lap_neu(&self, u: &[f64], out: &mut [f64]) {
        let m = self.m;
        let h2 = self.h * self.h;
        let lo = |i: usize| if i > 0 { i - 1 } else { 1 };
        let hi = |i: usize| if i + 1 < m { i + 1 } else { m - 2 };
        for j in 0..m {
            for i in 0..m {
                let k = i + j * m;
                out[k] =
                    (u[lo(i) + j * m] + u[hi(i) + j * m] + u[i + lo(j) * m] + u[i + hi(j) * m]
                        - 4.0 * u[k])
                        / h2;
            }
        }
    }

    /// 0.1 Δu − y D_x u + x D_y u (DIFKON operator).
    fn conv_diff(&self, u: &[f64], out: &mut [f64]) {
        let m = self.m;
        self.lap_dir(u, out);
        for j in 0..m {
            for i in 0..m {
                let k = i + j * m;
                let a = if i + 1 < m { u[k + 1] } else { 0.0 };
                let b = if i > 0 { u[k - 1] } else { 0.0 };
                let dx = (a - b) / (2.0 * self.h);
                let a = if j + 1 < m { u[k + m] } else { 0.0 };
                let b = if j > 0 { u[k - m] } else { 0.0 };
                let dy = (a - b) / (2.0 * self.h);
                out[k] = 0.1 * out[k] - self.xs[j] * dx + self.xs[i] * dy;
            }
        }
    }

    /// Source term and its time derivative (DIFKON, DIFFU2).
    fn source(&self, t: f64, g: &mut [f64], gt: &mut [f64]) {
        let m = self.m;
        let eps = 0.1;
        let e = (-2.0 * eps * PI * PI * t).exp();
        for j in 0..m {
            for i in 0..m {
                let (x, y) = (self.xs[i], self.xs[j]);
                let k = i + j * m;
                let s = (PI * x).sin() * (PI * y).sin();
                let cx = (PI * x).cos() * (PI * y).sin();
                let cy = (PI * x).sin() * (PI * y).cos();
                if self.kind == Kind::Difkon {
                    let g0 = 2.0 * eps * PI * PI * s + PI * (y * cx - x * cy);
                    let g1 = -8.0 * eps * PI * (y * cx + x * cy)
                        + 4.0 * PI * x * y * (y * cx - x * cy)
                        + 4.0 * s * (y * y - x * x);
                    g[k] = g0 + e * g1;
                    gt[k] = -2.0 * eps * PI * PI * e * g1;
                } else {
                    let w = y * cx + x * cy;
                    g[k] = 4.0 * x * y * s * t.cos()
                        + 2.0 * PI * PI * s * (1.0 + 4.0 * x * y * t.sin())
                        - 8.0 * PI * t.sin() * w;
                    gt[k] = -4.0 * x * y * s * t.sin() + 8.0 * PI * PI * x * y * s * t.cos()
                        - 8.0 * PI * t.cos() * w;
                }
            }
        }
    }
}

impl OdeSystem for Problem {
    fn dim(&self) -> usize {
        self.n()
    }

    fn autonomous(&self) -> bool {
        !matches!(self.kind, Kind::Difkon | Kind::Diffu2)
    }

    fn rhs(&self, t: f64, y: &[f64], f: &mut [f64]) {
        let np = self.np;
        match self.kind {
            Kind::Heat => self.lap_dir(y, f),
            Kind::Nilidi => {
                self.lap_dir(y, f);
                for k in 0..np {
                    let eu = y[k].exp();
                    f[k] = eu * f[k] + y[k] * (18.0 * eu - 1.0);
                }
            }
            Kind::Bruss | Kind::Food | Kind::Vdpl => {
                let (u, v) = y.split_at(np);
                let (fu, fv) = f.split_at_mut(np);
                self.lap_neu(u, fu);
                self.lap_neu(v, fv);
                for k in 0..np {
                    let (a, b) = (u[k], v[k]);
                    match self.kind {
                        Kind::Bruss => {
                            fu[k] = 0.02 * fu[k] + 1.0 + a * a * b - 4.0 * a;
                            fv[k] = 0.02 * fv[k] + 3.0 * a - a * a * b;
                        }
                        Kind::Food => {
                            fu[k] = 0.05 * fu[k] + a * (1.0 - 0.1 * b);
                            fv[k] += b * (-1000.0 + 100.0 * a);
                        }
                        _ => {
                            fu[k] = 0.05 * fu[k] + b;
                            fv[k] += 100.0 * (1.0 - a * a) * b - a;
                        }
                    }
                }
            }
            Kind::Difkon | Kind::Diffu2 => {
                if self.kind == Kind::Difkon {
                    self.conv_diff(y, f);
                } else {
                    self.lap_dir(y, f);
                }
                let mut g = vec![0.0; np];
                let mut gt = vec![0.0; np];
                self.source(t, &mut g, &mut gt);
                for k in 0..np {
                    f[k] += g[k];
                }
            }
        }
    }

    fn jac_vec(&self, _t: f64, y: &[f64], v: &[f64], z: &mut [f64]) {
        let np = self.np;
        match self.kind {
            Kind::Heat | Kind::Diffu2 => self.lap_dir(v, z),
            Kind::Difkon => self.conv_diff(v, z),
            Kind::Nilidi => {
                let mut ly = vec![0.0; np];
                self.lap_dir(y, &mut ly);
                self.lap_dir(v, z);
                for k in 0..np {
                    let eu = y[k].exp();
                    z[k] = eu * (z[k] + ly[k] * v[k]) + (18.0 * eu * (1.0 + y[k]) - 1.0) * v[k];
                }
            }
            Kind::Bruss | Kind::Food | Kind::Vdpl => {
                let (u, w) = y.split_at(np);
                let (va, vb) = v.split_at(np);
                let (za, zb) = z.split_at_mut(np);
                self.lap_neu(va, za);
                self.lap_neu(vb, zb);
                for k in 0..np {
                    let (a, b) = (u[k], w[k]);
                    match self.kind {
                        Kind::Bruss => {
                            za[k] = 0.02 * za[k] + (2.0 * a * b - 4.0) * va[k] + a * a * vb[k];
                            zb[k] = 0.02 * zb[k] + (3.0 - 2.0 * a * b) * va[k] - a * a * vb[k];
                        }
                        Kind::Food => {
                            za[k] = 0.05 * za[k] + (1.0 - 0.1 * b) * va[k] - 0.1 * a * vb[k];
                            zb[k] += 100.0 * b * va[k] + (-1000.0 + 100.0 * a) * vb[k];
                        }
                        _ => {
                            za[k] = 0.05 * za[k] + vb[k];
                            zb[k] += (-200.0 * a * b - 1.0) * va[k] + 100.0 * (1.0 - a * a) * vb[k];
                        }
                    }
                }
            }
        }
    }

    fn dfdt(&self, t: f64, _y: &[f64], ft: &mut [f64]) -> Result<(), Unsupported> {
        match self.kind {
            Kind::Difkon | Kind::Diffu2 => {
                let mut g = vec![0.0; self.np];
                self.source(t, &mut g, ft);
                Ok(())
            }
            _ => Err(Unsupported),
        }
    }
}

// ---------------------------------------------------------------- data files

pub fn data_dir() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

/// Read a whitespace table, skipping `#` lines.
pub fn read_table(path: &std::path::Path) -> Vec<Vec<f64>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.trim_start().starts_with('#'))
        .map(|l| l.split_whitespace().map(|x| x.parse().unwrap()).collect())
        .collect()
}

/// Reference solution at t_end (DOPRI5, 1e-10).
pub fn reference(kind: Kind, m: usize) -> Vec<f64> {
    read_table(&data_dir().join(format!("reference/{}_m{m}.txt", kind.name())))
        .into_iter()
        .map(|r| r[0])
        .collect()
}

/// Paper error measure ERR = sqrt(1/n Σ ((y − r)/(1 + |r|))²).
pub fn err_paper(y: &[f64], r: &[f64]) -> f64 {
    let s: f64 = y
        .iter()
        .zip(r)
        .map(|(a, b)| ((a - b) / (1.0 + b.abs())).powi(2))
        .sum();
    (s / y.len() as f64).sqrt()
}
