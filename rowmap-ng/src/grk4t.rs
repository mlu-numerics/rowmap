//! ROW method coefficients (single source).
//!
//! The data structure is method-generic (AGENTS.md §3.5). Defined: GRK4T
//! (Kaps–Rentrop, the ROWMAP method) and ROS34PW2 (Rang–Angermann, order 3;
//! phase 5b).
//!
//! Notation (Hairer–Wanner II, §IV.7; paper eq. (1.2)): stage `i`
//! ```text
//! (I − γhA) k̃_i = w_i,     k̃_i = k_i + Σ_{j<i} c_ij k_j,   c_ij = γ_ij/γ
//! w_i = f(t + α_i h, u + h Σ_{j<i} a_ij k_j) + Σ_{j<i} c_ij k_j [+ h d_i f_t]
//! u_new = u + h Σ b_i k_i,    err = h Σ (b̂_i − b_i) k_i
//! α_i = Σ_j a_ij,   d_i = γ + Σ_{j<i} γ_ij
//! ```

/// Maximal number of stages of any method table.
pub(crate) const MAXS: usize = 4;

/// Coefficients of an s-stage ROW method with embedded solution.
#[derive(Debug, Clone)]
pub(crate) struct RowMethod {
    /// Descriptive metadata (name and orders); kept for the order-3 method
    /// of phase 5b and for diagnostics.
    #[allow(dead_code)]
    pub name: &'static str,
    /// Number of stages.
    pub s: usize,
    /// Classical order of the main method and of the embedded method.
    /// (Also the minimal stage-1 Krylov dimension κ₁.)
    pub order: usize,
    #[allow(dead_code)]
    pub order_embedded: usize,
    /// Diagonal γ.
    pub gamma: f64,
    /// Strictly lower `a_ij` (stage arguments).
    pub a: [[f64; MAXS]; MAXS],
    /// Strictly lower `γ_ij`.
    pub g: [[f64; MAXS]; MAXS],
    /// Weights of the main and the embedded solution.
    pub b: [f64; MAXS],
    pub bhat: [f64; MAXS],
}

impl RowMethod {
    /// `c_ij = γ_ij / γ`.
    pub fn c(&self, i: usize, j: usize) -> f64 {
        self.g[i][j] / self.gamma
    }

    /// Time node `α_i = Σ_j a_ij`.
    pub fn alpha(&self, i: usize) -> f64 {
        self.a[i][..i].iter().sum()
    }

    /// Coefficient of `h f_t` in stage `i`: `d_i = γ + Σ_{j<i} γ_ij`.
    pub fn d(&self, i: usize) -> f64 {
        self.gamma + self.g[i][..i].iter().sum::<f64>()
    }

    /// Error weights `e_i = b̂_i − b_i`.
    pub fn e(&self, i: usize) -> f64 {
        self.bhat[i] - self.b[i]
    }

    /// If stage `i` evaluates f at exactly the same argument as an earlier
    /// stage `k` (identical rows of `a`, identical time node), returns `k`.
    /// For GRK4T: stage 4 reuses stage 3 (`f₄ = f₃`). The MAP relies on this
    /// (`w₄ ∈ K₃` in the autonomous case).
    pub fn reuses_f(&self, i: usize) -> Option<usize> {
        (1..i).find(|&k| {
            (0..i).all(|j| {
                let ak = if j < k { self.a[k][j] } else { 0.0 };
                self.a[i][j] == ak
            })
        })
    }
}

/// GRK4T (Kaps–Rentrop 1979), coefficients as in ROS4 / rowmap.F (`rowcoe`).
/// `bhat` are the embedded weights (rowmap.F stores `e = bhat − b`).
pub(crate) const GRK4T: RowMethod = RowMethod {
    name: "grk4t",
    s: 4,
    order: 4,
    order_embedded: 3,
    gamma: 0.231,
    a: [
        [0.0, 0.0, 0.0, 0.0],
        [0.462, 0.0, 0.0, 0.0],
        [-0.081_566_816_832_72, 0.961_775_150_166_06, 0.0, 0.0],
        [-0.081_566_816_832_72, 0.961_775_150_166_06, 0.0, 0.0],
    ],
    g: [
        [0.0, 0.0, 0.0, 0.0],
        [-0.270_629_667_752_44, 0.0, 0.0, 0.0],
        [0.311_254_483_294_09, 0.008_524_456_284_82, 0.0, 0.0],
        [
            0.282_816_832_043_53,
            -0.457_959_483_280_73,
            -0.111_208_333_333_33,
            0.0,
        ],
    ],
    b: [
        0.217_487_371_652_73,
        0.486_229_037_990_12,
        0.0,
        0.296_283_590_357_15,
    ],
    bhat: [
        -0.717_088_504_499_33,
        1.776_179_121_761_04,
        -0.059_090_617_261_71,
        0.0,
    ],
};

/// ROS34PW2 (Rang & Angermann, BIT 45 (2005) 761–787): 4 stages, order 3
/// with embedded order 2, stiffly accurate, L-stable, W-method of order 2.
/// Coefficients as in OrdinaryDiffEq.jl v6.66.0 (`ROS34PW2Tableau`, given
/// there in the original α, Γ, b, b̂ form). The order-3 option
/// (`Method::Ros34pw2`). ROS3P was examined first and dropped: its embedded
/// estimate vanishes identically on linear problems (docs/map.md §8).
// literals copied verbatim from the source (17 significant digits)
#[allow(clippy::excessive_precision)]
pub(crate) const ROS34PW2: RowMethod = RowMethod {
    name: "ros34pw2",
    s: 4,
    order: 3,
    order_embedded: 2,
    gamma: 4.358_665_215_084_590_0e-1,
    a: [
        [0.0, 0.0, 0.0, 0.0],
        [8.717_330_430_169_180_1e-1, 0.0, 0.0, 0.0],
        [
            8.445_706_001_536_942_3e-1,
            -1.129_906_423_648_418_5e-1,
            0.0,
            0.0,
        ],
        [0.0, 0.0, 1.0, 0.0],
    ],
    g: [
        [0.0, 0.0, 0.0, 0.0],
        [-8.717_330_430_169_180_1e-1, 0.0, 0.0, 0.0],
        [
            -9.033_805_701_304_408_2e-1,
            5.418_067_238_809_532_6e-2,
            0.0,
            0.0,
        ],
        [
            2.421_238_070_609_534_6e-1,
            -1.223_250_583_904_514_7,
            5.452_602_553_351_021_4e-1,
            0.0,
        ],
    ],
    b: [
        2.421_238_070_609_534_6e-1,
        -1.223_250_583_904_514_7,
        1.545_260_255_335_102_0,
        4.358_665_215_084_590_0e-1,
    ],
    bhat: [
        3.781_090_314_581_936_9e-1,
        -9.604_229_221_242_317_8e-2,
        0.5,
        2.179_332_607_542_295_0e-1,
    ],
};

#[cfg(test)]
mod tests {
    use super::*;

    /// β_ij = a_ij + γ_ij (j<i), β_ii = γ — the matrix of the ROW
    /// order conditions (Hairer–Wanner II, Table IV.7.1).
    fn beta(m: &RowMethod) -> [[f64; MAXS]; MAXS] {
        let mut bt = [[0.0; MAXS]; MAXS];
        for i in 0..m.s {
            for j in 0..i {
                bt[i][j] = m.a[i][j] + m.g[i][j];
            }
        }
        bt
    }

    /// Residuals of the ROW order conditions up to order 4 for weights `w`
    /// (8 trees). Uses β'_i = Σ_{j<i} β_ij and α_i.
    fn order_conditions(m: &RowMethod, w: &[f64; MAXS]) -> Vec<(usize, f64)> {
        let s = m.s;
        let g = m.gamma;
        let bt = beta(m);
        let al: Vec<f64> = (0..s).map(|i| m.alpha(i)).collect();
        let bp: Vec<f64> = (0..s).map(|i| bt[i][..i].iter().sum()).collect();
        let sum = |f: &dyn Fn(usize) -> f64| (0..s).map(f).sum::<f64>();
        let bb = |i: usize, v: &[f64]| (0..i).map(|k| bt[i][k] * v[k]).sum::<f64>();
        let aa = |i: usize, v: &[f64]| (0..i).map(|k| m.a[i][k] * v[k]).sum::<f64>();
        let al2: Vec<f64> = al.iter().map(|x| x * x).collect();
        let bbp: Vec<f64> = (0..s).map(|k| bb(k, &bp)).collect();
        vec![
            (1, sum(&|i| w[i]) - 1.0),
            (2, sum(&|i| w[i] * bp[i]) - (0.5 - g)),
            (3, sum(&|i| w[i] * al2[i]) - 1.0 / 3.0),
            (3, sum(&|i| w[i] * bb(i, &bp)) - (1.0 / 6.0 - g + g * g)),
            (4, sum(&|i| w[i] * al[i] * al2[i]) - 0.25),
            (
                4,
                sum(&|i| w[i] * al[i] * aa(i, &bp)) - (1.0 / 8.0 - g / 3.0),
            ),
            (4, sum(&|i| w[i] * bb(i, &al2)) - (1.0 / 12.0 - g / 3.0)),
            (
                4,
                sum(&|i| w[i] * bb(i, &bbp)) - (1.0 / 24.0 - g / 2.0 + 1.5 * g * g - g * g * g),
            ),
        ]
    }

    /// Stability function R(z) = 1 + z wᵀ (I − zB)⁻¹ 𝟙, B = (β_ij) incl. diagonal γ.
    fn stab(m: &RowMethod, w: &[f64; MAXS], z: (f64, f64)) -> (f64, f64) {
        // complex forward substitution, (I − zB) lower triangular
        let s = m.s;
        let bt = beta(m);
        let mul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
        let div = |a: (f64, f64), b: (f64, f64)| {
            let d = b.0 * b.0 + b.1 * b.1;
            ((a.0 * b.0 + a.1 * b.1) / d, (a.1 * b.0 - a.0 * b.1) / d)
        };
        let mut x = vec![(0.0, 0.0); s];
        for i in 0..s {
            let mut r = (1.0, 0.0);
            for (k, xk) in x.iter().enumerate().take(i) {
                let t = mul(mul(z, (bt[i][k], 0.0)), *xk);
                r = (r.0 + t.0, r.1 + t.1);
            }
            let diag = (1.0 - z.0 * m.gamma, -z.1 * m.gamma);
            x[i] = div(r, diag);
        }
        let mut acc = (0.0, 0.0);
        for i in 0..s {
            acc = (acc.0 + w[i] * x[i].0, acc.1 + w[i] * x[i].1);
        }
        let t = mul(z, acc);
        (1.0 + t.0, t.1)
    }

    fn abs(c: (f64, f64)) -> f64 {
        c.0.hypot(c.1)
    }

    #[test]
    fn grk4t_order_conditions() {
        let m = &GRK4T;
        for (p, r) in order_conditions(m, &m.b) {
            assert!(
                r.abs() < 1e-12,
                "main method, order-{p} condition residual {r:e}"
            );
        }
        let emb = order_conditions(m, &m.bhat);
        for &(p, r) in emb.iter().filter(|(p, _)| *p <= 3) {
            assert!(
                r.abs() < 1e-12,
                "embedded, order-{p} condition residual {r:e}"
            );
        }
        // embedded method is genuinely of order 3, not 4
        assert!(
            emb.iter()
                .filter(|(p, _)| *p == 4)
                .any(|(_, r)| r.abs() > 1e-3)
        );
    }

    #[test]
    fn grk4t_structure() {
        let m = &GRK4T;
        assert_eq!(m.b[2], 0.0, "b3 = 0 (Thm 2.1)");
        assert_eq!(m.reuses_f(3), Some(2), "f4 = f3");
        assert_eq!(m.reuses_f(1), None);
        assert_eq!(m.reuses_f(2), None);
        assert!((m.alpha(1) - 0.462).abs() < 1e-15);
        assert!((m.alpha(2) - 0.880_208_333_333_34).abs() < 1e-13);
        // d_i: f_t coefficients; d3 != d4 is why w4 is not in K3 in the
        // non-autonomous case (AGENTS.md §3.6)
        assert!((m.d(0) - 0.231).abs() < 1e-15);
        assert!((m.d(2) - 0.550_778_939_578_91).abs() < 1e-12);
        assert!((m.d(3) - (-0.055_350_984_570_53)).abs() < 1e-12);
    }

    #[test]
    fn grk4t_stability() {
        let m = &GRK4T;
        // values from AGENTS.md §3.1
        let rinf = abs(stab(m, &m.b, (-1e12, 0.0)));
        assert!((rinf - 0.4536).abs() < 1e-3, "|R(inf)| = {rinf}");
        let rhinf = abs(stab(m, &m.bhat, (-1e12, 0.0)));
        assert!((rhinf - 2.6023).abs() < 1e-3, "|Rhat(inf)| = {rhinf}");
        let (mut rmax, mut ymax) = (0.0, 0.0);
        for k in 0..20000 {
            let y = 10f64.powf(-3.0 + 9.0 * k as f64 / 20000.0);
            let r = abs(stab(m, &m.b, (0.0, y)));
            if r > rmax {
                (rmax, ymax) = (r, y);
            }
        }
        assert!((rmax - 1.028215).abs() < 1e-5, "max |R(iy)| = {rmax}");
        assert!((ymax - 4.43).abs() < 0.05, "at y = {ymax}");
        // A-stable on the negative real axis
        for k in 0..2000 {
            let x = -(10f64.powf(-3.0 + 12.0 * k as f64 / 2000.0));
            assert!(abs(stab(m, &m.b, (x, 0.0))) <= 1.0 + 1e-14);
        }
        // consistency: R(z) ≈ e^z for small z (order 4 ⇒ error O(z^5))
        let z = 1e-2;
        assert!((stab(m, &m.b, (z, 0.0)).0 - z.exp()).abs() < 1e-11);
    }

    #[test]
    fn ros34pw2_order_conditions() {
        let m = &ROS34PW2;
        let main = order_conditions(m, &m.b);
        for &(p, r) in main.iter().filter(|(p, _)| *p <= 3) {
            assert!(r.abs() < 1e-14, "ROS34PW2 order-{p} residual {r:e}");
        }
        assert!(
            main.iter()
                .filter(|(p, _)| *p == 4)
                .any(|(_, r)| r.abs() > 1e-4)
        );
        let emb = order_conditions(m, &m.bhat);
        for &(p, r) in emb.iter().filter(|(p, _)| *p <= 2) {
            assert!(
                r.abs() < 1e-14,
                "ROS34PW2 embedded order-{p} residual {r:e}"
            );
        }
    }

    #[test]
    fn ros34pw2_structure_and_stability() {
        let m = &ROS34PW2;
        // stiffly accurate: b = last row of (α + Γ), with b₄ = γ
        for j in 0..3 {
            assert!((m.b[j] - (m.a[3][j] + m.g[3][j])).abs() < 1e-15, "b{j}");
        }
        assert!((m.b[3] - m.gamma).abs() < 1e-15);
        assert_eq!(
            (1..4).filter_map(|i| m.reuses_f(i)).count(),
            0,
            "no stage reuse"
        );
        let rinf = abs(stab(m, &m.b, (-1e12, 0.0)));
        assert!(rinf < 1e-10, "L-stable: |R(inf)| = {rinf:e}");
        let rhinf = abs(stab(m, &m.bhat, (-1e12, 0.0)));
        assert!((rhinf - 0.48).abs() < 0.01, "|Rhat(inf)| = {rhinf}");
        let mut rmax: f64 = 0.0;
        for k in 0..20000 {
            let y = 10f64.powf(-3.0 + 9.0 * k as f64 / 20000.0);
            rmax = rmax.max(abs(stab(m, &m.b, (0.0, y))));
        }
        assert!(rmax <= 1.0 + 1e-12, "A-stable: max |R(iy)| = {rmax}");
        // the embedded estimate sees linear problems: R̂ ≠ R (ROS3P failed here)
        let d = (stab(m, &m.b, (-0.5, 0.0)).0 - stab(m, &m.bhat, (-0.5, 0.0)).0).abs();
        assert!(d > 1e-4, "R − R̂ at z = −0.5: {d:e}");
    }
}
