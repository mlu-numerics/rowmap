//! BLAS-1 style kernels on `&[f64]`.

/// Number of independent partial sums in reductions. Without fast-math the
/// compiler cannot reorder a single sum; independent accumulators let it use
/// SIMD (phase 7). This changes the summation order, not the accuracy class.
const LANES: usize = 8;

/// Σ g(x_i, y_i) with `LANES` independent accumulators.
#[inline(always)]
fn reduce2(x: &[f64], y: &[f64], g: impl Fn(f64, f64) -> f64) -> f64 {
    debug_assert_eq!(x.len(), y.len());
    let mut acc = [0.0; LANES];
    let (xc, xr) = x.as_chunks::<LANES>();
    let (yc, yr) = y.as_chunks::<LANES>();
    for (a, b) in xc.iter().zip(yc) {
        for l in 0..LANES {
            acc[l] += g(a[l], b[l]);
        }
    }
    let mut s: f64 = acc.iter().sum();
    for (a, b) in xr.iter().zip(yr) {
        s += g(*a, *b);
    }
    s
}

// ---------------------------------------------------------------- dispatch
//
// The hot kernels have a const-generic body: `FMA = false` is portable
// (mul + add), `FMA = true` uses `mul_add` and is compiled only inside
// functions with `#[target_feature(enable = "avx2,fma")]`, selected at runtime
// (as OpenBLAS does). Phase 7. FMA rounds once instead of twice, which changes
// results at rounding level only.

/// `a·b + c`.
#[inline(always)]
fn madd<const FMA: bool>(a: f64, b: f64, c: f64) -> f64 {
    if FMA { a.mul_add(b, c) } else { a * b + c }
}

/// True if the CPU supports AVX2 and FMA (the result is cached by std).
#[inline]
fn use_fma() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("fma")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// Generates a runtime-dispatched public kernel `$name` from the generic body
/// `$body::<FMA>`.
macro_rules! dispatch {
    ($(#[$doc:meta])* fn $name:ident / $fma:ident = $body:ident ($($arg:ident : $ty:ty),*) $(-> $ret:ty)?) => {
        $(#[$doc])*
        pub(crate) fn $name($($arg: $ty),*) $(-> $ret)? {
            #[cfg(target_arch = "x86_64")]
            if use_fma() {
                // SAFETY: the CPU supports avx2 and fma (checked above).
                return unsafe { $fma($($arg),*) };
            }
            $body::<false>($($arg),*)
        }
        #[cfg(target_arch = "x86_64")]
        #[target_feature(enable = "avx2,fma")]
        fn $fma($($arg: $ty),*) $(-> $ret)? {
            $body::<true>($($arg),*)
        }
    };
}

#[inline(always)]
fn dot_k<const FMA: bool>(x: &[f64], y: &[f64]) -> f64 {
    debug_assert_eq!(x.len(), y.len());
    let mut acc = [0.0; LANES];
    let (xc, xr) = x.as_chunks::<LANES>();
    let (yc, yr) = y.as_chunks::<LANES>();
    for (a, b) in xc.iter().zip(yc) {
        for l in 0..LANES {
            acc[l] = madd::<FMA>(a[l], b[l], acc[l]);
        }
    }
    let mut s: f64 = acc.iter().sum();
    for (a, b) in xr.iter().zip(yr) {
        s = madd::<FMA>(*a, *b, s);
    }
    s
}

#[inline(always)]
fn axpy_k<const FMA: bool>(a: f64, x: &[f64], y: &mut [f64]) {
    debug_assert_eq!(x.len(), y.len());
    for (yi, xi) in y.iter_mut().zip(x) {
        *yi = madd::<FMA>(a, *xi, *yi);
    }
}

dispatch! {
    /// Dot product `xᵀy`.
    fn dot / dot_fma = dot_k(x: &[f64], y: &[f64]) -> f64
}

dispatch! {
    /// `y += a·x`.
    fn axpy / axpy_fma = axpy_k(a: f64, x: &[f64], y: &mut [f64])
}

/// `x *= a`.
pub(crate) fn scale(a: f64, x: &mut [f64]) {
    for xi in x.iter_mut() {
        *xi *= a;
    }
}

/// Euclidean norm. Fast path: plain sum of squares. Falls back to the scaled
/// computation if that overflows or underflows (extreme entries).
pub(crate) fn nrm2(x: &[f64]) -> f64 {
    let s = reduce2(x, x, |a, _| a * a);
    // safe range: no overflow, and not so small that squares underflowed
    if s.is_finite() && s > 1e-280 {
        return s.sqrt();
    }
    // squares are non-negative: a NaN sum means NaN entries (the scaled path
    // below would ignore them, since f64::max skips NaN)
    if s.is_nan() {
        return f64::NAN;
    }
    if s == 0.0 && x.iter().all(|&v| v == 0.0) {
        return 0.0;
    }
    let amax = x.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    if amax == 0.0 || !amax.is_finite() {
        return amax;
    }
    let s: f64 = x.iter().map(|v| (v / amax) * (v / amax)).sum();
    amax * s.sqrt()
}

/// Weighted root-mean-square norm `sqrt(1/n Σ (x_i·w_i)²)`, where
/// `w_i = 1/(atol_i + rtol_i·|u_i|)` are precomputed inverse scales.
pub(crate) fn wrms(x: &[f64], w: &[f64]) -> f64 {
    debug_assert_eq!(x.len(), w.len());
    if x.is_empty() {
        return 0.0;
    }
    let s = reduce2(x, w, |a, b| (a * b) * (a * b));
    (s / x.len() as f64).sqrt()
}

/// Row-block length of the blocked Gram–Schmidt kernels: a block of `v`
/// (4 KB) stays in L1 while all basis columns pass over it, so `v` is read
/// once per sweep instead of once per column (phase 7).
const BLOCK: usize = 512;

#[inline(always)]
fn qt_v_k<const FMA: bool>(q: &[f64], n: usize, k: usize, v: &[f64], c: &mut [f64]) -> f64 {
    debug_assert!(q.len() >= n * k && v.len() == n && c.len() >= k);
    c[..k].fill(0.0);
    let mut vv = 0.0;
    let mut start = 0;
    while start < n {
        let end = (start + BLOCK).min(n);
        let vb = &v[start..end];
        for (j, cj) in c[..k].iter_mut().enumerate() {
            *cj += dot_k::<FMA>(&q[j * n + start..j * n + end], vb);
        }
        vv += dot_k::<FMA>(vb, vb);
        start = end;
    }
    vv
}

#[inline(always)]
fn v_minus_qc_k<const FMA: bool>(
    q: &[f64],
    n: usize,
    k: usize,
    c: &[f64],
    v: &mut [f64],
    mut next: Option<&mut [f64]>,
) -> f64 {
    debug_assert!(q.len() >= n * k && v.len() == n && c.len() >= k);
    if let Some(nx) = next.as_deref_mut() {
        nx[..k].fill(0.0);
    }
    let mut vv = 0.0;
    let mut start = 0;
    while start < n {
        let end = (start + BLOCK).min(n);
        let vb = &mut v[start..end];
        for (j, &cj) in c[..k].iter().enumerate() {
            axpy_k::<FMA>(-cj, &q[j * n + start..j * n + end], vb);
        }
        if let Some(nx) = next.as_deref_mut() {
            for (j, nj) in nx[..k].iter_mut().enumerate() {
                *nj += dot_k::<FMA>(&q[j * n + start..j * n + end], vb);
            }
        }
        vv += dot_k::<FMA>(vb, vb);
        start = end;
    }
    vv
}

#[inline(always)]
fn q_times_k<const FMA: bool>(q: &[f64], n: usize, k: usize, l: &[f64], out: &mut [f64]) {
    debug_assert!(q.len() >= n * k && out.len() == n && l.len() >= k);
    out.fill(0.0);
    let mut start = 0;
    while start < n {
        let end = (start + BLOCK).min(n);
        let ob = &mut out[start..end];
        for (j, &lj) in l[..k].iter().enumerate() {
            axpy_k::<FMA>(lj, &q[j * n + start..j * n + end], ob);
        }
        start = end;
    }
}

dispatch! {
    /// `c[j] = q_jᵀ v` for the `k` columns of `q` (column-major, `n × k`), in
    /// one blocked sweep. Returns `vᵀv` (fused, for the norm).
    fn qt_v / qt_v_fma = qt_v_k(q: &[f64], n: usize, k: usize, v: &[f64], c: &mut [f64]) -> f64
}

dispatch! {
    /// `v ← v − Σ_j c_j q_j` in one blocked sweep. If `next` is given, the
    /// coefficients of the following Gram–Schmidt pass, `next[j] = q_jᵀ v_new`,
    /// are accumulated in the same sweep (each block of `v` is final once all
    /// columns have been subtracted). Returns `v_newᵀ v_new`.
    fn v_minus_qc / v_minus_qc_fma = v_minus_qc_k(q: &[f64], n: usize, k: usize, c: &[f64], v: &mut [f64], next: Option<&mut [f64]>) -> f64
}

dispatch! {
    /// `out = Σ_j l_j q_j` (blocked).
    fn q_times / q_times_fma = q_times_k(q: &[f64], n: usize, k: usize, l: &[f64], out: &mut [f64])
}

/// True if all entries are finite.
pub(crate) fn all_finite(x: &[f64]) -> bool {
    x.iter().all(|v| v.is_finite())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kernels() {
        let x = [1.0, -2.0, 3.0];
        let mut y = [4.0, 5.0, 6.0];
        assert_eq!(dot(&x, &y), 4.0 - 10.0 + 18.0);
        axpy(2.0, &x, &mut y);
        assert_eq!(y, [6.0, 1.0, 12.0]);
        scale(0.5, &mut y);
        assert_eq!(y, [3.0, 0.5, 6.0]);
        assert!((nrm2(&[3.0, 4.0]) - 5.0).abs() < 1e-15);
        assert_eq!(nrm2(&[0.0, 0.0]), 0.0);
        assert!((wrms(&[1.0, 2.0], &[2.0, 1.0]) - 2.0).abs() < 1e-15);
    }

    #[test]
    fn nrm2_extreme_scales() {
        assert!((nrm2(&[3e200, 4e200]) / 5e200 - 1.0).abs() < 1e-15);
        assert!((nrm2(&[3e-200, 4e-200]) / 5e-200 - 1.0).abs() < 1e-15);
        assert!(nrm2(&[1.0, f64::NAN]).is_nan());
        assert!(
            nrm2(&[f64::NAN, f64::NAN]).is_nan(),
            "all-NaN must not look like zero"
        );
        assert!(!all_finite(&[1.0, f64::NAN]));
    }
}
