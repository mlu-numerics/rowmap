//! The integrator: step-size control, initial step, dense output, statistics
//! (AGENTS.md §2, §3.7).

use std::fmt;
use std::ops::ControlFlow;

use crate::grk4t::{GRK4T, MAXS, ROS34PW2, RowMethod};
use crate::krylov::Basis;
use crate::linalg::{all_finite, axpy, wrms};
use crate::map::{
    GRK4T_AUTONOMOUS, GRK4T_NONAUTONOMOUS, ROS34PW2_AUTONOMOUS, ROS34PW2_NONAUTONOMOUS, StageKind,
    StepConfig, StepFail, THETA, TestHook, step,
};
use crate::system::OdeSystem;

/// A tolerance: one value for all components, or one per component.
#[derive(Debug, Clone, PartialEq)]
pub enum Tol {
    Scalar(f64),
    Vector(Vec<f64>),
}

impl Tol {
    fn get(&self, i: usize) -> f64 {
        match self {
            Tol::Scalar(v) => *v,
            Tol::Vector(v) => v[i],
        }
    }
    fn check(&self, n: usize, name: &str) -> Result<(), RowmapError> {
        let ok = |v: f64| v.is_finite() && v >= 0.0;
        match self {
            Tol::Scalar(v) if ok(*v) => Ok(()),
            Tol::Vector(v) if v.len() == n && v.iter().all(|&x| ok(x)) => Ok(()),
            Tol::Vector(v) if v.len() != n => Err(RowmapError::InvalidInput(format!(
                "{name} has length {}, expected {n}",
                v.len()
            ))),
            _ => Err(RowmapError::InvalidInput(format!(
                "{name} must be finite and >= 0"
            ))),
        }
    }
}

impl From<f64> for Tol {
    fn from(v: f64) -> Self {
        Tol::Scalar(v)
    }
}

impl From<Vec<f64>> for Tol {
    fn from(v: Vec<f64>) -> Self {
        Tol::Vector(v)
    }
}

impl From<&[f64]> for Tol {
    fn from(v: &[f64]) -> Self {
        Tol::Vector(v.to_vec())
    }
}

/// Errors. On error, the state vector holds the last accepted solution at
/// the time given in the error.
#[derive(Debug, Clone, PartialEq)]
pub enum RowmapError {
    /// Inconsistent input (dimensions, tolerances, time interval).
    InvalidInput(String),
    /// The problem is not autonomous but does not implement `dfdt`.
    MissingDfdt,
    /// The step size fell below the resolution of `t`.
    StepSizeTooSmall { t: f64, h: f64 },
    /// More than `max_steps` steps (accepted + rejected).
    MaxStepsExceeded { t: f64 },
    /// The projected stage matrix was singular even after 3 step halvings.
    SingularStageMatrix { t: f64 },
    /// f or J·v returned NaN/Inf even after 3 step halvings.
    NonFinite { t: f64 },
}

impl fmt::Display for RowmapError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RowmapError::InvalidInput(s) => write!(f, "invalid input: {s}"),
            RowmapError::MissingDfdt => {
                write!(
                    f,
                    "non-autonomous problem without dfdt (implement dfdt or declare autonomous)"
                )
            }
            RowmapError::StepSizeTooSmall { t, h } => {
                write!(f, "step size too small (h = {h:e}) at t = {t}")
            }
            RowmapError::MaxStepsExceeded { t } => {
                write!(f, "maximum number of steps exceeded at t = {t}")
            }
            RowmapError::SingularStageMatrix { t } => {
                write!(f, "stage matrix repeatedly singular at t = {t}")
            }
            RowmapError::NonFinite { t } => write!(f, "f or J·v not finite at t = {t}"),
        }
    }
}

impl std::error::Error for RowmapError {}

/// Min/max/mean of the Krylov dimension at the end of one stage.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct DimStats {
    pub min: usize,
    pub max: usize,
    pub mean: f64,
}

/// Integration statistics.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Stats {
    /// Accepted + rejected steps (as Fortran `iwork(5)`).
    pub steps: usize,
    pub accepted: usize,
    pub rejected: usize,
    /// f evaluations (including the initial value and the initial step guess).
    pub f_evals: usize,
    /// Jacobian-vector products.
    pub jv: usize,
    /// ∂f/∂t evaluations.
    pub dfdt_evals: usize,
    /// Krylov dimensions per stage, over all steps.
    pub krylov: [DimStats; MAXS],
    /// Inserted stage vectors found dependent (e.g. `w₄ ∈ K₃`).
    pub dropped_insert: usize,
    /// Remainders found dependent (Krylov breakdown).
    pub dropped_promote: usize,
}

/// How the integration ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    Finished,
    /// Stopped by the callback.
    Interrupted,
}

/// Result of a successful integration.
#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    /// Final time (`t_end` if finished).
    pub t: f64,
    /// Predicted next step size.
    pub h_next: f64,
    pub status: Status,
    pub stats: Stats,
}

/// An accepted step, passed to the callback of [`Rowmap::solve_with`].
#[derive(Debug, Clone, Copy)]
pub struct Step<'a> {
    pub t_old: f64,
    pub t_new: f64,
    pub y_old: &'a [f64],
    pub y_new: &'a [f64],
    pub f_old: &'a [f64],
    pub f_new: &'a [f64],
}

impl Step<'_> {
    /// Cubic Hermite interpolation (order 3) at `s ∈ [t_old, t_new]`
    /// (rowmap.F `rowcon`).
    pub fn interpolate(&self, s: f64, out: &mut [f64]) {
        let h = self.t_new - self.t_old;
        let th = (s - self.t_old) / h;
        let th2 = th * th;
        let thm1 = th - 1.0;
        let g0 = th2 * (3.0 - 2.0 * th);
        let g1 = h * th * thm1 * thm1;
        let g2 = h * th2 * thm1;
        for i in 0..out.len() {
            out[i] = (1.0 - g0) * self.y_old[i]
                + g0 * self.y_new[i]
                + g1 * self.f_old[i]
                + g2 * self.f_new[i];
        }
    }
}

/// Step-size controller constants (paper §3, rowmap.F defaults).
const FAC_MIN: f64 = 0.25;
const FAC_MAX: f64 = 2.0;
const SAFETY: f64 = 0.8;
/// Consecutive step halvings after a singular stage matrix or a non-finite
/// value before giving up (rowmap.F: 3).
const MAX_FAILS: usize = 3;

/// The underlying ROW method.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Method {
    /// GRK4T (Kaps–Rentrop), order 4(3): the ROWMAP method of the paper.
    /// Minimal Krylov dimensions 4/7/8/11 (autonomous), 4/9/10/15
    /// (non-autonomous).
    #[default]
    Grk4t,
    /// ROS34PW2 (Rang–Angermann 2005), order 3(2), stiffly accurate,
    /// L-stable, 4 stages. Often cheaper at crude tolerances. Minimal Krylov
    /// dimensions 3/4/4/4 (autonomous), 3/4/8/8 (non-autonomous);
    /// docs/map.md §8.
    Ros34pw2,
}

/// The ROWMAP integrator. Create with [`Rowmap::new`]; the other settings
/// have defaults that fit most problems.
#[derive(Debug, Clone)]
pub struct Rowmap {
    rtol: Tol,
    atol: Tol,
    method: Method,
    max_krylov: usize,
    max_steps: usize,
    initial_step: Option<f64>,
    max_step: Option<f64>,
    test_sched: Option<Vec<StageKind>>,
    reorth: bool,
}

impl Rowmap {
    /// Integrator with relative and absolute tolerances (scalars or vectors).
    /// The error is measured in the weighted RMS norm with weights
    /// `atol_i + rtol_i |y_i|`.
    pub fn new(rtol: impl Into<Tol>, atol: impl Into<Tol>) -> Self {
        Rowmap {
            rtol: rtol.into(),
            atol: atol.into(),
            method: Method::Grk4t,
            max_krylov: 70,
            max_steps: 100_000,
            initial_step: None,
            max_step: None,
            test_sched: None,
            reorth: crate::krylov::REORTH_DEFAULT,
        }
    }

    /// Internal: switch the second Gram–Schmidt pass on or off (default off,
    /// `krylov::REORTH_DEFAULT`); for comparisons.
    #[doc(hidden)]
    pub fn __reorthogonalize(&mut self, on: bool) -> &mut Self {
        self.reorth = on;
        self
    }

    /// The ROW method. Default [`Method::Grk4t`].
    pub fn method(&mut self, m: Method) -> &mut Self {
        self.method = m;
        self
    }

    /// Maximal Krylov dimension (memory: n × max_krylov). Default 70.
    pub fn max_krylov(&mut self, m: usize) -> &mut Self {
        self.max_krylov = m;
        self
    }

    /// Maximal number of steps (accepted + rejected). Default 100 000.
    pub fn max_steps(&mut self, n: usize) -> &mut Self {
        self.max_steps = n;
        self
    }

    /// Initial step size. Default: estimated (Hairer–Wanner).
    pub fn initial_step(&mut self, h: f64) -> &mut Self {
        self.initial_step = Some(h);
        self
    }

    /// Maximal step size. Default: unlimited.
    pub fn max_step(&mut self, h: f64) -> &mut Self {
        self.max_step = Some(h);
        self
    }

    /// Internal: replace the MAP schedule (schedule comparisons, phase 4).
    #[doc(hidden)]
    pub fn __schedule_without_stage4_closure(&mut self) -> &mut Self {
        let mut s = GRK4T_AUTONOMOUS.to_vec();
        s[3] = StageKind::Extend {
            insert: true,
            rounds: 0,
        };
        self.test_sched = Some(s);
        self
    }

    /// Integrate from `t0` to `t_end`; `y` holds the initial value and is
    /// overwritten with the solution. Only forward integration is supported
    /// (`t_end >= t0`).
    pub fn solve<S: OdeSystem + ?Sized>(
        &self,
        sys: &S,
        t0: f64,
        y: &mut [f64],
        t_end: f64,
    ) -> Result<Solution, RowmapError> {
        self.solve_with(sys, t0, y, t_end, |_| ControlFlow::Continue(()))
    }

    /// Like [`solve`](Self::solve), calling `callback` after every accepted
    /// step; returning `ControlFlow::Break(())` stops the integration.
    pub fn solve_with<S, F>(
        &self,
        sys: &S,
        t0: f64,
        y: &mut [f64],
        t_end: f64,
        mut callback: F,
    ) -> Result<Solution, RowmapError>
    where
        S: OdeSystem + ?Sized,
        F: FnMut(&Step) -> ControlFlow<()>,
    {
        let n = sys.dim();
        if y.len() != n {
            return Err(RowmapError::InvalidInput(format!(
                "y has length {}, system dimension {n}",
                y.len()
            )));
        }
        self.rtol.check(n, "rtol")?;
        self.atol.check(n, "atol")?;
        if (0..n).any(|i| self.rtol.get(i) == 0.0 && self.atol.get(i) == 0.0) {
            return Err(RowmapError::InvalidInput(
                "rtol and atol must not both be zero".into(),
            ));
        }
        if !(t0.is_finite() && t_end.is_finite()) || t_end < t0 {
            return Err(RowmapError::InvalidInput("need finite t0 <= t_end".into()));
        }
        if self.max_krylov < 12 {
            return Err(RowmapError::InvalidInput("max_krylov must be >= 12".into()));
        }
        if !all_finite(y) {
            return Err(RowmapError::InvalidInput("initial value not finite".into()));
        }
        let autonomous = sys.autonomous();
        let mut ft = vec![0.0; n];
        if !autonomous && sys.dfdt(t0, y, &mut ft).is_err() {
            return Err(RowmapError::MissingDfdt);
        }
        let method: &RowMethod = match self.method {
            Method::Grk4t => &GRK4T,
            Method::Ros34pw2 => &ROS34PW2,
        };
        let sched: &[StageKind] = match (self.method, &self.test_sched, autonomous) {
            (Method::Grk4t, Some(s), true) => s,
            (Method::Grk4t, _, true) => &GRK4T_AUTONOMOUS,
            (Method::Grk4t, _, false) => &GRK4T_NONAUTONOMOUS,
            (Method::Ros34pw2, _, true) => &ROS34PW2_AUTONOMOUS,
            (Method::Ros34pw2, _, false) => &ROS34PW2_NONAUTONOMOUS,
        };
        let cfg = StepConfig {
            max_krylov: self.max_krylov,
            theta: THETA,
            hook: TestHook::default(),
        };

        let mut stats = Stats::default();
        for d in stats.krylov.iter_mut() {
            d.min = usize::MAX;
        }
        let mut dim_sum = [0usize; MAXS];
        let mut dim_count = 0usize;
        let mut t = t0;
        let mut f = vec![0.0; n];
        sys.rhs(t, y, &mut f);
        stats.f_evals += 1;
        if !all_finite(&f) {
            return Err(RowmapError::NonFinite { t });
        }
        if !autonomous {
            stats.dfdt_evals += 1;
            if !all_finite(&ft) {
                return Err(RowmapError::NonFinite { t });
            }
        }
        if t_end == t0 {
            return Ok(Solution {
                t,
                h_next: 0.0,
                status: Status::Finished,
                stats,
            });
        }
        let weights = |y: &[f64]| -> Vec<f64> {
            (0..n)
                .map(|i| 1.0 / (self.atol.get(i) + self.rtol.get(i) * y[i].abs()))
                .collect()
        };
        let hmax = self.max_step.unwrap_or(f64::INFINITY);
        let mut h = match self.initial_step {
            Some(h0) if h0 > 0.0 => h0,
            _ => initial_step(sys, t, y, &f, t_end, &weights(y), method.order, &mut stats),
        };
        h = h.min(hmax);
        let mut basis = Basis::new(n);
        basis.set_reorth(self.reorth);
        let mut fails = 0usize;
        let mut last_rejected = false;
        let mut f_new = vec![0.0; n];

        loop {
            if t >= t_end {
                break;
            }
            if stats.steps >= self.max_steps {
                return Err(RowmapError::MaxStepsExceeded { t });
            }
            let mut last = false;
            if t + 1.005 * h >= t_end {
                h = t_end - t;
                last = true;
            }
            let w = weights(y);
            stats.steps += 1;
            // f_t is evaluated once per accepted point (t, y) and reused after rejections
            let ft_arg = if autonomous {
                None
            } else {
                Some(ft.as_slice())
            };
            let res = {
                let _g = crate::prof::section("step (total)");
                step(
                    sys, method, sched, t, y, &f, ft_arg, h, &w, &cfg, &mut basis,
                )
            };
            stats.jv = basis.stats.jv;
            stats.dropped_insert = basis.stats.dropped_insert;
            stats.dropped_promote = basis.stats.dropped_promote;
            let out = match res {
                Ok(out) => out,
                Err(e) => {
                    fails += 1;
                    stats.rejected += 1;
                    if fails >= MAX_FAILS {
                        return Err(match e {
                            StepFail::Singular => RowmapError::SingularStageMatrix { t },
                            StepFail::NonFinite => RowmapError::NonFinite { t },
                        });
                    }
                    h *= 0.5;
                    last_rejected = true;
                    continue;
                }
            };
            stats.f_evals += out.nfev;
            for i in 0..MAXS {
                let d = out.dims[i];
                stats.krylov[i].min = stats.krylov[i].min.min(d);
                stats.krylov[i].max = stats.krylov[i].max.max(d);
                dim_sum[i] += d;
            }
            dim_count += 1;
            let err = wrms(&out.err, &w);
            let fac = if err == 0.0 {
                FAC_MAX
            } else {
                // exponent 1/(q+1), q = order of the embedded method
                // (GRK4T: 1/4 as in rowmap.F; ROS34PW2: 1/3)
                (SAFETY * err.powf(-1.0 / (method.order_embedded as f64 + 1.0)))
                    .clamp(FAC_MIN, FAC_MAX)
            };
            let mut h_new = h * fac;
            if 0.1 * h_new.abs() <= f64::EPSILON * t.abs() || !h_new.is_finite() {
                return Err(RowmapError::StepSizeTooSmall { t, h: h_new });
            }
            if err < 1.0 {
                let t_new = if last { t_end } else { t + h };
                sys.rhs(t_new, &out.u_new, &mut f_new);
                stats.f_evals += 1;
                if !all_finite(&f_new) {
                    fails += 1;
                    stats.rejected += 1;
                    if fails >= MAX_FAILS {
                        return Err(RowmapError::NonFinite { t });
                    }
                    h *= 0.5;
                    last_rejected = true;
                    continue;
                }
                stats.accepted += 1;
                fails = 0;
                let flow = callback(&Step {
                    t_old: t,
                    t_new,
                    y_old: y,
                    y_new: &out.u_new,
                    f_old: &f,
                    f_new: &f_new,
                });
                y.copy_from_slice(&out.u_new);
                std::mem::swap(&mut f, &mut f_new);
                t = t_new;
                if !autonomous && t < t_end && !flow.is_break() {
                    // dfdt succeeded at t0, so it is implemented
                    let _ = sys.dfdt(t, y, &mut ft);
                    stats.dfdt_evals += 1;
                    if !all_finite(&ft) {
                        return Err(RowmapError::NonFinite { t });
                    }
                }
                if last_rejected {
                    h_new = h_new.min(h);
                }
                last_rejected = false;
                h = h_new.min(hmax);
                if flow.is_break() {
                    finish_dims(&mut stats, &dim_sum, dim_count);
                    return Ok(Solution {
                        t,
                        h_next: h,
                        status: Status::Interrupted,
                        stats,
                    });
                }
            } else {
                stats.rejected += 1;
                last_rejected = true;
                h = h_new;
            }
        }
        finish_dims(&mut stats, &dim_sum, dim_count);
        Ok(Solution {
            t,
            h_next: h,
            status: Status::Finished,
            stats,
        })
    }
}

/// Finalize the Krylov dimension statistics; `count` = number of steps whose
/// stages all completed (accepted or rejected by the error test).
fn finish_dims(stats: &mut Stats, sum: &[usize; MAXS], count: usize) {
    for i in 0..MAXS {
        if count == 0 {
            stats.krylov[i] = DimStats::default();
        } else {
            stats.krylov[i].mean = sum[i] as f64 / count as f64;
        }
    }
}

/// Initial step size (Hairer–Wanner I, II.4) for a method of order `p`.
#[allow(clippy::too_many_arguments)]
fn initial_step<S: OdeSystem + ?Sized>(
    sys: &S,
    t: f64,
    y: &[f64],
    f: &[f64],
    t_end: f64,
    w: &[f64],
    p: usize,
    stats: &mut Stats,
) -> f64 {
    let span = t_end - t;
    let d0 = wrms(y, w);
    let d1 = wrms(f, w);
    let mut h0 = if d0 < 1e-5 || d1 < 1e-5 {
        1e-6
    } else {
        0.01 * d0 / d1
    };
    h0 = h0.min(span);
    let mut y1 = y.to_vec();
    axpy(h0, f, &mut y1);
    let mut f1 = vec![0.0; y.len()];
    sys.rhs(t + h0, &y1, &mut f1);
    stats.f_evals += 1;
    if !all_finite(&f1) {
        return h0 * 1e-3;
    }
    axpy(-1.0, f, &mut f1);
    let d2 = wrms(&f1, w) / h0;
    let dm = d1.max(d2);
    let h1 = if dm <= 1e-15 {
        (h0 * 1e-3).max(1e-6)
    } else {
        (0.01 / dm).powf(1.0 / (p as f64 + 1.0))
    };
    (100.0 * h0).min(h1).min(span)
}
