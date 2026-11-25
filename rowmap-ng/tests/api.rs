//! Public API behaviour: errors, dense output, callback, simple exact cases.

use std::ops::ControlFlow;

use rowmap::{Method, OdeSystem, Rowmap, RowmapError, Status};

/// y' = λ y (componentwise, different λ), exact solution exp(λ t) y0.
struct Decay {
    lam: Vec<f64>,
}

impl OdeSystem for Decay {
    fn dim(&self) -> usize {
        self.lam.len()
    }
    fn rhs(&self, _t: f64, y: &[f64], f: &mut [f64]) {
        for i in 0..y.len() {
            f[i] = self.lam[i] * y[i];
        }
    }
    fn jac_vec(&self, _t: f64, _y: &[f64], v: &[f64], out: &mut [f64]) {
        for i in 0..v.len() {
            out[i] = self.lam[i] * v[i];
        }
    }
    fn autonomous(&self) -> bool {
        true
    }
}

/// Non-autonomous without dfdt.
struct NoDfdt;

impl OdeSystem for NoDfdt {
    fn dim(&self) -> usize {
        1
    }
    fn rhs(&self, t: f64, _y: &[f64], f: &mut [f64]) {
        f[0] = t;
    }
    fn jac_vec(&self, _t: f64, _y: &[f64], _v: &[f64], out: &mut [f64]) {
        out[0] = 0.0;
    }
}

/// Non-autonomous, stiff: y' = −λ(y − sin t) + cos t, exact y = sin t + y0 e^{−λt}.
struct Forced {
    lam: Vec<f64>,
}

impl OdeSystem for Forced {
    fn dim(&self) -> usize {
        self.lam.len()
    }
    fn rhs(&self, t: f64, y: &[f64], f: &mut [f64]) {
        for i in 0..y.len() {
            f[i] = -self.lam[i] * (y[i] - t.sin()) + t.cos();
        }
    }
    fn jac_vec(&self, _t: f64, _y: &[f64], v: &[f64], out: &mut [f64]) {
        for i in 0..v.len() {
            out[i] = -self.lam[i] * v[i];
        }
    }
    fn dfdt(&self, t: f64, _y: &[f64], ft: &mut [f64]) -> Result<(), rowmap::Unsupported> {
        for (fi, l) in ft.iter_mut().zip(&self.lam) {
            *fi = l * t.cos() - t.sin();
        }
        Ok(())
    }
}

/// ROS34PW2 (public order-3 option): accuracy scales with the tolerance,
/// autonomous and non-autonomous, on the very stiff spectrum λ = −1 … −981
/// (L-stable). Stage-1 Krylov dimension ≥ 3 (= order).
#[test]
fn ros34pw2_method() {
    let (sys, y0) = decay();
    let forced = Forced {
        lam: sys.lam.iter().map(|l| -l).collect(),
    };
    let t_end = 1.0;
    for (label, nonaut) in [("decay", false), ("forced", true)] {
        let mut prev = f64::INFINITY;
        for tol in [1e-3, 1e-5, 1e-7] {
            let mut y = y0.clone();
            let mut solver = Rowmap::new(tol, tol);
            solver.method(Method::Ros34pw2);
            let sol = if nonaut {
                solver.solve(&forced, 0.0, &mut y, t_end).unwrap()
            } else {
                solver.solve(&sys, 0.0, &mut y, t_end).unwrap()
            };
            let err = y
                .iter()
                .zip(&sys.lam)
                .map(|(v, l)| {
                    // decay: y = e^{λt} (λ < 0); forced: y = sin t + e^{λt}
                    let exact = if nonaut {
                        t_end.sin() + (l * t_end).exp()
                    } else {
                        (l * t_end).exp()
                    };
                    (v - exact).abs()
                })
                .fold(0.0, f64::max);
            assert!(err < 50.0 * tol, "{label} tol {tol}: err {err:e}");
            assert!(err < prev, "{label}: no improvement at tol {tol}");
            prev = err;
            assert!(sol.stats.krylov[0].min >= 3);
        }
    }
}

/// Diagnostic: error vs tolerance for both public methods.
#[test]
#[ignore]
fn methods_diagnose() {
    let (sys, _) = decay();
    let lam = sys.lam.clone();
    for m in [Method::Grk4t, Method::Ros34pw2] {
        for tol in [1e-3, 1e-4, 1e-5, 1e-6, 1e-7] {
            let mut y = vec![1.0; 50];
            let sol = Rowmap::new(tol, tol)
                .method(m)
                .solve(&sys, 0.0, &mut y, 1.0)
                .unwrap();
            let (imax, err) = y
                .iter()
                .zip(&lam)
                .map(|(v, l)| (v - l.exp()).abs())
                .enumerate()
                .fold((0, 0.0), |a, (i, e)| if e > a.1 { (i, e) } else { a });
            eprintln!(
                "{m:?} tol {tol:.0e}: err {err:.2e} (comp {imax}), steps {} rej {}, k1 [{}, {}]",
                sol.stats.steps,
                sol.stats.rejected,
                sol.stats.krylov[0].min,
                sol.stats.krylov[0].max
            );
        }
    }
}

#[test]
fn nonautonomous_forced_problem() {
    let lam: Vec<f64> = (0..40).map(|i| 1.0 + 25.0 * i as f64).collect();
    let sys = Forced { lam: lam.clone() };
    let t_end = 2.0;
    let mut prev = f64::INFINITY;
    for tol in [1e-3, 1e-5, 1e-7] {
        let mut y = vec![1.0; 40];
        let sol = Rowmap::new(tol, tol)
            .solve(&sys, 0.0, &mut y, t_end)
            .unwrap();
        let err = y
            .iter()
            .zip(&lam)
            .map(|(v, l)| (v - (t_end.sin() + (-l * t_end).exp())).abs())
            .fold(0.0, f64::max);
        assert!(err < 50.0 * tol, "tol {tol}: err {err:e}");
        assert!(err < prev);
        prev = err;
        assert_eq!(sol.stats.dfdt_evals, sol.stats.accepted);
        // non-autonomous minimal dimensions: stage 2 adds up to 5, stage 4 up to 5
        assert!(sol.stats.krylov[0].min >= 4);
    }
}

fn decay() -> (Decay, Vec<f64>) {
    let lam: Vec<f64> = (0..50).map(|i| -(1.0 + i as f64 * 20.0)).collect();
    (Decay { lam }, vec![1.0; 50])
}

#[test]
fn stiff_linear_decay_accuracy_scales_with_tolerance() {
    let (sys, y0) = decay();
    let mut prev = f64::INFINITY;
    for tol in [1e-3, 1e-5, 1e-7] {
        let mut y = y0.clone();
        let sol = Rowmap::new(tol, tol).solve(&sys, 0.0, &mut y, 1.0).unwrap();
        assert_eq!(sol.status, Status::Finished);
        assert_eq!(sol.t, 1.0);
        let err = y
            .iter()
            .zip(&sys.lam)
            .map(|(a, l)| (a - l.exp()).abs())
            .fold(0.0, f64::max);
        assert!(err < 50.0 * tol, "tol {tol}: err {err:e}");
        assert!(err < prev);
        prev = err;
        assert!(sol.stats.accepted > 0 && sol.stats.steps >= sol.stats.accepted);
        assert!(sol.stats.krylov[0].min >= 4);
    }
}

#[test]
fn dense_output_and_interrupt() {
    let (sys, y0) = decay();
    let mut y = y0.clone();
    let mut outputs = Vec::new();
    let mut next = 0.1;
    let sol = Rowmap::new(1e-8, 1e-8)
        .solve_with(&sys, 0.0, &mut y, 1.0, |s| {
            let mut buf = vec![0.0; s.y_new.len()];
            while next <= s.t_new && next <= 0.5 + 1e-12 {
                s.interpolate(next, &mut buf);
                outputs.push((next, buf[0]));
                next += 0.1;
            }
            if s.t_new >= 0.5 {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        })
        .unwrap();
    assert_eq!(sol.status, Status::Interrupted);
    assert!(sol.t >= 0.5 && sol.t < 1.0);
    assert_eq!(outputs.len(), 5);
    for (t, v) in outputs {
        // component 0: λ = −1; Hermite interpolation is order 3
        assert!(
            (v - (-t).exp()).abs() < 1e-5,
            "t = {t}: {v} vs {}",
            (-t).exp()
        );
    }
}

#[test]
fn input_errors() {
    let (sys, _) = decay();
    let r = Rowmap::new(1e-6, 1e-6);
    let mut short = vec![1.0; 3];
    assert!(matches!(
        r.solve(&sys, 0.0, &mut short, 1.0),
        Err(RowmapError::InvalidInput(_))
    ));
    let mut y = vec![1.0; 50];
    assert!(matches!(
        r.solve(&sys, 1.0, &mut y, 0.0),
        Err(RowmapError::InvalidInput(_))
    ));
    let bad = Rowmap::new(vec![1e-6; 7], 1e-6);
    assert!(matches!(
        bad.solve(&sys, 0.0, &mut y, 1.0),
        Err(RowmapError::InvalidInput(_))
    ));
    let zero = Rowmap::new(0.0, 0.0);
    assert!(matches!(
        zero.solve(&sys, 0.0, &mut y, 1.0),
        Err(RowmapError::InvalidInput(_))
    ));
    let mut y1 = vec![0.0];
    assert_eq!(
        r.solve(&NoDfdt, 0.0, &mut y1, 1.0),
        Err(RowmapError::MissingDfdt)
    );
}

#[test]
fn max_steps_and_zero_interval() {
    let (sys, y0) = decay();
    let mut y = y0.clone();
    let err = Rowmap::new(1e-10, 1e-10)
        .max_steps(3)
        .solve(&sys, 0.0, &mut y, 1.0);
    assert!(matches!(err, Err(RowmapError::MaxStepsExceeded { .. })));
    let mut y = y0.clone();
    let sol = Rowmap::new(1e-6, 1e-6)
        .solve(&sys, 0.5, &mut y, 0.5)
        .unwrap();
    assert_eq!((sol.t, sol.stats.steps), (0.5, 0));
    assert_eq!(y, y0);
}

#[test]
fn max_step_is_respected() {
    let (sys, y0) = decay();
    let mut y = y0.clone();
    let mut hmax_seen: f64 = 0.0;
    Rowmap::new(1e-3, 1e-3)
        .max_step(0.05)
        .solve_with(&sys, 0.0, &mut y, 1.0, |s| {
            hmax_seen = hmax_seen.max(s.t_new - s.t_old);
            ControlFlow::Continue(())
        })
        .unwrap();
    assert!(hmax_seen <= 0.05 * 1.005 + 1e-15, "{hmax_seen}");
}
