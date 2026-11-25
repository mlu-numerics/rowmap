//! Internal timers for performance work (cargo feature `profile`). Without
//! the feature, `section` is a no-op that the compiler removes.

#[cfg(feature = "profile")]
mod imp {
    use std::cell::RefCell;
    use std::time::Instant;

    thread_local! {
        static T: RefCell<Vec<(&'static str, f64, u64)>> = const { RefCell::new(Vec::new()) };
    }

    pub struct Guard(&'static str, Instant);

    impl Drop for Guard {
        fn drop(&mut self) {
            let dt = self.1.elapsed().as_secs_f64();
            T.with(|t| {
                let mut t = t.borrow_mut();
                match t.iter_mut().find(|e| e.0 == self.0) {
                    Some(e) => {
                        e.1 += dt;
                        e.2 += 1;
                    }
                    None => t.push((self.0, dt, 1)),
                }
            });
        }
    }

    pub fn section(name: &'static str) -> Guard {
        Guard(name, Instant::now())
    }

    /// Take and reset the accumulated timers: (name, seconds, calls).
    pub fn take() -> Vec<(&'static str, f64, u64)> {
        T.with(|t| std::mem::take(&mut *t.borrow_mut()))
    }
}

#[cfg(feature = "profile")]
pub use imp::{section, take};

/// Zero-sized stand-in for the timer guard.
#[cfg(not(feature = "profile"))]
pub struct Guard;

#[cfg(not(feature = "profile"))]
#[inline(always)]
pub fn section(_name: &'static str) -> Guard {
    Guard
}
