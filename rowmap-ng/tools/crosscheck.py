# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "scipy"]
# ///
"""Cross-check the Python problem definitions against the Fortran ones.

Uses `driver dump` output (y0, f, J*v with v_k = sin k, f_t) and also checks
the sparse Jacobian (used by Radau) against jv. Usage:
    uv run tools/crosscheck.py [problem ...]     (default: all, small and paper grids)
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import problems as P  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
DRIVER = ROOT / "fortran/experiments/build/driver"


def rel(a, b, scale=0.0):
    """Max abs difference relative to max|b|, or to `scale` if larger.

    `scale` accounts for cancellation (f is small compared to the O(1/h^2)
    terms that cancel in it)."""
    return np.max(np.abs(a - b)) / max(np.max(np.abs(b)), scale, 1e-300)


def check(name, m):
    with tempfile.TemporaryDirectory() as d:
        out = Path(d) / "dump.txt"
        subprocess.run([DRIVER, "dump", name, str(m), out], check=True)
        data = np.loadtxt(out, comments="#")
    pr = P.make(name, m)
    v = np.sin(np.arange(1, pr.n + 1, dtype=float))
    y0 = pr.y0
    fscale = np.max(abs(pr.jac(pr.t0, y0)) @ np.abs(y0))
    res = {
        "y0": rel(y0, data[:, 0]),
        "f": rel(pr.f(pr.t0, y0), data[:, 1], fscale),
        "jv": rel(pr.jv(pr.t0, y0, v), data[:, 2]),
        "ft": rel(pr.ft(pr.t0, y0), data[:, 3]) if np.any(data[:, 3]) else 0.0,
        "jac": rel(pr.jac(pr.t0, y0) @ v, pr.jv(pr.t0, y0, v)),
    }
    ok = all(r < 1e-12 for r in res.values())
    print(f"{name:7s} m={m:3d} n={pr.n:5d}  " + "  ".join(f"{k}={r:.1e}" for k, r in res.items())
          + ("  OK" if ok else "  MISMATCH"))
    return ok


def main():
    names = sys.argv[1:] or list(P.PAPER_M)
    ok = True
    for name in names:
        for m in sorted({5 if name != "food3" else 4, P.PAPER_M[name]}):
            ok &= check(name, m)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
