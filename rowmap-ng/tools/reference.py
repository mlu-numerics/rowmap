# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "scipy"]
# ///
"""One-off reference solutions (docs/problems.md, AGENTS.md §5.2).

DOPRI5 (SciPy RK45), rtol = atol = 1e-10. FOOD3 is far too stiff for an
explicit method and uses Radau with the analytic sparse Jacobian.
Output: tests/data/reference/<problem>_m<m>.txt (header + one y_end value per line).

    uv run tools/reference.py <problem> [m]      (m defaults to the paper grid)
"""

import sys
import time
from datetime import date
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).parent))
import problems as P  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "tests/data/reference"
TOL = 1e-10


def main():
    name = sys.argv[1]
    m = int(sys.argv[2]) if len(sys.argv) > 2 else P.PAPER_M[name]
    pr = P.make(name, m)
    if name == "food3":
        method, kw = "Radau", {"jac": pr.jac}
    else:
        method, kw = "RK45", {}
    c0 = time.process_time()
    # t_eval=[tend]: do not keep every step in memory (VDPL needs ~1e6 steps)
    sol = solve_ivp(pr.f, (pr.t0, pr.tend), pr.y0, method=method, rtol=TOL, atol=TOL,
                    t_eval=[pr.tend], **kw)
    cpu = time.process_time() - c0
    if not sol.success:
        sys.exit(f"{name} m={m}: {sol.message}")
    y = sol.y[:, -1]
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{name}_m{m}.txt"
    solver = "Radau (SciPy)" if method == "Radau" else "DOPRI5 (SciPy RK45)"
    with open(path, "w") as fh:
        fh.write(f"# problem = {name}\n# m = {m}  n = {pr.n}\n")
        fh.write(f"# t0 = {pr.t0!r}  t_end = {pr.tend!r}\n")
        fh.write(f"# solver = {solver}  rtol = atol = {TOL:g}\n")
        fh.write(f"# nfev = {sol.nfev}  cpu = {cpu:.1f} s  date = {date.today()}\n")
        for v in y:
            fh.write(f"{v: .16e}\n")
    print(f"{name:7s} m={m:3d} n={pr.n:5d} {method}: nfev={sol.nfev} "
          f"cpu={cpu:.1f}s -> {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
