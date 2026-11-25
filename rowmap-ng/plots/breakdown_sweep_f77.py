# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "matplotlib"]
# ///
"""ε-sweep near an exact Krylov breakdown: rowmap.F vs Rust vs exact ROW.

Inputs (see AGENTS.md quick reference):
    cargo test sweep_export -- --ignored          -> results/sweep_input.txt, results/sweep_rust.txt
    fortran/experiments/build/sweep results/sweep_input.txt results/sweep_f77.txt
Output: plots/out/sweep_f77.(pdf|png) and a table on stdout.
"""

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"


def numbers(path, skip=0):
    lines = [l for l in path.read_text().splitlines() if l.strip() and not l.startswith("#")]
    return np.array([float(l) for l in lines[skip:]])


def main():
    head = (RES / "sweep_input.txt").read_text().splitlines()[1].split()
    n, neps = int(head[0]), int(head[1])
    eps = numbers(RES / "sweep_input.txt", skip=1)[:neps]
    rs = numbers(RES / "sweep_rust.txt").reshape(neps, 3, n)
    f77 = numbers(RES / "sweep_f77.txt").reshape(neps, n)
    row, k4, kdef = rs[:, 0], rs[:, 1], rs[:, 2]
    nrm = np.linalg.norm
    dev = {
        "rowmap.F": [nrm(f77[i] - row[i]) / nrm(row[i]) for i in range(neps)],
        "Rust κ₁=4": [nrm(k4[i] - row[i]) / nrm(row[i]) for i in range(neps)],
        "Rust default": [nrm(kdef[i] - row[i]) / nrm(row[i]) for i in range(neps)],
    }
    chg_row = [nrm(row[i] - row[0]) for i in range(neps)]
    chg = {
        "rowmap.F": [nrm(f77[i] - f77[0]) for i in range(neps)],
        "Rust κ₁=4": [nrm(k4[i] - k4[0]) for i in range(neps)],
    }
    print("| ε | dev rowmap.F | dev Rust κ₁=4 | dev Rust default | change ROW | change F77 | change Rust |")
    print("|---|---|---|---|---|---|---|")
    for i in range(neps):
        print(f"| {eps[i]:.1e} | {dev['rowmap.F'][i]:.1e} | {dev['Rust κ₁=4'][i]:.1e} | "
              f"{dev['Rust default'][i]:.1e} | {chg_row[i]:.1e} | {chg['rowmap.F'][i]:.1e} | {chg['Rust κ₁=4'][i]:.1e} |")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for k, v in dev.items():
        ax[0].loglog(eps, np.maximum(v, 1e-17), "o-", ms=3, label=k)
    ax[0].set_title("deviation from exact ROW, ‖u − u_ROW‖/‖u_ROW‖")
    ax[1].loglog(eps, np.maximum(chg_row, 1e-17), "k-", lw=2.5, alpha=0.3, label="exact ROW")
    for k, v in chg.items():
        ax[1].loglog(eps, np.maximum(v, 1e-17), "o-", ms=3, label=k)
    ax[1].set_title("change ‖u(ε) − u(0)‖ (continuity)")
    for a in ax:
        a.set_xlabel("ε")
        a.grid(True, alpha=0.3)
        a.legend(fontsize=8)
    fig.suptitle("y' = (diag(−1,−10,−100) + εG) y, one step h = 0.1, n = 30")
    fig.tight_layout()
    (ROOT / "plots/out").mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(ROOT / f"plots/out/sweep_f77.{ext}", dpi=130)
    print("wrote plots/out/sweep_f77.(pdf|png)")


if __name__ == "__main__":
    main()
