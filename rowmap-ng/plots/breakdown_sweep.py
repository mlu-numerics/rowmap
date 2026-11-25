# /// script
# requires-python = ">=3.11"
# dependencies = ["matplotlib"]
# ///
"""ε-continuity sweep near an exact Krylov breakdown (AGENTS.md §5.3 item 4).

Data: results/sweep_autonomous.toml, written by
    cargo test sweep_report -- --ignored --nocapture
Output: plots/out/sweep_autonomous.(pdf|png)
"""

import tomllib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "results/sweep_autonomous.toml"
OUT = ROOT / "plots/out"


def main():
    with open(DATA, "rb") as fh:
        d = tomllib.load(fh)
    fig, axes = plt.subplots(1, len(d), figsize=(11, 4.2), sharey=True)
    for ax, (label, s) in zip(axes, d.items()):
        ax.loglog(s["eps"], s["change_row"], "k-", lw=2.5, alpha=0.3, label="‖u_ROW(ε) − u_ROW(0)‖")
        ax.loglog(s["eps"], s["change_k"], "o-", ms=3, label="‖u_K(ε) − u_K(0)‖")
        ax.loglog(s["eps"], s["dev_row"], "s-", ms=3, label="‖u_K(ε) − u_ROW(ε)‖ / ‖u_ROW‖")
        ax.loglog(s["eps"], s["orth"], "^-", ms=3, label="max |QᵀQ − I|")
        ax2 = ax.twinx()
        ax2.step(s["eps"], s["dim4"], "C3:", where="mid", label="κ₄")
        ax2.set_ylabel("κ₄", color="C3")
        ax.set_title(f"diag(3 eigenvalues) + εG, h = 0.1 — {label}")
        ax.set_xlabel("ε")
        ax.grid(True, which="major", alpha=0.3)
    axes[0].legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"sweep_autonomous.{ext}", dpi=130)
    print("wrote plots/out/sweep_autonomous.(pdf|png)")


if __name__ == "__main__":
    main()
