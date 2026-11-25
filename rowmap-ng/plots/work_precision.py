# /// script
# requires-python = ">=3.11"
# dependencies = ["matplotlib", "numpy"]
# ///
"""Work-precision diagrams in the style of the paper (Figs. 3.1-3.6, 4.7-4.8).

Reads every run TOML (tests/data/fortran/*.toml, results/*.toml) and draws
lg(ERR) vs CPU time per problem and grid, one series per producer/method/sft.

    uv run plots/work_precision.py [--grid paper|coarse|all]
Output: plots/out/wp_<problem>_m<m>.pdf / .png
"""

import argparse
import math
import tomllib
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FormatStrFormatter, LogLocator, NullFormatter  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SOURCES = [ROOT / "tests/data/fortran", ROOT / "results"]
OUT = ROOT / "plots/out"
PAPER_M = {"heat": 50, "nilidi": 69, "bruss": 50, "food": 50, "food3": 10,
           "vdpl": 20, "difkon": 69, "diffu2": 69}
PRODUCER = {"rowmap-f77": "F77", "rowmap-rs": "Rust"}


def load():
    runs = []
    for d in SOURCES:
        for f in sorted(d.glob("*.toml")) if d.exists() else []:
            with open(f, "rb") as fh:
                t = tomllib.load(fh)
            if "run" in t and "stats" in t:  # skip other TOML files (e.g. sweeps)
                runs.append(t)
    return runs


def label(run, show_sft=True):
    r = run["run"]
    s = f"{PRODUCER.get(r['producer'], r['producer'])} {r['method'].upper()}"
    if show_sft and r.get("sft", 0):
        s += f" sft={r['sft']}"
    if r.get("schedule", "paper") != "paper":
        s += f" ({r['schedule']})"
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", choices=["paper", "coarse", "all"], default="paper")
    ap.add_argument("--native-only", action="store_true",
                    help="only sft 0/2 (no augmented-system strategies)")
    ap.add_argument("--out", type=Path, default=OUT, help="output directory")
    args = ap.parse_args()
    out = args.out

    groups = defaultdict(lambda: defaultdict(list))
    for run in load():
        r, st = run["run"], run["stats"]
        if args.grid == "paper" and r["m"] != PAPER_M.get(r["problem"]):
            continue
        if args.grid == "coarse" and r["m"] == PAPER_M.get(r["problem"]):
            continue
        err = st["err"]
        if r["idid"] < 0 or not math.isfinite(err) or err <= 0:
            continue   # failed runs are missing marks, as in the paper
        if args.native_only and r.get("sft", 0) not in (0, 2):
            continue
        if "single-pass" in r.get("schedule", ""):
            continue   # comparison-only variant (tools/jv_cost.py)
        groups[(r["problem"], r["m"], r["n"])][label(run, not args.native_only)].append(
            (r["rtol"], max(st["cpu_time"], 1e-4), err))

    out.mkdir(parents=True, exist_ok=True)
    markers = "os^vD<>p*"
    # one fixed style per series label across all figures
    all_labels = sorted({lab for s in groups.values() for lab in s})
    style = {lab: (f"C{i % 10}", markers[i % len(markers)]) for i, lab in enumerate(all_labels)}
    for (prob, m, n), series in sorted(groups.items()):
        fig, ax = plt.subplots(figsize=(5.2, 4.0))
        for lab, pts in sorted(series.items()):
            pts.sort(reverse=True)  # crude tolerance first
            color, marker = style[lab]
            ax.plot([p[1] for p in pts], [math.log10(p[2]) for p in pts],
                    color=color, marker=marker, label=lab, lw=1.2, ms=5)
        ax.set_xscale("log")
        # plain tick labels at 1, 2, 5 × 10^k (the default minor labels overlap)
        ax.xaxis.set_major_locator(LogLocator(subs=(1.0, 2.0, 5.0)))
        ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel("CPU time [s]")
        ax.set_ylabel("lg(ERR)")
        ax.set_title(f"{prob.upper()}, N={n}")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        for ext in ("pdf", "png"):
            fig.savefig(out / f"wp_{prob}_m{m}.{ext}", dpi=130)
        plt.close(fig)
        print(f"wrote {out}/wp_{prob}_m{m}.(pdf|png)")


if __name__ == "__main__":
    main()
