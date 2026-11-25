# /// script
# requires-python = ">=3.11"
# ///
"""CPU-time ratio Rust/Fortran at equal accuracy (paper grids).

For every Fortran run (GRK4T, sft 0/2) the Rust CPU time at the same ERR is
interpolated log-log from the Rust series; the table shows the geometric mean
of these ratios (and how many Fortran points lie inside the Rust ERR range).
Values < 1: Rust is faster at equal accuracy.

    uv run tools/wp_ratio.py [--tex FILE]
"""

import argparse
import math
import tomllib
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAPER_M = {"heat": 50, "nilidi": 69, "bruss": 50, "food": 50, "vdpl": 20,
           "difkon": 69, "diffu2": 69}
SERIES = [("grk4t", "paper", "GRK4T"),
          ("grk4t", "no-stage4-closure", "GRK4T reduced"),
          ("ros34pw2", "paper", "ROS34PW2")]


def load():
    s = defaultdict(list)
    for d in [ROOT / "tests/data/fortran", ROOT / "results"]:
        for f in sorted(d.glob("*.toml")):
            t = tomllib.loads(f.read_text())
            if "run" not in t or "stats" not in t:
                continue
            r, st = t["run"], t["stats"]
            if PAPER_M.get(r["problem"]) != r["m"] or not math.isfinite(st["err"]):
                continue
            if r.get("sft", 0) not in (0, 2) or r.get("idid", 1) < 0:
                continue
            prod = "f77" if r["producer"] == "rowmap-f77" else "rs"
            s[(r["problem"], prod, r["method"], r.get("schedule", "paper"))].append(
                (st["err"], st["cpu_time"]))
    return s


def cpu_at(pts, e):
    pts = sorted(pts)
    for (e0, c0), (e1, c1) in zip(pts, pts[1:]):
        if e0 <= e <= e1 and e1 > e0:
            w = (math.log(e) - math.log(e0)) / (math.log(e1) - math.log(e0))
            return math.exp(math.log(c0) + w * (math.log(c1) - math.log(c0)))
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tex", type=Path)
    args = ap.parse_args()
    s = load()
    rows = []
    for p in PAPER_M:
        ref = s[(p, "f77", "grk4t", "paper")]
        cells = []
        for meth, sched, _ in SERIES:
            rs = s.get((p, "rs", meth, sched))
            rat = [cpu_at(rs, e) / c for e, c in ref if rs and cpu_at(rs, e)]
            cells.append((math.exp(sum(map(math.log, rat)) / len(rat)), len(rat))
                         if rat else None)
        rows.append((p, cells))

    hdr = ["problem"] + [lab for *_, lab in SERIES]
    print(" | ".join(hdr))
    for p, cells in rows:
        print(" | ".join([p.upper()] + [f"{c[0]:.2f} ({c[1]})" if c else "--" for c in cells]))

    if args.tex:
        lines = [r"\begin{tabular}{l" + "r" * len(SERIES) + "}", r"\toprule",
                 " & ".join(["Problem"] + [lab for *_, lab in SERIES]) + r" \\", r"\midrule"]
        for p, cells in rows:
            lines.append(" & ".join([p.upper()] + [
                f"{c[0]:.2f} ({c[1]})" if c else "--" for c in cells]) + r" \\")
        lines += [r"\bottomrule", r"\end{tabular}"]
        args.tex.write_text("\n".join(lines) + "\n")
        print(f"wrote {args.tex}")


if __name__ == "__main__":
    main()
