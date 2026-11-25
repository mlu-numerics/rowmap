# /// script
# requires-python = ">=3.11"
# ///
"""Compiler/memory comparison: cost per Jacobian-vector product (paper grids).

The algorithms differ in how many Krylov vectors they build (rowmap.F stops
stages 2-4 early; Rust completes the proven schedule), so the CPU time is
normalized by the number of J·v products. Each J·v comes with one
orthogonalization against the current basis, so CPU/njv measures the cost
of one basis extension (user J·v + Gram-Schmidt + bookkeeping). With
`--single-pass` Rust does one classical Gram-Schmidt pass, as rowmap.F does
one MGS pass: the per-J·v work is then nearly identical.

Columns per problem (geometric means over tol = 1e-2 … 1e-7):
  us/Jv   CPU time per J·v product in microseconds
  Jv/step J·v products per step
  ratio   Rust/Fortran of us/Jv

    uv run tools/jv_cost.py [--tex FILE]
"""

import argparse
import math
import tomllib
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAPER_M = {"heat": 50, "nilidi": 69, "bruss": 50, "food": 50, "vdpl": 20,
           "difkon": 69, "diffu2": 69}
SERIES = [("f77", "grk4t", "paper", "F77"),
          ("rs", "grk4t", "single-pass", "Rust 1-pass"),
          ("rs", "grk4t", "paper", "Rust 2-pass")]


def gmean(xs):
    return math.exp(sum(map(math.log, xs)) / len(xs)) if xs else float("nan")


def load():
    s = defaultdict(dict)
    for d in [ROOT / "tests/data/fortran", ROOT / "results"]:
        for f in sorted(d.glob("*.toml")):
            t = tomllib.loads(f.read_text())
            if "run" not in t or "stats" not in t:
                continue
            r, st = t["run"], t["stats"]
            if PAPER_M.get(r["problem"]) != r["m"] or r.get("sft", 0) not in (0, 2):
                continue
            if r.get("idid", 1) < 0 or "njv" not in st or st["njv"] == 0:
                continue
            prod = "f77" if r["producer"] == "rowmap-f77" else "rs"
            key = (r["problem"], prod, r["method"], r.get("schedule", "paper"))
            s[key][r["rtol"]] = (st["cpu_time"], st["njv"], st["nstep"])
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tex", type=Path)
    args = ap.parse_args()
    s = load()
    rows = []
    for p in PAPER_M:
        cells = []
        ref = s.get((p, "f77", "grk4t", "paper"), {})
        for prod, meth, sched, _ in SERIES:
            runs = s.get((p, prod, meth, sched), {})
            tols = sorted(set(runs) & set(ref))
            if not tols:
                cells.append(None)
                continue
            us = gmean([runs[t][0] / runs[t][1] * 1e6 for t in tols])
            jvs = gmean([runs[t][1] / runs[t][2] for t in tols])
            ratio = gmean([(runs[t][0] / runs[t][1]) / (ref[t][0] / ref[t][1]) for t in tols])
            cells.append((us, jvs, ratio))
        rows.append((p, cells))

    head = ["problem"] + [f"{lab}: us/Jv  Jv/step  ratio" for *_, lab in SERIES]
    print(" | ".join(head))
    for p, cells in rows:
        out = [p.upper()]
        for c in cells:
            out.append("--" if c is None else f"{c[0]:7.1f} {c[1]:6.1f} {c[2]:5.2f}")
        print(" | ".join(out))

    if args.tex:
        lines = [r"\begin{tabular}{l" + "rrr" * len(SERIES) + "}", r"\toprule",
                 "& " + " & ".join(rf"\multicolumn{{3}}{{c}}{{{lab}}}" for *_, lab in SERIES)
                 + r" \\",
                 "Problem" + r" & $\mu$s/$Jv$ & $Jv$/step & ratio" * len(SERIES) + r" \\",
                 r"\midrule"]
        for p, cells in rows:
            parts = [p.upper()]
            for c in cells:
                parts += ["--"] * 3 if c is None else [f"{c[0]:.1f}", f"{c[1]:.1f}", f"{c[2]:.2f}"]
            lines.append(" & ".join(parts) + r" \\")
        lines += [r"\bottomrule", r"\end{tabular}"]
        args.tex.write_text("\n".join(lines) + "\n")
        print(f"wrote {args.tex}")


if __name__ == "__main__":
    main()
