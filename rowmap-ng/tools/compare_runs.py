# /// script
# requires-python = ">=3.11"
# ///
"""Side-by-side comparison of Rust (results/*.toml) and Fortran
(tests/data/fortran/*.toml) runs: ERR, nstep, nfeval, CPU time.

    uv run tools/compare_runs.py [--schedule paper|no-stage4-closure] [--markdown FILE]
"""

import argparse
import math
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def load(d, pattern):
    out = {}
    for f in sorted((ROOT / d).glob(pattern)):
        with open(f, "rb") as fh:
            t = tomllib.load(fh)
        r = t["run"]
        out[(r["problem"], r["m"], r.get("sft", 0), round(-math.log10(r["rtol"])))] = t
    return out


def reduced_vs_paper():
    """Rust runs: schedule without stage-4 closure vs paper schedule."""
    paper = load("results", "*_grk4t_sft0_tol*.toml")
    red = load("results", "*_grk4t_reduced_sft0_tol*.toml")
    lines = [
        "| problem | m | tol | ERR paper | ERR reduced | ratio | nstep paper/reduced | CPU paper | CPU reduced |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for key in sorted(red):
        if key not in paper:
            continue
        a, b = paper[key]["stats"], red[key]["stats"]
        p, m, _, k = key
        lines.append(
            f"| {p} | {m} | 1e-{k} | {a['err']:.2e} | {b['err']:.2e} | {b['err'] / a['err']:.2f} "
            f"| {a['nstep']}/{b['nstep']} | {a['cpu_time']:.4f} | {b['cpu_time']:.4f} |"
        )
    return "\n".join(lines)


def summary():
    """Per problem and grid: geometric means over tol = 1e-2..1e-7 of the
    Rust/Fortran ratios of ERR, nstep and CPU time."""
    rs = load("results", "*_grk4t_sft*_tol*.toml")
    f77 = load("tests/data/fortran", "*_grk4t_sft*_tol*.toml")
    groups = {}
    for key, a in rs.items():
        if key not in f77:
            continue
        p, m, sft, _ = key
        b = f77[key]["stats"]
        a = a["stats"]
        groups.setdefault((p, m, sft), []).append(
            (a["err"] / b["err"], a["nstep"] / b["nstep"], a["cpu_time"] / b["cpu_time"]))
    gm = lambda xs: math.exp(sum(math.log(x) for x in xs) / len(xs))
    lines = [
        "| problem | m | sft | ERR ratio (geo. mean) | ERR ratio range | nstep ratio | CPU ratio |",
        "|---|---|---|---|---|---|---|",
    ]
    order = ["bruss", "heat", "nilidi", "food", "vdpl", "difkon", "diffu2"]
    for (p, m, sft) in sorted(groups, key=lambda k: (order.index(k[0]), -k[1])):
        g = groups[(p, m, sft)]
        e = [x[0] for x in g]
        lines.append(
            f"| {p} | {m} | {sft} | {gm(e):.2f} | {min(e):.2f}–{max(e):.2f} | "
            f"{gm([x[1] for x in g]):.2f} | {gm([x[2] for x in g]):.1f} |")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--schedule", default="paper")
    ap.add_argument("--reduced-vs-paper", action="store_true")
    ap.add_argument("--summary", action="store_true")
    ap.add_argument("--markdown")
    args = ap.parse_args()
    if args.summary:
        text = summary()
        if args.markdown:
            Path(args.markdown).write_text(text + "\n")
            print(f"wrote {args.markdown}")
        else:
            print(text)
        return
    if args.reduced_vs_paper:
        text = reduced_vs_paper()
        if args.markdown:
            Path(args.markdown).write_text(text + "\n")
            print(f"wrote {args.markdown}")
        else:
            print(text)
        return
    tag = "grk4t" if args.schedule == "paper" else "grk4t_reduced"
    rs = load("results", f"*_{tag}_sft*_tol*.toml")
    f77 = load("tests/data/fortran", "*_grk4t_sft*_tol*.toml")
    lines = [
        f"| problem | m | tol | ERR rs | ERR f77 | ratio | nstep rs | nstep f77 | nfeval rs | nfeval f77 | CPU rs | CPU f77 |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for key in sorted(rs):
        if key not in f77:
            continue
        a, b = rs[key]["stats"], f77[key]["stats"]
        p, m, _, k = key
        ratio = a["err"] / b["err"] if b["err"] > 0 else float("nan")
        lines.append(
            f"| {p} | {m} | 1e-{k} | {a['err']:.2e} | {b['err']:.2e} | {ratio:.2f} | {a['nstep']} | {b['nstep']} "
            f"| {a['nfeval']} | {b['nfeval']} | {a['cpu_time']:.4f} | {b['cpu_time']:.4f} |"
        )
    text = "\n".join(lines)
    if args.markdown:
        Path(args.markdown).write_text(text + "\n")
        print(f"wrote {args.markdown}")
    else:
        print(text)


if __name__ == "__main__":
    main()
