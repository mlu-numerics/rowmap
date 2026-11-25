#!/usr/bin/env bash
# Runs all paper experiments with the original Fortran ROWMAP (AGENTS.md phase 1).
#
#   fortran/run_experiments.sh [paper|coarse|all]   (default: all)
#
# Output: tests/data/fortran/<problem>_m<m>_<method>_sft<s>_tol<k>.toml
# Needs the reference solutions in tests/data/reference/ (tools/reference.py)
# to compute ERR; without one, err = nan.
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
exp="$root/fortran/experiments"
out="$root/tests/data/fortran"
ref="$root/tests/data/reference"
which=${1:-all}

make -s -C "$exp"
mkdir -p "$out"

declare -A paper_m=([heat]=50 [nilidi]=69 [bruss]=50 [food]=50 [food3]=10 [vdpl]=20 [difkon]=69 [diffu2]=69)
declare -A coarse_m=([heat]=10 [nilidi]=10 [bruss]=10 [food]=10 [food3]=5 [vdpl]=8 [difkon]=10 [diffu2]=10)
methods=([2]=grk4t [5]=lstab)   # iwork(2) -> name; only GRK4T is run
maxkry=66   # as in the paper

run() { # problem m method sft tolexp
    local p=$1 m=$2 meth=$3 sft=$4 k=$5 bin=driver rf
    [[ $sft -ge 4 ]] && bin=driver_sft
    rf="$ref/${p}_m${m}.txt"
    [[ -f $rf ]] || rf=-
    local f="$out/${p}_m${m}_${methods[$meth]}_sft${sft}_tol${k}.toml"
    "$exp/build/$bin" run "$p" "$m" "1e-$k" "$meth" "$maxkry" "$sft" "$rf" "$f" > "$f.log" 2>&1 || true
    [[ -s $f.log ]] || rm -f "$f.log"
    printf '%-7s m=%-3s %-6s sft=%s tol=1e-%s  %s\n' "$p" "$m" "${methods[$meth]}" "$sft" "$k" \
        "$(rg -N '^(idid|err|nstep|nfeval|cpu_time) =' "$f" | tr -s ' ' | tr '\n' ' ')"
}

grids=()
[[ $which == paper || $which == all ]] && grids+=(paper)
[[ $which == coarse || $which == all ]] && grids+=(coarse)

for g in "${grids[@]}"; do
    for p in heat nilidi bruss food food3 vdpl difkon diffu2; do
        if [[ $g == paper ]]; then m=${paper_m[$p]}; else m=${coarse_m[$p]}; fi
        for k in 2 3 4 5 6 7; do
            case $p in
                difkon|diffu2) for s in 2 3 4 5; do run "$p" "$m" 2 "$s" "$k"; done ;;
                *)             run "$p" "$m" 2 0 "$k" ;;
            esac
        done
    done
done
