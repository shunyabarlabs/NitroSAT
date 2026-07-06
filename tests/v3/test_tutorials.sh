#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/../.." && pwd)
bin=${1:-"$root/src/c/v3/nitrosatv3"}
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

python3 "$root/examples/v3/sudoku.py" encode \
    "$root/examples/v3/puzzle.txt" "$tmp/sudoku.cnf"
"$bin" "$tmp/sudoku.cnf" --epochs 5 --indexed-finisher \
    --indexed-flips 5000 --indexed-checkpoint 250 --exact \
    --exact-memory-mb 256 --exact-max-conflicts 1000000 \
    --solution "$tmp/sudoku.sol" >"$tmp/sudoku.json"
python3 - "$tmp/sudoku.json" <<'PY'
import json, sys
result = json.load(open(sys.argv[1]))
assert result["solved"] is True and result["satisfied"] == 3270, result
PY
python3 "$root/examples/v3/sudoku.py" decode "$tmp/sudoku.sol" \
    >"$tmp/sudoku-grid.txt"
grep -q '| 5 3 4 | 6 7 8 | 9 1 2 |' "$tmp/sudoku-grid.txt"

python3 "$root/examples/v3/workforce.py" generate "$tmp/workforce.wcnf"
"$bin" "$tmp/workforce.wcnf" --epochs 100 --finisher-passes 1000 \
    --finisher-max-clauses 10000 --solution "$tmp/workforce.sol" \
    >"$tmp/workforce.json"
python3 - "$tmp/workforce.json" <<'PY'
import json, sys
result = json.load(open(sys.argv[1]))
assert result["feasible"] is True and result["hard_unsatisfied"] == 0, result
assert result["soft_cost"] == 4 and result["total_soft_weight"] == 29, result
PY
python3 "$root/examples/v3/workforce.py" decode "$tmp/workforce.sol" \
    >"$tmp/workforce-schedule.txt"
grep -q 'preference score: 25/29' "$tmp/workforce-schedule.txt"

models=(
    university-timetabling
    vehicle-assignment
    meeting-room-scheduling
    graph-coloring
    kubernetes-pod-placement
    product-configuration
    exam-scheduling
    nurse-rostering
)

for model in "${models[@]}"; do
    formula="$tmp/$model.formula"
    solution="$tmp/$model.sol"
    python3 "$root/examples/v3/tutorial_models.py" generate "$model" "$formula"
    if [[ "$model" == graph-coloring ]]; then
        "$bin" "$formula" --epochs 20 --indexed-finisher \
            --indexed-flips 10000 --indexed-checkpoint 100 --exact \
            --solution "$solution" >"$tmp/$model.json"
    else
        "$bin" "$formula" --epochs 100 --finisher-passes 1000 \
            --finisher-max-clauses 10000 --solution "$solution" \
            >"$tmp/$model.json"
    fi
    python3 "$root/examples/v3/tutorial_models.py" validate "$model" "$solution"
done

echo "all executable tutorials passed"
