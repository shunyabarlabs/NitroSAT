#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/../.." && pwd)
bin=${1:-"$root/src/c/v3/nitrosatv3"}
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

out=$("$bin" "$root/tests/v3/base.cnf" --add "$root/tests/v3/increment.cnf" \
    --store "$tmp/combined.nsv3" --epochs 40 --batch-clauses 2 \
    --active-clauses 16 --active-literals 64)
python3 - "$out" "$tmp/combined.nsv3" <<'PY'
import json, os, sys
d = json.loads(sys.argv[1])
assert d["solved"] is True, d
assert d["variables"] == 5, d
assert d["clauses"] == 6, d
assert d["literals"] == 9, d
assert d["increment_files"] == 1, d
# 40-byte header + six 16-byte records + nine 4-byte literals.
assert os.path.getsize(sys.argv[2]) == 40 + 6 * 16 + 9 * 4
PY

# Raw WCNF: hard clauses must hold and the lower-cost side of a conflicting
# soft pair must be selected. Full satisfaction is impossible by construction.
wout=$("$bin" "$root/tests/v3/weighted.wcnf" --epochs 40 --batch-clauses 1 \
    --active-clauses 16 --active-literals 64)
python3 - "$wout" <<'PY'
import json, sys
d = json.loads(sys.argv[1])
assert d["format"] == "wcnf", d
assert d["feasible"] is True and d["solved"] is False, d
assert d["hard_unsatisfied"] == 0, d
assert d["soft_unsatisfied"] == 1, d
assert d["soft_cost"] == 3, d
assert d["total_soft_weight"] == 8, d
PY

# Regression: soft repair must not overwrite a feasible lower-cost state. The
# known optimum sets every variable true: zero hard violations and cost 300.
awk 'BEGIN { n=10; reps=10; print "p wcnf 10 300 10"; for (r=0;r<reps;r++) for (v=1;v<=n;v++) { print "10 " v " 0"; print "5 " v " 0"; print "3 -" v " 0" } }' \
    >"$tmp/repair-guard.wcnf"
guard_out=$("$bin" "$tmp/repair-guard.wcnf" --epochs 10 --batch-clauses 32 \
    --active-clauses 512 --active-literals 1024)
python3 - "$guard_out" <<'PY'
import json, sys
d = json.loads(sys.argv[1])
assert d["feasible"] is True and d["hard_unsatisfied"] == 0, d
assert d["soft_cost"] == 300 and d["total_soft_weight"] == 800, d
PY

hout=$("$bin" "$root/tests/v3/hard_marker.wcnf" --epochs 20 --batch-clauses 1 \
    --active-clauses 8 --active-literals 16)
python3 - "$hout" <<'PY'
import json, sys
d = json.loads(sys.argv[1])
assert d["format"] == "wcnf", d
assert d["solved"] is True and d["hard_unsatisfied"] == 0, d
assert d["total_soft_weight"] == 4 and d["soft_cost"] == 0, d
PY

set +e
"$bin" "$root/tests/v3/hard_conflict.wcnf" --epochs 5 --batch-clauses 1 \
    >"$tmp/hard-conflict.json"
hard_conflict_rc=$?
set -e
test "$hard_conflict_rc" -eq 1
python3 - "$tmp/hard-conflict.json" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
assert d["feasible"] is False and d["hard_unsatisfied"] == 1, d
assert d["finisher_passes"] > 0, d
PY

# Mixed persistent stream: ordinary CNF records remain hard while WCNF records
# retain their explicit weights and top-derived hard status.
"$bin" "$root/tests/v3/base.cnf" --add "$root/tests/v3/weighted_increment.wcnf" \
    --store "$tmp/mixed.nsv3" --store-only >"$tmp/mixed.json"
"$bin" --load-store "$tmp/mixed.nsv3" --store-only >"$tmp/mixed-loaded.json"
python3 - "$tmp/mixed.json" "$tmp/mixed-loaded.json" <<'PY'
import json, sys
for path in sys.argv[1:]:
    d = json.load(open(path))
    assert d["format"] == "wcnf", d
    assert (d["clauses"], d["literals"]) == (5, 9), d
PY

# Persistent incremental path: build once, append without reparsing the base,
# solve from the store, and write a complete assignment.
"$bin" "$root/tests/v3/base.cnf" --store "$tmp/persistent.nsv3" --store-only >"$tmp/base-store.json"
"$bin" --load-store "$tmp/persistent.nsv3" --add "$root/tests/v3/increment.cnf" \
    --store-only >"$tmp/extended-store.json"
persist_out=$("$bin" --load-store "$tmp/persistent.nsv3" --epochs 40 \
    --batch-clauses 2 --active-clauses 16 --active-literals 64 \
    --solution "$tmp/assignment.sol")
python3 - "$tmp/base-store.json" "$tmp/extended-store.json" "$persist_out" "$tmp/assignment.sol" <<'PY'
import json, sys
base = json.load(open(sys.argv[1]))
extended = json.load(open(sys.argv[2]))
solved = json.loads(sys.argv[3])
assignment = open(sys.argv[4]).read().split()
assert (base["clauses"], base["literals"]) == (3, 6), base
assert (extended["clauses"], extended["literals"]) == (6, 9), extended
assert solved["solved"] is True, solved
assert len(assignment) == 6 and assignment[-1] == "0", assignment
PY

set +e
"$bin" "$root/tests/v3/base.cnf" --add "$root/tests/v3/contradiction.cnf" \
    --epochs 5 --batch-clauses 2 >"$tmp/unsat.json"
unsat_rc=$?
"$bin" "$root/tests/v3/malformed.cnf" --store "$tmp/bad.nsv3" --store-only \
    >"$tmp/bad.out" 2>"$tmp/bad.err"
bad_rc=$?
before_size=$(stat -c %s "$tmp/persistent.nsv3")
before_hash=$(sha256sum "$tmp/persistent.nsv3" | cut -d' ' -f1)
"$bin" --load-store "$tmp/persistent.nsv3" --add "$root/tests/v3/malformed.cnf" \
    --store-only >"$tmp/rollback.out" 2>"$tmp/rollback.err"
rollback_rc=$?
"$bin" "$root/tests/v3/malformed.wcnf" --store "$tmp/bad-wcnf.nsv3" --store-only \
    >"$tmp/bad-wcnf.out" 2>"$tmp/bad-wcnf.err"
bad_wcnf_rc=$?
after_size=$(stat -c %s "$tmp/persistent.nsv3")
after_hash=$(sha256sum "$tmp/persistent.nsv3" | cut -d' ' -f1)
set -e
test "$unsat_rc" -eq 1
test "$bad_rc" -eq 2
test "$rollback_rc" -eq 2
test "$bad_wcnf_rc" -eq 2
test ! -e "$tmp/bad.nsv3"
test ! -e "$tmp/bad-wcnf.nsv3"
test "$before_size" = "$after_size"
test "$before_hash" = "$after_hash"

cp "$tmp/persistent.nsv3" "$tmp/corrupt.nsv3"
printf '\001' >>"$tmp/corrupt.nsv3"
set +e
"$bin" --load-store "$tmp/corrupt.nsv3" --store-only >"$tmp/corrupt.out" 2>"$tmp/corrupt.err"
corrupt_rc=$?
set -e
test "$corrupt_rc" -eq 2
python3 - "$tmp/unsat.json" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
assert d["solved"] is False, d
assert d["clauses"] == 5, d
assert d["increment_files"] == 1, d
PY

echo "v3 end-to-end tests passed"
