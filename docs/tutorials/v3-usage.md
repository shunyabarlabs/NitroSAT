# NitroSAT V3 — Hands-On Tutorial

A practical guide to the bounded-memory SAT/MaxSAT solver. Every example is
runnable on the binaries and CNF/WCNF inputs already in this repository.

## 0. What V3 Is (and Is Not)

V3 is a **streaming SAT/MaxSAT solver**. It does **not** load the clause table,
clause state, or a variable-to-clause index into RAM. Instead it converts the
input into a binary stream and reads it sequentially on every global pass.

| Use V3 when... | Use V2 when... |
|---|---|
| Formula exceeds practical RAM | Formula fits in RAM and you need raw speed |
| You need WCNF / partial MaxSAT | You need DRAT proofs in the default path |
| You need incremental `--add` constraint files | The instance is dense but small |
| The solver must run inside a memory ceiling | You want the dense topological repair pipeline |

V3 is heuristic — it returns a satisfying assignment when it finds one, and a
DRAT proof only via `--exact`. There is no optimality proof for partial
MaxSAT cost; the reported `soft_cost` is the best found.

---

## 1. Build and Verify

```bash
make -C src/c v3
src/c/v3/nitrosatv3 --help
```

A successful build produces `src/c/v3/nitrosatv3` (~60 KB). Run a smoke test:

```bash
src/c/v3/nitrosatv3 tests/v3/base.cnf --cinematic
```

Expected:

```text
epoch 1: unsatisfied=0/3 cached=0
NitroSAT V3 (CNF): all clauses satisfied
Variables: 4  Clauses: 3  Literals: 6
```

---

## 2. The Three Input Formats

| Format | Header | Use |
|---|---|---|
| **CNF** | `p cnf VARS CLAUSES` | Pure Boolean satisfiability |
| **WCNF (with top)** | `p wcnf VARS CLAUSES TOP` | Weighted partial MaxSAT — weights ≥ TOP are hard, others soft |
| **WCNF (legacy)** | `p wcnf VARS CLAUSES` | Every clause is soft |
| **WCNF (modern, prefix)** | `p wcnf VARS CLAUSES` + `h ...` prefixes | Soft clause marked hard inline |

Tiny CNF example (`tests/v3/base.cnf`):

```text
c base formula
p cnf 4 3
1 2 0
-1 3 0
-2 -3 4 0
```

Tiny WCNF example (`tests/v3/weighted.wcnf`):

```text
p wcnf 2 3 10
10 1 0
10 -1 0
3 2 0
```

— `1` and `-1` are hard (weight ≥ top=10); clause 3 is soft (weight 3).

---

## 3. Understanding the JSON Output

By default V3 prints a single JSON object on stdout. The key fields:

```json
{
  "solver": "NitroSAT V3",
  "format": "cnf" | "wcnf",
  "status": "SATISFIED" | "PARTIAL" | "UNSAT",
  "solved":   true | false,   // every clause satisfied
  "feasible": true | false,   // no hard clause violated
  "variables": 729,
  "clauses":   3270,
  "literals":  8778,
  "satisfied": 3270,
  "unsatisfied": 0,
  "hard_unsatisfied":  0,    // WCNF only
  "soft_unsatisfied":  0,    // WCNF only
  "soft_cost":         0,    // WCNF only — sum of soft clause weights violated
  "total_soft_weight": 0,    // WCNF only — sum of all soft weights
  "exact_attempted":   false,
  "exact_result":      "NOT_RUN" | "SAT" | "UNSAT" | "UNKNOWN" | "OOM" | "TIMEOUT",
  "epochs":            30,
  "finisher_passes":   256,
  "indexed_flips":     0,
  "increment_files":   0,
  "bounded_memory_mb": 4.87,
  "store_build_ms":    0.84,
  "solve_ms":          567.6
}
```

Add `--cinematic` for plain-text output instead.

---

## 4. End-to-End Example: Solving a Sudoku

### 4.1 Encode the puzzle to CNF

```bash
python3 examples/v3/sudoku.py encode examples/v3/puzzle.txt /tmp/sudoku.cnf
```

Output: `wrote /tmp/sudoku.cnf: 729 variables, 3270 clauses`

### 4.2 Solve with the indexed finisher (recommended for Sudoku)

```bash
src/c/v3/nitrosatv3 /tmp/sudoku.cnf \
  --epochs 100 \
  --indexed-finisher --indexed-flips 200000 --indexed-checkpoint 500 \
  --solution /tmp/sudoku.sol > /tmp/sudoku.json
cat /tmp/sudoku.json | python3 -m json.tool
```

Key result fields: `"solved": true`, `"unsatisfied": 0`, `"solve_ms": ~2400`.

### 4.3 Decode and visualize

```bash
python3 examples/v3/sudoku.py decode /tmp/sudoku.sol
```

```text
+-------+-------+-------+
| 5 3 4 | 6 7 8 | 9 1 2 |
| 6 7 2 | 1 9 5 | 3 4 8 |
| 1 9 8 | 3 4 2 | 5 6 7 |
+-------+-------+-------+
| 8 5 9 | 7 6 1 | 4 2 3 |
| 4 2 6 | 8 5 3 | 7 9 1 |
| 7 1 3 | 9 2 4 | 8 5 6 |
+-------+-------+-------+
| 9 6 1 | 5 3 7 | 2 8 4 |
| 2 8 7 | 4 1 9 | 6 3 5 |
| 3 4 5 | 2 8 6 | 1 7 9 |
+-------+-------+-------+
```

### 4.4 If the heuristic stalls

Try the exact CDCL fallback:

```bash
src/c/v3/nitrosatv3 /tmp/sudoku.cnf --exact \
  --exact-max-clauses 200000 --exact-memory-mb 2048 \
  --solution /tmp/sudoku.sol
```

`--exact` invokes the embedded CDCL solver when the heuristic plateau is
detected or when the user requests it. The exit code is `10` for SAT, `20`
for UNSAT, `1` for failure.

---

## 5. End-to-End Example: Partial MaxSAT (WCNF)

The university-timetabling tutorial produces a WCNF file with hard
"exactly-one-slot" constraints and soft "prefer this slot" preferences.

```bash
python3 examples/v3/tutorial_models.py generate \
  university-timetabling /tmp/uni.wcnf

src/c/v3/nitrosatv3 /tmp/uni.wcnf \
  --epochs 100 --finisher-passes 1000 --finisher-max-clauses 10000 \
  --solution /tmp/uni.sol > /tmp/uni.json
```

Check the JSON:

```bash
python3 -c "import json; d=json.load(open('/tmp/uni.json')); \
  print('feasible:', d['feasible'], 'hard_unsat:', d['hard_unsatisfied'], \
        'soft_cost:', d['soft_cost'], '/', d['total_soft_weight'])"
```

Validate independently:

```bash
python3 examples/v3/tutorial_models.py validate \
  university-timetabling /tmp/uni.sol
```

---

## 6. Incremental Solving with `--add`

V3 lets you start with a base formula and append constraint sets later.
Variables may extend the base's declared range — V3 uses the maximum
across all inputs.

```bash
src/c/v3/nitrosatv3 tests/v3/base.cnf \
  --add tests/v3/increment.cnf \
  --store /tmp/combined.nsv3 --epochs 40 \
  --batch-clauses 2 --active-clauses 16 --active-literals 64
```

Output JSON will report `"increment_files": 1` and the merged clause count
(`base.clauses + increment.clauses = 6`).

If any `--add` file is malformed, V3 **rolls back the in-progress store to
its prior byte length** and returns exit code `2`. The previous store is
untouched.

---

## 7. Persistent Stores and Re-solving

Build a store once, solve or extend it many times without re-parsing
the base CNF:

```bash
# Build once
src/c/v3/nitrosatv3 tests/v3/base.cnf \
  --store /tmp/persistent.nsv3 --store-only > /tmp/base.json

# Re-solve from the store
src/c/v3/nitrosatv3 --load-store /tmp/persistent.nsv3 \
  --epochs 40 --solution /tmp/assignment.sol

# Append more clauses transactionally
src/c/v3/nitrosatv3 --load-store /tmp/persistent.nsv3 \
  --add tests/v3/increment.cnf --store-only
```

The on-disk store has a 40-byte header (magic, version, header size,
variable count, format flags, clause count, literal count) followed by
records of `uint32 literal_count | uint32 flags | uint64 weight | int32[] literals`.

---

## 8. Heuristic Tuning Cheat-Sheet

### 8.1 Knobs that matter

| Flag | Default | When to change |
|---|---|---|
| `--epochs N` | 20 | Increase (100–1000) on structured instances where the heuristic plateaus |
| `--finisher-passes N` | 256 | Local-search sweeps over the streamed store |
| `--finisher-batch-flips N` | 64 | Increase (128–512) for dense coloring; decrease (1) for tight parity |
| `--finisher-max-clauses N` | 100,000 | Default skips finisher above this ceiling — override for big formulas |
| `--batch-clauses N` | 65,536 | Smaller = more gradient updates per pass |
| `--learning-rate X` | 0.015 | Tune down if the optimizer oscillates |
| `--seed N` | time-based | Set explicitly for reproducible runs |

### 8.2 Recipes from the V3 README

**Tight parity / small difficult formulas** — WalkSAT-heavy local search:

```bash
nitrosatv3 parity.cnf --epochs 2000 \
  --finisher-passes 10000 --finisher-batch-flips 1 \
  --finisher-max-clauses 100000
```

**Planted graph coloring** — heat-batch heavy:

```bash
nitrosatv3 planted.cnf --epochs 200 --finisher-passes 500 \
  --finisher-batch-flips 64 --finisher-max-clauses 500000
```

**Difficult structured CNF** — disk-indexed WalkSAT:

```bash
nitrosatv3 hard.cnf --epochs 20 --indexed-finisher \
  --indexed-flips 200000 --indexed-checkpoint 1000
```

**Difficult WCNF with hard constraints** — hard-feasibility indexed repair,
then streamed soft-cost improvement:

```bash
nitrosatv3 model.wcnf --epochs 20 --indexed-finisher \
  --indexed-flips 200000 --indexed-checkpoint 1000 \
  --finisher-passes 1000
```

For WCNF, the indexed finisher is deliberately **hard-only**. It uses the
disk-backed index to drive `hard_unsatisfied` to zero, then hands the best
hard-feasible assignment to the normal streamed finisher for soft-cost
improvement. This means `feasible=true` is the primary correctness signal;
`soft_cost` is still heuristic and not an optimality proof.

---

## 9. UNSAT Detection

V3 detects UNSAT by GF(2) parity / topological probe ordering on small
instances. For formal proof, use `--exact`:

```bash
src/c/v3/nitrosatv3 tests/v3/contradiction.cnf --exact \
  --proof /tmp/unsat.drat
```

JSON `"exact_result": "UNSAT"` and `"proof_generated": true`.

For huge formulas where `--exact` would exhaust memory, fall back to the
heuristic alone; the `"status": "UNSAT"` line means V3 found structural
inconsistency without a formal proof.

---

## 10. Memory and Performance Guardrails

| Flag | Default | Purpose |
|---|---|---|
| `--active-clauses N` | 65,536 | Local repair clause ceiling |
| `--active-literals N` | 1,048,576 | Local repair literal ceiling |
| `--exact-max-clauses N` | 1,000,000 | Hard ceiling on `--exact` CDCL |
| `--exact-memory-mb N` | 1,024 | CDCL memory cap |
| `--exact-max-conflicts N` | 5,000,000 | CDCL conflict cap |

`"bounded_memory_mb"` in the JSON is the **peak solver allocation** reported
by the runtime. For 37M-clause enterprise timetabling on this machine the
peak was **9.62 MiB** with default settings.

---

## 11. Quick-Reference Command Index

```text
# Solve
nitrosatv3 base.cnf [--epochs N] [--cinematic] [--solution OUT]

# Partial MaxSAT
nitrosatv3 model.wcnf [--epochs N] [--finisher-passes N]
nitrosatv3 model.wcnf --indexed-finisher [--indexed-flips N] [--finisher-passes N]

# Incremental
nitrosatv3 base.cnf --add c1.cnf --add c2.cnf --store s.nsv3

# Persistent
nitrosatv3 base.cnf --store s.nsv3 --store-only
nitrosatv3 --load-store s.nsv3 [--add more.cnf] [--epochs N]

# Refinement
nitrosatv3 hard.cnf --indexed-finisher --indexed-flips N

# Exact / DRAT proof
nitrosatv3 formula.cnf --exact --proof out.drat

# Reproducibility
nitrosatv3 formula.cnf --seed 42
```

---

## 12. Common Pitfalls

1. **Assuming WCNF `--indexed-finisher` optimizes soft clauses directly** →
   it does not. The indexed phase repairs hard clauses only; soft-cost
   improvement happens afterward through the streamed finisher.
2. **Forgetting `--solution`** → no assignment file is written. The JSON
   alone is not enough to reconstruct the model.
3. **Trusting `solved` for WCNF** → check `feasible` first. `feasible=true`
   means all hard clauses hold; `solved=true` additionally means all soft
   clauses hold.
4. **Assuming `soft_cost` is optimal** → V3 reports the best found, not a
   proof of optimality.
5. **`--exact` on huge formulas** → almost always OOM. Use it as a
   fallback on formulas that fit the clause / memory ceilings.
6. **Running V3 instead of V2 on small instances** → V2 is 5–10× faster on
   anything that fits in RAM.

---

## 13. End-to-End Test Script

`tests/v3/test_v3.sh` is the canonical end-to-end test. It exercises:

- base + `--add` combination
- WCNF partial MaxSAT cost reporting
- soft-repair guard (does not overwrite a feasible lower-cost state)
- hard-conflict rejection (`exit 1`)
- malformed input rollback (`exit 2`, store unchanged)
- corrupt store rejection (`exit 2`)
- persistent store round-trip

Run it:

```bash
bash tests/v3/test_v3.sh
echo "exit=$?"
```

A clean run prints `v3 end-to-end tests passed`.

---

## 14. Further Reading

- `src/c/v3/README.md` — V3 architecture, scaling invariants, store format
- `docs/tutorials/README.md` — executable tutorial library (9 problems)
- `docs/tutorials/graph-coloring.md` — smallest end-to-end CNF example
- `docs/tutorials/university-timetabling.md` — canonical WCNF example
- `docs/TUTORIAL_SUDOKU_AND_WORKFORCE_V3.md` — large worked example
- `docs/MATH.md` — mathematical foundation (continuous relaxation, BAHA,
  persistent homology, prime weighting)
- `AGENTS.md` — contributor / agent guide
