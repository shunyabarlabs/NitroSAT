# NitroSAT Agent Guide

This file is the tactical operating contract for coding agents working in this
repository. Read it before editing source, tests, benchmarks, or documentation.

## What this repository is

NitroSAT is a SAT/MaxSAT research and engineering project with three distinct C
solver paths:

- **V1** is the original in-memory research solver retained for historical,
  regression, and cross-version comparison work.
- **V2** is the feature-rich in-memory heuristic solver. It keeps global clause
  and variable-occurrence structures in RAM and uses aggressive repair.
- **V3** is the large-instance path. It streams a binary clause store, supports
  incremental CNF/WCNF inputs, offers an optional disk-indexed CNF finisher,
  and can use a resource-limited exact CDCL fallback.

Do not treat heuristic failure as proof of UNSAT. Only exact mode may report a
formal UNSAT result, and generated proofs still require independent checking
for formal workflows.

## Repository map

| Path | Purpose |
|---|---|
| `src/c/v1/nitrosat.c` | Original V1 in-memory solver |
| `src/c/v2/nitrosatv2.c` | In-memory V2 solver |
| `src/c/v3/nitrosatv3.c` | Streaming/indexed V3 solver and CLI |
| `src/c/v3/cdcl_internal.h` | Internal exact CDCL implementation |
| `src/c/v3/README.md` | V3 architecture and operating guide |
| `tests/v3/` | V3 regression, exact, and tutorial tests |
| `examples/v3/` | Executable tutorial encoders and decoders |
| `docs/tutorials/` | Searchable tutorial knowledge base |
| `benchmarks/` | Benchmark tooling and results |

Large CNFs, generated solutions, compiled binaries, temporary stores, CSV
results, and personal notes may exist locally without being tracked. Do not add
them to Git unless the task explicitly requires it.

## Version lineage: V1, V2, and V3

The versions are different solver strategies, not interchangeable release
labels. Read the source for the version being changed.

| Property | V1 | V2 | V3 |
|---|---|---|---|
| Source | `src/c/v1/nitrosat.c` | `src/c/v2/nitrosatv2.c` | `src/c/v3/nitrosatv3.c` |
| Primary role | Original in-memory research solver | Refined, aggressive in-memory heuristic | Streaming, incremental, bounded-memory-oriented solver |
| Input | CNF | CNF | CNF, WCNF, persistent V3 stores, incremental `--add` files |
| Clause storage | Full formula in RAM | Full formula and global occurrence data in RAM | Binary stream on disk; bounded active state in RAM |
| Optimizer lineage | NADAM-style continuous relaxation | WAdam-style amplitude/phase momentum | Batched streaming Adam-style optimizer |
| Temperature/LR | Earlier fixed/ζ-gain behaviour | Quasi-periodic learning rate and evolving ζ-driven β | Streaming learning-rate decay by epoch |
| Unsatisfied tracking | Earlier full-scan-oriented tracking | Incremental unsatisfied list/positions plus stronger validation | Full streamed verification and bounded active cache |
| Local repair | Topological/adelic and stochastic fallback | Strong in-memory repair using variable-to-clause CSR | Streamed finisher or optional disk-indexed WalkSAT finisher |
| Exact mode | Historical proof machinery | Historical proof machinery | Resource-limited CDCL with atomic CNF DRAT output |
| Best fit | Historical comparisons and hard stochastic cases | Formulas that fit RAM and benefit from aggressive repair | Huge, incremental, WCNF, or memory-constrained workloads |

The detailed historical V1→V2 design narrative is in `CHANGELOG.md`. It
documents the intended optimizer, annealing, topology, and tracking changes.
However, historical prose can drift from current source. For example, the
changelog describes WalkSAT removal during the V2 transition, while the current
V2 source contains a `baha_walksat()` finisher. Treat current code and tests as
authoritative for present behaviour; treat the changelog as release history.

When changing a version:

- Do not copy V2 in-memory structures into V3 without an explicit memory-bound
  design.
- Do not “modernize” V1 merely because V2 has a newer implementation; V1 is a
  historical and comparative solver path.
- Do not silently move CLI or JSON semantics between versions.
- If a cross-version benchmark is claimed, rebuild every version from its own
  source and run the same input, hardware, compiler flags, timeout, and
  verifier.
- Record which binary produced every result. `nitrosat`, `nitrosatv2`, and
  `nitrosatv3` are not aliases for the same algorithm.

## Benchmark map and source-of-truth rules

Benchmark material has several layers:

| Path | Role | Agent treatment |
|---|---|---|
| `README.md` | Public headline results | Edit only with reproducible evidence and caveats |
| `benchmarks/README.md` | Historical benchmark heritage and consolidated narrative | Useful context; verify links and claims before repeating |
| `benchmarks/full_benchmark_results.csv` | Tracked historical raw telemetry | Preserve; do not rewrite without an explicit benchmark migration |
| `benchmarks/run_nitrosat_bench.py` | Current generic CNF runner | Review solver flags and output schema before use |
| `benchmarks/cnf_suite/` | Local/generated benchmark inputs | May be untracked; do not assume CI or redistribution availability |
| `benchmarks/*_results.csv` | Generated run outputs | Usually local artifacts; do not stage by default |
| `benchmarks/*.summary.txt` | Generated summaries | Derived artifacts, not primary evidence |
| `tests/v3/` | Correctness/regression fixtures | Tests correctness; not a performance benchmark suite |

Important caveats:

- `benchmarks/README.md` consolidates multiple historical eras. Some referenced
  legacy files may no longer exist at the repository root. A missing source
  link means the claim needs archaeological verification, not silent reuse.
- Historical tables may mix C, LuaJIT, V1, V2, API, different machines, and
  different dates. Never aggregate them as one controlled experiment.
- Generated CSV filenames such as `latest`, `expanded`, `heuristic`, or `exact`
  do not by themselves establish provenance. Capture the command, commit,
  solver hash, input hash, host, and timeout.
- A benchmark runner parsing JSON must understand the selected version's
  schema. V2 exposes `satisfaction_rate`; V3 commonly requires computing it
  from `satisfied / clauses`.
- Return code `1` may represent a valid partial or UNSAT result rather than a
  process failure. Judge the structured result and verifier output together.

### Reproducing the local CNF suite

Build the intended solver first, then use an explicit output path:

```bash
make -C src/c v2
python3 benchmarks/run_nitrosat_bench.py \
  --cnf-dir benchmarks/cnf_suite \
  --solver src/c/v2/nitrosatv2 \
  --timeout 300 \
  --output /tmp/nitrosat-v2-results.csv
```

For V3, confirm the runner supports the desired mode and JSON fields before
publishing results. Prefer a dedicated command manifest that records indexed,
exact, memory, and conflict limits. Never pass V2-only flags to V3 or assume
`--exact` has the same implementation across versions.

### Evidence hierarchy

Use this order when evaluating a claim:

1. Independently verified assignment or externally checked proof.
2. Raw per-instance result with input hash, solver commit, command, and host.
3. Reproducible generated CSV and summary.
4. Tracked benchmark table with documented provenance.
5. Historical prose or marketing copy.

Do not promote a lower-level claim to a higher confidence level without
reproduction.

## Build

Build all C solvers:

```bash
make -C src/c all
```

Build only V2 or V3:

```bash
make -C src/c v2
make -C src/c v3
```

Strict V3 diagnostic build:

```bash
gcc -O2 -std=c99 -Wall -Wextra -Werror -Wpedantic \
  -Wconversion -Wshadow src/c/v3/nitrosatv3.c -lm \
  -o /tmp/nitrosatv3-strict
```

Sanitizer build:

```bash
gcc -O1 -g -std=c99 -Wall -Wextra -Werror \
  -fsanitize=address,undefined src/c/v3/nitrosatv3.c -lm \
  -o /tmp/nitrosatv3-asan
```

Write diagnostic binaries to `/tmp` unless the repository build explicitly
expects an artifact under `src/c/`.

## Safe test commands

Run the V3 regression suite:

```bash
make -C src/c test-v3
```

Run exact-mode and proof tests:

```bash
python3 tests/v3/test_exact.py
```

Execute every tutorial end to end:

```bash
bash tests/v3/test_tutorials.sh src/c/v3/nitrosatv3
```

Check Python syntax:

```bash
python3 -m py_compile examples/v3/*.py
```

Remove generated `__pycache__` directories before staging files. Do not treat a
timeout as a passing test. Do not reuse a stale solver binary when validating
source changes.

## Files and areas requiring explicit scope

- Do not modify `src/c/v2/nitrosatv2.c` while working on V3 unless the user
  explicitly requests a V2 change.
- Do not commit compiled binaries such as `src/c/v2/nitrosatv2` or
  `src/c/v3/nitrosatv3`.
- Do not commit large generated CNF/WCNF files, `.sol` files, temporary stores,
  benchmark CSVs, or local notes without explicit approval.
- Do not rewrite benchmark history or published result tables merely to make a
  new implementation look better.
- Preserve unrelated local and untracked files. This workspace may be dirty.

Before committing, inspect the exact staged set:

```bash
git diff --cached --check
git diff --cached --stat
git status --short
```

## V3 architectural invariants

### Streaming memory model

- The default V3 path must not load the complete clause set or a global
  variable-to-clause index into RAM.
- Memory should remain proportional to variables, configured batch limits, and
  configured active-cache limits—not total clause count.
- The optional indexed finisher may create a temporary disk-backed occurrence
  index. Report its disk/mapping contribution honestly.
- Any new unbounded allocation requires a documented limit and an explicit
  failure or `UNKNOWN` result.

### Store format and transactions

- Preserve the 40-byte V3 store header and version checks.
- Preserve little-endian record encoding.
- Clause records contain literal count, flags, weight, and signed literals.
- `--add` must remain transactional. If any appended file is malformed, restore
  both the original byte length and original header.
- Never leave a partially appended store that looks valid.
- Validate record flags, literal ranges, clause counts, and trailing bytes.
- Explicit store paths are persistent artifacts; temporary stores must be
  cleaned up on success and failure.

### SAT, WCNF, and result semantics

- Verify every final assignment by streaming the complete store.
- Preserve the best globally verified assignment across phases.
- For WCNF, hard-constraint feasibility dominates soft-cost improvement.
- `solved=true` means every clause is satisfied.
- `feasible=true` means every hard clause is satisfied.
- A feasible WCNF result may legitimately have `solved=false`.
- Do not claim WCNF soft-cost optimality unless an independent exact optimizer
  proves it.
- Use `UNKNOWN`/`LIMIT` when exact resource limits are exhausted.

### Exact solving and proofs

- Exact mode must obey clause, memory, and conflict limits.
- Distinguish `SAT`, `UNSAT`, `LIMIT`, `ERROR`, and `NOT_RUN` internally.
- Independently verify a CDCL SAT assignment before accepting it.
- DRAT output must be atomic: write a temporary file and rename only after a
  successful UNSAT derivation.
- Never overwrite a valid existing proof on SAT, input failure, or limit exit.
- WCNF proof requests are currently rejected because the hard-clause projection
  is not directly checkable against the original WCNF input.
- Proof tests must validate derivations, not merely check that a file ends in
  an empty clause.

### Indexed finisher

- The indexed finisher currently supports CNF only.
- Its temporary occurrence index must be cleaned up on every return path.
- Make/break flips are heuristic; only complete streamed checkpoints establish
  improvement.
- Restore the best verified assignment after local search.
- Keep flip and checkpoint budgets explicit and report actual flips performed.

## Benchmark discipline

- Do not claim universal SAT or MaxSAT superiority.
- Identify planted instances as planted and state how they were generated.
- Separate heuristic satisfaction from exact SAT/UNSAT conclusions.
- Preserve visible failures and partial results; do not report only successful
  seeds.
- For every benchmark, report when available:
  - variables, clauses, and literals;
  - SAT/WCNF type and generation method;
  - seed and solver options;
  - satisfaction or hard/soft result;
  - indexed flips or exact conflicts;
  - solve time and store-build time;
  - peak RSS and temporary index size;
  - independent assignment/proof verification result.
- Reproduce important claims from a clean build.
- Compare against established SAT/MaxSAT solvers before making competitive
  claims.

## Adding executable tutorials

Every tutorial must follow this structure:

```text
Problem
  ↓
Constraint Model
  ↓
CNF / WCNF Generation
  ↓
NitroSAT Solve
  ↓
Decode Solution
  ↓
Visualize Result
  ↓
Complexity Discussion
```

Requirements:

- Provide a runnable encoder or model generator under `examples/v3/`.
- Provide a decoder that translates literals into domain language.
- Include exact commands and expected output in the tutorial.
- Verify hard constraints independently.
- For small WCNF examples, independently enumerate or otherwise prove the
  expected preference optimum.
- Add the tutorial to `docs/tutorials/README.md`.
- Add it to `tests/v3/test_tutorials.sh` so CI executes it.
- Keep examples small enough for CI and clearly distinguish educational models
  from production-ready formulations.

## Documentation style

- Prefer commands that have been executed successfully in this repository.
- State prerequisites and working directory assumptions.
- Explain variable mappings and why each clause family is correct.
- Show result-status checks before decoding an assignment.
- Include at least one decoded or visualized result.
- Discuss variable/clause growth and production encoding alternatives.
- Use precise terms: heuristic, feasible, solved, exact, verified, planted, and
  optimal are not interchangeable.

## Change workflow

1. Inspect current status and relevant source before editing.
2. Preserve unrelated changes and untracked files.
3. Make the smallest coherent change.
4. Build from source with strict warnings.
5. Run targeted tests, then the relevant regression suites.
6. Independently verify important assignments or proofs.
7. Run formatting and whitespace checks.
8. Stage explicit paths, never broad generated directories.
9. Review the staged diff before committing.
10. Push only when the user explicitly requests it.
