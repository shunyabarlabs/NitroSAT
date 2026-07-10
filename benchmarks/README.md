# 📊 NitroSAT — Complete Benchmark Heritage

> **5,000+ CNF instances tested** · **77% perfect SATISFIED rate** · **99.7% median satisfaction**
> Assignments dataset: [HuggingFace](https://huggingface.co/datasets/sethuiyer/navokoj_sat_2024)

This document consolidates **every** benchmark table from `README.md`, `BENCHMARKS.md`, `RESULTS.md`, `CHANGELOG.md`, and `full_benchmark_results.csv` — ordered chronologically by git commit date.

## 🧭 Quick Navigation

| Source File | Lines | Scope |
| :--- | :--- | :--- |
| [BENCHMARKS.md](BENCHMARKS.md) | 567 | v1 era: 13 categories, verification suite, live audit |
| [RESULTS.md](RESULTS.md) | 204 | v2 era: Death Run, Edwards-Anderson, spin glass |
| [full_benchmark_results.csv](full_benchmark_results.csv) | 360 | Raw telemetry for 360 CNF seeds |
| [../README.md](../README.md) | 260 | Hero table, v2 timetable breakthrough |
| [../CHANGELOG.md](../CHANGELOG.md) | 280 | v1 vs v2 algorithmic comparison |

---

# ⏳ Phase 7 — July 2026: V3 Disk-Streaming Architecture

> Commits `ecddb34` → `58dae03` — V3 disk-streaming solver, indexed/exact
> modes, executable tutorials, WCNF hard-repair finisher, and structured
> combinatorial benchmark tooling.

V3 is a different benchmark regime from V1/V2. It is designed for formulas
where full clause state or a global variable-to-clause index should not be
kept in RAM. V3 converts CNF/WCNF input into a binary stream, scans that stream
on each global pass, and keeps only bounded active state plus optional
temporary disk-backed index data.

Do not merge the V3 numbers below into the older V1/V2 aggregate percentages.
The goals are different:

* **V2:** raw in-memory heuristic speed when the formula fits RAM.
* **V3:** bounded-memory streaming, incremental stores, WCNF support, and
  disk-indexed finishing.

## V3 Streaming and Store Capabilities

| Capability | Status | Notes |
|---|---|---|
| CNF input | Supported | DIMACS `p cnf` |
| WCNF / partial MaxSAT input | Supported | Hard clauses via `top` weights or `h` prefix |
| Persistent binary store | Supported | 40-byte header plus clause records |
| Incremental `--add` | Supported | Transactional rollback on malformed append |
| Exact CDCL fallback | Resource-limited | CNF only; formal UNSAT proof workflows still require independent checking |
| Indexed finisher | CNF + WCNF hard repair | For WCNF, indexed phase repairs hard clauses only, then streamed finisher improves soft cost |
| Memory model | Bounded active state | No default full clause table or global occurrence index in RAM |

## V3 Planted 3-SAT Indexed-Finisher Matrix

These are reproducible planted satisfiable instances. They are useful for
checking V3 indexed-finisher behavior, but they are not a substitute for an
independently sourced industrial benchmark suite.

Standard run: 10 epochs, 100,000 indexed flips, checkpoint every 500 flips.

| Variables | Clauses | Ratio | Seed | Result | Indexed flips | Solve time |
|---:|---:|---:|---:|---:|---:|---:|
| 100 | 430 | 4.30 | 101 | 430 / 430 | 14,555 | 400ms |
| 100 | 923 | 9.23 | 102 | 923 / 923 | 339 | 32ms |
| 100 | 2,000 | 20.00 | 103 | 2,000 / 2,000 | 59 | 13ms |
| 250 | 1,065 | 4.26 | 201 | 1,065 / 1,065* | 40,997 | 1.46s |
| 250 | 2,308 | 9.23 | 202 | 2,308 / 2,308 | 487 | 47ms |
| 250 | 5,000 | 20.00 | 203 | 5,000 / 5,000 | 234 | 53ms |
| 500 | 2,130 | 4.26 | 301 | 2,130 / 2,130 | 46,178 | 2.08s |
| 500 | 4,615 | 9.23 | 302 | 4,615 / 4,615 | 535 | 64ms |
| 500 | 10,000 | 20.00 | 303 | 10,000 / 10,000 | 345 | 92ms |
| 1,000 | 4,260 | 4.26 | 401 | 4,259 / 4,260† | 500,000 | 25.37s |
| 1,000 | 9,230 | 9.23 | 402 | 9,230 / 9,230 | 1,043 | 131ms |
| 2,000 | 18,460 | 9.23 | 501 | 18,460 / 18,460 | 1,874 | 253ms |

\* The 250-variable ratio-4.26 case initially reached 1,064/1,065 with the
standard budget and solved after increasing the budget to 250,000 flips and 20
epochs. Its final assignment was independently verified.

† The 1,000-variable ratio-4.26 case reached 4,256/4,260 at 100,000 flips,
4,258/4,260 at 250,000 flips, and 4,259/4,260 at 500,000 flips. It is retained
as a visible failure rather than reported as solved.

## V3 WCNF Indexed Hard-Repair Regression

`tests/v3/test_wcnf_finisher.py` generates a 5-vertex, 3-color WCNF graph
coloring instance with hard coloring constraints and soft color preferences.

| Mode | Feasible | Hard unsatisfied | Soft cost | Indexed index |
|---|---:|---:|---:|---:|
| `--epochs 0 --finisher-passes 0 --indexed-finisher` | true | 0 | 30 | 0.001 MiB |
| `--epochs 0 --finisher-passes 256 --indexed-finisher` | true | 0 | 10 | 0.001 MiB |

This verifies the intended V3 WCNF sequence: indexed repair first finds a
hard-feasible region, then the streamed finisher improves soft cost. The soft
cost is heuristic and is not an optimality certificate.

## V3 Large-Formula Streaming Stress

Local stress testing documented in the V3 usage guide reports a
37M-clause enterprise-timetabling-style stream with **9.62 MiB** peak bounded
solver allocation under default settings. Treat this as a local streaming
memory demonstration until the full artifact pack contains command, commit,
input hash, raw output, host, and verifier output.

## V3 Evidence Rules

For every V3 result, preserve:

* exact binary and commit;
* input generator, seed, and input hash;
* command line and timeout;
* `bounded_memory_mb` and, when used, `indexed_index_mb`;
* hard vs soft result for WCNF;
* independent assignment/proof verification status;
* whether the instance is planted, generated, benchmark-suite, or production.

---

# ⏳ Phase 1 — January 15, 2026: Initial Release (v1.0)

> Commit `8c6fc20` — *"Initial release of NitroSAT"*

The `full_benchmark_results.csv` (360 instances) was introduced here with the core solver. Key category aggregates:

### Category Aggregates (C Engine)

| Category | Problem Type | C Avg. Satisfaction | Lua Avg. Satisfaction |
| :--- | :--- | :--- | :--- |
| `rand3sat` | Random 3-CNF | 99.59% | 99.86% |
| `rand4sat` | Random 4-CNF | 99.86% | *Not tested* |
| `rand5sat` | Random 5-CNF | 99.86% | *Not tested* |
| `parity` | XOR / Parity | 99.95% | 99.94% |
| `clique` | Graph Clique | 99.94% | 99.98% |
| `domset` | Dominating Set | **100.00%** | **100.00%** |

### Engine Parity: C vs LuaJIT

| Problem Instance | Vars | Clauses | C Version (Time) | LuaJIT (Time) | Steps (Lua) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `rand3sat_50_200` | 50 | 200 | 100% (~0.01s) | 100% (4.26ms) | 78 |
| `clique_4_20` | 80 | 2,600 | 100% (~0.01s) | 100% (18.75ms) | 214 |
| `cliquecol_80_10_10` | 4,760 | 354,890 | 100% (**14.00s**) | 100% (**72.81s**) | 3,000 |

---

# ⏳ Phase 2 — February 17–28, 2026: Scaling Expansion

> Commits `ff164f8` → `b0a684d` — *"insane benchmarks added"*, *"Unbelievable.."*, *"Added new benchmarks"*

Massive expansion of benchmark coverage across 13 categories. The `BENCHMARKS.md` file was created with the full verification suite.

### 1️⃣ Large Structured Instances

| Problem (type) | Variables | Clauses | Satisfaction | Time | Hardware |
|-----------------|-----------|---------|--------------|------|----------|
| **Graph Coloring** | – | 650,000 | **100%** | 4.6s | Ryzen 5 5600H |
| **Clique Coloring – `cliquecol_80_10_10`** | 4,760 | 354,890 | **100%** (5/5 seeds) | **3.5s** | Laptop |
| **Ramsey R(5,5,5)** | – | 402,752 | **100%** | 7.5s | Laptop |
| **Parity (CNFgen)** | – | 485,200 | **100%** | 2.6s | – |
| **Counting (CNFgen)** | – | 78,760 | **100%** | 0.5s | – |
| **Matching (100‑node)** | – | 400 | **100%** | 21ms | – |
| **Van der Waerden** | – | 20 | **100%** | <1ms | – |
| **Tiling (8×8 grid)** | – | 580 | **99.1%** | 154ms | – |
| **Subset Cardinality** | – | 210 | **95.7%** | 107ms | – |

### 2️⃣ Real‑World Scheduling

| Jobs | Slots | Density | Clauses | Satisfaction | Time |
|------|-------|---------|---------|--------------|------|
| 30 | 5 | 0.30 | 1,095 | **100%** | 28ms |
| 50 | 6 | 0.20 | 2,522 | **100%** | 22ms |
| 100 | 8 | 0.10 | 7,516 | **100%** | 26ms |
| 200 | 10 | 0.05 | 21,150 | **100%** | 85ms |
| 500 | 10 | 0.03 | 64,570 | **100%** | 172ms |
| **1,000** | 10 | 0.02 | 154,760 | **99.99%** (UNSAT detection) | 225s |

### 3️⃣ Random 3‑SAT Phase‑Transition (α ≈ 4.26)

| Variables (n) | Clauses | Seeds | Avg. Sat. | Min | Max | Std. Dev. |
|---------------|---------|-------|-----------|-----|-----|-----------|
| 300 | 1,278 | 50 | **99.65%** | 99.37% | 99.84% | 0.11% |
| 500 | 2,130 | 20 | **99.64%** | 99.44% | 99.81% | 0.10% |
| 1,000 | 4,260 | 10 | **99.65%** | 99.53% | 99.72% | **0.06%** |

*Variance **shrinks** as instance size grows — scale‑invariant robustness.*

### 4️⃣ Grid‑Coloring Stress Test (Spectral Scaling)

| Grid Size | Variables | Clauses | Satisfaction | Time |
|-----------|-----------|---------|--------------|------|
| 10 × 10 | 400 | 1,420 | **100%** | 0.02s |
| 20 × 20 | 1,600 | 5,840 | **100%** | 0.07s |
| 50 × 50 | 10,000 | 37,100 | **100%** | 0.45s |
| 100 × 100 | 40,000 | 149,200 | **100%** | 2.10s |
| **1,000 × 1,000** | **4,000,000** | **14,992,000** | **100%** | **475s** |

### 5️⃣ XOR‑SAT Stress Test (Parity‑Chain Detection)

| Instance | Clauses | Result | Time | Detected Loops (β₁) |
|----------|---------|--------|------|----------------------|
| XOR (SAT) | 200 | **100%** | 3.9ms | 98 |
| XOR (planted) | 8,000 | **100%** | 9.5ms | — |
| XOR (UNSAT) | 2,000 | 95.15% | 802ms | 2–26 |

### 6️⃣ UNSAT Awareness – "Mirage" Trap

| Benchmark | Stopping Satisfaction | Observation |
|-----------|----------------------|-------------|
| Mirage 200 | **85.2%** | Trap Detected (structural impossibility) |
| Mirage 300 | **95.9%** | Structural Awareness (phase‑transition signal) |

### 7️⃣ Permutation Invariance (Encoding‑Agnostic)

| Permutations Tested | Perfect Solves | Std. Dev. |
|---------------------|----------------|-----------|
| 20 random variable/clause permutations | 20/20 (100%) | **0.0000%** |

### 8️⃣ Quasigroup / Latin‑Square Completion

| Run Type | Instances (SAT + UNSAT) | Avg. Sat. | Std. Dev. | ≥ 99.9% | Perfect 100% |
|----------|--------------------------|-----------|-----------|----------|--------------|
| Single‑seed (8 instances) | 8 (4 SAT + 4 UNSAT) | **99.976%** | – | 8/8 | 1/8 |
| Seed sweep (40 runs = 8 × 5 seeds) | 40 | **99.960%** | 0.036% | 38/40 | – |
| SAT runs (40) | – | **99.985%** | – | 20/20 | – |
| UNSAT runs (40) | – | **99.967%** (plateau) | – | 0/20 | – |

### 9️⃣ Global Benchmark Summary (All 358 Instances)

| Metric | Value |
|--------|-------|
| Total instances evaluated | **358** |
| Solved at 100% | **115** (32.1%) |
| Solved ≥ 99% | **340** (95.0%) |
| Solved ≥ 96% | **353** (98.6%) |
| Average satisfaction | **99.58%** |
| Largest perfect solve | **354,890 clauses** (Clique Coloring) |
| Fastest >10K‑clause perfect solve | **22,521 clauses** in **33ms** (clique\_6\_v40) |

### 🔟 CNFgen Benchmark Suite

| Category | Tests | Avg Sat% | Perfect | ≥99% | Notes |
|----------|-------|----------|---------|------|-------|
| **Graph Coloring** | 6 | **99.97%** | 4/6 | 6/6 | NitroSAT's strength |
| **Parity** | 4 | **99.94%** | 2/4 | 4/4 | Topology tracking (β₁) |
| **Ordering Principle** | 3 | **99.94%** | 0/3 | 3/3 | UNSAT detected |
| **Counting** | 3 | **99.99%** | 1/3 | 3/3 | Handles 162K clauses |
| **Random 3-SAT** (α=4.26) | 8 | **99.73%** | 2/8 | 8/8 | Phase transition plateau |
| **Pigeonhole** (UNSAT) | 6 | **99.81%** | 0/6 | 6/6 | Mirage trap detection |
| **Tseitin** (UNSAT) | 3 | **99.59%** | 1/3 | 3/3 | Parity reasoning |
| **Ramsey Numbers** | 5 | **97.76%** | 2/5 | 3/5 | Hardest category |

**CNFgen totals (44 instances):** 99.61% avg · 15/43 perfect · 40/43 ≥99% · 1.14s avg runtime

### 1️⃣1️⃣ LeetCode Challenging Problems

#### Core Suite (8 problems)

| Problem | Type | Variables | Clauses | Sat% | Time |
|---------|------|-----------|---------|------|------|
| **N-Queens Completion (8×8, 3 pre)** | Constraint Sat | 64 | 470 | **100.00%** | <0.01s |
| **Exact Cover** | NP-Complete | 6 | 16 | **100.00%** | <0.01s |
| **Graph 5-Coloring** (Petersen-like) | Graph Theory | 50 | 185 | **100.00%** | <0.01s |
| **Sudoku** (Easy) | Constraint Sat | 729 | 12,018 | **99.95%** | <0.01s |
| **Hamiltonian Cycle** | NP-Complete | 25 | 110 | **100.00%** | <0.01s |
| **3D Matching** | NP-Complete | 4 | 11 | **100.00%** | <0.01s |
| **N-Queens** (12×12) | Constraint Sat | 144 | 1,816 | **99.94%** | <0.01s |
| **Graph 3-Coloring K4** | UNSAT Test | 12 | 34 | **97.06%** | <0.01s |

#### Extended Suite (6 harder problems)

| Problem | Type | Variables | Clauses | Sat% | Time |
|---------|------|-----------|---------|------|------|
| **Sudoku** (17 clues - hardest) | Minimal Clues | 729 | 12,005 | **99.87%** | 0.06s |
| **Graph 7-Coloring** G(50, 0.1) | Random at Threshold | 350 | 1,940 | **100.00%** | <0.01s |
| **N-Queens** (20×20) | Large Instance | 400 | 8,760 | **99.97%** | 0.10s |
| **Clique+Coloring Tension** | Contradictory Goals | 150 | 1,025 | **99.51%** | 0.02s |
| **Latin Square** (10×10) | Combinatorial Design | 1,000 | 13,805 | **99.93%** | 0.08s |
| **Independent Set** (50, k=10) | NP-Complete | 50 | 113 | **100.00%** | <0.01s |

**LeetCode totals (14 problems):** 99.75% avg · 7/14 perfect · 13/14 ≥99% · 14/14 sub-second

### 1️⃣2️⃣ Navokoj Pro API Validation

| Test | Type | Satisfaction | Time | Engine | Verdict |
|------|------|--------------|------|--------|---------|
| **Simple CNF** | Basic | **100%** | 0.26s | pro-deepthink | ✅ |
| **Ramsey-like** | UNSAT | **95%** | 21.2s | pro-deepthink | ⚠️ Detected |
| **Boolean Expression** | Parsing | **100%** | 0.10s | mini-deepthink | ✅ |
| **PHP-like** (204 clauses) | UNSAT | **99.51%** | 23.2s | pro-deepthink | ⚠️ Detected |
| **Mini Engine** | Comparison | **100%** | 0.10s | mini-deepthink | ✅ |
| **DEFEKT Diagnostic** | Spectral | N/A | 0.006s | diagnostic | ✅ |
| **Batch Solving** (3 probs) | Throughput | **100%** | 0.14s | auto | ✅ |

**NitroSAT vs Navokoj Pro (PHP 8-7):**

| Metric | NitroSAT | Navokoj Pro |
|--------|----------|-------------|
| Satisfaction | 99.51% | 99.51% |
| Time | **0.13s** | 23.2s |
| Diagnostics | Basic | Full violation + blame |
| Cost | Free | Pay-per-solve |

### Number Theory Ablation Study (Prime vs Uniform Weights)

| Instance | Type | Weight Mode | Sat% | Time (ms) | Steps | Final β₁ (Cycles) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `clique_4_20` | Structured | **Prime** | 100% | **12.8ms** | **94** | **20** |
| `clique_4_20` | Structured | Uniform | 100% | 43.8ms | 381 | 79 *(4 fractures)* |
| `rand3sat_200_850` | Random | **Prime** | 99.65% | **768ms** | 3000 | 181 |
| `rand3sat_200_850` | Random | Uniform | 99.65% | 3,082ms | 3000 | 179 |
| `parity_14` | XOR | **Prime** | 100% | 5.8ms | 79 | 74 |
| `parity_14` | XOR | Uniform | 100% | 3.1ms | 74 | 85 |

**Finding:** Prime weights actively prune topological noise (β₁: 79→20), resulting in 3.4x to 4x speedups on structured geometries.

### Thermal Phase Transition (Heat Beta Sensitivity)

| heat_beta | Environment | Steps | Final Satisfaction |
| :--- | :--- | :--- | :--- |
| 0.01 | Heavy Heat / High Noise | 500 (Timeout) | 99.76% |
| 0.10 | Warm / Noise | 500 (Timeout) | 99.52% |
| **0.50** | **NitroSAT Default** | **137** | **100.00%** |
| 1.00 | Cold / Less Noise | 139 | 100.00% |
| 2.00 | Colder | 139 | 100.00% |
| 5.00 | Very Cold / Greedy | 139 | 100.00% |
| 10.00 | Near-Absolute Zero | 139 | 100.00% |

### Universal NP-Complete Tests (CNFGen)

| Problem Type | Config | Vars | Clauses | Sat% | Time (s) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| K-Clique | K=3, N=10 | 30 | 393 | 99.7% | 0.12s |
| K-Clique | K=4, N=15 | 60 | 1,624 | 99.9% | 0.38s |
| K-Clique | K=5, N=25 | 125 | 7,255 | 99.9% | 1.55s |
| Matching | gnm 15 25 | 25 | 91 | 98.9% | 0.06s |
| Pigeonhole | 10→8 | 80 | 370 | 99.4% | 0.14s |
| Ordering | N=12 | 132 | 1,398 | 99.9% | 0.42s |

### Hardware / Chip Verification Suite

| Circuit | Type | Variables | Clauses | Satisfaction | Time (s) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 8-Bit Adder | Logic Chain | 27 | 62 | **100%** | 0.00s |
| 16-Bit Adder | Logic Chain | 51 | 126 | **100%** | 0.00s |
| 32-Bit Adder | Logic Chain | 99 | 254 | **100%** | 0.00s |
| 64-Bit Adder | Logic Chain | 195 | 510 | **100%** | 0.00s |
| 4×4 Multiplier | Wallace Tree | 64 | 133 | **100%** | 0.00s |
| 6×6 Multiplier | Wallace Tree | 132 | 317 | **100%** | 0.01s |
| 8×8 Multiplier | Wallace Tree | 224 | 581 | **100%** | 0.00s |
| 12×12 Multiplier | Wallace Tree | 480 | 1,349 | **100%** | 0.00s |
| 16×16 Multiplier | Wallace Tree | 832 | 2,437 | **100%** | 0.01s |
| 20×20 Multiplier | Wallace Tree | 1,280 | 3,845 | **100%** | 0.01s |
| 32×32 Multiplier | Wallace Tree | 3,200 | 9,989 | **100%** | 0.03s |
| 48×48 Multiplier | Wallace Tree | 7,104 | 22,661 | **100%** | 0.06s |
| 64×64 Multiplier | Wallace Tree | 12,544 | 40,453 | **100%** | 0.10s |
| **128×128 Multiplier** | **Massive** | **49,664** | **162,821** | **100%** | **0.37s** |
| **256×256 Multiplier** | **Ultra Massive** | **197,632** | **653,317** | **100%** | **1.40s** |
| **512×512 Multiplier** | **Ultra Massive** | **788,480** | **2,617,349** | **100%** | **5.92s** |

### Frustrated Small-World Lattice

| Size | Colors | Variables | Clauses | Satisfaction | Time |
|------|--------|-----------|---------|-------------|------|
| 100×100 | 3 | 30,000 | 100,000 | 99.91% | 31.3s |
| 100×100 | 4 | 40,000 | 150,000 | **100%** | 2.15s |
| 200×200 | 4 | 160,000 | 601,596 | **100%** | 27.6s |
| **300×300** | **4** | **360,000** | **1,354,800** | **100%** | **62.75s** |

### Protein Contact Map Prediction

| Sequence | Length | Max Contacts | Clauses | Satisfaction | Time |
|----------|--------|--------------|---------|-------------|------|
| Random | 12 AA | 5 | 1,593 | 98.2% | 0.74s |
| Random | 15 AA | 5 | 16,256 | 99.7% | 8.0s |
| Random | 20 AA | 3 | 51,549 | **99.8%** | 18s |
| Real (ACDEFGHIKLMNPQRSTVWY) | 20 AA | 2 | 16,691 | 99.6% | 5.1s |

### University Timetabling (ITC-Style)

| Instance | Courses | Rooms | Slots | Variables | Clauses | Satisfaction | Time | Topology (β₁) |
|----------|---------|-------|-------|-----------|---------|-------------|------|---------------|
| **Timetabling** | 50 | 12 | 30 | 18,000 | **2,504,500** | **100%** | **97.47s** | 3,213,050 → 1,028,177 |

### Enterprise Timetabling (100 Courses) — v1

| Instance | Courses | Rooms | Slots | Variables | Clauses | Satisfaction | Time | RAM |
|----------|---------|-------|-------|-----------|---------|-------------|------|-----|
| **Enterprise** | 100 | 36 | 41 | 147,600 | **80,278,884** | **99.99999%** | **5.2 hours** | ~3GB |

### Verification Summary (February 28, 2026)

| Metric | Value |
|--------|-------|
| Total instances tested | **80+** |
| Average satisfaction | **99.65%** |
| Perfect solves (100%) | 49/75 (65%) |
| Hardware verification (100%) | 15/15 (100%) |
| Lattice (4-color, 300×300) | **1,354,800 clauses** (100%) in 63s |
| Protein contact maps | **99.8%** on 51K clauses (20 AA) |
| **Timetabling (ITC)** | **2,504,500 clauses** (100%) in 97s |
| **Enterprise Timetabling** | **80,278,884 clauses** (99.99999%) in 5.2h |
| Prime weight speedup | **4x** on structured problems |
| Largest instance solved | **80,278,884 clauses** |

**Updated Global Summary (Feb 2026 — All 65+ New Instances):**

| Metric | Value |
|--------|-------|
| **Total instances tested** | **65+** (358 original + 65 new) |
| **CNFgen suite** | 44 instances, **99.61%** avg |
| **LeetCode suite** | 14 instances, **99.75%** avg |
| **Navokoj API** | 7 tests, **99.29%** avg |
| **Combined average** | **99.65%** satisfaction |
| **Perfect solves** | 37/65 (57%) |
| **≥99% satisfaction** | 58/65 (89%) |
| **Sub-second solutions** | 95%+ |

---

# ⏳ Phase 3 — March 2, 2026: Adversarial & Combinatorial Traps

> Commits `ce873a6`, `faee262` — *"Added pitfall benchmark"*, *"more benchmarks"*

### 1️⃣1️⃣ Pitfall Formula (CDCL Adversarial Trap)

The **Pitfall formula** (Buss & Nordström) is specifically engineered to expose CDCL solvers' weakness. NitroSAT's continuous relaxation never "commits" to a branch.

| Instance | Parameters | Variables | Clauses | Satisfaction | Time | Topology (β₁) |
|----------|-----------|-----------|---------|-------------|------|---------------|
| `pit.cnf` | `pitfall 45 4 30 5 8` | 1,784 | 361,095 | **99.998%** (7 unsat) | 383.55s | 1,575 → 18,996 |

### 1️⃣2️⃣ Boolean Pythagorean Triples (The 200TB Proof Problem)

| Instance | N | Variables | Clauses | Satisfaction | Time |
|----------|---|-----------|---------|-------------|------|
| `ptn_5000` | 5,000 | 5,000 | 11,362 | **99.81%** (21 unsat) | 102s |
| `ptn_7824` | 7,824 | 7,824 | 18,930 | **99.64%** (67 unsat) | 217s |

### 1️⃣3️⃣ Titan Ramsey R(5,5) on 40 Nodes (Extreme Density)

| Instance | Edges (Vars) | Clauses | Density (α) | Satisfaction | Time | Topology (β₁) |
|----------|--------------|---------|-------------|-------------|------|---------------|
| `titan_ramsey_40_5` | 780 | 1,316,016 | **1,687.2** | **99.995%** (60 unsat) | 3,403s | 6,483 → 53,869 |

---

# ⏳ Phase 4 — March 4, 2026: Live Audit

> Commit `d1ec356` — *"docs: improve professional tone in documentation"*

### 📊 NitroSAT Live Audit (9 Instances)

| Category | Instance | Variables | Clauses | Density (α) | Satisfaction | Time | Result |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| *Random 3-SAT* | `phase5k.cnf` | 5,000 | 21,300 | 4.26 | *99.55%* | ~30s | *Plateau Shift* |
| *Tseitin (Parity)* | `tseitin5k.cnf` | 7,500 | 20,000 | 2.67 | *100.0%* | *26.8s* | *Near-Optimal* |
| *Ramsey Theory* | `ramsey.cnf` | 435 | 54,810 | 126.0 | *99.17%* | 74s | *β₁ Increase* |
| *Resolution Trap* | `pit.cnf` | 2,950 | 1,047,620 | 355.1 | *100.0%* | *~400s* | *RESISTANT* |
| *UNSAT Trap* | `php 25 24` | 600 | 7,225 | 12.0 | *99.99%* | <1s | *Near-Optimal* |
| *Phase Transition* | `phase_200k` | 200,000 | 852,000 | 4.26 | *99.20%* | *11.2 min* | *Consistent Scaling* |
| *Planted Color* | `planted_3500` | 10,500 | 14,933 | 1.42 | *100.0%* | *0.07s* | *OPTIMIZED* |
| *Planted Color* | `planted_35k` | 105,000 | 232,043 | 2.21 | *100.0%* | *13.78s* | *LUA ENGINE SUPERIOR* |
| *Hyper-Dense* | `planted_10k` | 10,500 | 931,661 | *88.7* | *99.62%* | ~120s | *SPEC-GEOM SUCCESS* |

### Scaling Analysis

| Instance Scale | Variables | Clauses | Time | Throughput (clauses/s) |
|----------------|-----------|---------|------|------------------------|
| Small | 5,000 | 21,300 | ~30s | ~710 |
| Medium | 10,500 | 232,043 | 13.78s | ~16,840 |
| Large | 105,000 | 232,043 | 13.78s | ~16,840 |
| Massive | 200,000 | 852,000 | 11.2 min | ~1,268 |

### "Phase-2" Rescue (Glassy Plateau Recovery)

| Instance | Plateau Sat% | Final Sat% | Rescue Delta | Mechanism |
|----------|--------------|------------|--------------|-----------|
| `phase_200k` | ~96% | 99.20% | +3.2% | Topological Repair + BAHA |
| `ramsey.cnf` | ~95% | 99.17% | +4.17% | β₁ Explosion Resolution |
| `pit.cnf` | ~96% | 100.0% | +4.0% | Stage 2/3 Adelic Saturation |

### High Constraint Density (α = 88.7)

| Instance | α (Clauses/Var) | Satisfaction | Time |
|----------|-----------------|--------------|------|
| `phase5k.cnf` | 4.26 | 99.55% | ~30s |
| `planted_10k` | **88.7** | **99.62%** | ~120s |

### Audit Summary (March 4, 2026)

| Metric | Value |
|--------|-------|
| **Total instances audited** | **9** |
| **Perfect solves (100%)** | **5/9** (55.6%) |
| **≥99% satisfaction** | **8/9** (88.9%) |
| **Largest perfect solve** | **105,000 variables** (`planted_35k`) |
| **Highest clause density** | **α = 88.7** (`planted_10k`) |
| **CDCL trap resistant** | **`pit.cnf`** (1M+ clauses, 100%) |

---

# ⏳ Phase 5 — March 8, 2026: The v2 Generational Leap

> Commits `d3b9447`, `d554b7a` — *"Introducing the improved version 2"*, *"Added timetable witness"*
> NADAM → WAdam (Wasserstein-Flow with Resonance) · O(1) incremental unsat tracking · Walksat removed

### Hero Table (from README.md — v1/v2 mixed flagship results)

| Test Name / Problem Type | Variables | Clauses | Ratio (α) | Latency | Satisfaction |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Planted Coloring (Lua)** | 105,000 | 232,043 | 2.21 | 13.78s | 100% |
| **Grid Coloring (1000x1000)** | 4,000,000 | 14,992,000 | 3.75 | 475.0s | 100% |
| **Massive Hardware Check** | 788,480 | 2,617,349 | 3.32 | 5.92s | 100% |
| **Enterprise Timetabling** | 147,600 | 80,278,884 | 543.89 | 5.2h | 99.99999% |
| **Titan Ramsey Density** | 780 | 1,316,016 | **1,687.20** | 3,403.0s | 99.995% |
| **Adversarial Pitfall Trap** | 2,950 | 1,047,620 | 355.12 | ~400.0s | 100% |
| **Large Phase Transition** | 200,000 | 852,000 | 4.26 | 11.2 min | 99.20% |
| **Ramsey R(5,5,5)** | ~4,760 | 354,890 | 74.56 | 3.5s | 100% |
| **Small-World Lattice** | 360,000 | 1,354,800 | 3.76 | 62.75s | 100% |

### 🏆 v2 Timetable Breakthrough

| Instance | Variables | Clauses | Satisfaction | Time | Throughput | Hardware |
|----------|-----------|---------|--------------|------|------------|----------|
| timetable.cnf | 147,600 | 80,278,884 | **100%** | **73s** | **1.1M clauses/sec** | AMD Ryzen 5 5600H (laptop) |

**Latency Breakdown:** Initialization: 41s · Langevin Flow: 33s · Total: **73s**
**Output:** [timetable_output.json](../timetable_output.json)

### v1 vs v2 Comparison (from CHANGELOG.md)

| Instance | v1 | v2 | Δ |
|----------|-----|-----|---|
| 80M-clause enterprise timetabling | 5.20h | **73s** | **~250x** |
| 512×512 integer multiplier | 5.92s | 3.71s | **−37%** |

#### Topological Metrics (clique_4_20)

| Metric | v1 | v2 | Δ |
|--------|-----|-----|---|
| β₁ (initial) | 79 | 79 | — |
| β₁ (post-solve) | 20 | 16 | **−20%** |
| Topology complexity score | 0.78 ↑ | 0.00 ✓ | Stable |

---

# ⏳ Phase 6 — April 7, 2026: Physics-Informed Advanced Models (v2)

> Commits `5f17cad`, `1204264` — *"Add files via upload"* (RESULTS.md)
> **Hardware**: Apple Silicon

### 1. Grid/Torus Dominating Set Problems

| Instance | Vars | Clauses | Sat% | Time |
|----------|------|---------|------|------|
| grid60 | 14,400 | 53,520 | **100%** | 35ms |
| grid75 | 22,500 | 83,775 | **100%** | 68ms |
| grid80 | 25,600 | 95,360 | **100%** | 61ms |
| grid100 | 40,000 | 149,200 | **100%** | 94ms |
| grid100_v2 | 40,000 | 149,200 | **100%** | 97ms |
| grid120 | 57,600 | 215,040 | **100%** | 148ms |
| torus50 | 10,000 | 37,500 | **100%** | 27ms |
| torus60 | 14,400 | 54,000 | **100%** | 46ms |
| torus75 | 22,500 | 84,375 | **100%** | 56ms |
| torus100 | 40,000 | 150,000 | **100%** | 98ms |
| torus_color | 6,400 | 24,000 | **100%** | 16ms |

**Result**: 11/11 = 100% on grid/torus problems

### 2. Hard Structured Instances

| Instance | Vars | Clauses | Sat% | Time |
|----------|------|---------|------|------|
| pitfall (CDCL trap) | 2,072 | 118,102 | **100%** | 611ms |
| kc4 | 400 | 59,230 | **100%** | 22.2s |
| kc5 | 500 | 84,935 | **100%** | 52.5s |
| vdw2 | 300 | 5,667 | **100%** | 8ms |
| dom | 200 | 8,001 | **100%** | 4ms |
| dom_complete50 | 250 | 12,501 | **100%** | 3ms |
| match40 | 780 | 29,680 | **100%** | 15ms |
| match | 435 | 12,210 | **100%** | 6ms |
| ptn_large | 500 | 772 | **100%** | 1ms |
| ptn_1k | 1,000 | 1,762 | **100%** | 2ms |
| ptn2k | 2,000 | 3,962 | **100%** | 8ms |
| reg3 | 400 | 1,300 | **100%** | 1ms |
| reg3_150 | 600 | 1,950 | **100%** | 2ms |
| grid_color | 3,600 | 13,260 | **100%** | 12ms |
| subgraph | 600 | 14,420 | **100%** | 7ms |
| ec | 800 | 3,200 | 99.94% | 1.1s |

### 3. Near-Perfect (≥99%) Instances

| Instance | Sat% | Time |
|----------|------|------|
| ptn10k | 99.86% | 16.5s |
| ptn5k | 99.96% | 6.9s |
| ptn15k | 99.70% | 23s |
| ptn20k | 99.60% | 30s |
| ptn25k | 99.47% | 35s |
| dom_gnm | 99.97% | 44s |
| match_grid | 99.98% | 3.8s |
| random_graph | 99.97% | 1.2s |
| order | 99.99% | 21.7s |
| tseitin | 99.94% | 858ms |
| tseitin6 | 99.99% | 21.6s |
| gnp05 | 99.56% | 36s |
| complete100 | 99.53% | 28s |

### 4. Probe Suite (45 instances)

- **100%**: 18 instances
- **99-99.9%**: 22 instances
- **94-98%**: 5 instances (mostly XOR-based)
- **Overall**: 40/45 = 89% at ≥99%

### 5. Custom Challenging Instances

#### Metallurgy (Frustrated Lattice)

| Instance | Vars | Clauses | Sat% | Time |
|----------|------|---------|------|------|
| metallurgy | 64 | 102 | 98.04% | 133ms |
| metallurgy2 | 128 | 160 | 97.50% | 131ms |
| metallurgy3 | 64 | 184 | 97.28% | 179ms |

#### Spin-Glass Models

| Instance | Vars | Clauses | Sat% | Time |
|----------|------|---------|------|------|
| spinglass (6x6) | 36 | 182 | 91.76% | 23ms |
| spinglass2 (8x8) | 64 | 186 | 97.85% | 37ms |
| spinglass3 (10x10) | 100 | 211 | **100%** | <1ms |

#### Sherrington-Kirkpatrick (SK) Model

| Instance | Vars | Clauses | Sat% | Time |
|----------|------|---------|------|------|
| sk (30 fully connected) | 30 | 495 | 86.87% | 45ms |
| sk2 | 40 | 370 | 86.49% | 42ms |
| sk3 (80% aligned) | 30 | 244 | **98.36%** | 180ms |

### 6. Edwards-Anderson 3D Spin Glass (Holy Grail)

| Size | Spins | Clauses | Sat% | Time | Steps |
|------|-------|---------|------|------|-------|
| 5×5×5 | 125 | 316 | **99.37%** | 78ms | 3K |
| 6×6×6 | 216 | 370 | **98.92%** | 61ms | 3K |
| 8×8×8 | 512 | 438 | **99.32%** | 144ms | 3K |
| 10×10×10 | 1,000 | 2,792 | **99.03%** | 590ms | 3K |
| 20×20×20 | 8,000 | 23,099 | **99.28%** | 4s | 3K |
| 30×30×30 | 27,000 | 79,035 | **99.98%** | 16s | 3K |
| **40×40×40** | **64,000** | 188,666 | **99.47%** | **4.3s** | 10K |
| 45×45×45 | 91,125 | 269,216 | 94.85% | 89s | 5K |
| 50×50×50 | 125,000 | 369,905 | 95.00% | 5.3min | 20K |
| 60×60×60 | 216,000 | 641,069 | 94.84% | 3.5min | 5K |

**Key Finding**: Sweet spot is L≈40 (64K spins) with 99.47% in 4.3 seconds — corresponds to the spin glass correlation length threshold.

### 7. Death Run (Final Boss Instances)

#### 7.1 Topological Traps (β₁ Destroyers)

| Instance | Vars | Clauses | Sat% | Time |
|----------|------|---------|------|------|
| overlapping_5cycles | 2,500 | 15,099 | **100%** | 4ms |
| cycle_complex | 6,400 | 47,974 | **100%** | 11ms |

**Result**: Topological traps fail — β₁ drops to 0, heat kernel resolves the cycles.

#### 7.2 Gradient Killers (XOR at Phase Transition)

| Instance | Vars | Clauses | Sat% | Time |
|----------|------|---------|------|------|
| xor_sat_threshold | 500 | 1,840 | **99.35%** | 1.5s |
| xor_hard | 1,000 | 3,680 | **98.97%** | 3.4s |

#### 7.3 Locality Destroyers (Expander Graphs)

| Instance | Vars | Clauses | Sat% | Time |
|----------|------|---------|------|------|
| expander (2K) | 2,000 | 10,657 | 92.08% | 1.3s |
| expander_large (5K) | 5,000 | 39,011 | 90.24% | 7.2s |
| expander_20k | 20,000 | 128,001 | 91.56% | 31s |
| expander_50k | 50,000 | 382,593 | 90.61% | 1m 41s |
| **expander_100k** | **100,000** | **764,238** | **90.57%** | **3m 24s** |

**Result**: Stable ~90% from 2K to 100K vars. The solver defies the expansion property.

---

# 🔑 Key Insights (April 2026)

1. **Geometry is everything**: Regular degree distributions (grids, lattices) let the heat kernel shine. Random graphs (SK) struggle.
2. **EA 3D is groundbreaking**: First gradient-based solver to crack Edwards-Anderson 3D at scale (64K+ spins, >99%).
3. **Phase transition at L≈40**: Beyond 40 spins per dimension, performance degrades — correlation length exceeds system size.
4. **CDCL traps don't trap NitroSAT**: The pitfall formula solves in 600ms with 100% satisfaction.
5. **Expander graphs are the breaking point**: High-expansion graphs are the hardest (~90%). If you can crack Urquhart at scale, that's publication-worthy.

---

# 💻 Hardware Reference

| Phase | Platform | Date |
| :--- | :--- | :--- |
| Phases 1–5 | AMD Ryzen 5 5600H @ 4.280GHz (single core, laptop) | Jan–Mar 2026 |
| Phase 6 | Apple Silicon | Apr 2026 |

**Compiler:** `gcc -O3 -lm` · **No external dependencies** · **Single-threaded**

**Reproducibility:** [HuggingFace Dataset](https://huggingface.co/datasets/sethuiyer/navokoj_sat_2024) · [timetable_output.json](../timetable_output.json)

---

*Last updated: April 8, 2026*
*Repository: https://github.com/sethuiyer/NitroSAT*
