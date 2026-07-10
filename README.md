<div align="center">

# NitroSAT

![NitroSAT Logo](img/logo.png)

**A High-Performance O(M) MaxSAT approximator with unusually good scaling behavior and high satisfaction rates on massive instances. Using Physics-Informed Continuous Relaxation**. Made in India by [Shunyabar Labs](https://shunyabar.foo/).

[![License: Apache 2.0](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.18753235-blue)](https://doi.org/10.5281/zenodo.18753235)
[![Codeberg](https://img.shields.io/badge/Codeberg-Lua%20Suite-2185d0.svg?logo=codeberg)](https://codeberg.org/sethuiyer/NitroSAT)
[![Verified Math](https://img.shields.io/badge/Math-Empirically_Verified-success.svg)](MATH.md#empirical-verification-2026-independent-audit)

</div>

---

## Repository Status

This repository is the **open-source NitroSAT research release**. It includes
the public solver lineage currently present in this tree:

| Solver path | Source | Status |
|---|---|---|
| V1 research solver | [`src/c/v1/nitrosat.c`](src/c/v1/nitrosat.c) | Open source / historical reference |
| V2 in-memory heuristic solver | [`src/c/v2/nitrosatv2.c`](src/c/v2/nitrosatv2.c) | Open source / feature-rich local solver |
| V3 disk-streaming solver | [`src/c/v3/nitrosatv3.c`](src/c/v3/nitrosatv3.c) | Open source / bounded-memory research release |

Further production development of NitroSAT is **not open source**. The
commercial version is available as the `nitro` engine through the
[Navokoj API](https://navokoj.shunyabar.foo/), where ShunyaBar Labs continues
work on hosted execution, production hardening, telemetry, billing, enterprise
deployment, and commercial support.

The goal of this split is simple: keep the research release inspectable and
reproducible, while keeping the production infrastructure and future commercial
engine development inside Navokoj.

---

## Overview

NitroSAT is a high-performance MaxSAT approximator with unusually good scaling
behavior and high satisfaction rates on massive instances. Unlike traditional
CDCL-based solvers, NitroSAT treats Boolean satisfiability as a
physics-informed dynamical system on a Riemannian manifold, using continuous
relaxation, spectral methods, and topological analysis.

We threw thousands of problems at NitroSAT—everything from domset_4 to 7-million-clause planted coloring monsters to random 3SAT, 5,000+ CNF instances and solver came out with a 77% perfect SATISFIED rate and a 99.7% median satisfaction. Assignments: https://huggingface.co/datasets/sethuiyer/navokoj_sat_2024 and benchmarks/README.md.

> **Note from GPT-5.5:** NitroSAT v2 looks very strong as a fast approximate/heuristic satisfiability engine, especially for high-satisfaction or MaxSAT-style use, but the formal proof and benchmark hygiene need tightening before making stronger solver-theory claims.

Some notable results:

| Test Name / Problem Type | Variables | Clauses | Ratio ($\alpha$) | Latency | Satisfaction |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Planted Coloring (Lua)** | 105,000 | 232,043 | 2.21 | 13.78s | 100% |
| **Grid Coloring (1000x1000)** | 4,000,000 | 14,992,000 | 3.75 | 475.0s | 100% |
| **Massive Hardware Check** | 788,480 | 2,617,349 | 3.32 | 5.92s | 100% |
| **Enterprise Timetabling** | 147,600 | 80,278,884 | 543.89 |73sec | 100% |
| **Titan Ramsey Density** | 780 | 1,316,016 | **1,687.20** | 3,403.0s | 99.995% |
| **Adversarial Pitfall Trap** | 2,950 | 1,047,620 | 355.12 | ~400.0s | 100% |
| **Large Phase Transition** | 200,000 | 852,000 | 4.26 | 11.2 min | 99.20% |
| **Ramsey R(5,5,5)** | ~4,760 | 354,890 | 74.56 | 3.5s | 100% |
| **Small-World Lattice** | 360,000 | 1,354,800 | 3.76 | 62.75s | 100% |
| **Planted Random 3-SAT (V3 indexed)** | 100 | 923 | 9.23 | 24ms | 100% |

---

### V3 indexed-finisher sanity check

A reproducible planted random 3-SAT instance (`seed=20260706`) with 100
variables and 923 clauses was solved completely by V3 after 209 indexed flips:

| Metric | Result |
|---|---:|
| Clauses satisfied | 923 / 923 |
| Indexed flips | 209 |
| Solve time | 24ms |
| Peak RSS | 2.4 MiB |
| Independent assignment verification | Passed |
| Exact fallback required | No |

The instance was generated around a planted assignment, so satisfiability was
guaranteed by construction. The returned assignment was then checked against
all 923 clauses independently. This is a focused correctness and integration
test, not evidence that every instance at the same clause-to-variable ratio is
equally easy.

```bash
src/c/v3/nitrosatv3 planted-100v-923c.cnf \
  --epochs 20 \
  --indexed-finisher \
  --indexed-flips 50000 \
  --indexed-checkpoint 500 \
  --solution assignment.sol
```

### Planted-SAT matrix

V3 was also tested on twelve reproducible planted 3-SAT instances across
multiple sizes and clause densities. The standard run used 10 epochs, a
100,000-flip indexed budget, and checkpoints every 500 flips.

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

Eleven of the twelve instances were solved completely, and their returned
assignments were checked independently against every clause. The difficult
cases cluster near the random 3-SAT phase-transition ratio around 4.26;
higher-density planted instances expose a stronger planted signal and were
easier for this search configuration. These are planted satisfiable instances,
not a substitute for an independently sourced industrial benchmark suite.

---

## Key Features

- **Linear Scaling** — O(M) time complexity relative to clause count
- **Structural Awareness** — Detects UNSAT via thermodynamic phase transitions
- **Zero Configuration** — Works out-of-the-box on Scheduling, Ramsey, Coloring, N-Queens
- **Multi-Phase Architecture** — Langevin flow, topological repair, and adelic saturation phases
- **UNSAT Detection** — GF(2) parity checking and topological probe ordering
- **No External Dependencies** — Standalone C99 implementation
- **Lua Implementation** — LuaJIT version available via [Codeberg](https://codeberg.org/sethuiyer/NitroSAT)

---


## Algorithm Overview

NitroSAT maps the Boolean satisfiability problem to an energy landscape, using spectral initialization and Branch-Aware Holonomy Annealing (BAHA) to navigate complex basins.

1. **Continuous Relaxation** — Maps Boolean variables to continuous values in [0,1]
2. **Prime Weighting** — Number-theoretic clause weights based on the Prime Number Theorem
3. **Heat Kernel Diffusion** — Spectral smoothing via the heat kernel `exp(t * Δ)`
4. **BAHA (Branch-Aware Holonomy Annealing)** — Uses Lambert-W function for phase transitions
5. **Persistent Homology** — Tracks topological holes (Betti numbers) to guide repair
6. **NADAM Optimization** — Nesterov-accelerated adaptive moment estimation

See [MATH.md](MATH.md) for the complete mathematical theory.

---

![Satisfaction Distribution](img/satisfaction_hist.png)
![Scaling Analysis](img/scaling_plot.png)

---

## Quick Start

For inputs too large to hold clause state in RAM, use the V3 disk-streaming
solver. V2 remains the feature-complete in-memory solver.

```bash
make -C src/c v3
src/c/v3/nitrosatv3 base.cnf --add incremental.cnf \
  --store formula.nsv3 --solution assignment.sol
```

V3 accepts raw CNF and WCNF, uses 64-bit clause/literal totals, supports
transactional incremental constraint files, and bounds memory by variable
count plus configurable batch and active-clause limits. See
[the V3 documentation](src/c/v3/README.md).

### Executable tutorial library

The [NitroSAT tutorial library](docs/tutorials/README.md) contains tested,
runnable examples for Sudoku, workforce scheduling, university timetabling,
vehicle assignment, meeting rooms, graph coloring, Kubernetes pod placement,
product configuration, exam scheduling, and nurse rostering. Every tutorial
follows the same model-to-result workflow and is executed in CI to prevent
documentation drift.

### 1. Compile

```bash
gcc -O3 -march=native -std=c99 nitrosat.c -o nitrosat -lm
```

**Requirements:** GCC or any C99-compatible compiler + standard math library (`-lm`)

### 2. Run

```bash
./nitrosat problem.cnf
```

JSON is the **default output format**. Progress milestones are printed to stderr so stdout stays clean for piping:

```
NitroSAT  problem.cnf  | 625 vars  24825 clauses
[pass 1/5]   [stagnation] 99.99% plateau → jumping to finisher
  [phase-2] topological repair... done → 24823/24825 (99.99%)
  [phase-3] adelic saturation... done → 24820/24825 (99.99%)
```

### 3. Command-Line Options

| Option | Description |
|--------|-------------|
| `<cnf-file>` | Input file in DIMACS CNF format |
| `[max-steps]` | Maximum optimization steps (default: adaptive) |
| `--no-dcw` | Disable prime-weighted clause learning |
| `--no-topo` | Disable topological repair phase |
| `--json` | Output JSON format (default) |
| `--cinematic` | Verbose per-step logging |
| `--proof <path>` | Generate DRAT proof for UNSAT instances |
| `--proof-format drat` | Proof format (default: drat) |

### 4. Piping & Integration

```bash
# Extract the assignment
./nitrosat problem.cnf | jq '.assignment'

# Check if solved
./nitrosat problem.cnf | jq '.status'

# Save result to file (progress still prints to terminal)
./nitrosat problem.cnf > result.json

# Parse in Python
./nitrosat problem.cnf | python3 -c "import json,sys; d=json.load(sys.stdin); print(d['satisfaction_rate'])"
```

### 5. Cinematic Mode (Verbose)

```bash
./nitrosat problem.cnf --cinematic
```

### 6. UNSAT Certificate Mode (DRAT)

```bash
./nitrosat problem.cnf --proof proof.drat --proof-format drat --json
```

The proof backend is **solver-aware** — it uses NitroSAT's own UNSAT detection signals to guide proof generation:

1. **UNSAT core extraction** — extracts the irreducible core after solving
2. **GF(2) parity check** — detects odd-cycle contradictions in the binary clause implication graph
3. **Topological probe ordering** — failed-literal probing ordered by prime-weighted irreducibility pressure
4. **Solver diagnostics** — reports UNSAT awareness signals: unsatisfied count, fractures, beta1, convexity

Proof status values in JSON output:
- `generated_top_level_up_conflict` — contradiction found by unit propagation
- `generated_gf2_parity_refutation` — contradiction from odd parity cycle
- `generated_topo_failed_literal_refutation` — contradiction via math-guided probing
- `inconclusive` — no proof derived; signals explain why

> **Note:** The backend is exact (sound) but incomplete. Pigeonhole-style instances require extended resolution for short proofs.

---

## JSON Output Format

```json
{
  "status": "SATISFIED",
  "satisfaction_rate": 1.0,
  "satisfied": 425,
  "unsatisfied": 0,
  "variables": 150,
  "clauses": 425,
  "assignment": [-1, -2, 3, -4, -5, 6],
  "latency": {
    "total_ms": 54.12,
    "throughput": "7.9K clauses/sec"
  }
}
```

**Key fields:**

| Field | Description |
|-------|-------------|
| `confidence` | Raw continuous values (0.0–1.0) before rounding to Boolean |
| `critical_beta` | Predicted phase transition point from the Prime Number Theorem |
| `convexity_status` | `STABLE` (convex regime) or `NON_CONVEX` (glassy phase) |
| `betti_0` / `betti_1` | Topological invariants of the unsatisfied clause complex |
| `blame` | Top-20 unsatisfied clauses ranked by prime weight (UNSAT debugging) |
| `proof` | UNSAT certificate request/result metadata |

---

## Citation

If you use NitroSAT in your research, please cite:

```bibtex
@software{sethurathienam_iyer_2026_18753235,
  author       = {Sethurathienam Iyer},
  title        = {NitroSAT: A Physics-Informed MaxSAT Solver},
  year         = 2026,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.18753235},
  url          = {https://doi.org/10.5281/zenodo.18753235}
}
```

---

## Reproducibility

The CNF files used to test NitroSAT are available at:
- [HuggingFace Dataset](https://huggingface.co/datasets/sethuiyer/navokoj_sat_2024)

---

### Implementations

| Implementation | Location |
|----------------|----------|
| NitroSAT (C) | [GitHub](https://github.com/sethuiyer/NitroSAT) · [Zenodo](https://zenodo.org/records/18753235) |
| NitroSAT (Lua) | [Codeberg](https://codeberg.org/sethuiyer/NitroSAT) |
| Navokoj API | [navokoj.shunyabar.foo](https://navokoj.shunyabar.foo/) |

---

## License

NitroSAT is licensed under the Apache License 2.0. See the [LICENSE](LICENSE) file for details.

[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2Fsethuiyer%2FNitroSAT.svg?type=large\&issueType=license)](https://app.fossa.com/projects/git%2Bgithub.com%2Fsethuiyer%2FNitroSAT?ref=badge_large&issueType=license)

---

## Support & Contributing

- **Issues**: Report bugs and feature requests via GitHub Issues
- **Email**: Contact the maintainer at shunyabarlabs@zohomail.com
- **ORCID**: [0009-0008-5446-2856](https://orcid.org/0009-0008-5446-2856)
- **Sponsor**: Support development at [github.com/sponsors/sethuiyer](https://github.com/sponsors/sethuiyer/)

Bug reports, reproducibility notes, documentation fixes, and benchmark
artifacts for this open-source research release are welcome. Production
`nitro` engine development happens inside Navokoj and is not developed through
this public repository.

Please follow the [Code of Conduct](CODE_OF_CONDUCT.md).
