<div align="center">

# NitroSAT

![NitroSAT Logo](img/logo.png)

**A High-Performance O(M) MaxSAT Solver Using Physics-Informed Continuous Relaxation**. Made in India by [Shunyabar Labs](https://shunyabar.foo/).

[![License: Apache 2.0](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.18753235-blue)](https://doi.org/10.5281/zenodo.18753235)
[![Codeberg](https://img.shields.io/badge/Codeberg-Lua%20Suite-2185d0.svg?logo=codeberg)](https://codeberg.org/sethuiyer/NitroSAT)
[![Verified Math](https://img.shields.io/badge/Math-Empirically_Verified-success.svg)](MATH.md#empirical-verification-2026-independent-audit)

</div>

---

## Overview

NitroSAT is a high-performance MaxSAT solver that achieves **O(M) linear time complexity** relative to the number of clauses. Unlike traditional CDCL-based solvers, NitroSAT treats Boolean satisfiability as a physics-informed dynamical system on a Riemannian manifold, using continuous relaxation, spectral methods, and topological analysis.


| Test Name / Problem Type | Variables | Clauses | Ratio ($\alpha$) | Latency | Satisfaction |
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

The solver consistently achieves **99.5%+ satisfaction** on million-clause instances and has been verified even on problems with over **80 million clauses** on a single laptop core. The version 2 of the solver
features improved convergence rate to perfect satisfaction on the same problem.

| Instance | Variables | Clauses | Satisfaction | Time | Throughput | Hardware |
|----------|-----------|---------|--------------|------|------------|----------|
| timetable.cnf | 147,600 | 80,278,884 | **100%** | **73s** | **1.1M clauses/sec** | AMD Ryzen 5 5600H (laptop) |

**Latency Breakdown:**
- Initialization: 41s
- Langevin Flow: 33s
- Total: 73s

*Run command:* `./nitrosatv2 timetable.cnf > output.json`

**Output:** [timetable_output.json](benchmarks/timetable_output.json)

CNF instance: [timetable.cnf](https://huggingface.co/datasets/sethuiyer/navokoj_sat_2024/blob/main/tests_cnf.zip) from HuggingFace

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

![Benchmarks Summary](img/benchmarks.png)


## Performance Highlights
For complete benchmark results, see [benchmarks/BENCHMARKS.md](benchmarks/BENCHMARKS.md).

![Satisfaction Distribution](img/satisfaction_hist.png)
![Scaling Analysis](img/scaling_plot.png)

---

## Quick Start

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

---

## Support & Contributing

- **Issues**: Report bugs and feature requests via GitHub Issues
- **Email**: Contact the maintainer at shunyabarlabs@zohomail.com
- **ORCID**: [0009-0008-5446-2856](https://orcid.org/0009-0008-5446-2856)
- **Sponsor**: Support development at [github.com/sponsors/sethuiyer](https://github.com/sponsors/sethuiyer/)

Contributions are welcome. Please follow the [Code of Conduct](CODE_OF_CONDUCT.md).

