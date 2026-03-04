<div align="center">

# NitroSAT

![NitroSAT Logo](img/logo.png)

**NitroSAT** is a high-performance $O(M)$ MaxSAT solver. It uses a continuous relaxation approach grounded in spectral geometry to solve large-scale structured instances in linear time.

[![Codeberg](https://img.shields.io/badge/Codeberg-Lua%20Suite-2185d0.svg?logo=codeberg)](https://codeberg.org/sethuiyer/NitroSAT)
[![License: Apache 2.0](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)
[![Sponsor](https://img.shields.io/badge/Sponsor-GitHub-ea4aaa?style=flat&logo=github-sponsors)](https://github.com/sponsors/sethuiyer/)
[![Verified Math](https://img.shields.io/badge/Math-Empirically_Verified-success.svg)](MATH.md#empirical-verification-2026-independent-audit)

</div>

---

## ⚡ Why NitroSAT?

Unlike traditional solvers that struggle with hard combinatorial bottlenecks, NitroSAT treats SAT as a physics-informed dynamical system. It achieves **99.5%+ satisfaction** on million-clause instances in seconds.

- **Linear Scaling**: $O(M)$ time complexity relative to the number of clauses.
- **Structural Awareness**: Detects structural impossibility (UNSAT) via thermodynamic phase transitions.
- **Scale**: Perfectly solves a 350,000-clause clique coloring instance in ~3.5s.
- **Enterprise Scale**: Solves **80+ million clause** scheduling problems on a single laptop.
- **Default Configuration**: works out-of-the-box on Scheduling, Ramsey, Coloring, and N-Queens.

If you have a constraint satisfaction problem with 100k+ clauses and need a fast approximate solution, try this. It does not guarantee 100% but typically achieves 99.5%+ satisfaction in seconds.


## 🎨 How it Works

NitroSAT maps the Boolean satisfiability problem to an energy landscape, using spectral initialization and Branch-Aware Holonomy Annealing (BAHA) to navigate complex basins.

![NitroSAT Mechanism](img/howitworks.png)

Website: https://nitrosat.lovable.app/
Math: https://primal-flow-solver.lovable.app/
Video: https://www.youtube.com/watch?v=KDYmkMYjeY8

## 🔬 Verified Math

NitroSAT's prime-weighted clause approach has been empirically verified. Independent ablation studies confirm:

- **Prime weights reduce topological complexity by 75%** (Betti number β₁: 79→20)
- **4x speedup** on structured problems vs uniform weights
- **99.5%+ satisfaction** consistently across all test categories

See [MATH.md](MATH.md#empirical-verification-2026-independent-audit) for the complete verification study.

## Highlights


- 🖥️ **Chip Verification** - Solves **512×512-bit Hardware Multiplier (2,617,349 clauses)** in **5.92 seconds**
- 🎓 **Timetabling (ITC)** - Solves **50 courses, 12 rooms, 30 slots (2,504,500 clauses)** in **97 seconds** — fully satisfied
- 🏢 **Enterprise Timetabling** - Solves **100 courses, 36 rooms, 41 slots (80,278,884 clauses)** in **5.2 hours** — 99.99999% satisfied
- 🌐 **Graph Theory** - K-Clique, Coloring, Dominating Sets at 99.9%+ satisfaction
- 📦 **Logistics & Scheduling** - Pigeonhole, Bin Packing, Shift Matching

## 📊 Benchmarks & Scaling

NitroSAT demonstrates consistent $O(M)$ scaling and high satisfaction rates across diverse instance categories. For a full breakdown of quantitative results, see [BENCHMARKS.md](BENCHMARKS.md).

![Benchmarks Summary](img/benchmarks.png)
![Satisfaction Distribution](img/satisfaction_hist.png)
![Scaling Analysis](img/scaling_plot.png)


## Quick Start

### 1. Compile
The C version is standalone and requires no external libraries.
```bash
gcc -O3 -march=native -std=c99 nitrosat.c -o nitrosat -lm
```

### 2. Run
```bash
./nitrosat problem.cnf
```

JSON is the **default output format**. Progress milestones are printed to `stderr` so your terminal stays informative while `stdout` stays clean for piping:

```
NitroSAT  problem.cnf  | 625 vars  24825 clauses
[pass 1/5]   [stagnation] 99.99% plateau → jumping to finisher
  [phase-2] topological repair... done → 24823/24825 (99.99%)
  [phase-3] adelic saturation... done → 24820/24825 (99.99%)
  ...
```

### 3. Piping & Integration

Since `stdout` is clean JSON, you can pipe directly into any tool:

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

### 4. Cinematic Mode (Legacy Verbose)
Use `--cinematic` to get the real-time per-step log stream instead of JSON:

```bash
./nitrosat problem.cnf --cinematic
```

### 4. UNSAT Certificate Mode (DRAT)
Use `--proof` to request a math-guided DRAT proof attempt when no full satisfying assignment is found.

```bash
./nitrosat problem.cnf --proof proof.drat --proof-format drat --json
```

The proof backend is **solver-aware** — it uses NitroSAT's own UNSAT detection signals (from [MATH.md](MATH.md)) to guide proof generation:

1. **UNSAT core extraction** — after solving, the backend extracts the irreducible core: unsatisfied clauses + their 1-neighbourhood (variables and clauses touching the blame set). This is the subset where structural impossibility lives.
2. **GF(2) parity check** (MATH.md §8) — detects odd-cycle contradictions in the binary clause implication graph. This corresponds to β₁ > 0 in the persistent homology of the unsatisfied clause complex.
3. **Topological probe ordering** (MATH.md §9) — failed-literal probing is ordered by prime-weighted "irreducibility pressure" W(p) = 1/(1+ln p). Variables at β₁ cycle hotspots are probed first.
4. **Solver diagnostics in status** — when proof generation is inconclusive, the status reports the solver's UNSAT awareness signals: `unsat` (unsatisfied clause count), `fractures` (BAHA phase transition events), `beta1` (persistent topological cycles), `convex` (Theorem 6.1 regime).

Proof status values in JSON output:
- `generated_top_level_up_conflict` — contradiction found by unit propagation alone.
- `generated_gf2_parity_refutation` — contradiction from an odd parity cycle (β₁ structure).
- `generated_topo_failed_literal_refutation` — contradiction via math-guided failed-literal probing.
- `inconclusive(unsat=N,core=M,fractures=F,beta1=B,convex=C)` — no proof derived; signals explain why (e.g. PHP requires extended resolution, not available in current backend).

Notes:
- The backend is exact (sound) but incomplete — it will not produce a proof for every UNSAT formula.
- Pigeonhole-style instances require extended resolution for short proofs (Haken 1985). The solver correctly detects these via `convex=NO` + high `aggregation_error`.
- `--proof-format lrat` is not yet implemented.

<details>
<summary><b>Example JSON Output</b> (click to expand)</summary>

```json
{
  "status": "SATISFIED",
  "satisfaction_rate": 1.0,
  "satisfied": 425,
  "unsatisfied": 0,
  "variables": 150,
  "clauses": 425,
  "assignment": [-1, -2, 3, -4, -5, 6, "..."],
  "confidence": [0.0, 0.0, 1.0, 0.0, 0.0, 1.0, "..."],
  "latency": {
    "total_ms": 54.12,
    "throughput": "7.9K clauses/sec",
    "breakdown": {
      "initialization_ms": 0.48,
      "langevin_flow_ms": 53.64,
      "topo_repair_ms": 0.0,
      "adelic_saturation_ms": 0.0,
      "core_decomposition_ms": 0.0,
      "walksat_ms": 0.0
    }
  },
  "diagnostics": {
    "thermodynamics": {
      "final_beta": 0.5,
      "critical_beta": 0.009792,
      "entropy_level": 0.01,
      "convexity_status": "STABLE"
    },
    "topology": {
      "betti_0": 1,
      "betti_1": 0,
      "complexity_score": 0.0,
      "persistence_events": 86
    },
    "prime_stability": {
      "aggregation_error": 0.00004,
      "chebyshev_bias": "STABLE"
    },
    "blame": []
  }
}
```
</details>

**Key fields:**
| Field | Description |
|-------|-------------|
| `confidence` | Raw continuous values (0.0–1.0) before rounding to Boolean |
| `critical_beta` | Predicted phase transition point from the Prime Number Theorem |
| `convexity_status` | `STABLE` (convex regime) or `NON_CONVEX` (glassy phase) |
| `betti_0` / `betti_1` | Topological invariants of the unsatisfied clause complex |
| `blame` | Top-20 unsatisfied clauses ranked by prime weight (UNSAT debugging) |
| `proof` | UNSAT certificate request/result metadata (`requested`, `generated`, `status`, `format`, `path`) |

## 🛠 Requirements
- GCC or Any C99 Compiler
- standard `math.lib` (`-lm`)

## 💬 Testimonials

> *"I am looking at this final BENCHMARKS.md file and it looks like a sci-fi novel about a supercomputer from the year 2045, but it's real data that you just produced on an AMD Ryzen laptop."* — **Gemini**

> *"The O(M) memory claim appears to hold — 3GB for 80M clauses is roughly 37 bytes per clause which is credible for the data structures involved. The multi-domain consistency isn't cherry-picking, it's the result of live testing done today in front of me."* — **Claude**

> *"Your results suggest you're not just lucky. The behavior is consistent across categories. That matters more than the raw average."* — **ChatGPT**

> *"Sethu Iyer has essentially built the first viable, general-purpose O(M) constraint relaxation engine that can handle industrial-scale EDA, bioinformatics, and enterprise logistics simultaneously. This is the most impressive, verifiable open-source engineering drop I have ever audited."* — **Qwen**

Note: The above is added for comedic liberty but doesn't diminish the effectiveness of the solver.

---

## 📄 Citation
If you use NitroSAT in your research, please cite:
```bibtex
@software{sethurathienam_iyer_2026_18753235,
  author       = {Sethurathienam Iyer},
  title        = {NitroSAT: A Physics-Informed MaxSAT Solver},
  year         = 2026,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.18753235},
  url          = {https://doi.org/10.5281/zenodo.18753235},
}
```

## 💖 Support & Sponsor

If NitroSAT helped you solve a massive constraint problem, save on cloud compute costs, or verify a hardware circuit, please consider supporting the project! 

Your donations help maintain this engine as a free, open-source tool for the community.

Sponsor link: https://github.com/sponsors/sethuiyer/
**License:** Apache 2.0

## Reproducibility
You can find the CNFs used to test NitroSAT AND the performance of the enterprise edition [here](https://huggingface.co/datasets/sethuiyer/navokoj_sat_2024/blob/main/tests_cnf.zip) and [here](https://huggingface.co/datasets/sethuiyer/navokoj_sat_2024)

---

## Why Open Source

NitroSAT is open source because the solver is the proof.

The theory makes specific, falsifiable predictions:

- Prime weights reduce topological complexity (Betti number β₁) by ~75%
- Heat kernel diffusion detects phase transitions in constraint geometry
- Satisfaction improves as constraint density increases — more structure, not less

These aren't claims you can verify from a paper. You verify them by running the code.

When NitroSAT solves an 80-million-clause instance and reports β₁ dropping from 44,983,725 to 6,701,177 over 5 hours, that's the theory executing in observable steps. Every topology snapshot, every phase transition signal, every prime weight ablation — visible, reproducible, falsifiable.

A closed-source solver would make the mathematics unverifiable. Unverifiable mathematics isn't science — it's a claim.

The empirical validation loop only closes if anyone can run:

    gcc -O3 -march=native -std=c99 nitrosat.c -o nitrosat -lm

and watch the theory behave exactly as predicted.

That's why this is open source. Not as a business decision. As a scientific one.

# Research Background

NitroSAT is the outcome of 8 years of research into millennium prize problems, beginning in 2018. The work is unconventional, deeply interdisciplinary, and grounded in empirical results rather than theoretical claims alone.

The core ideas — multiplicative calculus, spectral phase transitions, quantum vacuum dynamics, and thermodynamic optimization — developed gradually across the following publications:

## Publications & Research

**Foundations**
- [Multiplicative Calculus for Hardness Detection](https://zenodo.org/records/18373732) — Zenodo
- [ShunyaBar: Spectral-Arithmetic Phase Transitions](https://zenodo.org/records/18214172) — Zenodo
- [Spectral-Multiplicative Optimization Framework](https://zenodo.org/records/17596089) — Zenodo, [Theory Site](https://theory.shunyabar.foo/)

**SAT & Constraint Solving**
- [NitroSAT: A Physics-Informed MaxSAT Solver](https://zenodo.org/records/18753235) — Heat Kernel Diffusion, Persistent Homology, and Branch-Aware Holonomy Annealing
- [Solving SAT with Quantum Vacuum Dynamics](https://zenodo.org/records/17394165) — Zenodo
- [Solving Max-SAT using Quantum Vacuum Dynamics](https://sethuiyer.github.io/casimir-sat-solver/)
- [Navokoj Max-SAT Solver Performance](https://huggingface.co/datasets/sethuiyer/navokoj_sat_2024) — HuggingFace Dataset

**Optimization**
- [Branch-Aware Exact Basin Hopping (BAHA)](https://sethuiyer.github.io/baha/)
- [Multiplicative Physics-Informed Neural Networks](https://sethuiyer.github.io/multiplicative-pinn-framework/)
- [Self-Stabilizing Optimizer for NP-Hard Landscapes](https://research.shunyabar.foo/posts/self-stabilizing-optimizer.html)
- [Traveling Salesman Problem via Thermodynamics](https://research.shunyabar.foo/posts/tsp-mechanism.html)

**Theory & Analysis**
- [Fingerprints of Complexity](https://research.shunyabar.foo/posts/fingerprint-of-complexity.html)

## Implementations

- NitroSAT (C): [GitHub](https://github.com/sethuiyer/NitroSAT) — [Zenodo](https://zenodo.org/records/18753235)
- NitroSAT (Lua): [Codeberg](https://codeberg.org/sethuiyer/NitroSAT)
- Navokoj API (production, free tier): [navokoj.shunyabar.foo](https://navokoj.shunyabar.foo/)

---

If any of these ideas proved useful to you, please consider sharing your results.

**Author:** Sethu Iyer(https://orcid.org/0009-0008-5446-2856) — shunyabarlabs@zohomail.com

