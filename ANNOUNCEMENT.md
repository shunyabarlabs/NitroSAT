# Open Sourcing NitroSAT: The Thermodynamics of NP-Complete Problems


![Prime-Weighted Hyperbolic Manifold](img/announcement_1.png)
*Figure 1: Visualization of the Prime-Weighted Hyperbolic Manifold and Spectral Heat Kernels.*

Imagine submitting a constraint problem so large, so tangled with dependencies and resource conflicts, that traditional solvers simply give up.

You have likely been there. You run a shift scheduling optimization for 1,000 nurses across 15 departments. You configure a microservices dependency validation for a Kubernetes cluster. Or perhaps you are trying to find a Ramsey number for mathematical research. You hit "solve," and you watch the progress bar stall. The CPU spikes, memory usage climbs, and minutes turn into hours. Eventually, the process times out, crashes, or returns a solution that violates your most critical constraints.

Now, imagine a different reality.

Imagine submitting that same massive, messy problem and getting a near-perfect assignment back in milliseconds. Imagine the solver not only providing the result but also telling you exactly why it couldn't reach 100% ; whether it's a structural parity contradiction in your logic or a genuine resource conflict that requires human intervention.

That is the reality that **NitroSAT** delivers.

**Stop searching. Start evolving.**

Today, I am open-sourcing **NitroSAT**, the physics-informed MaxSAT engine that powers ShunyaBar Labs. NitroSAT represents a fundamental departure from traditional CSP/SAT theory. We have abandoned the search tree entirely in favor of a **continuous physics-based simulation**.


![Phase Transition](img/announcement_2.png)
*Figure 2: Real-time detection of thermodynamic phase transitions (The Mirage Trap).*


## Why am I open-sourcing this?

Because NitroSAT isn't just a solver anymore. It has become a **mathematical instrument**.

When I saw the benchmarks hitting **$O(N)$ linear scaling** on problems that mathematically "should" scale exponentially, I realized this needed to be audited by the community. When I saw the solver achieving **0.0000% permutation variance**; being completely blind to how you label your variables because it only sees the underlying structure of the problem, I knew we have something special here which demands scrutiny from the larger community.


## The Mind-Blowing Realities of NitroSAT

### 1. $O(N)$ Scaling is the New Reality
NitroSAT solves a **million-clause** grid coloring problem in ~12 seconds. It doesn't branch. It maintains a flat state array $x \in [0, 1]^V$ and updates it simultaneously across all variables using a highly optimized Sparse CSR structure. The physics scale linearly, even when the logic is NP-hard.

### 2. Zero Permutation Variance (Label Independence)
In a standard solver, shuffling your CNF file can change the runtime from 1ms to 1 hour. In NitroSAT, it changes **nothing**. 
Each clause is weighted by a unique prime number. This makes the solver **Label Independent**. It doesn't care about your variable names; it only cares about the core shape of the problem. It is the first solver that treats a SAT instance as a deterministic physical law.

### 3. Smart Escapes from Local Minima
Hidden inside the **BAHA (Branch-Aware Holonomy Annealing)** module is an advanced jump mechanism. 
Optimization often hits walls where progress stalls. Most solvers get stuck there. NitroSAT detects these boundaries and **teleports** the system into new stable energy basins. It doesn't just climb hills; it tunnels through them.

## From "Search" to "Physics"

For decades, the dominant approach to solving SAT (Boolean Satisfiability) and MaxSAT problems has been combinatorial search. Solvers like CDCL (Conflict-Driven Clause Learning) operate by navigating a massive binary search tree, encountering conflicts, and learning new clauses to avoid those conflicts in the future.

This approach works well for many problems, but it has a fatal flaw: **Structure**.

When problems are "structured" ; when they exhibit deep symmetries, tight variable couplings, or phase-transition boundaries ; search trees explode. NitroSAT evolves a solution by treating the constraint hypergraph as a physical system governed by energy, heat, and topology.

### Core Modules

1.  **Continuous Heat Diffusion**: Smoothing the energy landscape, preventing the solver from getting stuck in shallow local valleys.
2.  **Topology Tracking**: Detecting "holes" or loops in the problem structure (e.g., parity chains in XOR-SAT) and breaking them explicitly.
3.  **Prime Resonance**: Injecting prime-frequency harmonies to shake the system and help the solver tunnel through energy barriers.
4.  **Phase-Transition Jumps (BAHA)**: Detecting when the system hits a wall and calculating jumps to new energy basins.


![Heat Kernel Diffusion](img/announcement_3.png)
*Figure 3: Multi-scale relaxation via heat kernel diffusion operators.*


### The Math is the Code
Inside `nitrosat.c`, you won't find complex branching logic or massive heuristic tables. You will find:
- **Heat Diffusion**: Smoothing out the problem landscape.
- **Resonance Injections**: Shaking the system to prevent it from getting trapped.
- **Unique Prime Weights**: Ensuring every logical constraint has a unique, mathematically pure mass.

![Metric Geometry](img/math_disk.png)
*Figure 4: The Inverted Poincaré Disk;where constraints live on the boundary of the infinite vacuum.*

## The Zero-Tuning Approach

NitroSAT relies on fundamental mathematical constants (primes, the golden ratio, eigenvalues of the graph Laplacian), allowing the same parameters to work across wildly different domains without manual tweaking.



![Diffusion Flow](img/math_diffusion.png)
*Figure 6: Heat diffusion flow traversing the constraint hypergraph.*


### Performance Highlights

-   **N-Queens (25x25)**: 100% in 3.7 seconds.
-   **Exact Cover**: 100% in 4ms.
-   **Planted 3-SAT (α=4.26)**: 100% in 5ms.
-   **Hamiltonian Cycle**: 99.99% (missing 1 clause).

## UNSAT Awareness: The Mirage Trap

By monitoring "heat capacity" (energy function variance) and spectral properties, NitroSAT can detect when a solution is physically impossible (a "Mirage Trap"). This transforms the solver from a simple optimizer into a diagnostic engine that quantifies structural impossibility.

## Global Benchmarks

Across 358 benchmark instances:
-   **Average Clause Satisfaction**: 99.58%
-   **Instances Solved at 100%**: 115 (32.1%)
-   **Solved ≥ 99%**: 340 (95.0%)
-   **Speed**: Fastest 10K+ clause solve in 33ms.

## How to Get Started

NitroSAT is a bare-metal, dependency-free C implementation designed for extreme throughput.

### Compile (C)
```bash
gcc -O3 -march=native -std=c99 nitrosat.c -lm -o nitrosat
./nitrosat <cnf-file> [max-steps]
```

### Script (Lua)
```lua
local nitrosat = require("nitrosat")
local solver = nitrosat.NitroSat.new(instance, { seed = 42 })
local success, steps, sat_count = solver:solve()
```

## An Invitation to Audit

Optimization is filled with black-box solvers. I am releasing this under the **Apache 2.0 license** because the implications for complex systems and physics-based computation are too significant to keep behind closed doors. Whether you are an engineer looking for scalable optimization or someone interested in entirely new ways to solve logic problems, the code is now yours.

![Prime Resonance](img/math_zeta.png)
*Figure 5: Spectral perturbations derived from prime numbers.*

> *A note to the community: NitroSAT does not claim to magically bypass NP-hardness. Rather, it introduces a physical framework where prime numbers and thermodynamics guide the solver, offering a new empirical way to study problem difficulty.*

---
**Sethu Iyer**  
Founder, ShunyaBar Labs

---
**Codeberg (Full Suite):** [codeberg.org/sethuiyer/NitroSAT](https://codeberg.org/sethuiyer/NitroSAT)  

**GitHub (C Engine):** [github.com/sethuiyer/NitroSAT](https://github.com/sethuiyer/NitroSAT)

**Paper**: [codeberg.org/sethuiyer/NitroSAT/paper.pdf](https://codeberg.org/sethuiyer/NitroSAT/src/branch/master/nitrosat_paper.pdf)

**Navokoj**: [navokoj.shunyabar.foo](https://navokoj.shunyabar.foo/)

---

# NitroSAT vs Navokoj: What Separates Them

NitroSAT is the **open-source physics core**. Navokoj is the **production intelligence platform** built on top of it. They share the same underlying logic, but Navokoj adds everything needed to turn a research engine into industrial infrastructure.

---

## The One-Line Summary

> **NitroSAT** is the equation. **Navokoj** is the machine that runs it at scale, tells you *why* it failed, and charges you per solve.

---

## Architecture Comparison

| Dimension | NitroSAT (Open Source) | Navokoj (Production API) |
|---|---|---|
| **Language** | Pure C (`nitrosat.c`) + Lua wrapper | Python + C + CUDA backend |
| **Interface** | CLI / Lua scripting | REST API (`/v1/solve`, `/v1/diagnose`, `/v1/schedule`) |
| **Input Format** | DIMACS CNF files | CNF, Boolean expressions, XOR+CNF hybrid, Q-SAT, scheduling DSL |
| **Output** | Assignment + satisfaction % | Assignment + satisfaction % + **violation diagnostics** + **variable blame** + billing |
| **Hardware** | Single CPU core | CPU → L4 GPU → H100 GPU (tiered) |
| **Scale Limit** | ~1M clauses (LuaJIT 2GB RAM) | 8M+ clauses (H100), 1M+ variables |
| **Concurrency** | Single-threaded | 2-4 concurrent NP-complete solves per tier |

---

## What Navokoj Adds Beyond NitroSAT

### 1. Diagnostic Intelligence (DEFEKT)

NitroSAT tells you *what* it solved. Navokoj tells you *whether you should even try*, *why it failed*, and *whose fault it is*.

- **Pre-solve diagnostics** (`/v1/diagnose`): Predicts solvability score (0–100) before running the solver ; avoids wasting compute on likely-UNSAT instances.
- **Violation summary**: Identifies which constraints were violated and *why*.
- **Variable blame attribution**: Pinpoints the specific variables causing conflicts ; turns "UNSAT" from a dead end into a debugging tool.

NitroSAT returns `Unsatisfied: 3`. Navokoj returns `Variable 1 appears in 3 violated clauses; constraint [-1, -2] conflicts with clauses [1] and [2]`.

### 2. ⚡ Multi-Engine Selection

NitroSAT has one engine (the Langevin flow with BAHA). Navokoj exposes **four engines** optimized for different regimes:

| Engine | Role | Sweet Spot |
|---|---|---|
| **Nano** | Ultra-fast, real-time | Massive scale (N=100k+), sub-100ms |
| **Mini** | Balanced general-purpose | Medium problems, good quality |
| **Pro** | Deep optimization | Complex problems, 100% target |
| **Q-SAT** | N-ary satisfaction | Multi-valued logic (Sudoku, graph coloring with k states) |

Plus `auto` routing that selects the best engine based on problem structure.

### 3. 🔧 Production Features NitroSAT Doesn't Have

| Feature | Details |
|---|---|
| **Anytime algorithm** | Returns best-so-far on timeout ; never returns empty |
| **Timeout budget control** | `timeout_budget_seconds` for real-time SLAs |
| **Batch solving** | Multiple problems in one API call (30+ solves/second) |
| **Constraint weighting** | Soft vs hard constraints with priority weights |
| **Boolean expression parser** | Natural syntax (`a -> b`, `x <-> y`, `!z`) ; no CNF conversion needed |
| **XOR-native hybrid mode** | First-class XOR constraints without CNF explosion |
| **QBF / PSPACE support** | Quantified Boolean Formulas (`forall x exists y`) |
| **Schedule API** | Domain-specific endpoint for workforce scheduling |
| **Billing per solve** | Pay-per-use: $0.01/min (CPU) to $1/min (H100) |

### 4. GPU Acceleration

NitroSAT runs on a single CPU core. Navokoj deploys on:

- **Standard (CPU)**: 5k vars, 35k clauses ; free tier
- **L4 GPU**: 100k vars, 300k clauses ; $0.25 + $0.10/min
- **H100 GPU**: 1M vars, 8M clauses ; $1.50 + $1.00/min, ensemble of 500 parallel jobs

The H100 tier is what produces the headline results: R(5,5,5) N=52 in 17 minutes, R(6,6) N=35 in 24 minutes, Kubernetes 2M-clause placement at 100%.

### 5. Industrial Track Record

NitroSAT's benchmarks are research-grade (358 instances, ~300 variables typical). Navokoj has been tested on production workloads:

| Benchmark | NitroSAT | Navokoj |
|---|---|---|
| **SAT 2024 Industrial Track** | ; | 92.57% (4,199 problems) |
| **Ramsey R(5,5,5) N=52** | ; | 100% (7.8M clauses, 17 min) |
| **K8s Placement (2M clauses)** | ; | 100% (17 min, ~$2) |
| **129-SAT (1M clauses)** | ; | 100% (9-10 min) |
| **Random 3-SAT 1M vars** | ; | 92.15% (4.26M clauses) |
| **PSPACE (Sokoban, Pebbling)** | ; | 21/21 tests passed |
| **Reversible Pebbling 2.4M clauses** | ; | 97.55% |

### 6. Production Hardening

Things NitroSAT doesn't worry about but Navokoj must:

- **30-minute timeout guarantee** with best partial result on timeout
- **Request tracing** and billing infrastructure
- **Graceful degradation** ; never returns empty, always returns best-effort
- **Rate limiting** and concurrency control (NP-complete solves are expensive)
- **Enterprise SLA roadmap** ($100 base + $10/hour for dedicated 3-5 hour runs)

---

## What They Share

Both NitroSAT and Navokoj share the **same core foundation**:

- **Continuous physics-based solver**
- **Prime-weighted clauses**
- **Heat diffusion** for landscape smoothing
- **Entropy regularization** for symmetry breaking
- **Topology tracking** for structural repair
- **Smart phase-transition jumps (BAHA)**
- **Label Independence** (0.0000% permutation variance)

The core physics is identical. Navokoj adds the engineering layer: GPU parallelism, multi-engine routing, diagnostic intelligence, domain-specific APIs, and production guarantees.

---

## The Business Model

| | NitroSAT | Navokoj |
|---|---|---|
| **License** | Apache 2.0 (open) | Commercial API |
| **Pricing** | Free | Pay-per-solve ($0.01 – $1/min) |
| **Target** | Researchers, mathematicians, auditors | Engineers, enterprises, production systems |
| **Why it exists** | Audit the math, advance the theory | Run the math at scale, make money |

**NitroSAT is how you prove the physics works. Navokoj is how you ship it.**






