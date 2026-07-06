# NitroSAT Business Use Cases and Commercialization Guide

## Executive summary

NitroSAT can be commercialized as an optimization engine, a managed solving
API, or an enterprise integration product. The strongest initial positioning
is not “another general-purpose SAT solver.” It is:

> A fast constraint-optimization service that produces useful assignments for
> large scheduling, allocation, configuration, and graph problems, including
> workloads that exceed normal in-memory limits.

NitroSAT V2 and V3 serve different operating requirements:

| Product | Best use | Main tradeoff |
|---|---|---|
| V2 | Maximum heuristic strength on machines with sufficient RAM | Holds global clause and occurrence data in memory |
| V3 streaming | Very large CNF/WCNF workloads with predictable memory usage | Pure streaming repair is less aggressive |
| V3 indexed finisher | Large CNF workloads needing stronger local search | Uses a temporary disk index and more I/O |
| V3 exact mode | Smaller CNF workloads requiring SAT/UNSAT certainty | CDCL memory and runtime can grow rapidly |

The Apache 2.0 license permits commercial use. Revenue therefore comes from
hosting, reliability, integrations, support, proprietary workflow software,
and domain expertise—not from assuming customers must pay for source access.

## Where customers obtain value

### Workforce and production scheduling

Represent shifts, skills, availability, legal constraints, machine capacity,
and preferences as hard and soft clauses.

- Hard clauses: safety rules, certifications, capacity, required coverage.
- Soft clauses: employee preferences, overtime reduction, changeover cost.
- Commercial output: a feasible schedule plus its remaining soft-constraint
  cost and an explanation of unresolved constraints.

Likely buyers include manufacturers, hospitals, logistics operators, call
centres, universities, and scheduling-software vendors.

### Timetabling and resource allocation

Use NitroSAT for classrooms, examinations, meeting rooms, compute clusters,
vehicle assignments, warehouse slots, and maintenance windows.

The product should accept business objects through JSON or CSV and perform the
CNF/WCNF translation internally. Most customers should never need to author
DIMACS files.

### Product and infrastructure configuration

Encode compatibility, dependency, entitlement, regional, and capacity rules.
Examples include telecom plans, cloud deployments, industrial bills of
materials, software package selection, and configurable products.

The sellable feature is a configuration API that either returns a valid
configuration or identifies conflicting rule groups.

### Graph optimization

Graph coloring and related encodings apply to frequency assignment, register
allocation, network segmentation, test planning, and conflict-free resource
allocation. The V3 indexed finisher is particularly relevant when a graph
encoding contains millions of binary clauses.

### EDA and verification assistance

Potential applications include bounded verification, test-vector generation,
constraint debugging, and approximate MaxSAT preprocessing. Formal workflows
must use exact mode or an independently validated downstream solver before
claiming proof of correctness.

### Constraint-health analytics

Even when no complete solution is found, NitroSAT can expose satisfaction,
hard violations, soft cost, variable importance, and difficult clauses. This
supports a “constraint observability” product for diagnosing over-constrained
business rules.

## Choosing V2 or V3

### Use V2 when

- The formula fits comfortably in RAM.
- Fast global repair is more important than predictable memory.
- The objective is a high-quality assignment or approximate MaxSAT result.
- Multiple restarts and aggressive local search are valuable.

Build and run:

```bash
make -C src/c v2
src/c/v2/nitrosatv2 problem.cnf --json > result.json
```

Useful operating modes:

```bash
# Human-readable progress
src/c/v2/nitrosatv2 problem.cnf --cinematic

# Disable optional repair behavior for controlled comparisons
src/c/v2/nitrosatv2 problem.cnf --no-dcw --no-topo --json
```

Treat V2 as a heuristic unless its assignment is independently checked. A
partial result can still have commercial value for MaxSAT-style objectives,
but it must not be described as a proof that the input is UNSAT.

### Use V3 streaming when

- Clause data may exceed available RAM.
- Inputs arrive incrementally.
- Predictable resident memory matters.
- CNF or weighted CNF must be supported.

```bash
make -C src/c v3

src/c/v3/nitrosatv3 base.cnf \
  --add incremental.cnf \
  --store workload.nsv3 \
  --solution assignment.sol
```

Persistent stores allow a service to ingest a base model once and append new
constraints transactionally:

```bash
src/c/v3/nitrosatv3 base.cnf --store workload.nsv3 --store-only

src/c/v3/nitrosatv3 --load-store workload.nsv3 \
  --add daily-constraints.cnf --store-only

src/c/v3/nitrosatv3 --load-store workload.nsv3 \
  --solution assignment.sol
```

### Use the V3 indexed finisher when

- The input is CNF rather than WCNF.
- Streaming optimization reaches a high-quality plateau.
- Repeated full-store finisher scans are too expensive.
- Temporary disk usage of roughly eight bytes per literal is acceptable.

```bash
src/c/v3/nitrosatv3 large.cnf \
  --epochs 5 \
  --indexed-finisher \
  --indexed-flips 200000 \
  --indexed-checkpoint 1000 \
  --solution assignment.sol
```

`--indexed-checkpoint` controls how often the complete formula is streamed to
verify the current state. Only globally verified improvements are retained.
Monitor `indexed_flips`, `indexed_index_mb`, `bounded_memory_mb`, satisfaction,
and total solve time.

### Use V3 exact mode when

- A smaller CNF instance needs a definitive SAT or UNSAT result.
- Resource limits can be defined explicitly.
- An UNSAT proof is required.

```bash
src/c/v3/nitrosatv3 problem.cnf \
  --exact \
  --exact-max-clauses 1000000 \
  --exact-memory-mb 1024 \
  --exact-max-conflicts 5000000 \
  --proof proof.drat \
  --solution assignment.sol
```

Interpret exact results explicitly:

- `SAT`: the returned assignment was verified against the complete formula.
- `UNSAT`: the exact solver derived a contradiction.
- `LIMIT` / `UNKNOWN`: a configured resource limit was reached.
- `NOT_RUN`: the heuristic already solved the formula or exact mode was not
  requested.

DRAT proof output currently supports CNF only. Validate certificates with an
independent proof checker before using them in formal or regulated workflows.

## Commercial product models

### Managed optimization API

Offer endpoints such as:

```text
POST /v1/models
POST /v1/models/{id}/constraints
POST /v1/models/{id}/solve
GET  /v1/jobs/{id}
GET  /v1/jobs/{id}/assignment
GET  /v1/jobs/{id}/diagnostics
```

Charge for compute time, clause volume, retained models, and service level.
The service can route jobs automatically:

1. V2 for in-memory heuristic workloads.
2. V3 streaming for large or incremental workloads.
3. V3 indexed finishing for large CNF plateaus.
4. Exact mode when requested and within configured limits.

### Vertical SaaS

Build a complete application for one domain, such as employee scheduling or
university timetabling. This usually supports stronger pricing than selling a
generic solver because customers pay for imports, rule modelling, workflows,
approvals, explanations, and integrations.

### Enterprise deployment

Sell an on-premises or private-cloud package containing:

- containerized solver workers;
- job queue and API gateway;
- usage and performance dashboard;
- SSO, audit logs, and role-based access;
- model/version management;
- support and upgrade commitments.

### Integration and optimization consulting

Charge for translating operational rules into CNF/WCNF, benchmarking against
the customer’s incumbent process, tuning solve policies, and integrating the
results into existing ERP, MES, EDA, or scheduling systems.

### Support and reliability subscriptions

Apache-licensed software can still support paid offerings for response-time
SLAs, tested releases, long-term support, security reviews, architecture help,
and production incident assistance.

## Illustrative pricing

These figures are starting hypotheses, not validated market prices:

| Offering | Possible price |
|---|---:|
| Developer API | $49–$299 per month plus usage |
| Team optimization service | $500–$2,500 per month |
| Vertical scheduling SaaS | $10–$50 per scheduled employee/resource monthly |
| Enterprise private deployment | $25,000–$150,000 annually |
| Proof/verification worker tier | Premium compute pricing per job |
| Integration pilot | $10,000–$50,000 fixed scope |
| Custom optimization engagement | $1,000–$2,500 per engineering day |

Usage pricing can combine:

- clauses processed;
- solver CPU seconds;
- peak memory or indexed-disk allocation;
- exact-mode conflicts;
- retained model storage;
- priority and concurrency.

Avoid pricing solely per SAT call because workloads vary by several orders of
magnitude.

## Recommended go-to-market plan

### Phase 1: prove one economic outcome

Choose one narrow problem where improvement is measurable—for example,
reducing scheduling violations or producing timetables faster. Collect 10–20
representative customer instances and compare:

- solution quality;
- hard violations;
- soft cost;
- runtime and memory;
- operational savings versus the existing process.

### Phase 2: package the workflow

Build a domain adapter that accepts normal business data, generates constraints,
runs the appropriate NitroSAT mode, verifies the assignment, and converts it
back into business objects. Add reproducible job records and downloadable
diagnostics.

### Phase 3: paid pilot

Offer a fixed-scope pilot with a clear acceptance test, such as:

- all mandatory constraints satisfied;
- at least 20% fewer preference violations;
- solve completed within a defined window;
- deployment within the customer’s security boundary.

### Phase 4: recurring product

Convert successful pilots into annual contracts covering hosted usage,
private deployment, support, or a vertical SaaS subscription.

## Production requirements before selling

- Maintain a benchmark suite split by domain and difficulty.
- Independently verify every returned assignment.
- Verify DRAT proofs with an external checker in CI.
- Report `UNKNOWN` honestly when limits are reached.
- Add cancellation, deadlines, quotas, and worker isolation.
- Protect uploaded formulas and generated assignments as customer-confidential
  data.
- Record solver version, seed, options, input hash, result, and verification.
- Compare against established SAT/MaxSAT solvers on the customer’s workload.
- Avoid universal “linear-time SAT solver” or guaranteed-solve claims; runtime
  and completion remain instance-dependent.
- Obtain legal advice for commercial contracts, privacy obligations, patent
  questions, export requirements, and Apache 2.0 notice compliance.

## Practical first product

The lowest-risk commercial starting point is a managed “constraint-solving
worker” with a small SDK:

1. Customer uploads CNF/WCNF or uses a domain adapter.
2. The service estimates size and selects V2 or V3.
3. A heuristic result is produced under a deadline.
4. The assignment is independently verified.
5. Optional exact solving runs within customer-defined limits.
6. The customer receives the assignment, quality metrics, diagnostics, and a
   reproducible job manifest.

This product monetizes operational reliability and ease of use while keeping
the open-source solver available as the adoption funnel.
