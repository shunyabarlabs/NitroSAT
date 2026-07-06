# NitroSAT V3

V3 is the bounded-memory path for CNF and weighted CNF (WCNF) instances whose
clause set is larger than RAM.
It does not load clause offsets, clause state, or a variable-to-clause index.
Instead, DIMACS input is converted into an appendable binary stream and every
global solver pass reads that stream sequentially.

Build:

```sh
make -C src/c v3
```

Solve a base formula with incremental constraint files:

```sh
src/c/v3/nitrosatv3 base.cnf \
  --add constraints-1.cnf --add constraints-2.cnf \
  --store formula.nsv3 --epochs 20
```

Each `--add` file may be DIMACS `p cnf` or `p wcnf`. Its clauses are appended in
argument order, so CNF and WCNF inputs may share one store. Variable identifiers
may extend the base formula's declared range; V3 uses the maximum declared
variable count across all inputs. Parsing is transactional: a malformed input
removes the incomplete store.

For `p wcnf variables clauses top`, weights greater than or equal to `top` are
stored as hard clauses; lower weights are soft costs. The older WCNF form
without `top` is accepted and treats every clause as soft. Weights must be
positive signed-64-bit integers. JSON reports `hard_unsatisfied`,
`soft_unsatisfied`, `soft_cost`, and `total_soft_weight`. Optimization retains
the best assignment lexicographically: hard violations first, then soft cost.
Modern top-less WCNF clauses prefixed with `h` are also supported. A WCNF run
returns success when all hard constraints are satisfied even when soft clauses
remain; `solved` specifically means every hard and soft clause is satisfied,
while `feasible` means no hard clause is violated. V3 is heuristic and does not
claim an optimality proof for the reported soft cost.

## Streamed make/break finisher

For small difficult formulas, V3 performs globally informed local search
without constructing a variable-to-clause index. Each finisher pass streams the
entire store, accumulates make/break scores and degree-normalized heat forces in
O(variables) memory, and applies a conflict-screened batch of positive-gain
flips. At local minima it uses a reservoir-sampled unsatisfied clause and a
bounded WalkSAT noise move.

```sh
nitrosatv3 parity.cnf --epochs 2000 \
  --finisher-passes 10000 --finisher-batch-flips 1 \
  --finisher-max-clauses 100000
```

Planted graph-coloring instances benefit from the default 64-flip heat batch:

```sh
nitrosatv3 planted.cnf --epochs 200 --finisher-passes 500 \
  --finisher-batch-flips 64 --finisher-max-clauses 500000
```

The default finisher ceiling is 100,000 clauses. Larger instances skip this
phase unless explicitly overridden, preventing accidental repeated scans of a
billion-clause store.

## Disk-indexed finisher

Large CNF instances can opt into a temporary disk-backed variable occurrence
index. This supports WalkSAT make/break flips without rescanning the entire
formula for every move:

```sh
nitrosatv3 planted.cnf --epochs 20 --indexed-finisher \
  --indexed-flips 200000 --indexed-checkpoint 1000
```

The index uses eight bytes per literal on disk and an eight-byte offset table
per variable. Index construction temporarily maps the index, so
`bounded_memory_mb` includes that peak bound and `indexed_index_mb` reports its
contribution. Full streamed verification runs at every checkpoint and only
verified improvements are retained. This path currently supports CNF only;
WCNF requests are rejected rather than applying incorrect weighted semantics.

An existing store can be solved again or extended transactionally without
re-reading the base CNF:

```sh
nitrosatv3 --load-store formula.nsv3 --add late-constraints.cnf \
  --solution assignment.sol
```

If an added file is malformed, V3 truncates the store back to its original byte
length and restores its original header before returning an error.

Use `--store-only` to build and validate a store without allocating solver
state. If `--store` is omitted, V3 creates and removes a temporary store.

## Scaling invariants

- Clause and literal totals are unsigned 64-bit values.
- Clause records contain a 32-bit literal count, flags, a 64-bit weight, and
  signed 32-bit literals. There is no per-clause offset array.
- Prime-inspired clause weights are computed from the clause number; there is
  no prime table.
- Gradients use a sparse touched-variable batch, avoiding a full variable-array
  clear for each batch.
- Local repair retains at most `--active-clauses` and `--active-literals`.
- Final satisfaction is always established by a complete streamed pass.
- Best assignments are retained across gradient, repair, and finisher phases.

The binary store is little-endian and currently internal to V3. Its 40-byte
header contains magic, version, header size, variable count, format flags,
clause count, and literal count. A version-2 record is `uint32 literal_count`,
`uint32 flags`, `uint64 weight`, followed by that many `int32` literals. V1
stores are intentionally rejected because they did not preserve WCNF data.
