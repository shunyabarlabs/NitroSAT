# Graph Coloring with NitroSAT V3

## Problem

Color a five-vertex cycle with Red, Green, or Blue so adjacent vertices use
different colors.

## Constraint Model

`X(vertex,color)` means a vertex has a color. Exactly-one clauses choose one
color per vertex. For every edge `(u,v)` and color `c`, the clause
`NOT X(u,c) OR NOT X(v,c)` prevents a collision.

## CNF / WCNF Generation

```bash
python3 examples/v3/tutorial_models.py generate \
  graph-coloring /tmp/coloring.cnf
```

## NitroSAT Solve

```bash
src/c/v3/nitrosatv3 /tmp/coloring.cnf --epochs 20 \
  --indexed-finisher --indexed-flips 10000 --indexed-checkpoint 100 \
  --exact --solution /tmp/coloring.sol > /tmp/coloring.json
```

## Decode Solution

```bash
python3 examples/v3/tutorial_models.py validate \
  graph-coloring /tmp/coloring.sol
```

## Visualize Result

```text
Vertex 1 — Green — Vertex 2 — Red — Vertex 3 — Blue
    |                                      |
Vertex 5 — Red — Vertex 4 — Green --------+
```

## Complexity Discussion

For `N` vertices, `E` edges, and `K` colors, this encoding uses `N×K`
variables, `O(N×K²)` exactly-one clauses, and `E×K` edge clauses. Sparse graph
encodings scale linearly with edges once `K` is fixed.
