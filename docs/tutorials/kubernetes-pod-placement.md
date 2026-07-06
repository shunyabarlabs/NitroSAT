# Kubernetes Pod Placement with NitroSAT V3

## Problem

Place API, Worker, Database, and Inference pods on `edge-a`, `edge-b`, or
`gpu`. Database and API require anti-affinity; Database cannot use `edge-a`;
Inference requires the GPU node.

## Constraint Model

`X(pod,node)` chooses a node. Exactly-one clauses place every pod. Unit clauses
represent node selectors and hardware requirements. Pairwise clauses implement
anti-affinity. Soft clauses reward preferred locality.

## CNF / WCNF Generation

```bash
python3 examples/v3/tutorial_models.py generate \
  kubernetes-pod-placement /tmp/pods.wcnf
```

## NitroSAT Solve

```bash
src/c/v3/nitrosatv3 /tmp/pods.wcnf --epochs 100 \
  --finisher-passes 1000 --solution /tmp/pods.sol > /tmp/pods.json
```

## Decode Solution

```bash
python3 examples/v3/tutorial_models.py validate \
  kubernetes-pod-placement /tmp/pods.sol
```

## Visualize Result

```text
edge-a  [api]
edge-b  [worker, database]
gpu     [inference]
```

## Complexity Discussion

For `P` pods and `N` nodes, placement uses `P×N` variables. Exactly-one clauses
cost `O(P×N²)`. Affinity rules add clauses per affected pod pair and node.
Production models also need cardinality encodings for CPU, memory, and GPU
capacity.
