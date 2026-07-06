# University Timetabling with NitroSAT V3

## Problem

Place Math, Physics, and History into three teaching slots. Math and Physics
share students and cannot occur together. Departments also have preferred
slots.

## Constraint Model

`X(course,slot)` means a course uses a slot. Each course chooses exactly one
slot. Pairwise negative clauses prevent Math and Physics from sharing a slot.
Soft unit clauses reward Math on Monday morning, Physics on Monday afternoon,
and History on Tuesday morning.

## CNF / WCNF Generation

```bash
python3 examples/v3/tutorial_models.py generate \
  university-timetabling /tmp/university.wcnf
```

## NitroSAT Solve

```bash
src/c/v3/nitrosatv3 /tmp/university.wcnf --epochs 100 \
  --finisher-passes 1000 --finisher-max-clauses 10000 \
  --solution /tmp/university.sol > /tmp/university.json
```

Require `feasible=true` and `hard_unsatisfied=0`.

## Decode Solution

```bash
python3 examples/v3/tutorial_models.py validate \
  university-timetabling /tmp/university.sol
```

## Visualize Result

```text
Mon-AM  Math
Mon-PM  Physics
Tue-AM  History
```

## Complexity Discussion

For `C` courses and `T` slots, the direct encoding uses `C×T` variables and
`O(C×T²)` exactly-one clauses, plus one binary clause per conflict and slot.
Real deployments should add rooms, teachers, capacity, and multi-period events.
