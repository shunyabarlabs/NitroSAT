# Nurse Rostering with NitroSAT V3

## Problem

Cover day and night shifts over two days using Asha, Bea, and Carla. Every
shift needs one nurse, adjacent shifts must respect rest rules, and Asha is not
night-qualified. Preferences reward desirable assignments.

## Constraint Model

`X(shift,nurse)` selects a nurse. Exactly-one clauses cover shifts. Pairwise
clauses prevent a nurse working prohibited adjacent shifts. Unit clauses encode
night qualifications. Soft clauses represent preferences and fairness goals.

## CNF / WCNF Generation

```bash
python3 examples/v3/tutorial_models.py generate \
  nurse-rostering /tmp/nurses.wcnf
```

## NitroSAT Solve

```bash
src/c/v3/nitrosatv3 /tmp/nurses.wcnf --epochs 100 \
  --finisher-passes 1000 --solution /tmp/nurses.sol > /tmp/nurses.json
```

## Decode Solution

```bash
python3 examples/v3/tutorial_models.py validate \
  nurse-rostering /tmp/nurses.sol
```

## Visualize Result

```text
Mon-Day    Asha
Mon-Night  Bea
Tue-Day    Carla
Tue-Night  Bea
```

## Complexity Discussion

For `S` shifts and `N` nurses, direct assignment uses `S×N` variables and
`O(S×N²)` coverage clauses. Rest, qualification, and fairness rules dominate
real models. Weekly hour bounds should use compact cardinality encodings rather
than enumerating invalid schedules.
