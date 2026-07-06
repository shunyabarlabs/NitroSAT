# NitroSAT Executable Tutorial Library

Every tutorial follows the same path:

```text
Problem → Constraint Model → CNF/WCNF Generation → NitroSAT Solve
        → Decode Solution → Visualize Result → Complexity Discussion
```

| Tutorial | Encoding | Executable model |
|---|---|---|
| [Sudoku and Workforce Scheduling](../TUTORIAL_SUDOKU_AND_WORKFORCE_V3.md) | CNF + WCNF | `sudoku.py`, `workforce.py` |
| [University Timetabling](university-timetabling.md) | WCNF | `university-timetabling` |
| [Vehicle Assignment](vehicle-assignment.md) | WCNF | `vehicle-assignment` |
| [Meeting Room Scheduling](meeting-room-scheduling.md) | WCNF | `meeting-room-scheduling` |
| [Graph Coloring](graph-coloring.md) | CNF | `graph-coloring` |
| [Kubernetes Pod Placement](kubernetes-pod-placement.md) | WCNF | `kubernetes-pod-placement` |
| [Product Configuration](product-configuration.md) | WCNF | `product-configuration` |
| [Exam Scheduling](exam-scheduling.md) | WCNF | `exam-scheduling` |
| [Nurse Rostering](nurse-rostering.md) | WCNF | `nurse-rostering` |

Generate any shared model with:

```bash
python3 examples/v3/tutorial_models.py generate MODEL /tmp/model.wcnf
```

Decode and independently validate a solution with:

```bash
python3 examples/v3/tutorial_models.py validate MODEL /tmp/model.sol
```

The validator enumerates every assignment in these intentionally small
examples. It checks all hard clauses and confirms that WCNF preference scores
match the independently computed optimum.

CI executes the complete library through `tests/v3/test_tutorials.sh`.
