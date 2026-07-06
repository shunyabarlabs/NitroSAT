# Exam Scheduling with NitroSAT V3

## Problem

Schedule Math, Physics, Chemistry, and History across Monday, Tuesday, and
Wednesday. Exams sharing students cannot use the same day, while departments
have preferred days.

## Constraint Model

`X(exam,day)` assigns one day to an exam. Exactly-one clauses schedule every
exam. Each student-conflict edge produces one binary exclusion clause per day.
Soft clauses reward preferred dates.

## CNF / WCNF Generation

```bash
python3 examples/v3/tutorial_models.py generate \
  exam-scheduling /tmp/exams.wcnf
```

## NitroSAT Solve

```bash
src/c/v3/nitrosatv3 /tmp/exams.wcnf --epochs 100 \
  --finisher-passes 1000 --solution /tmp/exams.sol > /tmp/exams.json
```

## Decode Solution

```bash
python3 examples/v3/tutorial_models.py validate \
  exam-scheduling /tmp/exams.sol
```

## Visualize Result

```text
Monday     Math, History
Tuesday    Physics
Wednesday  Chemistry
```

## Complexity Discussion

With `E` exams, `T` slots, and `C` conflict edges, the encoding uses `E×T`
variables, `O(E×T²)` exactly-one clauses, and `C×T` collision clauses. Room
capacity and multi-session exams add assignment dimensions or cardinality
constraints.
