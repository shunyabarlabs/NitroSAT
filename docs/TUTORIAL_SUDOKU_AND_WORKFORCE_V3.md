# Tutorial: Solve Sudoku and a Workforce Problem with NitroSAT V3

This tutorial shows the complete path from a familiar problem to a NitroSAT V3
result:

1. Model decisions as Boolean variables.
2. Translate rules into CNF or weighted CNF.
3. Run NitroSAT V3.
4. Check the JSON result.
5. Convert the Boolean assignment back into a useful answer.

The repository includes runnable helpers under `examples/v3/`.

## Prerequisites

Build V3 from the repository root:

```bash
make -C src/c v3
```

The examples require Python 3. `jq` is convenient for reading JSON but is not
required.

## Part 1: 9×9 Sudoku

### Model the decisions

For every row `r`, column `c`, and digit `d`, create one Boolean variable:

```text
X(r,c,d) = true when cell (r,c) contains digit d
```

There are `9 × 9 × 9 = 729` variables. The encoder maps them to DIMACS numbers:

```text
variable(r,c,d) = 81 × (r - 1) + 9 × (c - 1) + d
```

For example, variable `1` means row 1, column 1 contains digit 1; variable `9`
means the same cell contains digit 9.

### Translate Sudoku rules into clauses

#### Every cell has at least one digit

For cell `(1,1)`:

```text
X(1,1,1) OR X(1,1,2) OR ... OR X(1,1,9)
```

This becomes one nine-literal clause.

#### Every cell has at most one digit

For every pair of different digits in a cell:

```text
NOT X(1,1,1) OR NOT X(1,1,2)
```

If digit 1 is selected, digit 2 cannot also be selected.

#### Every row, column, and box contains every digit

For digit 7 in row 1:

```text
X(1,1,7) OR X(1,2,7) OR ... OR X(1,9,7)
```

Equivalent clauses cover each column and each 3×3 box. Because every cell has
exactly one digit, requiring every digit at least once also forces every digit
to occur exactly once in each nine-cell region.

#### Given cells are unit clauses

If the puzzle says row 1, column 1 is 5:

```text
X(1,1,5)
```

### Encode the included puzzle

The input uses `0` for empty cells:

```text
530070000
600195000
098000060
800060003
400803001
700020006
060000280
000419005
000080079
```

Generate DIMACS CNF:

```bash
python3 examples/v3/sudoku.py encode \
  examples/v3/puzzle.txt /tmp/sudoku.cnf
```

The generated formula contains 729 variables and 3,270 clauses for this
puzzle: 3,240 general Sudoku clauses plus 30 clues.

### Solve it with V3

Sudoku is small enough for exact fallback, while the indexed finisher can often
find a solution before exact solving is needed:

```bash
src/c/v3/nitrosatv3 /tmp/sudoku.cnf \
  --epochs 20 \
  --indexed-finisher \
  --indexed-flips 100000 \
  --indexed-checkpoint 500 \
  --exact \
  --exact-memory-mb 256 \
  --exact-max-conflicts 1000000 \
  --solution /tmp/sudoku.sol \
  > /tmp/sudoku-result.json
```

Inspect the result before decoding:

```bash
python3 -c 'import json; print(json.load(open("/tmp/sudoku-result.json")))'
```

Proceed only when `status` is `SATISFIED` and `solved` is `true`. A partial
assignment file is not a valid Sudoku solution.

### Decode the grid

```bash
python3 examples/v3/sudoku.py decode /tmp/sudoku.sol
```

Expected solution:

```text
+-------+-------+-------+
| 5 3 4 | 6 7 8 | 9 1 2 |
| 6 7 2 | 1 9 5 | 3 4 8 |
| 1 9 8 | 3 4 2 | 5 6 7 |
+-------+-------+-------+
| 8 5 9 | 7 6 1 | 4 2 3 |
| 4 2 6 | 8 5 3 | 7 9 1 |
| 7 1 3 | 9 2 4 | 8 5 6 |
+-------+-------+-------+
| 9 6 1 | 5 3 7 | 2 8 4 |
| 2 8 7 | 4 1 9 | 6 3 5 |
| 3 4 5 | 2 8 6 | 1 7 9 |
+-------+-------+-------+
```

### Detect an impossible puzzle

If two fixed cells in one row are both 5, the formula is UNSAT. With `--exact`,
V3 can distinguish an impossible puzzle from a heuristic that merely failed to
find a solution:

```bash
src/c/v3/nitrosatv3 impossible.cnf \
  --exact --proof /tmp/impossible.drat
```

For formal use, validate the DRAT certificate with an independent proof
checker.

## Part 2: Workforce shift assignment

Sudoku is pure SAT: every rule is mandatory. Business optimization usually has
both mandatory rules and preferences, so weighted CNF is a better fit.

The included example assigns four employees to three shifts:

```text
Employees: Alice, Bob, Chen, Devi
Shifts:    Morning, Afternoon, Night
```

### Model the decisions

Create one variable for every employee-shift combination:

```text
X(employee, shift) = true when the employee works that shift
```

Four employees and three shifts produce twelve variables.

### Hard business rules

The example requires:

- Every shift has exactly one employee.
- An employee works at most one shift.
- Alice cannot work Night.
- Bob cannot work Morning.
- Alice is not qualified for Afternoon.

Hard clauses receive weight `100`, the WCNF `top` value. Violating any of these
rules makes an assignment infeasible.

For Morning coverage:

```text
X(Alice,Morning) OR X(Bob,Morning) OR
X(Chen,Morning) OR X(Devi,Morning)
```

For Alice working at most one of Morning and Night:

```text
NOT X(Alice,Morning) OR NOT X(Alice,Night)
```

### Soft preferences

Preferences receive weights below `top`:

| Preference | Weight |
|---|---:|
| Alice prefers Morning | 8 |
| Bob prefers Afternoon | 7 |
| Chen prefers Night | 10 |
| Devi prefers Morning | 4 |

NitroSAT minimizes the total weight of unsatisfied soft clauses. Equivalently,
this example tries to maximize the preference score while satisfying every
hard rule.

### Generate the WCNF

```bash
python3 examples/v3/workforce.py generate /tmp/workforce.wcnf
```

### Solve the scheduling problem

```bash
src/c/v3/nitrosatv3 /tmp/workforce.wcnf \
  --epochs 100 \
  --finisher-passes 1000 \
  --finisher-max-clauses 10000 \
  --solution /tmp/workforce.sol \
  > /tmp/workforce-result.json
```

The indexed finisher is currently CNF-only, so this WCNF example uses the
streaming weighted optimizer and streamed make/break finisher.

Inspect the business metrics:

```bash
python3 - <<'PY'
import json
result = json.load(open('/tmp/workforce-result.json'))
for field in ('feasible', 'hard_unsatisfied', 'soft_unsatisfied',
              'soft_cost', 'total_soft_weight'):
    print(f'{field}: {result[field]}')
PY
```

For an operational schedule, require:

```text
feasible = true
hard_unsatisfied = 0
```

`solved=false` can still be a valid WCNF business result: it may mean every
hard rule is satisfied while mutually incompatible preferences remain. V3 does
not currently prove that a reported soft cost is globally optimal.

### Decode the schedule

```bash
python3 examples/v3/workforce.py decode /tmp/workforce.sol
```

The expected high-preference assignment is:

```text
Morning    -> Alice
Afternoon  -> Bob
Night      -> Chen
preference score: 25/29
```

Devi's Morning preference cannot also be satisfied because each shift accepts
exactly one employee. Its weight of 4 becomes the unavoidable soft cost for
this assignment.

## Turning this into a real business application

A production workforce service would replace the fixed Python lists with data
from CSV, an HR system, or an API. Typical extensions include:

- required staffing levels per shift;
- employee skills and certifications;
- maximum weekly hours;
- rest periods between shifts;
- labour-law constraints;
- overtime cost;
- fairness and preference weights;
- absence and last-minute incremental constraints.

The application layer should retain a mapping from every business rule to the
clauses it generated. That allows the UI to translate solver output such as
“hard clause 184 is unsatisfied” back into “Night shift lacks a certified
operator.”

For large recurring schedules, build a persistent V3 store and append daily
constraints transactionally:

```bash
src/c/v3/nitrosatv3 base-model.wcnf \
  --store schedule.nsv3 --store-only

src/c/v3/nitrosatv3 --load-store schedule.nsv3 \
  --add todays-availability.wcnf --store-only

src/c/v3/nitrosatv3 --load-store schedule.nsv3 \
  --solution schedule.sol
```

## Production checklist

- Verify every returned assignment against the original rules.
- Treat `PARTIAL`, `UNKNOWN`, and `LIMIT` as distinct outcomes.
- Do not interpret heuristic failure as proof of UNSAT.
- Use exact mode only within explicit clause, memory, and conflict limits.
- Expose hard violations and soft cost separately to business users.
- Record the input hash, solver version, seed, options, and result JSON.
- Benchmark against established SAT/MaxSAT solvers on the actual customer
  workload.
