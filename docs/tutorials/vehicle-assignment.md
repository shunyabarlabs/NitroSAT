# Vehicle Assignment with NitroSAT V3

## Problem

Assign Medical, Groceries, and Furniture deliveries to a Van or Bike. Furniture
cannot use the Bike, and Medical plus Furniture cannot both consume Van
capacity. Preferences reward efficient vehicle choices.

## Constraint Model

`X(delivery,vehicle)` selects one vehicle for a delivery. Exactly-one clauses
cover every delivery. Unit and pairwise clauses encode capability and capacity.
Soft clauses reward preferred assignments.

## CNF / WCNF Generation

```bash
python3 examples/v3/tutorial_models.py generate \
  vehicle-assignment /tmp/vehicles.wcnf
```

## NitroSAT Solve

```bash
src/c/v3/nitrosatv3 /tmp/vehicles.wcnf --epochs 100 \
  --finisher-passes 1000 --solution /tmp/vehicles.sol > /tmp/vehicles.json
```

## Decode Solution

```bash
python3 examples/v3/tutorial_models.py validate \
  vehicle-assignment /tmp/vehicles.sol
```

## Visualize Result

```text
Medical    -> Bike
Groceries  -> Bike
Furniture  -> Van
```

## Complexity Discussion

With `D` deliveries and `V` vehicles, the direct model has `D×V` variables and
`O(D×V²)` exactly-one clauses. Capacity greater than one is commonly encoded
with cardinality networks or sequential counters instead of enumerating every
invalid combination.
