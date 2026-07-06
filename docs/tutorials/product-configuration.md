# Product Configuration with NitroSAT V3

## Problem

Configure a product by choosing a database, cloud, and service tier. PostgreSQL
requires Pro, while GCP is incompatible with MySQL. Preferences reward the
recommended configuration.

## Constraint Model

Each option is Boolean. Exactly-one clauses choose one option from each family.
Binary incompatibility clauses encode cross-feature rules. Weighted soft unit
clauses encode commercial recommendations or customer preferences.

## CNF / WCNF Generation

```bash
python3 examples/v3/tutorial_models.py generate \
  product-configuration /tmp/product.wcnf
```

## NitroSAT Solve

```bash
src/c/v3/nitrosatv3 /tmp/product.wcnf --epochs 100 \
  --finisher-passes 1000 --solution /tmp/product.sol > /tmp/product.json
```

## Decode Solution

```bash
python3 examples/v3/tutorial_models.py validate \
  product-configuration /tmp/product.sol
```

## Visualize Result

```text
Database  PostgreSQL
Cloud     GCP
Tier      Pro
```

## Complexity Discussion

The variable count is the total number of feature options. Exactly-one costs
are quadratic within each option family; cross-feature compatibility is
usually sparse. Large configurators should preserve rule-to-clause mappings so
conflicts can be explained in business language.
