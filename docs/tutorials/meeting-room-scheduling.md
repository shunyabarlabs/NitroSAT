# Meeting Room Scheduling with NitroSAT V3

## Problem

Assign Planning, Sales, and AllHands meetings to Atlas or Birch. Planning and
Sales overlap, so they cannot share a room. AllHands is too large for Birch.

## Constraint Model

Each `X(meeting,room)` variable represents one room assignment. Exactly-one
clauses assign every meeting. Overlap clauses prohibit the same room for
simultaneous meetings. A unit clause prohibits AllHands in Birch. Soft clauses
express room preferences.

## CNF / WCNF Generation

```bash
python3 examples/v3/tutorial_models.py generate \
  meeting-room-scheduling /tmp/rooms.wcnf
```

## NitroSAT Solve

```bash
src/c/v3/nitrosatv3 /tmp/rooms.wcnf --epochs 100 \
  --finisher-passes 1000 --solution /tmp/rooms.sol > /tmp/rooms.json
```

## Decode Solution

```bash
python3 examples/v3/tutorial_models.py validate \
  meeting-room-scheduling /tmp/rooms.sol
```

## Visualize Result

```text
Planning  -> Atlas
Sales     -> Birch
AllHands  -> Atlas
```

## Complexity Discussion

For `M` meetings and `R` rooms, the assignment layer uses `M×R` variables.
Every overlap edge contributes `R` binary clauses. Large systems should first
compute the interval-overlap graph rather than compare every meeting pair.
