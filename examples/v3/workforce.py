#!/usr/bin/env python3
"""Generate and decode a small workforce-assignment WCNF example."""

import argparse


EMPLOYEES = ["Alice", "Bob", "Chen", "Devi"]
SHIFTS = ["Morning", "Afternoon", "Night"]
TOP = 100
PREFERENCES = {
    ("Alice", "Morning"): 8,
    ("Bob", "Afternoon"): 7,
    ("Chen", "Night"): 10,
    ("Devi", "Morning"): 4,
}


def variable(employee, shift):
    return EMPLOYEES.index(employee) * len(SHIFTS) + SHIFTS.index(shift) + 1


def generate(path):
    clauses = []

    # Every shift must have exactly one employee.
    for shift in SHIFTS:
        candidates = [variable(employee, shift) for employee in EMPLOYEES]
        clauses.append((TOP, candidates, "hard: cover shift"))
        for left in range(len(candidates)):
            for right in range(left + 1, len(candidates)):
                clauses.append((TOP, [-candidates[left], -candidates[right]],
                                "hard: one employee per shift"))

    # No employee may work more than one shift.
    for employee in EMPLOYEES:
        for left in range(len(SHIFTS)):
            for right in range(left + 1, len(SHIFTS)):
                clauses.append((TOP, [-variable(employee, SHIFTS[left]),
                                      -variable(employee, SHIFTS[right])],
                                "hard: maximum one shift"))

    # Availability and qualification rules.
    clauses.append((TOP, [-variable("Alice", "Night")], "hard: Alice unavailable at night"))
    clauses.append((TOP, [-variable("Bob", "Morning")], "hard: Bob unavailable in morning"))
    clauses.append((TOP, [-variable("Alice", "Afternoon")], "hard: Alice not afternoon-qualified"))

    # A satisfied soft unit clause earns its preference weight.
    for (employee, shift), weight in PREFERENCES.items():
        clauses.append((weight, [variable(employee, shift)], "soft preference"))

    with open(path, "w", encoding="ascii") as output:
        output.write("c Workforce assignment example\n")
        output.write(f"p wcnf {len(EMPLOYEES) * len(SHIFTS)} {len(clauses)} {TOP}\n")
        for weight, literals, comment in clauses:
            output.write(f"c {comment}\n")
            output.write(f"{weight} " + " ".join(map(str, literals)) + " 0\n")
    print(f"wrote {path}: 12 variables, {len(clauses)} weighted clauses")


def decode(path):
    assignment = {}
    with open(path, encoding="ascii") as source:
        for token in source.read().split():
            literal = int(token)
            if literal == 0:
                break
            assignment[abs(literal)] = literal > 0
    chosen = []
    score = 0
    for employee in EMPLOYEES:
        for shift in SHIFTS:
            if assignment.get(variable(employee, shift), False):
                chosen.append((shift, employee))
                score += PREFERENCES.get((employee, shift), 0)
    for shift, employee in sorted(chosen, key=lambda item: SHIFTS.index(item[0])):
        print(f"{shift:10} -> {employee}")
    print(f"preference score: {score}/{sum(PREFERENCES.values())}")


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    generate_parser = subparsers.add_parser("generate")
    generate_parser.add_argument("wcnf")
    decode_parser = subparsers.add_parser("decode")
    decode_parser.add_argument("solution")
    arguments = parser.parse_args()
    if arguments.command == "generate":
        generate(arguments.wcnf)
    else:
        decode(arguments.solution)


if __name__ == "__main__":
    main()
