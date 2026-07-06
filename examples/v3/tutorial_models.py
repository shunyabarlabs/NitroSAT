#!/usr/bin/env python3
"""Executable CNF/WCNF models for the NitroSAT V3 tutorial library."""

import argparse
import itertools
from dataclasses import dataclass


TOP = 1000


@dataclass
class Model:
    title: str
    variables: list
    groups: list
    hard: list
    preferences: dict

    @property
    def weighted(self):
        return bool(self.preferences)


def model(title, groups, forbidden=(), prohibited=(), implications=(), preferences=None):
    variables = [item for group in groups for item in group]
    number = {item: index + 1 for index, item in enumerate(variables)}
    hard = []
    for group in groups:
        hard.append([number[item] for item in group])
        for left, right in itertools.combinations(group, 2):
            hard.append([-number[left], -number[right]])
    for left, right in forbidden:
        hard.append([-number[left], -number[right]])
    for item in prohibited:
        hard.append([-number[item]])
    for premise, consequence in implications:
        hard.append([-number[premise], number[consequence]])
    weighted = {number[item]: weight for item, weight in (preferences or {}).items()}
    return Model(title, variables, groups, hard, weighted)


def build_models():
    university_groups = [
        [f"{course} @ {slot}" for slot in ("Mon-AM", "Mon-PM", "Tue-AM")]
        for course in ("Math", "Physics", "History")
    ]
    university_forbidden = [
        (f"Math @ {slot}", f"Physics @ {slot}")
        for slot in ("Mon-AM", "Mon-PM", "Tue-AM")
    ]

    vehicle_groups = [
        [f"{delivery} -> {vehicle}" for vehicle in ("Van", "Bike")]
        for delivery in ("Medical", "Groceries", "Furniture")
    ]

    meeting_groups = [
        [f"{meeting} -> {room}" for room in ("Atlas", "Birch")]
        for meeting in ("Planning", "Sales", "AllHands")
    ]

    graph_groups = [
        [f"Vertex {vertex} = {color}" for color in ("Red", "Green", "Blue")]
        for vertex in range(1, 6)
    ]
    graph_edges = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)]
    graph_forbidden = [
        (f"Vertex {left} = {color}", f"Vertex {right} = {color}")
        for left, right in graph_edges for color in ("Red", "Green", "Blue")
    ]

    pod_groups = [
        [f"{pod} -> {node}" for node in ("edge-a", "edge-b", "gpu")]
        for pod in ("api", "worker", "database", "inference")
    ]

    product_groups = [
        ["Database = PostgreSQL", "Database = MySQL"],
        ["Cloud = AWS", "Cloud = GCP"],
        ["Tier = Basic", "Tier = Pro"],
    ]

    exam_groups = [
        [f"{exam} @ {slot}" for slot in ("Monday", "Tuesday", "Wednesday")]
        for exam in ("Math", "Physics", "Chemistry", "History")
    ]
    exam_conflicts = (("Math", "Physics"), ("Math", "Chemistry"),
                      ("Physics", "Chemistry"), ("Chemistry", "History"))

    nurse_groups = [
        [f"{shift} -> {nurse}" for nurse in ("Asha", "Bea", "Carla")]
        for shift in ("Mon-Day", "Mon-Night", "Tue-Day", "Tue-Night")
    ]

    return {
        "university-timetabling": model(
            "University Timetabling", university_groups, university_forbidden,
            preferences={"Math @ Mon-AM": 8, "Physics @ Mon-PM": 7,
                         "History @ Tue-AM": 5}),
        "vehicle-assignment": model(
            "Vehicle Assignment", vehicle_groups,
            forbidden=(("Medical -> Van", "Furniture -> Van"),),
            prohibited=("Furniture -> Bike",),
            preferences={"Medical -> Bike": 3, "Groceries -> Bike": 6,
                         "Furniture -> Van": 10}),
        "meeting-room-scheduling": model(
            "Meeting Room Scheduling", meeting_groups,
            forbidden=(("Planning -> Atlas", "Sales -> Atlas"),
                       ("Planning -> Birch", "Sales -> Birch")),
            prohibited=("AllHands -> Birch",),
            preferences={"Planning -> Atlas": 5, "Sales -> Birch": 7,
                         "AllHands -> Atlas": 10}),
        "graph-coloring": model("Graph Coloring", graph_groups, graph_forbidden),
        "kubernetes-pod-placement": model(
            "Kubernetes Pod Placement", pod_groups,
            forbidden=(("api -> edge-a", "database -> edge-a"),
                       ("api -> edge-b", "database -> edge-b"),
                       ("api -> gpu", "database -> gpu")),
            prohibited=("database -> edge-a", "inference -> edge-a",
                        "inference -> edge-b"),
            preferences={"api -> edge-a": 7, "worker -> edge-b": 6,
                         "database -> edge-b": 9, "inference -> gpu": 10}),
        "product-configuration": model(
            "Product Configuration", product_groups,
            forbidden=(("Database = PostgreSQL", "Tier = Basic"),
                       ("Cloud = GCP", "Database = MySQL")),
            preferences={"Database = PostgreSQL": 8, "Cloud = GCP": 5,
                         "Tier = Pro": 7}),
        "exam-scheduling": model(
            "Exam Scheduling", exam_groups,
            forbidden=[(f"{left} @ {slot}", f"{right} @ {slot}")
                       for left, right in exam_conflicts
                       for slot in ("Monday", "Tuesday", "Wednesday")],
            preferences={"Math @ Monday": 8, "Physics @ Tuesday": 6,
                         "Chemistry @ Wednesday": 7, "History @ Monday": 4}),
        "nurse-rostering": model(
            "Nurse Rostering", nurse_groups,
            forbidden=[(f"{first} -> {nurse}", f"{second} -> {nurse}")
                       for first, second in (("Mon-Day", "Mon-Night"),
                                             ("Mon-Night", "Tue-Day"),
                                             ("Tue-Day", "Tue-Night"))
                       for nurse in ("Asha", "Bea", "Carla")],
            prohibited=("Mon-Night -> Asha", "Tue-Night -> Asha"),
            preferences={"Mon-Day -> Asha": 7, "Mon-Night -> Bea": 9,
                         "Tue-Day -> Carla": 6, "Tue-Night -> Bea": 5}),
    }


MODELS = build_models()


def write_formula(item, path):
    soft = [[variable] for variable in item.preferences]
    clause_count = len(item.hard) + len(soft)
    with open(path, "w", encoding="ascii") as output:
        output.write(f"c {item.title} tutorial model\n")
        if item.weighted:
            output.write(f"p wcnf {len(item.variables)} {clause_count} {TOP}\n")
            for clause in item.hard:
                output.write(f"{TOP} " + " ".join(map(str, clause)) + " 0\n")
            for variable, weight in item.preferences.items():
                output.write(f"{weight} {variable} 0\n")
        else:
            output.write(f"p cnf {len(item.variables)} {clause_count}\n")
            for clause in item.hard:
                output.write(" ".join(map(str, clause)) + " 0\n")
    print(f"wrote {path}: {len(item.variables)} variables, {clause_count} clauses")


def read_assignment(path):
    assignment = {}
    with open(path, encoding="ascii") as source:
        for token in source.read().split():
            literal = int(token)
            if literal == 0:
                break
            assignment[abs(literal)] = literal > 0
    return assignment


def clause_satisfied(clause, assignment):
    return any(assignment.get(abs(literal), False) == (literal > 0)
               for literal in clause)


def score(item, assignment):
    hard_ok = all(clause_satisfied(clause, assignment) for clause in item.hard)
    earned = sum(weight for variable, weight in item.preferences.items()
                 if assignment.get(variable, False))
    return hard_ok, earned


def optimum(item):
    best = -1
    for values in itertools.product((False, True), repeat=len(item.variables)):
        assignment = {index + 1: value for index, value in enumerate(values)}
        hard_ok, earned = score(item, assignment)
        if hard_ok:
            best = max(best, earned)
    return best


def decode(item, path, validate=False):
    assignment = read_assignment(path)
    hard_ok, earned = score(item, assignment)
    print(item.title)
    print("=" * len(item.title))
    for group in item.groups:
        selected = [label for label in group
                    if assignment.get(item.variables.index(label) + 1, False)]
        print(selected[0] if selected else "UNASSIGNED")
    print(f"hard constraints: {'satisfied' if hard_ok else 'VIOLATED'}")
    if item.weighted:
        maximum = sum(item.preferences.values())
        print(f"preference score: {earned}/{maximum}")
    if validate:
        expected = optimum(item)
        optimal = hard_ok and (not item.weighted or earned == expected)
        print(f"independent optimum: {expected if item.weighted else 'n/a'}")
        print(f"validation: {'passed' if optimal else 'failed'}")
        if not optimal:
            raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("list")
    for command in ("generate", "decode", "validate"):
        child = subparsers.add_parser(command)
        child.add_argument("model", choices=sorted(MODELS))
        child.add_argument("path")
    arguments = parser.parse_args()
    if arguments.command == "list":
        for name in sorted(MODELS):
            print(f"{name:28} {MODELS[name].title}")
    elif arguments.command == "generate":
        write_formula(MODELS[arguments.model], arguments.path)
    else:
        decode(MODELS[arguments.model], arguments.path,
               validate=arguments.command == "validate")


if __name__ == "__main__":
    main()
