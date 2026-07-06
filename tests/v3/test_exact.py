#!/usr/bin/env python3
"""Deterministic integration tests for NitroSAT V3 exact mode."""

import itertools
import json
import os
import random
import subprocess
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
SOURCE = os.path.join(ROOT, "src/c/v3/nitrosatv3.c")
BINARY = os.path.join(tempfile.gettempdir(), "nitrosatv3-exact-test")


def build_solver():
    subprocess.run(
        ["gcc", "-O2", "-std=c99", "-Wall", "-Wextra", "-Werror", SOURCE, "-lm", "-o", BINARY],
        check=True,
    )


def write_cnf(path, num_vars, clauses):
    with open(path, "w", encoding="ascii") as fp:
        fp.write(f"p cnf {num_vars} {len(clauses)}\n")
        for clause in clauses:
            fp.write(" ".join(map(str, clause)) + " 0\n")


def brute_force_sat(num_vars, clauses):
    for values in itertools.product((False, True), repeat=num_vars):
        if all(any(values[abs(lit) - 1] == (lit > 0) for lit in clause) for clause in clauses):
            return True
    return False


def unit_conflict(clauses, assumptions=()):
    assignment = {}
    pending = list(assumptions)
    while True:
        while pending:
            lit = pending.pop()
            var, value = abs(lit), lit > 0
            if var in assignment and assignment[var] != value:
                return True
            assignment[var] = value
        changed = False
        for clause in clauses:
            unresolved = []
            satisfied = False
            for lit in clause:
                value = assignment.get(abs(lit))
                if value is None:
                    unresolved.append(lit)
                elif value == (lit > 0):
                    satisfied = True
                    break
            if satisfied:
                continue
            if not unresolved:
                return True
            if len(unresolved) == 1:
                pending.append(unresolved[0])
                changed = True
        if not changed:
            return False


def verify_rup_proof(original, proof_path):
    clauses = []
    for clause in original:
        normalized = set(clause)
        if any(-lit in normalized for lit in normalized):
            continue
        clauses.append(list(normalized))
    with open(proof_path, encoding="ascii") as fp:
        lines = [line.split() for line in fp if line.strip()]
    assert all(line[-1] == "0" for line in lines), "every proof clause must end in 0"
    proof = [[int(x) for x in line[:-1]] for line in lines]
    assert proof and proof[-1] == [], "proof must derive the empty clause"
    for learned in proof:
        assumptions = [-lit for lit in learned]
        assert unit_conflict(clauses, assumptions), f"non-RUP proof clause: {learned}"
        clauses.append(learned)


def run_exact(num_vars, clauses, extra=()):
    with tempfile.TemporaryDirectory() as tmp:
        cnf = os.path.join(tmp, "case.cnf")
        proof = os.path.join(tmp, "case.drat")
        solution = os.path.join(tmp, "case.sol")
        write_cnf(cnf, num_vars, clauses)
        result = subprocess.run(
            [BINARY, cnf, "--exact", "--epochs", "0", "--finisher-passes", "0",
             "--proof", proof, "--solution", solution, *extra],
            text=True, capture_output=True, timeout=10, check=False,
        )
        assert result.returncode in (0, 1), result.stderr
        data = json.loads(result.stdout)
        expected_sat = brute_force_sat(num_vars, clauses)
        if data["exact_attempted"]:
            assert data["exact_result"] == ("SAT" if expected_sat else "UNSAT")
        else:
            assert expected_sat and data["exact_result"] == "NOT_RUN"
        assert data["status"] == ("SATISFIED" if expected_sat else "UNSATISFIABLE")
        if expected_sat:
            assert data["solved"] is True
            assert os.path.exists(solution)
            values = [int(x) for x in open(solution, encoding="ascii").read().split() if x != "0"]
            assignment = {abs(lit): lit > 0 for lit in values}
            assert all(any(assignment[abs(lit)] == (lit > 0) for lit in clause) for clause in clauses)
            assert data["proof_generated"] is False
        else:
            assert data["solved"] is False
            assert not os.path.exists(solution)
            assert data["proof_generated"] is True
            verify_rup_proof(clauses, proof)


def test_limits_and_atomic_proofs():
    with tempfile.TemporaryDirectory() as tmp:
        cnf = os.path.join(tmp, "large-vars.cnf")
        write_cnf(cnf, 100000, [[1], [-1]])
        result = subprocess.run(
            [BINARY, cnf, "--exact", "--epochs", "0", "--finisher-passes", "0",
             "--exact-memory-mb", "1"], text=True, capture_output=True, check=False,
        )
        data = json.loads(result.stdout)
        assert data["status"] == "UNKNOWN" and data["exact_result"] == "LIMIT"

        # Pigeonhole(5,4) requires search, so a one-conflict ceiling must return UNKNOWN.
        clauses = []
        var = lambda pigeon, hole: pigeon * 4 + hole + 1
        for pigeon in range(5):
            clauses.append([var(pigeon, hole) for hole in range(4)])
        for hole in range(4):
            for left in range(5):
                for right in range(left + 1, 5):
                    clauses.append([-var(left, hole), -var(right, hole)])
        write_cnf(cnf, 20, clauses)
        result = subprocess.run(
            [BINARY, cnf, "--exact", "--epochs", "0", "--finisher-passes", "0",
             "--exact-max-conflicts", "1"], text=True, capture_output=True, check=False,
        )
        data = json.loads(result.stdout)
        assert data["status"] == "UNKNOWN" and data["exact_result"] == "LIMIT"

        proof = os.path.join(tmp, "existing.drat")
        with open(proof, "w", encoding="ascii") as fp:
            fp.write("preserve-me\n")
        write_cnf(cnf, 1, [[1]])
        subprocess.run([BINARY, cnf, "--proof", proof, "--epochs", "0"],
                       text=True, capture_output=True, check=False)
        assert open(proof, encoding="ascii").read() == "preserve-me\n"

        wcnf = os.path.join(tmp, "hard-conflict.wcnf")
        with open(wcnf, "w", encoding="ascii") as fp:
            fp.write("p wcnf 1 2 10\n10 1 0\n10 -1 0\n")
        result = subprocess.run([BINARY, wcnf, "--proof", proof],
                                text=True, capture_output=True, check=False)
        assert result.returncode == 2 and "supported only for CNF" in result.stderr


def test_indexed_finisher():
    with tempfile.TemporaryDirectory() as tmp:
        cnf = os.path.join(tmp, "units.cnf")
        clauses = [[var if var % 2 else -var] for var in range(1, 101)]
        write_cnf(cnf, 100, clauses)
        result = subprocess.run(
            [BINARY, cnf, "--epochs", "0", "--indexed-finisher",
             "--indexed-flips", "1000", "--indexed-checkpoint", "10"],
            text=True, capture_output=True, timeout=10, check=False,
        )
        data = json.loads(result.stdout)
        assert result.returncode == 0 and data["solved"] is True, data
        assert 0 < data["indexed_flips"] <= 1000, data
        assert data["indexed_index_mb"] > 0, data


def main():
    build_solver()
    run_exact(1, [[1], [-1]])
    run_exact(1, [[1, 1], [-1]])  # duplicated literal; clause length exceeds variable count
    run_exact(2, [[1, 2], [1, -2], [-1, 2], [-1, -2]])  # requires learned unit
    run_exact(3, [[1, -1, 2], [3]])  # tautology normalization
    run_exact(1, [[1], [-1]], ("--finisher-max-clauses", "1"))  # exact is not finisher-gated

    rng = random.Random(12345)
    for _ in range(40):
        num_vars = rng.randint(1, 7)
        clauses = []
        for _ in range(rng.randint(1, 30)):
            clauses.append([rng.choice((-1, 1)) * rng.randint(1, num_vars)
                            for _ in range(rng.randint(1, min(6, num_vars + 2)))])
        run_exact(num_vars, clauses)

    test_limits_and_atomic_proofs()
    test_indexed_finisher()

    print("all exact-mode tests passed")


if __name__ == "__main__":
    main()
