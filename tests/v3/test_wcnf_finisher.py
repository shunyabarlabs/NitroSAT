import os
import subprocess
import json
import tempfile
import sys

def create_wcnf(path):
    # 5 vertices, 3 colors. Vars = v * 3 + c (v from 0 to 4, c from 1 to 3)
    def var(v, c):
        return v * 3 + c

    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]

    hard_clauses = []
    # At least one color
    for v in range(5):
        hard_clauses.append([var(v, 1), var(v, 2), var(v, 3)])
        # At most one color
        hard_clauses.append([-var(v, 1), -var(v, 2)])
        hard_clauses.append([-var(v, 1), -var(v, 3)])
        hard_clauses.append([-var(v, 2), -var(v, 3)])

    # Adjacent different colors
    for u, v in edges:
        for c in range(1, 4):
            hard_clauses.append([-var(u, c), -var(v, c)])

    soft_clauses = [
        (10, [var(0, 1)]),
        (10, [var(1, 2)]),
        (10, [var(2, 1)]),
        (10, [var(3, 2)]),
        (10, [var(4, 1)])
    ]

    with open(path, "w") as f:
        top = 1000000
        f.write(f"p wcnf 15 {len(hard_clauses) + len(soft_clauses)} {top}\n")
        for c in hard_clauses:
            f.write(f"{top} " + " ".join(map(str, c)) + " 0\n")
        for w, c in soft_clauses:
            f.write(f"{w} " + " ".join(map(str, c)) + " 0\n")

def run_solver(bin_path, wcnf_path, extra_args):
    cmd = [bin_path, wcnf_path] + extra_args
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode not in (0, 1):
        print(f"Solver failed with code {result.returncode}")
        print("stderr:", result.stderr)
        sys.exit(1)

    try:
        data = json.loads(result.stdout)
        return data
    except json.JSONDecodeError:
        print("Failed to parse JSON output:")
        print(result.stdout)
        sys.exit(1)

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <path_to_nitrosatv3>")
        sys.exit(1)

    bin_path = sys.argv[1]

    with tempfile.TemporaryDirectory() as tmpdir:
        wcnf_path = os.path.join(tmpdir, "graph_coloring.wcnf")
        create_wcnf(wcnf_path)

        # Test 1: Only epochs=0, finisher-passes=0, indexed-finisher
        # This tests that the hard-only finisher can reach hard_unsatisfied == 0
        # but leaves soft clauses unoptimized.
        print("\n--- Test 1: Hard-only finisher (no soft optimization) ---")
        res1 = run_solver(bin_path, wcnf_path, ["--epochs", "0", "--finisher-passes", "0", "--indexed-finisher"])

        assert res1["feasible"] is True, "Expected feasible=true for WCNF with indexed finisher"
        assert res1["hard_unsatisfied"] == 0, "Expected hard_unsatisfied to be 0"

        print("Status:", res1["status"])
        print("Feasible:", res1["feasible"])
        print("Soft cost:", res1["soft_cost"])
        print("Indexed Index MB:", res1["indexed_index_mb"])

        assert res1["indexed_index_mb"] > 0.0, "Expected indexed_index_mb to be reported correctly"

        # Test 2: epochs=0, finisher-passes=256, indexed-finisher
        # This tests that after hard feasibility, the soft clauses are optimized sequentially.
        print("\n--- Test 2: Hard-only finisher + sequential soft optimization ---")
        res2 = run_solver(bin_path, wcnf_path, ["--epochs", "0", "--finisher-passes", "256", "--indexed-finisher"])

        assert res2["feasible"] is True, "Expected feasible=true"
        assert res2["hard_unsatisfied"] == 0, "Expected hard_unsatisfied to be 0"

        print("Status:", res2["status"])
        print("Feasible:", res2["feasible"])
        print("Soft cost:", res2["soft_cost"])

        # Typically soft cost should be equal or better after soft optimization
        assert res2["soft_cost"] <= res1["soft_cost"], f"Expected soft cost {res2['soft_cost']} <= {res1['soft_cost']}"

        print("\nAll tests passed successfully.")

if __name__ == "__main__":
    main()
