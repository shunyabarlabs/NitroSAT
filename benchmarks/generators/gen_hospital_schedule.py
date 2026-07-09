"""
Hard benchmark: Hospital Shift Scheduling as WCNF

This is a real-world NP-hard optimization problem commonly found in MaxSAT
evaluations. It generates dense, structured constraint interactions with
both hard feasibility requirements and soft preference optimization.

Problem structure:
- N nurses, D days, S shifts per day
- Hard constraints:
  * Each shift covered by exactly K nurses (via at-least-K + at-most-K)
  * No nurse works more than 1 shift per day
  * No nurse works more than max_consecutive consecutive days
  * Incompatible nurse pairs cannot share a shift
- Soft constraints (weighted):
  * Nurse shift preferences (higher weight = stronger preference)
  * Balanced workload: penalize deviation from target shifts
  * Weekend off preferences
"""

import random
import itertools


def generate_hospital_schedule(nurses, days, shifts, cover_k,
                               max_consecutive, num_conflicts,
                               filename, seed=42):
    random.seed(seed)
    top = 10000000  # hard clause weight

    # Variable: nurse n works shift s on day d
    # var(n, d, s) = n * days * shifts + d * shifts + s + 1
    def var(n, d, s):
        return n * days * shifts + d * shifts + s + 1

    num_vars = nurses * days * shifts
    hard_clauses = []
    soft_clauses = []

    # =========================================================
    # HARD CONSTRAINT 1: Each shift on each day covered by >= K nurses
    # At-least-K encoding using sequential counter (Sinz encoding)
    # For simplicity, use pairwise for at-most-K and direct for at-least-K
    # =========================================================
    for d in range(days):
        for s in range(shifts):
            shift_vars = [var(n, d, s) for n in range(nurses)]

            # At-least-K: enumerate all subsets of size (nurses - K + 1)
            # that must contain at least one selected nurse
            # More efficient: for each combination of (nurses - K + 1) nurses
            # NOT assigned, at least one must actually be assigned
            # Equivalent: every (nurses - K + 1)-subset must have >= 1 true
            # This is: for every set of (nurses - K) unselected, 
            # the remaining K must have >= 1 true
            # Simpler: at-least-K = negation of at-most-(K-1)
            # Use the clause approach: for every (N-K+1)-sized subset of 
            # negations, assert the clause
            if nurses <= 15:
                # Direct encoding for small N
                for combo in itertools.combinations(range(nurses), nurses - cover_k + 1):
                    clause = [var(n, d, s) for n in combo]
                    hard_clauses.append(clause)

            # At-most-K: no (K+1) nurses all assigned
            if nurses <= 15:
                for combo in itertools.combinations(range(nurses), cover_k + 1):
                    clause = [-var(n, d, s) for n in combo]
                    hard_clauses.append(clause)

    # =========================================================
    # HARD CONSTRAINT 2: No nurse works > 1 shift per day
    # =========================================================
    for n in range(nurses):
        for d in range(days):
            for s1 in range(shifts):
                for s2 in range(s1 + 1, shifts):
                    hard_clauses.append([-var(n, d, s1), -var(n, d, s2)])

    # =========================================================
    # HARD CONSTRAINT 3: No nurse works > max_consecutive consecutive days
    # If nurse works day d, d+1, ..., d+max_consecutive, that's a violation
    # For each window of (max_consecutive + 1) days, at least one must be off
    # =========================================================
    for n in range(nurses):
        for d_start in range(days - max_consecutive):
            window_days = range(d_start, d_start + max_consecutive + 1)
            # "at least one day off" = at least one day where no shift is worked
            # Encode: for each day in window, create aux var "works_day"
            # Then: not all works_day can be true
            # Simpler: just say not all shifts in the window can be true
            # For each assignment of one shift per day in the window:
            # Use the direct approach: for each combo of shifts across the window
            for shift_combo in itertools.product(range(shifts), repeat=max_consecutive + 1):
                clause = [-var(n, day, shift_combo[i])
                          for i, day in enumerate(window_days)]
                hard_clauses.append(clause)

    # =========================================================
    # HARD CONSTRAINT 4: Incompatible nurse pairs
    # =========================================================
    conflict_pairs = set()
    while len(conflict_pairs) < num_conflicts:
        n1 = random.randint(0, nurses - 1)
        n2 = random.randint(0, nurses - 1)
        if n1 != n2:
            conflict_pairs.add((min(n1, n2), max(n1, n2)))

    for n1, n2 in conflict_pairs:
        for d in range(days):
            for s in range(shifts):
                hard_clauses.append([-var(n1, d, s), -var(n2, d, s)])

    # =========================================================
    # SOFT CONSTRAINT 1: Nurse shift preferences
    # =========================================================
    for n in range(nurses):
        # Each nurse has 2-3 preferred shifts
        num_prefs = random.randint(2, 3)
        for _ in range(num_prefs):
            d = random.randint(0, days - 1)
            s = random.randint(0, shifts - 1)
            weight = random.randint(5, 20)
            soft_clauses.append((weight, [var(n, d, s)]))

    # =========================================================
    # SOFT CONSTRAINT 2: Weekend off preference
    # Weekends are days 5,6 in each week
    # =========================================================
    for n in range(nurses):
        for week in range(days // 7):
            for weekend_day in [5, 6]:
                d = week * 7 + weekend_day
                if d < days:
                    for s in range(shifts):
                        # Prefer NOT working weekends
                        soft_clauses.append((8, [-var(n, d, s)]))

    # =========================================================
    # SOFT CONSTRAINT 3: Workload balance
    # Target: each nurse works roughly (days * cover_k * shifts) / nurses shifts
    # Penalize if a nurse works too many shifts total
    # Approximate: penalize each shift assignment slightly
    # =========================================================
    target_shifts = (days * cover_k * shifts) // nurses
    for n in range(nurses):
        # Light penalty for each extra shift beyond target
        assigned_count = 0
        for d in range(days):
            for s in range(shifts):
                # Small cost for every shift worked (encourages balance)
                soft_clauses.append((1, [-var(n, d, s)]))

    num_clauses = len(hard_clauses) + len(soft_clauses)
    total_hard = len(hard_clauses)
    total_soft = len(soft_clauses)

    print(f"Hospital Shift Scheduling Benchmark")
    print(f"  Nurses: {nurses}, Days: {days}, Shifts/day: {shifts}")
    print(f"  Coverage: {cover_k} nurses per shift")
    print(f"  Max consecutive days: {max_consecutive}")
    print(f"  Incompatible pairs: {num_conflicts}")
    print(f"  Variables: {num_vars}")
    print(f"  Hard clauses: {total_hard}")
    print(f"  Soft clauses: {total_soft}")
    print(f"  Total clauses: {num_clauses}")
    print(f"  Constraint density (clauses/vars): {num_clauses / num_vars:.1f}")

    with open(filename, 'w') as f:
        f.write(f"p wcnf {num_vars} {num_clauses} {top}\n")
        for c in hard_clauses:
            f.write(f"{top} " + " ".join(map(str, c)) + " 0\n")
        for w, c in soft_clauses:
            f.write(f"{w} " + " ".join(map(str, c)) + " 0\n")


if __name__ == "__main__":
    # Configuration: 12 nurses, 14 days (2 weeks), 3 shifts/day
    # Coverage: 2 nurses per shift, max 5 consecutive days, 8 conflict pairs
    generate_hospital_schedule(
        nurses=12, days=14, shifts=3, cover_k=2,
        max_consecutive=5, num_conflicts=8,
        filename="/tmp/hospital_schedule.wcnf"
    )
    print(f"\nGenerated /tmp/hospital_schedule.wcnf")
