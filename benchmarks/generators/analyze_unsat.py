import sys
from collections import Counter

if len(sys.argv) != 4:
    print("Usage: python3 analyze_unsat.py <cnf_file> <map_file> <sol_file>")
    sys.exit(1)

cnf_file = sys.argv[1]
map_file = sys.argv[2]
sol_file = sys.argv[3]

# 1. Read assignment
assignment = {}
try:
    with open(sol_file, "r") as f:
        tokens = f.read().strip().split()
        for tok in tokens:
            val = int(tok)
            if val == 0: break
            assignment[abs(val)] = (val > 0)
except Exception as e:
    print(f"Error reading solution: {e}")
    sys.exit(1)

# 2. Read map file
clause_meta = {}
try:
    with open(map_file, "r") as f:
        for line in f:
            parts = line.strip().split(" ", 1)
            if len(parts) == 2:
                clause_meta[int(parts[0])] = parts[1]
except Exception as e:
    print(f"Error reading map: {e}")
    sys.exit(1)

# 3. Read CNF and evaluate clauses
unsat_indices = []
clause_idx = 1
try:
    with open(cnf_file, "r") as f:
        for line in f:
            if line.startswith('p') or line.startswith('c'):
                continue
            lits = [int(x) for x in line.strip().split() if int(x) != 0]
            if not lits: continue
            
            # Evaluate
            is_sat = False
            for lit in lits:
                var = abs(lit)
                val = assignment.get(var, False)
                if (lit > 0 and val) or (lit < 0 and not val):
                    is_sat = True
                    break
            
            if not is_sat:
                unsat_indices.append(clause_idx)
            clause_idx += 1
except Exception as e:
    print(f"Error reading CNF: {e}")
    sys.exit(1)

print(f"Evaluated {clause_idx - 1} clauses. Found {len(unsat_indices)} unsatisfied.")

if not unsat_indices:
    print("All clauses satisfied!")
    sys.exit(0)

# 4. Analyze
meta_counts = Counter()
category_counts = Counter()

for idx in unsat_indices:
    meta = clause_meta.get(idx, "UNKNOWN")
    meta_counts[meta] += 1
    
    # Extract high-level category
    if "OUTPUT_BIT" in meta:
        category_counts["Output Constraints"] += 1
    elif "PP_" in meta:
        category_counts["Partial Products (AND gates)"] += 1
    elif "ACCUM_" in meta:
        # Accumulator row
        row = meta.split("_row")[1].split("_")[0]
        category_counts[f"Adder Grid Row {row}"] += 1
    elif "CARRY_PROP_" in meta:
        category_counts["Carry Propagation"] += 1
    elif "INPUT" in meta:
        category_counts["Input Constraints"] += 1
    else:
        category_counts["Other"] += 1

print("\n--- UNSATISFIED CLAUSES STRUCTURAL BREAKDOWN ---")
for cat, count in category_counts.most_common():
    print(f"{count:5d} : {cat}")

print("\n--- DETAILED FAULT LOCATIONS (Top 15) ---")
for meta, count in meta_counts.most_common(15):
    print(f"{count:5d} : {meta}")
