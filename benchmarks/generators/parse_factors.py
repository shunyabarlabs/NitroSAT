import sys

if len(sys.argv) != 4:
    print("Usage: python3 parse_factors.py <bits_a> <bits_b> <solution_file>")
    sys.exit(1)

bits_a = int(sys.argv[1])
bits_b = int(sys.argv[2])
sol_file = sys.argv[3]

try:
    with open(sol_file, "r") as f:
        tokens = f.read().strip().split()
except FileNotFoundError:
    print(f"Could not find solution file {sol_file}")
    sys.exit(1)

sol = set()
for tok in tokens:
    val = int(tok)
    if val == 0:
        break
    if val > 0:
        sol.add(val)

A_val = 0
for i in range(bits_a):
    if (i + 1) in sol:
        A_val |= (1 << i)

B_val = 0
for i in range(bits_b):
    if (bits_a + i + 1) in sol:
        B_val |= (1 << i)

print(f"Factor A = {A_val}")
print(f"Factor B = {B_val}")
print(f"Product = {A_val * B_val}")
