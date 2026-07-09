import sys

class CNFGenerator:
    def __init__(self):
        self.clauses = []
        self.clause_meta = []
        self.num_vars = 0

    def new_var(self):
        self.num_vars += 1
        return self.num_vars

    def add_clause(self, lits, meta=""):
        self.clauses.append(lits)
        self.clause_meta.append(meta)

    def add_and(self, a, b, meta="AND"):
        c = self.new_var()
        self.add_clause([-a, -b, c], meta)
        self.add_clause([a, -c], meta)
        self.add_clause([b, -c], meta)
        return c

    def add_xor(self, a, b, meta="XOR"):
        c = self.new_var()
        self.add_clause([-a, -b, -c], meta)
        self.add_clause([a, b, -c], meta)
        self.add_clause([a, -b, c], meta)
        self.add_clause([-a, b, c], meta)
        return c

    def add_or(self, a, b, meta="OR"):
        c = self.new_var()
        self.add_clause([a, b, -c], meta)
        self.add_clause([-a, c], meta)
        self.add_clause([-b, c], meta)
        return c

    def half_adder(self, x, y, meta="HA"):
        s = self.add_xor(x, y, meta + "_SUM")
        c_out = self.add_and(x, y, meta + "_CARRY")
        return s, c_out

    def full_adder(self, x, y, c_in, meta="FA"):
        s1, c1 = self.half_adder(x, y, meta + "_HA1")
        s2, c2 = self.half_adder(s1, c_in, meta + "_HA2")
        c_out = self.add_or(c1, c2, meta + "_OR")
        return s2, c_out

def generate_multiplier(bits_a, bits_b, target_product, filename):
    cnf = CNFGenerator()
    
    A = [cnf.new_var() for _ in range(bits_a)]
    B = [cnf.new_var() for _ in range(bits_b)]
    
    # Partial products
    PP = [[0] * bits_a for _ in range(bits_b)]
    for j in range(bits_b):
        for i in range(bits_a):
            PP[j][i] = cnf.add_and(A[i], B[j], meta=f"PP_row{j}_col{i}")
            
    # Array multiplier accumulation
    accum = list(PP[0])
    
    for j in range(1, bits_b):
        c_in = None
        for i in range(bits_a):
            idx = j + i
            meta_prefix = f"ACCUM_row{j}_col{idx}"
            if idx >= len(accum):
                if c_in is None:
                    accum.append(PP[j][i])
                else:
                    s, c = cnf.half_adder(PP[j][i], c_in, meta=meta_prefix)
                    accum.append(s)
                    c_in = c
            else:
                if c_in is None:
                    s, c = cnf.half_adder(accum[idx], PP[j][i], meta=meta_prefix)
                    accum[idx] = s
                    c_in = c
                else:
                    s, c = cnf.full_adder(accum[idx], PP[j][i], c_in, meta=meta_prefix)
                    accum[idx] = s
                    c_in = c
        
        idx = j + bits_a
        while c_in is not None:
            meta_prefix = f"CARRY_PROP_row{j}_col{idx}"
            if idx >= len(accum):
                accum.append(c_in)
                c_in = None
            else:
                s, c = cnf.half_adder(accum[idx], c_in, meta=meta_prefix)
                accum[idx] = s
                c_in = c
            idx += 1
            
    O = accum
            
    # Constrain output to target_product
    product_bits = bin(target_product)[2:][::-1] # least significant bit first
    # Pad with zeros if necessary
    product_bits = product_bits + '0' * (len(O) - len(product_bits))
    
    if len(product_bits) > len(O):
        print("Error: Target product is too large for the multiplier.")
        sys.exit(1)
        
    for i in range(len(O)):
        bit_val = int(product_bits[i])
        if bit_val == 1:
            cnf.add_clause([O[i]], meta=f"OUTPUT_BIT_{i}_TARGET_1")
        else:
            cnf.add_clause([-O[i]], meta=f"OUTPUT_BIT_{i}_TARGET_0")
            
    # Constrain inputs to not be 1
    cnf.add_clause(A[1:], meta="INPUT_A_GT_1")
    cnf.add_clause(B[1:], meta="INPUT_B_GT_1")
    
    cnf.add_clause([A[0]], meta="INPUT_A_LSB_1")
    cnf.add_clause([B[0]], meta="INPUT_B_LSB_1")

    print(f"Generated {bits_a}x{bits_b} multiplier for target {target_product}")
    print(f"Variables: {cnf.num_vars}, Clauses: {len(cnf.clauses)}")
    
    with open(filename, 'w') as f:
        f.write(f"p cnf {cnf.num_vars} {len(cnf.clauses)}\n")
        for clause in cnf.clauses:
            f.write(" ".join(map(str, clause)) + " 0\n")

    map_filename = filename + ".map"
    with open(map_filename, 'w') as f:
        for i, meta in enumerate(cnf.clause_meta):
            f.write(f"{i+1} {meta}\n")
    print(f"Saved clause map to {map_filename}")

if __name__ == "__main__":
    if len(sys.argv) == 4:
        bits_a = int(sys.argv[1])
        bits_b = int(sys.argv[2])
        target = int(sys.argv[3])
        filename = f"/tmp/factor_{bits_a}x{bits_b}.cnf"
    else:
        # Default: factor 4292870399 = 65519 * 65521 (16x16)
        bits_a = 16
        bits_b = 16
        target = 4292870399
        filename = "/tmp/factor_16x16.cnf"
        
    generate_multiplier(bits_a, bits_b, target, filename)
