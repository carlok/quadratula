# This generates the list of 990 equations considered in 1890 by Schröder
# https://archive.org/details/bub_gb_P95LAAAAYAAJ
# Namely all equations of the form x?y?z=a?b?c with ? replaced by
# one of the three quasigroup operations, and a,b,c replaced by some
# permutation of x,y,z, with some choice of how the LHS and RHS are
# associated.

op_symbols = ["*", "//", "\\\\"]
op_declare = 'op(400, infix, "*"). op(400, infix, "//"). op(400, infix, "\\\\").\n'
quasigroup_eqs = "(y // x) * y = x.  y // (x * y) = x.  x * (y \\\\ x) = y.  (x * y) \\\\ x = y.  x \\\\ (y // x) = y.  (x \\\\ y) // x = y."

op_pairs = [(op1, op2) for op1 in op_symbols for op2 in op_symbols]
var_perms = ["xyz", "xzy", "yxz", "yzx", "zxy", "zyx"]
var_perms_restricted = ["xzy", "yxz", "yzx", "zyx"]
eq_list = [""]
for i, pair_i in enumerate(op_pairs):
    for j, pair_j in enumerate(op_pairs):
        if i > j:
            continue
        for perm in (var_perms_restricted if i == j else var_perms):
            eq_list.append(f"x {pair_i[0]} (y {pair_i[1]} z) = {perm[0]} {pair_j[0]} ({perm[1]} {pair_j[1]} {perm[2]})")
for pair_i in op_pairs:
    for pair_j in op_pairs:
        for perm in var_perms:
            eq_list.append(f"x {pair_i[0]} (y {pair_i[1]} z) = ({perm[0]} {pair_j[0]} {perm[1]}) {pair_j[1]} {perm[2]}")
for i, pair_i in enumerate(op_pairs):
    for j, pair_j in enumerate(op_pairs):
        if i > j:
            continue
        for perm in (var_perms_restricted if i == j else var_perms):
            eq_list.append(f"(x {pair_i[0]} y) {pair_i[1]} z = ({perm[0]} {pair_j[0]} {perm[1]}) {pair_j[1]} {perm[2]}")
assert len(eq_list) - 1 == 990

if __name__ == "__main__":
    with open("quasigroup-equations.txt", "w") as f:
        for i, eq in enumerate(eq_list):
            if i > 0:
                print(f"Sch-{i}: {eq}", file=f)
