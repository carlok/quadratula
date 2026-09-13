# This considers conjunctions of the 47 inequivalent quasigroup laws
# and finds which conjunctions imply which laws and especially when
# the conjunctions are equivalent.
#
# It would be computationally faster to keep track of countermodels,
# but this would require implementing equation evaluation, hence
# parsing of equations expressed in terms of the three quasigroup
# operations.  This is not a priority.

from ATP_utils_quasigroup import models, prove
from Schroeder_equations import eq_list, op_declare, quasigroup_eqs
from Schroeder_implications import representatives, implications, class_list
from collections import defaultdict

claimed_multi_rep = [(), (1,), (2,), (3,), (4,), (5,), (6,), (7,), (8,), (9,), (10,), (12,), (14,), (15,), (18,), (20,), (23,), (24,), (25,), (29,), (34,), (38,), (54,), (56,), (60,), (65,), (79,), (88,), (93,), (100,), (101,), (102,), (103,), (105,), (123,), (127,), (140,), (156,), (161,), (183,), (202,), (224,), (264,), (336,), (408,), (432,), (504,), (654,), (2, 18), (2, 54), (2, 56), (4, 29), (4, 56), (4, 100), (4, 101), (4, 102), (4, 183), (4, 224), (4, 408), (25, 54), (29, 56), (29, 93), (34, 54), (34, 56), (34, 93), (34, 102), (34, 140), (34, 432), (54, 100), (54, 102), (54, 127), (54, 202), (54, 224), (54, 264), (54, 408), (54, 504), (54, 654), (56, 102), (56, 202), (56, 654), (102, 202), (102, 224), (102, 336), (103, 127), (103, 202), (103, 654), (127, 140), (127, 432), (161, 224), (161, 408), (202, 224), (202, 408), (202, 432), (224, 432), (224, 654), (408, 432), (408, 654), (432, 654), (4, 29, 56), (4, 102, 224), (34, 54, 102), (54, 102, 202), (54, 102, 224), (54, 202, 224), (54, 202, 408), (54, 224, 654), (54, 408, 654), (56, 102, 202), (102, 202, 224), (202, 224, 432), (202, 408, 432), (224, 432, 654), (408, 432, 654), (54, 102, 202, 224)]
# len(claimed_multi_rep) == 113

consequences = defaultdict(set, {(i,): impl for i, impl in implications.items()}) # maps (ordered tuples of eq id) to the (set of eq id) they imply
consequences[()] = set()
unknown_consequences = defaultdict(set) # same structure as consequences

def multi_implication(eqs_id, eqB_id, size1, size2=32):
    if not eqs_id:
        raise ValueError(eqs_id)
    if eqB_id in consequences[eqs_id]:
        return True
    extra1 = quasigroup_eqs + "\n" + ". ".join(eq_list[i] for i in eqs_id) + "."
    extra2 = eq_list[eqB_id] + "."
    eqs_id = tuple(eqs_id)
    implies = None
    if models(extra0=op_declare+f"assign(iterate_up_to, {size1}). assign(selection_measure, 3). assign(max_models, 1).", extra1=extra1, extra2=extra2):
        implies = False
    elif prove(extra0=op_declare, extra1=extra1, extra2=extra2):
        implies = True
    else:
        if models(extra0=op_declare+f"assign(iterate_up_to, {size2}). assign(selection_measure, 1). assign(max_models, 1).", extra1=extra1, extra2=extra2, max_seconds=3):
            implies = False
    if implies is True:
        consequences[eqs_id].add(eqB_id)
        consequences[eqs_id].update(consequences[(eqB_id,)])
    elif implies is None:
        print(f"Unknown implication {eqs_id} ({extra1}) to {eqB_id} ({extra2})")
        unknown_consequences[eqs_id].add(eqB_id)
    return implies

multi_rep = set((i,) for i in representatives) # set of non-diminishable (ordered tuples of eq id)
multi_to_rep = {} # maps (ordered tuples of eq id) (obtained by appending a eq_list to a multi_rep) to their representative
pending = [(i,) for i in representatives]
while pending:
    eqsA = pending.pop(0)
    for moreA in representatives:
        if moreA <= eqsA[-1]:
            continue
        fullA = eqsA + (moreA,)
        print(f"current {fullA}; pending {len(pending)}; found {len(multi_rep)} (latest {pending[-1] if pending else ""})", end="  \r")
        for a, i in enumerate(fullA):
            partA = fullA[:a] + fullA[a + 1:]
            consequences_of_partA = consequences.get(partA, set())
            consequences[fullA].update(consequences_of_partA)
            if i in consequences_of_partA:
                multi_to_rep[fullA] = partA
        if fullA not in multi_to_rep:
            extra1A = quasigroup_eqs + "\n" + ". ".join(eq_list[i] for i in fullA) + "."
            size1 = len(models(extra0=op_declare+"assign(max_models, 15). assign(iterate_up_to, 6).", extra1=extra1A)[-1]) # adapt search size to equation hardness
            [multi_implication(fullA, j, size1) for j in [3, 7, 24, 25]] # updates (unknown_)consequences[fullA]
            [multi_implication(fullA, j, size1) for j in representatives] # updates (unknown_)consequences[fullA]
            for eqsB in multi_rep:
                if consequences[fullA].issuperset(eqsB) and consequences[eqsB].issuperset(fullA):
                    multi_to_rep[fullA] = eqsB
                    break
                if (all(j in consequences[fullA] or j in unknown_consequences[fullA] for j in eqsB)
                    and all(i in consequences[eqsB] or i in unknown_consequences[eqsB] for i in fullA)):
                    print(f"Unknown equivalence {fullA} vs {eqsB}")
        if fullA in multi_to_rep:
            # print(fullA, "equivalent to", multi_to_rep[fullA])
            # We found an equivalent, delete useless data
            consequences.pop(fullA, None)
            unknown_consequences.pop(fullA, None)
        else:
            # The new equivalence class may lead to further conjunctions
            multi_to_rep[fullA] = fullA
            multi_rep.add(fullA)
            pending.append(fullA)
            #if (len(multi_rep) % 100) == 0:
            #    print(len(multi_rep), sorted(multi_rep, key=(lambda x:(len(x), x))))

multi_rep.add(())

print(len(multi_rep), sorted(multi_rep, key=(lambda x:(len(x), x))))

if sorted(consequences.keys()) != sorted(multi_rep):
    print("Surprisingly")
    print(sorted(multi_rep))
    print("differs from")
    print(sorted(consequences.keys()))

with open("quasigroup_outcomes.py", "w") as f:
    print(f"multi_rep = {multi_rep}", file=f)
    print(f"multi_to_rep = {multi_to_rep}", file=f)
    print(f"consequences = {dict(consequences)}", file=f)
    print(f"unknown_consequences = {dict(unknown_consequences)}", file=f)

class_dict = {i:t for i, t in class_list}
all_consequences = {k: tuple(sorted(j for i in v for j in class_dict[i])) for k, v in consequences.items()}

long_classes = sorted(all_consequences.values(), key=(lambda x:(len(x), x)))
eq_to_long_class = {k[0]: v for k, v in all_consequences.items() if len(k) == 1}

with open("quasigroup_long_classes.py", "w") as f:
    print(f"long_classes = {long_classes}", file=f)
    print(f"eq_to_long_class = {eq_to_long_class}", file=f)

with open("quasigroup_long_classes.json", "w") as f:
    f.write('[\n    ' + ',\n    '.join(str(list(c)) for c in long_classes) + '\n]\n')

with open("quasigroup_eq_to_long_class.json", "w") as f:
    f.write('{\n    ' + ',\n    '.join(f'"{k}": {list(v)}' for k, v in eq_to_long_class.items()) + '\n}\n')
