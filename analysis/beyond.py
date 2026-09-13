"""Beyond order 6, exactly, without enumerating all quasigroups of orders 7..9.

Every separation still open after the order-6 floor involves a law or variety
V whose models are enumerated exhaustively by `qg models` at orders 7..9
(data/beyond/models_o{n}_{generators}.tsv.gz). Whether a class pair (a, b) is
witnessed at order <= k depends only on models of `a`; whether a variety V is
realised at order <= k depends only on models of V. So the separation curves
extend exactly to the largest order at which every needed run completed.

Writes data/derived/beyond.json.
"""

from __future__ import annotations

import gzip
import json
import sys
from collections import Counter
from pathlib import Path

import highspy
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gates  # noqa: E402
import lefloch  # noqa: E402
from compute import solve_cover, variety_name  # noqa: E402

ROOT = gates.ROOT
BEYOND = ROOT / "data" / "beyond"
N = lefloch.N_LAWS
FLOOR_MAX = 6


def read_models(path: Path) -> list[tuple[str, int, frozenset[int]]]:
    out = []
    with gzip.open(path, "rt") as f:
        next(f)
        for line in f:
            _o, table, aut, sat = line.rstrip("\n").split("\t")
            out.append((table, int(aut), frozenset() if sat == "-" else frozenset(map(int, sat.split(",")))))
    return out


def runs() -> dict[tuple[int, str], dict]:
    """(order, generators) -> run log fields; generators as '103,127'."""
    out = {}
    for line in (BEYOND / "runs.tsv").read_text().splitlines():
        f = line.split("\t")
        d = dict(zip(f[0::2], f[1::2]))
        out[(int(d["order"]), d["laws"])] = d
    return out


def main() -> None:
    gt = lefloch.load()
    res = gates.results()
    reps = gt.reps
    long_set = set(gt.long_classes)
    log = runs()

    def closure(laws: frozenset[int]) -> frozenset[int]:
        out = frozenset(range(1, N + 1))
        for c in gt.long_classes:
            if laws <= c:
                out &= c
        return out

    # floor signatures with first order
    floor = [(frozenset(s["satisfied"]), s["first_order"], s["least_table"]) for s in res["signatures"]]
    missing6 = res["per_order"][FLOOR_MAX - 1]["missing_varieties"]
    gens = {name: name.strip("()").replace(" ", "") for name in missing6}

    # models per (generator set, order); None if the run did not complete
    models: dict[tuple[str, int], list | None] = {}
    orders_run = sorted({n for n, _ in log})
    for name, g in gens.items():
        for n in orders_run:
            d = log.get((n, g))
            p = BEYOND / f"models_o{n}_{g.replace(',', '-')}.tsv.gz"
            models[(g, n)] = read_models(p) if d and "iso_models" in d and p.exists() else None

    def variety_of(g: str) -> frozenset[int]:
        return closure(frozenset(x for r in map(int, g.split(",")) for x in gt.classes[r]))

    # Order k is exact once the runs at order k finished for every variety still
    # unrealised at order k-1 (realisation is monotone, so settled varieties need
    # no further runs). Computed incrementally below, after `realised_at` exists.
    complete_upto = FLOOR_MAX

    # sanity: every model satisfies its generators, and every signature is a variety
    bad = [(g, n) for (g, n), ms in models.items() if ms for _, _, s in ms if not (variety_of(g) <= s and s in long_set)]
    assert not bad, bad[:5]

    # signatures available at order <= k (floor plus models of the enumerated generators)
    def sigs_upto(k: int) -> list[tuple[frozenset[int], int, str]]:
        out = [x for x in floor if x[1] <= k]
        for (g, n), ms in models.items():
            if ms and n <= k:
                out += [(s, n, t) for t, _, s in ms]
        return out

    def realised_names(k: int) -> set[str]:
        sigset = {s for s, _, _ in sigs_upto(k)}
        done = set()
        for name, g in gens.items():
            V = variety_of(g)
            inter = frozenset(range(1, N + 1))
            for s in sigset:
                if V <= s:
                    inter &= s
            if inter == V:
                done.add(name)
        return done

    # `beyond.py --open-before k`: print the generator sets still unrealised at order k-1
    if "--open-before" in sys.argv:
        k = int(sys.argv[sys.argv.index("--open-before") + 1])
        settled = realised_names(k - 1)
        print(" ".join(gens[name] for name in gens if name not in settled))
        return

    for n in orders_run:
        if complete_upto != n - 1:
            break
        still_open = [name for name in gens if name not in realised_names(n - 1)]
        if all(models[(gens[name], n)] is not None for name in still_open):
            complete_upto = n

    ps = res["pair_space"]
    rep_pos = {r: i for i, r in enumerate(reps)}
    curve = []
    least_variety_order: dict[str, int | None] = {}
    for k in range(FLOOR_MAX, complete_upto + 1):
        S = sigs_upto(k)
        sigset = {s for s, _, _ in S}
        # class pairs witnessed
        covered = {(a, b) for s in sigset for a in reps if a in s for b in reps if b not in s}
        nonimp = {(a, b) for a in reps for b in reps if a != b and b not in gt.rep_implies[a]}
        assert covered <= nonimp
        law_cov = sum(len(gt.classes[a]) * len(gt.classes[b]) for a, b in covered)
        # varieties realised: V equals the intersection of the signatures containing it
        realised = set()
        for V in long_set:
            inter = frozenset(range(1, N + 1))
            for s in sigset:
                if V <= s:
                    inter &= s
            if inter == V:
                realised.add(V)
        for name in gens:
            V = variety_of(gens[name])
            if V in realised and name not in least_variety_order:
                least_variety_order[name] = k
        cols = {r: tuple(sorted(i for i, s in enumerate(sorted(sigset, key=sorted)) if r in s)) for r in reps}
        curve.append({
            "order": k,
            "class_non_implications_witnessed": len(covered),
            "class_saturation": len(covered) / ps["class_non_implications"],
            "law_non_implications_witnessed": law_cov,
            "law_saturation": law_cov / ps["law_non_implications"],
            "class_blocks": len(set(cols.values())),
            "classes_separated": sum(1 for c in Counter(cols.values()).values() if c == 1),
            "varieties_realised": len(realised),
            "distinct_known_signatures": len(sigset),
            "missing_varieties": sorted(variety_name(gt, V) for V in long_set - realised),
        })
    for name in gens:
        least_variety_order.setdefault(name, None)

    # least witness order for the residual class pairs, exhaustively
    residual = [tuple(p) for p in res["unrefuted_class_non_implications_at_order_6"]]
    mace = {(r["a"], r["b"]): r["order"] for r in json.loads((gates.DERIVED / "residual.json").read_text()).get("rows", [])}
    least_pair = []
    for a, b in residual:
        g = str(a)
        # A pair unwitnessed at order 6 has an unrealised hypothesis closure, generated by (a) alone.
        assert g in gens.values(), f"no generator run for residual hypothesis Sch-{a}"
        found = None
        for n in range(FLOOR_MAX + 1, complete_upto + 1):
            ms = models.get((g, n))
            if ms is None:
                break
            if any(b not in s for _, _, s in ms):
                found = n
                break
        least_pair.append({"a": a, "b": b, "exhaustive_least_order": found, "mace4_order": mace.get((a, b))})

    # model counts by order for every class and variety (orders 1..6 from the floor, 7.. where known exactly)
    def counts_for(V: frozenset[int]) -> dict[str, int | None]:
        out: dict[str, int | None] = {}
        for n in range(1, FLOOR_MAX + 1):
            out[str(n)] = sum(s_["iso_classes_by_order"][str(n)] for s_ in res["signatures"] if V <= frozenset(s_["satisfied"]))
        for n in range(FLOOR_MAX + 1, complete_upto + 1):
            # exact if some enumerated generator set defines a variety contained in V
            src = [g for g in gens.values() if variety_of(g) <= V and models.get((g, n)) is not None]
            out[str(n)] = sum(1 for _, _, s in models[(src[0], n)] if V <= s) if src else None
        return out

    class_counts = {f"Sch-{r}": counts_for(gt.closure_of_rep[r]) for r in reps}
    variety_counts = {variety_name(gt, V): counts_for(V) for V in gt.long_classes}

    # cover: all class non-implications, using known signatures up to complete_upto;
    # lower bound: the same ILP with every one of the 114 varieties as a hypothetical signature
    universe_pairs = sorted({(a, b) for a in reps for b in reps if a != b and b not in gt.rep_implies[a]})
    uidx = {p: i for i, p in enumerate(universe_pairs)}
    S = sigs_upto(complete_upto)
    best: dict[frozenset[int], tuple[int, str]] = {}
    for s, n, t in S:
        if s not in best or (n, t) < best[s]:
            best[s] = (n, t)
    cand = sorted(best, key=lambda s: best[s])
    sets = [{uidx[(a, b)] for a in reps if a in s for b in reps if b not in s and (a, b) in uidx} for s in cand]
    covered_all = set().union(*sets)
    cover = None
    if len(covered_all) == len(universe_pairs):
        obj, chosen, status, _ = solve_cover(sets, list(range(len(universe_pairs))))
        obj2, chosen2, status2, _ = solve_cover(sets, list(range(len(universe_pairs))),
                                                weights=[float(best[s][0]) for s in cand], fix_count=round(obj))
        vsets = [{uidx[(a, b)] for a in reps if a in V for b in reps if b not in V and (a, b) in uidx} for V in gt.long_classes]
        lb, _, lb_status, _ = solve_cover(vsets, list(range(len(universe_pairs))))
        members = []
        for i in sorted(chosen2, key=lambda i: best[cand[i]]):
            n, t = best[cand[i]]
            members.append({"order": n, "table": t, "table_rows": [" ".join(t[r * n:(r + 1) * n]) for r in range(n)],
                            "variety": variety_name(gt, cand[i]), "class_pairs_covered": len(sets[i])})
        cover = {"universe_class_pairs": len(universe_pairs), "candidates": len(cand), "optimum_known_signatures": round(obj),
                 "status": status, "order_weighted_status": status2, "lower_bound_any_quasigroups": round(lb),
                 "lower_bound_status": lb_status, "members": members}

    # Independent control at orders 7..9: Sch-5 implies Sch-23, so the models of
    # Sch-23 that satisfy Sch-5 are exactly the semisymmetric quasigroups (OEIS A076017).
    seq = gates.oeis()["semisymmetric"]
    semisym_control = []
    for n in orders_run:
        ms = models.get(("23", n))
        if ms is not None and "23" in gens.values():
            got = sum(1 for _, _, s in ms if 5 in s)
            semisym_control.append({"order": n, "sch23_models": len(ms), "satisfying_sch5": got,
                                    "oeis": seq["oeis"], "expected": seq["values"].get(n), "ok": got == seq["values"].get(n)})

    out = {
        "semisymmetric_control": semisym_control,
        "complete_upto": complete_upto,
        "generators": gens,
        "runs": [log[k] for k in sorted(log)],
        "curve": curve,
        "least_variety_order": least_variety_order,
        "least_pair_order": least_pair,
        "class_model_counts": class_counts,
        "variety_model_counts": variety_counts,
        "full_cover": cover,
    }
    (gates.DERIVED / "beyond.json").write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    with open(ROOT / "data" / "signatures" / "model_counts.tsv", "w") as f:
        orders = [str(n) for n in range(1, complete_upto + 1)]
        f.write("kind\tname\t" + "\t".join(f"order{n}" for n in orders) + "\n")
        for kind, table in (("class", class_counts), ("variety", variety_counts)):
            for name, c in table.items():
                f.write(f"{kind}\t{name}\t" + "\t".join("" if c[n] is None else str(c[n]) for n in orders) + "\n")
    print(json.dumps({k: out[k] for k in ("complete_upto", "least_variety_order")}, indent=1))
    for row in curve:
        print({k: v for k, v in row.items() if k != "missing_varieties"}, row["missing_varieties"])
    print("pairs exhaustive==mace4:", all(p["exhaustive_least_order"] == p["mace4_order"] for p in least_pair))
    print("cover", {k: v for k, v in (cover or {}).items() if k != "members"})


if __name__ == "__main__":
    main()
