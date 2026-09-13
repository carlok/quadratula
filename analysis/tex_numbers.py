"""Write TeX macros for every result quoted in paper/ and didactic/ (stdlib only).

Usage: python3 analysis/tex_numbers.py <out.tex>

Reads data/derived/{results,beyond,residual}.json, the run logs, the order-3
floor and data/controls/oeis_fingerprints.tsv. Values not available yet become
an italic "pending".
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gates  # noqa: E402

ROOT = gates.ROOT
DER = gates.DERIVED
W = ["Zero", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine"]
UNKNOWN = r"\textit{pending}"  # not "?": the paper uses ? as the operation placeholder


def n(x) -> str:
    return UNKNOWN if x is None else f"{x:,}".replace(",", "{,}")


def p(x) -> str:
    return UNKNOWN if x is None else f"{100 * x:.1f}\\%"


def load(name: str) -> dict:
    path = DER / name
    return json.loads(path.read_text()) if path.exists() else {}


def macros() -> dict[str, str]:
    res = load("results.json")
    bey = load("beyond.json")
    m: dict[str, str] = {}

    c = res["lefloch_consistency"]
    ps = res["pair_space"]
    m["NLaws"] = n(c["n_laws"])
    m["NClasses"] = n(c["n_classes"])
    m["NVarieties"] = n(c["n_long_classes"])
    m["LawPairs"] = n(ps["ordered_law_pairs_offdiagonal"])
    m["LawImp"] = n(ps["law_implications_offdiagonal"])
    m["LawNonImp"] = n(ps["law_non_implications"])
    m["ClassPairs"] = n(ps["ordered_class_pairs_offdiagonal"])
    m["ClassImp"] = n(ps["class_implications_offdiagonal"])
    m["ClassNonImp"] = n(ps["class_non_implications"])

    g1 = {r["order"]: r for r in gates.g1()}
    rows = {r["order"]: r for r in res["per_order"]}
    cum = 0
    for k in range(1, 7):
        r, w = rows[k], W[k]
        cum += r["iso_classes_at_order"]
        m[f"Iso{w}"] = n(r["iso_classes_at_order"])
        m[f"IsoUpTo{w}"] = n(cum)
        m[f"Latin{w}"] = n(g1[k]["latin_count"])
        m[f"ClassSat{w}"] = p(r["class_saturation"])
        m[f"LawSat{w}"] = p(r["law_saturation"])
        m[f"ClassWit{w}"] = n(r["class_refuted_non_implication"])
        m[f"Blocks{w}"] = n(r["law_blocks"])
        m[f"ClassesSep{w}"] = n(r["classes_separated"])
        m[f"Var{w}"] = n(r["varieties_realized"])
        m[f"Sigs{w}"] = n(r["distinct_signatures_cumulative"])
        m[f"FalseRef{w}"] = n(r["law_refuted_but_implication"])
        m[f"ClassNew{w}"] = n(r["class_pairs_new"])
    m["UnrefSix"] = n(rows[6]["class_unrefuted_non_implication"])
    m["EnumSecsSix"] = f"{g1[6]['enumerate_seconds']:.1f}"
    m["SigSecsSix"] = f"{g1[6]['signature_seconds']:.1f}"
    m["Threads"] = n(g1[6]["threads"])
    machine = dict(line.split("\t", 1) for line in (DER / "machine.tsv").read_text().splitlines())
    m["CPU"] = machine["cpu"]

    auts = []
    with gzip.open(ROOT / "data" / "floor" / "order3.tsv.gz", "rt") as f:
        next(f)
        auts = [int(line.split("\t")[2]) for line in f]
    m["AutsThree"] = ", ".join(map(str, auts))
    m["OrbitsThree"] = " + ".join(str(6 // a) for a in auts)

    g3 = gates.g3()
    m["GThreeCells"] = n(len(g3["rows"]))
    m["GThreeAllMatch"] = "all" if all(r["ok"] for r in g3["rows"]) else "not all"
    g4 = gates.g4()
    m["GFourChecks"] = n(sum(r["bit_checks"] for r in g4))
    m["GFourMismatches"] = n(sum(r["mismatches"] for r in g4))
    m["GFourMaxOrder"] = n(max(r["order"] for r in g4))

    cov = res["cover"]
    m["CoverILP"] = n(cov["ilp_optimum"])
    m["CoverGreedy"] = n(cov["greedy_size"])
    m["CoverLP"] = f"{cov['lp_bound']:.0f}"
    m["CoverUniverse"] = n(cov["universe_class_pairs"])
    prof: dict[int, int] = {}
    for q in cov["optimal"]:
        prof[q["order"]] = prof.get(q["order"], 0) + 1
    m["CoverProfile"] = ", ".join(f"{v} of order {k}" for k, v in sorted(prof.items()))

    pa = res["parastrophe"]
    m["LawOrbits"] = n(pa["law_orbits"])
    m["ClassOrbits"] = n(pa["class_orbits"])
    m["VarOrbits"] = n(pa["variety_orbits"])
    m["ClassPairOrbits"] = n(pa["class_pair_orbits"])
    m["NonImpOrbits"] = n(pa["non_implication_orbits"])

    tiny = {t["variety"]: t for t in res["tiny_quasigroups"]}
    m["ZThreeLaws"] = n(tiny["(3)"]["laws_satisfied"])
    m["ZThreeShare"] = p(tiny["(3)"]["law_share_of_non_implications"])
    m["BestTinyLaws"] = n(tiny["(24)"]["laws_satisfied"])
    m["BestTinyShare"] = p(tiny["(24)"]["law_share_of_non_implications"])

    curve = {r["order"]: r for r in bey.get("curve", [])}
    m["CompleteUpTo"] = n(bey.get("complete_upto"))
    for k in (7, 8, 9):
        r, w = curve.get(k), W[k]
        m[f"ClassSat{w}"] = p(r["class_saturation"]) if r else UNKNOWN
        m[f"LawSat{w}"] = p(r["law_saturation"]) if r else UNKNOWN
        m[f"ClassWit{w}"] = n(r["class_non_implications_witnessed"]) if r else UNKNOWN
        m[f"Blocks{w}"] = n(r["class_blocks"]) if r else UNKNOWN
        m[f"ClassesSep{w}"] = n(r["classes_separated"]) if r else UNKNOWN
        m[f"Var{w}"] = n(r["varieties_realised"]) if r else UNKNOWN
    last = curve[max(curve)] if curve else None
    m["LastOrder"] = n(max(curve)) if curve else UNKNOWN
    m["VarMissingLast"] = n(len(last["missing_varieties"])) if last else UNKNOWN
    m["MissingLastList"] = (", ".join(f"\\Sch{{({v.strip('()')})}}" for v in last["missing_varieties"]) or "none") if last else UNKNOWN
    m["NGenerators"] = n(len(bey.get("generators", {}))) if bey else UNKNOWN

    lvo = bey.get("least_variety_order", {})
    hist: dict[str, int] = {}
    for v in lvo.values():
        key = str(v) if v is not None else f">{bey.get('complete_upto')}"
        hist[key] = hist.get(key, 0) + 1
    m["VarLeastHist"] = "; ".join(f"order {k}: {v}" for k, v in sorted(hist.items())) if hist else UNKNOWN

    def least(a: int):
        orders = [q["exhaustive_least_order"] for q in bey.get("least_pair_order", []) if q["a"] == a]
        return None if not orders or None in orders else max(orders)

    m["WitTwentyThree"] = n(least(23))
    m["WitThirtyEight"] = n(least(38))
    m["WitSeventyNine"] = n(least(79))
    lpo = bey.get("least_pair_order", [])
    m["ResidualPairs"] = n(len(lpo)) if lpo else UNKNOWN
    m["MaceAgree"] = n(sum(1 for q in lpo if q["exhaustive_least_order"] == q["mace4_order"])) if lpo else UNKNOWN

    semi = {r["order"]: r for r in bey.get("semisymmetric_control", [])}
    for k in (7, 8, 9):
        m[f"Semi{W[k]}"] = n(semi[k]["satisfying_sch5"]) if k in semi else UNKNOWN

    counts = bey.get("class_model_counts", {})

    def seq(name: str) -> str:
        """Counts for orders 1, 2, ... up to the last order known exactly."""
        cc = counts.get(name)
        if not cc:
            return UNKNOWN
        known = []
        for _, v in sorted(cc.items(), key=lambda kv: int(kv[0])):
            if v is None:
                break
            known.append(str(v))
        return ", ".join(known)

    for name, key in (("Sch-7", "Seven"), ("Sch-38", "ThirtyEight"), ("Sch-23", "TwentyThree"), ("Sch-5", "Five")):
        m[f"Counts{key}"] = seq(name)

    # Sch-7 implies Sch-38, so equal counts at an order mean equal sets of models there.
    a, b = counts.get("Sch-7"), counts.get("Sch-38")
    same_upto, first_diff = None, None
    if a and b:
        for k in sorted(map(int, a)):
            x, y = a[str(k)], b[str(k)]
            if x is None or y is None:
                break
            if x != y:
                first_diff = k
                break
            same_upto = k
    m["SevenThirtyEightSameUpTo"] = n(same_upto)
    m["SevenThirtyEightFirstDiff"] = n(first_diff)

    fc = bey.get("full_cover")
    m["FullCover"] = n(fc["optimum_known_signatures"]) if fc else UNKNOWN
    m["FullCoverLB"] = n(fc["lower_bound_any_quasigroups"]) if fc else UNKNOWN

    ic = load("independent_checks.json")
    sat = ic.get("sat", {})
    if sat.get("status") == "ran":
        m["SatQueries"] = n(len(sat["rows"]))
        m["SatAllOk"] = "All" if sat["all_ok"] else "Not all"
        m["SatMaxSeconds"] = f"{max(r['seconds'] for r in sat['rows']):.0f}"
    else:
        m["SatQueries"] = m["SatAllOk"] = m["SatMaxSeconds"] = UNKNOWN
    covs = ic.get("cover", {})
    m["IndepCoversOk"] = ("both" if covs.get("cover_order6", {}).get("ok") and covs.get("full_cover", {}).get("ok")
                          else UNKNOWN)
    mo = ic.get("models", {})
    m["IndepModelPairs"] = n(mo["pairs"]) if mo.get("all_found") else UNKNOWN

    fp = ROOT / "data" / "controls" / "oeis_fingerprints.tsv"
    if fp.exists():
        rows_fp = [line.split("\t") for line in fp.read_text().splitlines()[1:]]
        seqs = {r[0] for r in rows_fp}
        matched = sorted({r[2] for r in rows_fp if r[2] != "-"})
        m["OeisSeqs"] = n(len(seqs))
        m["OeisMatched"] = ", ".join(matched)
        m["OeisNone"] = n(len({r[0] for r in rows_fp if r[2] == "-"}))
    return m


def main() -> None:
    out = Path(sys.argv[1])
    m = macros()
    lines = ["% Generated by analysis/tex_numbers.py. Do not edit."]
    lines += [f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in sorted(m.items())]
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out} ({len(m)} macros)")


if __name__ == "__main__":
    main()
