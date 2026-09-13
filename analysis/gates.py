"""Gate evaluation from the run logs in data/derived and the OEIS terms in
data/controls/oeis.tsv. Used by both the report generator and pytest."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DERIVED = ROOT / "data" / "derived"
ORDERS = list(range(1, 7))


def kv_lines(path: Path) -> list[dict[str, str]]:
    """Lines of alternating key<TAB>value fields."""
    out = []
    for line in path.read_text().splitlines():
        f = line.split("\t")
        out.append(dict(zip(f[0::2], f[1::2])))
    return out


def oeis() -> dict[str, dict]:
    rows = (ROOT / "data" / "controls" / "oeis.tsv").read_text().splitlines()
    out = {}
    for line in rows[1:]:
        col, aid, offset, terms, name = line.split("\t")
        t = [int(x) for x in terms.split(",")]
        start = 0 if offset == "diagonal" else int(offset)
        out[col] = {"oeis": aid, "name": name, "values": {start + i: v for i, v in enumerate(t)}}
    return out


def g1() -> list[dict]:
    seq = oeis()
    enum = {int(d["order"]): d for d in kv_lines(DERIVED / "enumerate.tsv")}
    latin = {int(d["order"]): d for d in kv_lines(DERIVED / "count_latin.tsv")}
    rows = []
    for n in ORDERS:
        r = {
            "order": n,
            "iso_classes": int(enum[n]["iso_classes"]),
            "oeis_iso": seq["all"]["values"][n],
            "orbit_sum": int(enum[n]["sum_n!/aut"]),
            "latin_count": int(latin[n]["latin_squares"]),
            "oeis_latin": seq["latin_squares"]["values"][n],
            "enumerate_seconds": float(enum[n]["enumerate_seconds"]),
            "signature_seconds": float(enum[n]["signature_seconds"]),
            "latin_seconds": float(latin[n]["seconds"]),
            "threads": int(enum[n]["threads"]),
        }
        r["ok"] = r["iso_classes"] == r["oeis_iso"] and r["orbit_sum"] == r["latin_count"] == r["oeis_latin"]
        rows.append(r)
    return rows


CONTROL_COLUMNS = ["groups", "abelian_groups", "commutative", "semisymmetric", "totally_symmetric", "idempotent"]
SCHROEDER_TWINS = {"sch_assoc": "groups", "sch1_commutative": "commutative", "sch5_semisymmetric": "semisymmetric"}


def g3() -> dict:
    lines = (DERIVED / "controls.tsv").read_text().splitlines()
    assoc_id = lines[0].split("\t")[1]
    header = lines[1].split("\t")
    table = {}
    for line in lines[2:]:
        vals = dict(zip(header, line.split("\t")))
        table[int(vals["order"])] = {k: int(v) for k, v in vals.items() if k != "order"}
    seq = oeis()
    rows = []
    for col in CONTROL_COLUMNS:
        for n in ORDERS:
            rows.append({"control": col, "oeis": seq[col]["oeis"], "order": n, "computed": table[n][col],
                         "expected": seq[col]["values"][n], "ok": table[n][col] == seq[col]["values"][n]})
    twins = [{"schroeder": k, "direct": v, "order": n, "schroeder_count": table[n][k], "direct_count": table[n][v],
              "ok": table[n][k] == table[n][v]} for k, v in SCHROEDER_TWINS.items() for n in ORDERS]
    return {"associativity_id": assoc_id, "rows": rows, "twins": twins, "table": table}


def g4() -> list[dict]:
    return [{k: int(v) for k, v in d.items()} for d in kv_lines(DERIVED / "g4.tsv")]


def parastrophe_summary() -> list[dict]:
    rows = []
    for line in (DERIVED / "parastrophe.tsv").read_text().splitlines()[1:]:
        sigma, mapped, bij, _ = line.split("\t")
        rows.append({"sigma": sigma, "mapped": int(mapped), "bijection": bij == "true"})
    return rows


def results() -> dict:
    return json.loads((DERIVED / "results.json").read_text())
