"""Verification gates over the full-size runs (run `make reproduce` first)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "analysis"))

import gates  # noqa: E402
import lefloch  # noqa: E402


def test_lefloch_files_consistent():
    c = lefloch.consistency(lefloch.load())
    assert all(v is True for k, v in c.items() if isinstance(v, bool)), c
    assert (c["n_laws"], c["n_classes"], c["n_long_classes"]) == (990, 47, 114)


def test_g1_enumeration_matches_oeis():
    for r in gates.g1():
        assert r["ok"], r


def test_g2_no_false_refutation():
    res = gates.results()
    for row in res["per_order"]:
        assert row["law_refuted_but_implication"] == 0, row["order"]
        assert row["class_refuted_but_implication"] == 0, row["order"]
        assert row["signatures_not_among_long_classes"] == 0, row["order"]
        assert row["closed_sets_not_long_classes"] == 0, row["order"]
        assert row["classes_split_by_floor"] == 0, row["order"]


def test_g3_classical_counts_match_oeis():
    g = gates.g3()
    assert g["associativity_id"] != "absent"
    for r in g["rows"] + g["twins"]:
        assert r["ok"], r


def test_residual_witnesses_verified():
    import json

    path = gates.DERIVED / "residual.json"
    data = json.loads(path.read_text())
    if data.get("status", "ran") != "ran":
        return  # mace4 absent: nothing to check
    for r in data["rows"]:
        if r["table"] is None:
            continue
        assert r["order"] >= 7, r
        assert r["satisfies_hypothesis"] and r["violates_goal"], r
        assert r["signature_is_variety"], r


def test_g5_beyond_order_6():
    import json

    b = json.loads((gates.DERIVED / "beyond.json").read_text())
    assert b["complete_upto"] >= 7
    assert b["semisymmetric_control"] and all(r["ok"] for r in b["semisymmetric_control"])
    # exhaustive least witness orders agree with Mace4 wherever both exist
    for p in b["least_pair_order"]:
        if p["exhaustive_least_order"] is not None and p["mace4_order"] is not None:
            assert p["exhaustive_least_order"] == p["mace4_order"], p
    curve = b["curve"]
    for x, y in zip(curve, curve[1:]):
        assert y["class_non_implications_witnessed"] >= x["class_non_implications_witnessed"]
        assert y["varieties_realised"] >= x["varieties_realised"]


def test_g6_independent_checks():
    import json

    path = gates.DERIVED / "independent_checks.json"
    if not path.exists():
        return  # `make independent-checks` not run
    ic = json.loads(path.read_text())
    if "cover" in ic:
        assert ic["cover"]["cover_order6"]["ok"], ic["cover"]["cover_order6"]
        if "full_cover" in ic["cover"]:
            assert ic["cover"]["full_cover"]["ok"], ic["cover"]["full_cover"]
    if "models" in ic:
        assert ic["models"]["all_found"]
    if ic.get("sat", {}).get("status") == "ran":
        bad = [r for r in ic["sat"]["rows"] if not r["ok"]]
        assert not bad, bad


def test_g4_parastrophe_equivariance():
    assert all(r["bijection"] and r["mapped"] == 990 for r in gates.parastrophe_summary())
    rows = gates.g4()
    assert [r["order"] for r in rows] == [1, 2, 3, 4, 5]
    assert all(r["mismatches"] == 0 and r["bit_checks"] > 0 for r in rows)
    res = gates.results()
    assert all(res["parastrophe"]["signature_sets_closed_under_action"].values())
