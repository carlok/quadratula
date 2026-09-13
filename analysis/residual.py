"""Beyond the floor: smallest finite witnesses for the residual non-implications.

For every class-level non-implication with no witness of order <= 6, run Mace4
(quasigroup axioms, hypothesis law, negated goal) over domain sizes 7, 8, ...
Mace4 moves to size s+1 only after finishing size s without a model, so the
first model it reports, together with the exhaustive floor through order 6,
gives the smallest order of a witness. That minimality relies on Mace4's
search being complete. Every model is re-evaluated with our own engine
(`qg check`): it must satisfy the hypothesis, violate the goal, and have a
signature that is one of Le Floch's closed sets.

Writes data/derived/residual.json. Skips cleanly (status "not run") without mace4.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gates  # noqa: E402
import lefloch  # noqa: E402

ROOT = gates.ROOT
OUT = gates.DERIVED / "residual.json"
START, END, MAX_SECONDS = 7, 14, 600
OP_DECLARE = 'op(400, infix, "*"). op(400, infix, "//"). op(400, infix, "\\\\").'
AXIOMS = "(y // x) * y = x.  y // (x * y) = x.  x * (y \\\\ x) = y.  (x * y) \\\\ x = y.  x \\\\ (y // x) = y.  (x \\\\ y) // x = y."


def mace4_input(hyp: str, goal: str) -> str:
    return (
        f"{OP_DECLARE}\nassign(domain_size, {START}). assign(end_size, {END}). assign(max_seconds, {MAX_SECONDS}).\n"
        f"formulas(assumptions).\n{AXIOMS}\n{hyp}.\nend_of_list.\n"
        f"formulas(goals).\n{goal}.\nend_of_list.\n"
    )


def run(mace4: str, a: int, b: int, laws: dict[int, str]) -> dict:
    out = subprocess.run([mace4], input=mace4_input(laws[a], laws[b]), capture_output=True, text=True).stdout
    sizes = [int(s) for s in re.findall(r"DOMAIN SIZE (\d+)", out)]
    exit_reason = (re.findall(r"exit \((\w+)\)", out) or ["unknown"])[-1]
    m = re.search(r"interpretation\( (\d+),.*?function\(\*\(_,_\), \[(.*?)\]\)", out, re.S)
    row = {"a": a, "b": b, "sizes_searched": sizes, "exit": exit_reason, "order": None, "table": None}
    if m:
        n = int(m.group(1))
        entries = [int(x) for x in m.group(2).replace("\n", "").split(",")]
        assert len(entries) == n * n
        row["order"], row["table"] = n, entries
        # minimality needs every size from START to n-1 to have been searched to the end
        row["searched_all_smaller_sizes"] = sizes == list(range(START, n + 1))
    return row


def main() -> None:
    mace4 = shutil.which("mace4")
    if mace4 is None:
        OUT.write_text(json.dumps({"status": "not run: mace4 not on PATH"}, indent=1) + "\n")
        print("mace4 not found; residual step skipped")
        return
    gt = lefloch.load()
    res = gates.results()
    pairs = [tuple(p) for p in res["unrefuted_class_non_implications_at_order_6"]]
    with ThreadPoolExecutor(max_workers=10) as ex:
        rows = list(ex.map(lambda p: run(mace4, p[0], p[1], gt.laws), pairs))

    # independent re-evaluation with the Rust engine
    found = [r for r in rows if r["table"] is not None]
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
        for r in found:
            f.write(f"{r['a']}-{r['b']}\t{r['order']}\t{','.join(map(str, r['table']))}\n")
        tmp = f.name
    qg = ROOT / "target" / "release" / "qg"
    checked = subprocess.run([str(qg), "check", tmp], capture_output=True, text=True, check=True).stdout
    Path(tmp).unlink()
    sat = {}
    for line in checked.splitlines():
        key, s = line.split("\t")
        assert not s.startswith("NOT_LATIN"), line
        sat[key] = frozenset() if s == "-" else frozenset(map(int, s.split(",")))
    long_set = set(gt.long_classes)
    for r in found:
        s = sat[f"{r['a']}-{r['b']}"]
        r["satisfies_hypothesis"] = r["a"] in s
        r["violates_goal"] = r["b"] not in s
        r["signature_is_variety"] = s in long_set
        r["laws_satisfied"] = len(s)
        r["satisfied"] = sorted(s)

    # what the floor plus these witnesses realises (a lower bound for the floor at the largest witness order)
    reps = gt.reps
    sigs = {frozenset(x["satisfied"]) for x in res["signatures"]} | {frozenset(r["satisfied"]) for r in found}
    closed = {frozenset(range(1, lefloch.N_LAWS + 1))} | sigs
    while True:
        new = {x & y for x in closed for y in sigs} - closed
        if not new:
            break
        closed |= new
    columns = {i: tuple(sorted(k for k, s in enumerate(sorted(sigs, key=sorted)) if i in s)) for i in reps}
    witnessed_all = all(r["table"] is not None for r in rows)
    orders = [r["order"] for r in found]
    summary = {
        "status": "ran",
        # stdin must be closed: mace4 reads its input from stdin and would wait forever.
        "mace4": (re.findall(r"Mace4 \(\d+\) version [^,\n]+", subprocess.run([mace4, "--version"], capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=30).stdout) or ["mace4 (version unknown)"])[0],
        "search": {"start_size": START, "end_size": END, "max_seconds": MAX_SECONDS},
        "pairs": len(rows),
        "witnessed": len(found),
        "all_witnessed": witnessed_all,
        "all_verified": all(r["satisfies_hypothesis"] and r["violates_goal"] for r in found),
        "all_signatures_are_varieties": all(r["signature_is_variety"] for r in found),
        "all_minimal": all(r.get("searched_all_smaller_sizes") for r in found),
        "order_histogram": {str(o): orders.count(o) for o in sorted(set(orders))},
        "max_witness_order": max(orders) if orders else None,
        "distinct_new_signatures": len({frozenset(r["satisfied"]) for r in found} - {frozenset(x["satisfied"]) for x in res["signatures"]}),
        "with_floor_varieties_realised": len(closed & long_set),
        "with_floor_closed_sets_not_varieties": len(closed - long_set),
        "with_floor_class_blocks": len(set(columns.values())),
    }
    OUT.write_text(json.dumps({"summary": summary, "rows": rows}, indent=1, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
