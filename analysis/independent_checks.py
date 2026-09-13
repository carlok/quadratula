"""Independent re-checks of the headline claims (stdlib only, plus kissat for --sat).

Shares no code with the Rust engine or with compute.py/beyond.py: it has its own
law parser and evaluator, and reads Le Floch's raw files directly.

  --cover   Both covers witness what they claim, recomputed from the tables.
  --models  Every residual least witness order has a witness in the stored model
            files, re-evaluated here.
  --sat     SAT (kissat): encoding controls, nonexistence of witnesses below the
            least orders, existence at the least orders (models decoded and
            re-evaluated here).

Writes data/derived/independent_checks.json (merging with earlier runs).
"""

from __future__ import annotations

import gzip
import itertools
import json
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LF = ROOT / "data" / "lefloch"
DER = ROOT / "data" / "derived"
OUT = DER / "independent_checks.json"
OPS = {"*": 0, "//": 1, "\\\\": 2}
TOKEN = re.compile(r"\*|//|\\\\|\(|\)|=|[a-z]")


# ---------- parsing and evaluation -------------------------------------------------

def parse_law(text: str):
    toks = TOKEN.findall(text)
    assert "".join(toks) == text.replace(" ", ""), f"unparsed characters in {text!r}"
    pos = 0

    def atom():
        nonlocal pos
        t = toks[pos]
        pos += 1
        if t == "(":
            e = expr()
            assert toks[pos] == ")"
            pos += 1
            return e
        assert t.isalpha(), t
        return ("v", t)

    def expr():
        nonlocal pos
        left = atom()
        if pos < len(toks) and toks[pos] in OPS:
            o = OPS[toks[pos]]
            pos += 1
            return ("o", o, left, atom())
        return left

    lhs = expr()
    assert toks[pos] == "="
    pos += 1
    rhs = expr()
    assert pos == len(toks)
    return lhs, rhs


def load_laws() -> dict[int, tuple]:
    laws = {}
    for line in (LF / "quasigroup-equations.txt").read_text().splitlines():
        label, text = line.split(": ", 1)
        laws[int(label[4:])] = parse_law(text)
    assert sorted(laws) == list(range(1, 991))
    return laws


def variables(t, acc=None):
    acc = [] if acc is None else acc
    if t[0] == "v":
        if t[1] not in acc:
            acc.append(t[1])
    else:
        variables(t[2], acc)
        variables(t[3], acc)
    return acc


def compile_term(t, names):
    if t[0] == "v":
        i = names.index(t[1])
        return lambda e, T: e[i]
    o, left, right = t[1], compile_term(t[2], names), compile_term(t[3], names)
    return lambda e, T: T[o][left(e, T)][right(e, T)]


class Law:
    def __init__(self, lhs, rhs):
        self.names = variables(lhs, variables(rhs, []))
        self.left = compile_term(lhs, self.names)
        self.right = compile_term(rhs, self.names)

    def holds(self, T, n) -> bool:
        return all(self.left(e, T) == self.right(e, T) for e in itertools.product(range(n), repeat=len(self.names)))


def tables(star: list[list[int]]):
    n = len(star)
    for row in star:
        assert sorted(row) == list(range(n)), "row not a permutation"
    for col in zip(*star):
        assert sorted(col) == list(range(n)), "column not a permutation"
    slash = [[0] * n for _ in range(n)]
    back = [[0] * n for _ in range(n)]
    for c in range(n):
        for a in range(n):
            slash[a][star[c][a]] = c  # a // b = c  iff  c * a = b
    for b in range(n):
        for c in range(n):
            back[star[b][c]][b] = c  # a \\ b = c  iff  b * c = a
    return (star, slash, back), n


def from_digits(s: str, n: int):
    v = [int(ch, 36) for ch in s]
    return [v[i * n:(i + 1) * n] for i in range(n)]


def signature(laws: dict[int, Law], star) -> frozenset[int]:
    T, n = tables(star)
    return frozenset(i for i, law in laws.items() if law.holds(T, n))


# ---------- ground truth, read directly --------------------------------------------

def ground_truth():
    classes = {}
    for line in (LF / "quasigroup-classes.txt").read_text().splitlines():
        m = re.fullmatch(r"Sch-(\d+) class \[([\d, ]+)\]", line)
        classes[int(m.group(1))] = [int(x) for x in m.group(2).split(",")]
    closure = {int(k): frozenset(v) for k, v in json.loads((LF / "quasigroup_eq_to_long_class.json").read_text()).items()}
    closed_sets = {frozenset(c) for c in json.loads((LF / "quasigroup_long_classes.json").read_text())}
    reps = sorted(classes)
    nonimp = {(a, b) for a in reps for b in reps if a != b and b not in closure[a]}
    return reps, closure, closed_sets, nonimp


# ---------- checks -----------------------------------------------------------------

def check_cover(laws) -> dict:
    reps, closure, closed_sets, nonimp = ground_truth()
    res = json.loads((DER / "results.json").read_text())
    bey = json.loads((DER / "beyond.json").read_text())
    out = {"class_non_implications": len(nonimp)}

    def run(members):
        covered, bad = set(), []
        for q in members:
            n = q["order"]
            s = signature(laws, from_digits(q["table"], n))
            if s not in closed_sets:
                bad.append(q["table"])
            covered |= {(a, b) for a in reps if a in s for b in reps if b not in s}
        return covered, bad

    cov6, bad6 = run(res["cover"]["optimal"])
    out["cover_order6"] = {
        "members": len(res["cover"]["optimal"]),
        "claimed_universe": res["cover"]["universe_class_pairs"],
        "covered_non_implications": len(cov6 & nonimp),
        "covers_an_implication": len(cov6 - nonimp),
        "signatures_not_closed_sets": len(bad6),
        "ok": len(cov6 & nonimp) == res["cover"]["universe_class_pairs"] and not (cov6 - nonimp) and not bad6,
    }
    fc = bey.get("full_cover")
    if fc:
        covall, badall = run(fc["members"])
        out["full_cover"] = {
            "members": len(fc["members"]),
            "orders": sorted({q["order"] for q in fc["members"]}),
            "covered_non_implications": len(covall & nonimp),
            "covers_an_implication": len(covall - nonimp),
            "signatures_not_closed_sets": len(badall),
            "ok": covall == nonimp and not badall,
        }
    return out


def check_models(laws) -> dict:
    bey = json.loads((DER / "beyond.json").read_text())
    rows = []
    cache: dict[tuple[int, int], list] = {}
    for p in bey["least_pair_order"]:
        a, b, o = p["a"], p["b"], p["exhaustive_least_order"]
        key = (a, o)
        if key not in cache:
            path = ROOT / "data" / "beyond" / f"models_o{o}_{a}.tsv.gz"
            with gzip.open(path, "rt") as f:
                next(f)
                cache[key] = [from_digits(line.split("\t")[1], o) for line in f]
        found = False
        for star in cache[key]:
            T, n = tables(star)
            if not laws[b].holds(T, n) and laws[a].holds(T, n):
                found = True
                break
        rows.append({"a": a, "b": b, "order": o, "witness_in_stored_models": found})
    return {"pairs": len(rows), "all_found": all(r["witness_in_stored_models"] for r in rows), "rows": rows}


def cnf(n: int, parsed: dict[int, tuple], hyps: list[int], goal: int | None) -> tuple[int, list[list[int]]]:
    def m(a, b, c):
        return 1 + (a * n + b) * n + c

    def rel(o, x, y, z):  # the triple fact "x op y = z" as a literal
        return m(x, y, z) if o == 0 else (m(z, x, y) if o == 1 else m(y, z, x))

    cl: list[list[int]] = []
    rng = range(n)
    for a, b in itertools.product(rng, rng):
        for pick in (lambda c: m(a, b, c), lambda c: m(a, c, b), lambda c: m(c, a, b)):
            lits = [pick(c) for c in rng]
            cl.append(lits)
            cl += [[-x, -y] for x, y in itertools.combinations(lits, 2)]
    if n >= 2:
        cl.append([m(0, 0, 0), m(0, 0, 1)])  # every quasigroup is isomorphic to one with T[0][0] <= 1

    def inner(t, env):
        """(left value, right value, literals) for the top operation of side t."""
        assert t[0] == "o"
        opts = []
        for lv, ll in values(t[2], env):
            for rv, rl in values(t[3], env):
                opts.append((lv, rv, ll + rl))
        return opts

    def values(t, env):
        if t[0] == "v":
            return [(env[t[1]], [])]
        out = []
        for lv, rv, lits in inner(t, env):
            for v in rng:
                out.append((v, lits + [rel(t[1], lv, rv, v)]))
        return out

    for h in hyps:
        lhs, rhs = parsed[h]
        names = variables(lhs, variables(rhs, []))
        for vals in itertools.product(rng, repeat=len(names)):
            env = dict(zip(names, vals))
            for w, l1 in values(lhs, env):
                for lv, rv, l2 in inner(rhs, env):
                    cl.append([-x for x in l1 + l2] + [rel(rhs[1], lv, rv, w)])
    nvars = n ** 3
    if goal is not None:
        lhs, rhs = parsed[goal]
        names = variables(lhs, variables(rhs, []))
        selectors = []
        for vals in itertools.product(rng, repeat=len(names)):
            env = dict(zip(names, vals))
            nvars += 1
            s = nvars
            selectors.append(s)
            for w, l1 in values(lhs, env):
                for lv, rv, l2 in inner(rhs, env):
                    cl.append([-s] + [-x for x in l1 + l2] + [-rel(rhs[1], lv, rv, w)])
        cl.append(selectors)
    return nvars, cl


def solve(n, parsed, laws, hyps, goal, limit) -> dict:
    nvars, clauses = cnf(n, parsed, hyps, goal)
    with tempfile.NamedTemporaryFile("w", suffix=".cnf", delete=False) as f:
        f.write(f"p cnf {nvars} {len(clauses)}\n")
        for c in clauses:
            f.write(" ".join(map(str, c)) + " 0\n")
        path = f.name
    t0 = time.time()
    try:
        p = subprocess.run(["kissat", "-q", f"--time={limit}", path], capture_output=True, text=True, timeout=limit + 60)
        code, stdout = p.returncode, p.stdout
    except subprocess.TimeoutExpired:
        code, stdout = 0, ""
    Path(path).unlink()
    result = {10: "SAT", 20: "UNSAT"}.get(code, "UNKNOWN")
    row = {"order": n, "hypotheses": hyps, "goal_violated": goal, "result": result,
           "seconds": round(time.time() - t0, 1), "variables": nvars, "clauses": len(clauses)}
    if result == "SAT":
        true = {int(x) for line in stdout.splitlines() if line.startswith("v") for x in line[1:].split() if int(x) > 0}
        star = [[next(c for c in range(n) if 1 + (a * n + b) * n + c in true) for b in range(n)] for a in range(n)]
        T, _ = tables(star)
        row["model_rechecked"] = all(laws[h].holds(T, n) for h in hyps) and (goal is None or not laws[goal].holds(T, n))
        row["table"] = "".join(str(v) for r in star for v in r)
    return row


def check_sat(laws, parsed, limit: int) -> dict:
    if shutil.which("kissat") is None:
        return {"status": "not run: kissat not on PATH"}
    queries = [
        # encoding controls (Example: associativity without commutativity needs a non-abelian group)
        ("control", 5, [253], 1, "UNSAT"),
        ("control", 6, [253], 1, "SAT"),
        # nonexistence below the least witness orders
        ("claim", 7, [38], None, "UNSAT"),
        ("claim", 7, [79], None, "UNSAT"),
        ("claim", 7, [23], 5, "UNSAT"),
        ("claim", 8, [23], 5, "UNSAT"),
        # existence at the least witness orders
        ("claim", 8, [38], 1, "SAT"),
        ("claim", 8, [79], 1, "SAT"),
        ("claim", 9, [23], 5, "SAT"),
    ]
    rows = []
    for kind, n, hyps, goal, expect in queries:
        row = solve(n, parsed, laws, hyps, goal, limit)
        row.update({"kind": kind, "expected": expect, "ok": row["result"] == expect and row.get("model_rechecked", True)})
        rows.append(row)
        print(json.dumps({k: v for k, v in row.items() if k != "table"}), flush=True)
    return {"status": "ran", "time_limit_seconds": limit, "all_ok": all(r["ok"] for r in rows), "rows": rows,
            "note": "Sch-5 implies Sch-29, Sch-93 and Sch-183, so Sch-23 with not Sch-5 unsatisfiable covers all four "
                    "Sch-23 residual pairs; no model of Sch-38 or Sch-79 at all covers their pairs."}


def main() -> None:
    parsed = load_laws()
    laws = {i: Law(*t) for i, t in parsed.items()}
    out = json.loads(OUT.read_text()) if OUT.exists() else {}
    if "--cover" in sys.argv:
        out["cover"] = check_cover(laws)
        print(json.dumps(out["cover"], indent=1))
    if "--models" in sys.argv:
        out["models"] = check_models(laws)
        print(json.dumps({k: v for k, v in out["models"].items() if k != "rows"}, indent=1))
    if "--sat" in sys.argv:
        limit = int(sys.argv[sys.argv.index("--limit") + 1]) if "--limit" in sys.argv else 7200
        out["sat"] = check_sat(laws, parsed, limit)
    OUT.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
