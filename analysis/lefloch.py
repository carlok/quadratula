"""Load Le Floch's ancillary data (arXiv:2603.29909v1, CC BY 4.0) as ground truth.

Two different objects appear in the paper and must not be confused:
  * 47 equivalence classes of single laws (quasigroup-classes.txt), and
  * 114 logically closed subsets of the 990 laws, i.e. the varieties defined by
    conjunctions of laws (quasigroup_long_classes.json), including the empty
    set and the full set.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "lefloch"
N_LAWS = 990


@dataclass(frozen=True)
class GroundTruth:
    laws: dict[int, str]                      # Sch id -> law text
    classes: dict[int, tuple[int, ...]]       # representative -> members
    rep_of: dict[int, int]                    # law -> representative
    rep_implies: dict[int, frozenset[int]]    # representative -> implied representatives (incl. itself)
    long_classes: tuple[frozenset[int], ...]  # the 114 closed sets
    closure_of_rep: dict[int, frozenset[int]] # representative -> its closed set of laws

    @property
    def reps(self) -> list[int]:
        return sorted(self.classes)

    def implies(self, i: int, j: int) -> bool:
        """Law i implies law j (as quasigroup laws)."""
        return j in self.closure_of_rep[self.rep_of[i]]


def load() -> GroundTruth:
    laws = {}
    for line in (DATA / "quasigroup-equations.txt").read_text().splitlines():
        m = re.fullmatch(r"Sch-(\d+): (.*)", line)
        assert m, line
        laws[int(m.group(1))] = m.group(2)
    assert sorted(laws) == list(range(1, N_LAWS + 1))

    classes = {}
    for line in (DATA / "quasigroup-classes.txt").read_text().splitlines():
        m = re.fullmatch(r"Sch-(\d+) class \[([\d, ]+)\]", line)
        assert m, line
        classes[int(m.group(1))] = tuple(int(x) for x in m.group(2).split(","))
    rep_of = {law: rep for rep, members in classes.items() for law in members}
    assert sorted(rep_of) == list(range(1, N_LAWS + 1)), "classes must partition the laws"

    rep_implies = {}
    for line in (DATA / "quasigroup-implications.txt").read_text().splitlines():
        m = re.fullmatch(r"Sch-(\d+): .* ⇒ \{([\d, ]+)\}", line)
        assert m, line
        rep_implies[int(m.group(1))] = frozenset(int(x) for x in m.group(2).split(","))
    assert set(rep_implies) == set(classes)

    long_classes = tuple(frozenset(c) for c in json.loads((DATA / "quasigroup_long_classes.json").read_text()))
    e2l = json.loads((DATA / "quasigroup_eq_to_long_class.json").read_text())
    closure_of_rep = {int(k): frozenset(v) for k, v in e2l.items()}
    assert set(closure_of_rep) == set(classes)

    return GroundTruth(laws, classes, rep_of, rep_implies, long_classes, closure_of_rep)


def consistency(gt: GroundTruth) -> dict:
    """Internal checks on Le Floch's files; every value is reported in CONTROLS.md."""
    reps = gt.reps
    reflexive = all(r in gt.rep_implies[r] for r in reps)
    transitive = all(gt.rep_implies[b] <= gt.rep_implies[a] for a in reps for b in gt.rep_implies[a])
    # The closed set of a representative is the union of the classes it implies.
    closure_matches = all(
        gt.closure_of_rep[r] == frozenset(x for s in gt.rep_implies[r] for x in gt.classes[s]) for r in reps
    )
    sets = set(gt.long_classes)
    antisymmetric = all(
        not (b in gt.rep_implies[a] and a in gt.rep_implies[b]) for a in reps for b in reps if a != b
    )
    return {
        "n_laws": len(gt.laws),
        "n_classes": len(gt.classes),
        "n_long_classes": len(gt.long_classes),
        "n_long_classes_distinct": len(sets),
        "long_classes_contain_empty": frozenset() in sets,
        "long_classes_contain_full": frozenset(range(1, N_LAWS + 1)) in sets,
        "long_classes_intersection_closed": all((a & b) in sets for a in sets for b in sets),
        "long_classes_are_unions_of_classes": all(
            all(gt.rep_of[x] in {gt.rep_of[y] for y in c} and set(gt.classes[gt.rep_of[x]]) <= c for x in c)
            for c in sets
        ),
        "single_law_closures_among_long_classes": all(c in sets for c in gt.closure_of_rep.values()),
        "rep_implications_reflexive": reflexive,
        "rep_implications_transitive": transitive,
        "rep_implications_antisymmetric": antisymmetric,
        "closure_matches_implication_file": closure_matches,
    }
