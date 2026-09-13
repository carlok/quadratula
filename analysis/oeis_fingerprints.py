"""One-off OEIS lookup of model counts by order (needs network; not part of `make reproduce`).

For each of Le Floch's 47 classes, the number of quasigroups of order 1..6 (up to
isomorphism) satisfying the class is a short integer sequence. Sequences with
at least three terms greater than 1 are searched on OEIS and the matches are
pinned in data/controls/oeis_fingerprints.tsv with the retrieval date. A match
of six initial terms is a lead for characterising a class, not a proof.
"""

from __future__ import annotations

import datetime
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gates  # noqa: E402
import lefloch  # noqa: E402

OUT = gates.ROOT / "data" / "controls" / "oeis_fingerprints.tsv"
STRUCTURAL = ("quasigroup", "latin", "loop", "systems with", "binary operation", "groupoid", "magma")


def fetch_page(q: str, start: int, attempts: int = 5) -> list[dict] | None:
    """One page of OEIS results; [] when OEIS reports no (more) matches, None if it never answered with JSON."""
    for i in range(attempts):
        raw = subprocess.run(
            ["curl", "-sL", "-A", "Mozilla/5.0", "--get", "--data-urlencode", f"q={q}", "--data-urlencode", "fmt=json",
             "--data-urlencode", f"start={start}", "https://oeis.org/search"],
            capture_output=True, text=True,
        ).stdout
        try:
            d = json.loads(raw) if raw.strip() else None  # OEIS answers `null` when nothing matches
        except json.JSONDecodeError:
            time.sleep(2.0 * (i + 1))  # rate limit or error page: back off and retry
            continue
        if d is None:
            return []
        return d if isinstance(d, list) else (d.get("results") or [])
    return None


KEYWORDS = ("quasigroup", "latin", "loop")


def search(terms: list[int]) -> tuple[list[dict], bool]:
    """OEIS entries containing the terms and one of KEYWORDS. Returns (hits, every scan complete).

    Searching terms alone returns hundreds of unrelated entries for short sequences
    and OEIS ranks structural ones low, so the keyword is part of the query.
    """
    seen: dict[int, dict] = {}
    complete = True
    for kw in KEYWORDS:
        hits, ok = search_query(f"{','.join(map(str, terms))} {kw}")
        complete &= ok
        for h in hits:
            seen.setdefault(h["number"], h)
        time.sleep(1.0)
    return list(seen.values()), complete


def search_query(q: str, max_hits: int = 200) -> tuple[list[dict], bool]:
    """OEIS hits for a query, following pagination. Returns (hits, scan_complete)."""
    hits: list[dict] = []
    while len(hits) < max_hits:
        page = fetch_page(q, len(hits))
        if page is None:
            return hits, False
        hits += page
        if len(page) < 10:
            return hits, True
        time.sleep(1.0)
    return hits, False


def main() -> None:
    gt = lefloch.load()
    res = gates.results()
    seqs: dict[tuple[int, ...], list[int]] = {}
    for r in gt.reps:
        c = gt.closure_of_rep[r]
        counts = tuple(sum(s["iso_classes_by_order"][str(n)] for s in res["signatures"] if r in s["satisfied"]) for n in range(1, 7))
        seqs.setdefault(counts, []).append(r)
    today = datetime.date.today().isoformat()
    lines = ["counts_order1_6\tclasses\toeis\tname\tretrieved"]
    for counts, reps in sorted(seqs.items()):
        if sum(1 for x in counts if x > 1) < 3:
            continue
        hits, complete = search(list(counts))
        time.sleep(1.0)
        cls = ",".join(f"Sch-{r}" for r in reps)
        scanned = f"{len(hits)} hits scanned" + ("" if complete else ", scan incomplete")
        # Short sequences match many unrelated entries; keep only structural ones.
        keep = [h for h in hits if any(w in h["name"].lower() for w in STRUCTURAL)]
        if not keep:
            verdict = "no quasigroup-related entry" if complete else "none found (incomplete scan)"
            lines.append(f"{','.join(map(str, counts))}\t{cls}\t-\t{verdict} ({scanned})\t{today}")
        for h in keep:
            lines.append(f"{','.join(map(str, counts))}\t{cls}\tA{h['number']:06d}\t{h['name']} ({scanned})\t{today}")
    OUT.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
