"""Quality screening for candidate seed genomes.

A seed founds a whole transmission chain, so one bad genome is inherited by
every infection descending from it. On the Virginia Delta wave the second
importation was `USA/VA-CDC-LC0041349/2021`, a mixed Alpha/Delta assembly
(the Alpha HV69/70 deletion, three Alpha markers, 6.9% ambiguous bases, 47%
of the 3' 2kb missing). It founded the giant transmission component, and
Nextclade's diagnostics in the ncov build then excluded 10,180 of 14,400
sampled genomes -- 71% of the state's sample -- for reversions and putative
contamination inherited from that one seed.

Catching that by counting differences from a consensus does not work: its
backbone mismatch count was 9 against a seed-set median of 8, because a year
of real evolution produces just as many differences. What separates it is
*which* sites differ -- they form a known Alpha haplotype. Deciding that
needs a lineage-aware tool, so the primary screen runs Nextclade, applying
the same criteria as ncov's own `scripts/diagnostic.py` plus a check that the
called lineage sits inside the requested one. A seed that passes therefore
cannot found descendants that the build's diagnostics would exclude for
inherited reasons.

The lineage check rolls Nextclade's call up with the same pango_aliasor
machinery `prepare_clusters` uses to select the clusters, rather than a table
of clades per variant: one definition of "belongs to this lineage", and
nothing to extend when the next variant arrives.

`basic` mode is for when Nextclade cannot be reached. It screens on ambiguity
and terminal coverage, and flags consensus outliers. It would have caught the
Virginia seed on ambiguity alone, but that was luck rather than diagnosis,
so it is explicitly best-effort.
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, Iterable

import numpy as np

try:
    from .. import nextstrain_runtime as _runtime
except ImportError:                       # run as a plain script
    import nextstrain_runtime as _runtime

# ncov scripts/diagnostic.py defaults, restated so a seed that passes here
# cannot be excluded later for a reason it already carried.
MAX_CONTAMINATION = 5          # reversion + potential contaminant mutations
MAX_CLOCK_DEVIATION = 20
MAX_SNP_CLUSTERS = 1
# Screens that only apply to the seed itself.
MAX_AMBIGUOUS = 0.05           # fraction of non-ACGT bases
MAX_TERMINAL_MISSING = 0.50    # fraction of non-ACGT in either terminal 2kb
TERMINAL_WINDOW = 2000

# A called lineage the tool could not resolve. Never a rejection: the screen
# declines to judge rather than discarding a seed on the tool's uncertainty.
_UNCALLED = {"", "?", "na", "unassigned", "unknown", "recombinant"}


def lineage_matcher(pango_lineages, base_map_fn=None, recombinant: bool = False):
    """A predicate: is this called Pango lineage inside the requested one?

    Rolls the call up to its requested ancestor with the same pango_aliasor
    machinery `prepare_clusters` uses to choose the clusters, so "belongs to
    this lineage" means exactly what it meant upstream -- no table of clades
    to extend for each new variant, and recombinants and future lineages work
    for free.

    Returns a callable giving True (inside), False (outside) or None (cannot
    tell -- the tool did not call it, or no aliasor is available). Only False
    rejects a seed.
    """
    requested = {str(p).strip().rstrip("*") for p in pango_lineages if str(p).strip()}
    if not requested:
        return None

    base_map = {}
    if base_map_fn is not None:
        try:
            base_map = base_map_fn(sorted(requested), recombinant=recombinant) or {}
        except TypeError:                     # older signature without the flag
            try:
                base_map = base_map_fn(sorted(requested)) or {}
            except Exception:
                base_map = {}
        except Exception:                     # aliasor missing or unusable
            base_map = {}

    def inside(called) -> "bool | None":
        name = str(called or "").strip()
        if name.lower() in _UNCALLED:
            return None
        if name in requested:
            return True
        base = base_map.get(name)
        if base is not None:
            return base in requested
        if not base_map:
            return None                       # no aliasor: cannot decide
        # The aliasor partitions the lineage space it knows, so a call absent
        # from the map sits outside the requested lineages. Rejecting here is
        # the safe direction: a wrong reject costs one fallback to the
        # cluster's next sample, a wrong accept poisons a whole chain.
        return False

    return inside


# --------------------------------------------------------------------------
# basic screen -- no external tool
# --------------------------------------------------------------------------
def _as_array(seq: str) -> np.ndarray:
    return np.frombuffer(seq.upper().encode("ascii", "replace"), dtype="S1")


def ambiguity(seq: str) -> float:
    a = _as_array(seq)
    if a.size == 0:
        return 1.0
    return float(1.0 - np.isin(a, [b"A", b"C", b"G", b"T"]).mean())


def terminal_missing(seq: str, window: int = TERMINAL_WINDOW) -> float:
    """Worst ambiguity across the two terminal windows.

    A genome can look acceptable overall while an entire end is missing; the
    Virginia seed was 47% ambiguous across its 3' 2kb.
    """
    a = _as_array(seq)
    if a.size < 2 * window:
        return ambiguity(seq)
    acgt = np.isin(a, [b"A", b"C", b"G", b"T"])
    return float(max(1.0 - acgt[:window].mean(), 1.0 - acgt[-window:].mean()))


def basic_failures(records: Dict[str, str],
                   max_ambiguous: float = MAX_AMBIGUOUS,
                   max_terminal: float = MAX_TERMINAL_MISSING) -> Dict[str, str]:
    """{strain: reason} for sequences failing the tool-free screen."""
    out = {}
    for strain, seq in records.items():
        amb = ambiguity(seq)
        if amb > max_ambiguous:
            out[strain] = f"ambiguous bases {amb:.1%} > {max_ambiguous:.0%}"
            continue
        term = terminal_missing(seq)
        if term > max_terminal:
            out[strain] = f"terminal 2kb {term:.0%} missing > {max_terminal:.0%}"
    return out


def consensus_outliers(records: Dict[str, str], n_mad: float = 5.0,
                       floor: int = 10) -> Dict[str, int]:
    """Mismatch count at near-fixed sites, for sequences that look unusual.

    Advisory only. On the Virginia seed set the contaminated genome scored 9
    against a median of 8, so this cannot be a rejection criterion on its own;
    it is reported so a reviewer can see what the screen could not decide.
    """
    strains = [s for s in records if records[s]]
    if len(strains) < 50:
        return {}
    lengths = {len(records[s]) for s in strains}
    if len(lengths) != 1:                      # unaligned input: not comparable
        return {}
    arr = np.stack([_as_array(records[s]) for s in strains])
    acgt = np.isin(arr, [b"A", b"C", b"G", b"T"])
    clean = arr[(1.0 - acgt.mean(axis=1)) <= MAX_AMBIGUOUS]
    if len(clean) < 50:
        return {}
    cons = np.empty(arr.shape[1], dtype="S1")
    best = np.zeros(arr.shape[1])
    for base in (b"A", b"C", b"G", b"T"):
        frac = (clean == base).mean(axis=0)
        hit = frac > best
        cons[hit], best[hit] = base, frac[hit]
    fixed = best >= 0.95
    mism = ((arr != cons) & acgt & fixed).sum(axis=1)
    med = float(np.median(mism))
    mad = float(np.median(np.abs(mism - med)))
    thr = max(floor, med + n_mad * 1.4826 * mad)
    return {s: int(m) for s, m in zip(strains, mism) if m > thr}


# --------------------------------------------------------------------------
# nextclade screen
# --------------------------------------------------------------------------
def nextclade_command(ncov_dir: "Path | None" = None) -> "list | None":
    """The argv prefix that invokes nextclade directly, if there is one.

    `ncov_dir` is unused and kept for callers. Returns None when the tool has
    no reachable path of its own, which is the normal case for the docker and
    singularity runtimes; `nextclade_failures` then goes through the CLI
    instead (see nextstrain_runtime.run).
    """
    cmd, _source = _runtime.tool_command("nextclade")
    return cmd


def _dataset(ncov_dir: "Path | None") -> "Path | None":
    if not ncov_dir:
        return None
    zipped = Path(ncov_dir) / "data" / "sars-cov-2-nextclade-defaults.zip"
    return zipped if zipped.is_file() else None


def _rows(tsv: Path) -> Iterable[dict]:
    with open(tsv, newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def _num(row: dict, *names: str) -> float:
    for n in names:
        v = row.get(n, "")
        if v not in ("", None):
            try:
                return float(v)
            except ValueError:
                pass
    return 0.0


def nextclade_failures(records: Dict[str, str], pango_lineages: Iterable[str],
                       ncov_dir: "Path | None" = None,
                       max_contamination: int = MAX_CONTAMINATION,
                       workdir: "Path | None" = None,
                       base_map_fn=None) -> "tuple[Dict[str, str], str]":
    """{strain: reason} for sequences Nextclade rejects, plus a status line.

    Criteria mirror ncov's scripts/diagnostic.py, with the lineage check
    added. A missing tool or dataset is not an error: the caller falls back.
    """
    dataset = _dataset(ncov_dir)
    if dataset is None:
        return {}, (f"no Nextclade dataset at {ncov_dir}/data/"
                    f"sars-cov-2-nextclade-defaults.zip (run `phylogas ncov-checkout`)")

    inside = lineage_matcher(pango_lineages, base_map_fn=base_map_fn)
    requested = ", ".join(sorted({str(p).strip().rstrip("*") for p in pango_lineages}))

    # Scratch goes inside the checkout rather than /tmp: the docker and
    # singularity runtimes mount only the directory handed to the CLI, so a
    # system temp directory would be invisible to nextclade there.
    tmp = Path(workdir) if workdir else Path(ncov_dir) / "seed_qc"
    tmp.mkdir(parents=True, exist_ok=True)
    fasta, out_tsv = tmp / "candidates.fasta", tmp / "nextclade.tsv"
    if out_tsv.exists():
        out_tsv.unlink()
    with open(fasta, "w") as fh:
        for strain, seq in records.items():
            fh.write(f">{strain}\n{seq}\n")

    args = ["run", "--input-dataset", str(dataset),
            "--output-tsv", str(out_tsv), "--jobs", "4", str(fasta)]
    proc, source = _runtime.run("nextclade", args, ncov_dir=ncov_dir)
    if proc is None:
        return {}, source
    if proc.returncode != 0 or not out_tsv.is_file():
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-3:]
        return {}, f"nextclade failed via {source}: " + " / ".join(tail)

    failures, seen = {}, 0
    for row in _rows(out_tsv):
        strain = row.get("seqName", "").split()[0] if row.get("seqName") else ""
        if not strain:
            continue
        seen += 1
        rev = _num(row, "privateNucMutations.reversionSubstitutions",
                   "reversion_mutations")
        con = _num(row, "privateNucMutations.labeledSubstitutions",
                   "potential_contaminants")
        if rev + con > max_contamination:
            failures[strain] = (f"reversions+contaminants {rev + con:.0f} > "
                                f"{max_contamination}")
            continue
        if abs(_num(row, "clock_deviation", "totalMissing")) > MAX_CLOCK_DEVIATION \
                and "clock_deviation" in row:
            failures[strain] = "clock deviation > %d" % MAX_CLOCK_DEVIATION
            continue
        if _num(row, "qc.snpClusters.totalSNPs", "snp_clusters") > MAX_SNP_CLUSTERS:
            failures[strain] = "SNP clusters > %d" % MAX_SNP_CLUSTERS
            continue
        if str(row.get("qc.mixedSites.status", "")).lower() == "bad":
            failures[strain] = "mixed sites: bad"
            continue
        if inside is not None:
            called = next((row[c] for c in ("Nextclade_pango", "lineage", "pango_lineage")
                           if row.get(c)), "")
            if inside(called) is False:
                failures[strain] = f"lineage {called} is outside {requested}"
    status = f"nextclade screened {seen:,} candidate(s) via {source}"
    return failures, status


# --------------------------------------------------------------------------
# entry point used by the fetch loop
# --------------------------------------------------------------------------
def screen(records: Dict[str, str], pango_lineages: Iterable[str], mode: str = "nextclade",
           ncov_dir: "Path | None" = None, max_ambiguous: float = MAX_AMBIGUOUS,
           max_contamination: int = MAX_CONTAMINATION,
           base_map_fn=None) -> "tuple[Dict[str, str], str]":
    """Return ({strain: reason} for rejects, a one-line status).

    `mode` is "nextclade" (falls back to "basic" with a warning when the tool
    cannot be reached), "basic", or "off".
    """
    if mode == "off" or not records:
        return {}, "seed QC disabled" if mode == "off" else "nothing to screen"

    if mode == "nextclade":
        failures, status = nextclade_failures(
            records, pango_lineages, ncov_dir=ncov_dir,
            max_contamination=max_contamination, base_map_fn=base_map_fn)
        if not status.startswith("nextclade screened"):
            basic = basic_failures(records, max_ambiguous=max_ambiguous)
            return basic, (f"WARNING: {status}; fell back to the basic screen, which "
                           f"cannot detect mixed-lineage assemblies")
        # The ambiguity and terminal-coverage screens still apply: they are
        # about whether the genome is usable, not about which clade it is.
        failures.update({s: r for s, r in
                         basic_failures(records, max_ambiguous=max_ambiguous).items()
                         if s not in failures})
        return failures, status

    return basic_failures(records, max_ambiguous=max_ambiguous), (
        "basic screen (ambiguity and terminal coverage only; cannot detect "
        "mixed-lineage assemblies)")


def write_report(path: Path, rows: Iterable[dict]) -> None:
    """One row per candidate considered, in the order they were tried."""
    rows = list(rows)
    if not rows:
        return
    cols = ["cluster_id", "attempt", "strain", "verdict", "reason", "intro_date", "order"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
