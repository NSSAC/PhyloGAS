"""Molecular clock rate: measured, inferred, and assumed.

Three quantities, deliberately kept apart. `docs/clock_modes.md` is the
decision record for which answers which question; this module computes them.

mu_truth
    The simulation's own root-to-tip rate, from the painted sequences. A
    regression of Hamming distance from the reference on sampling date -- no
    tree and no phylogenetic inference, so it carries no inference error.

mu_sim / mu_real
    TreeTime's estimate, from `augur refine` with no --clock-rate, over ncov's
    own intermediates. The same procedure on either kind of data.

mu_operational
    The rate ncov assumed (its 0.0008 default), read back out of the build's
    branch_lengths.json. Recorded so a tree is never mistaken for a
    measurement.

Why root-to-tip rather than counting the substitutions the painter placed:
root-to-tip is what TempEst and TreeTime regress and what the ~8e-4
subs/site/year literature value derives from, so it is the only form
comparable to real data. Counting placed substitutions measures events,
including back-mutations and repeat hits that no real root-to-tip plot can
see. Both are computed -- their ratio is a saturation diagnostic -- but only
root-to-tip is comparable.
"""

from __future__ import annotations

import datetime as _dt
import gzip
import lzma
import math
from pathlib import Path

GENOME_LENGTH_FALLBACK = 29903


def _open_text(path):
    p = str(path)
    if p.endswith(".xz"):
        return lzma.open(p, "rt")
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    return open(p, "r")


def _iter_fasta(path):
    """Yield (id, sequence) without holding the file in memory.

    The painted FASTA runs to millions of records, so this streams.
    """
    name, chunks = None, []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].strip().split(None, 1)[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def read_reference(path) -> str:
    for _name, seq in _iter_fasta(path):
        return seq.upper()
    raise ValueError(f"no sequence found in {path}")


def hamming_to_reference(seq: str, ref: str) -> "int | None":
    """Observable differences from the reference.

    Ambiguity and gaps are skipped rather than counted: a real root-to-tip
    plot cannot see through an N either, and counting them would inflate
    divergence in proportion to sequence quality rather than to time.

    Returns None if the lengths differ, since that means the record is not
    aligned to the reference and silently comparing it would be wrong.
    """
    if len(seq) != len(ref):
        return None
    valid = set("ACGT")
    return sum(1 for a, b in zip(seq.upper(), ref) if a != b and a in valid and b in valid)


def read_dates(metadata, date_col="date", id_col="strain") -> dict:
    """Map sequence id -> date, from a metadata table."""
    out = {}
    with _open_text(metadata) as fh:
        first = fh.readline().rstrip("\n").rstrip("\r")
        delim = "\t" if first.count("\t") >= first.count(",") else ","
        cols = first.split(delim)
        if id_col not in cols or date_col not in cols:
            raise ValueError(
                f"{metadata} needs '{id_col}' and '{date_col}' columns; has: "
                f"{', '.join(cols[:12])}{' ...' if len(cols) > 12 else ''}")
        i_id, i_date = cols.index(id_col), cols.index(date_col)
        for line in fh:
            parts = line.rstrip("\n").rstrip("\r").split(delim)
            if len(parts) <= max(i_id, i_date):
                continue
            sid, raw = parts[i_id].strip(), parts[i_date].strip()[:10]
            if not sid or not raw:
                continue
            try:
                out[sid] = _dt.date.fromisoformat(raw)
            except ValueError:
                continue          # ambiguous or partial dates: unusable here
    return out


def linregress(xs, ys):
    """Least-squares slope, intercept and R^2. Avoids a scipy dependency."""
    n = len(xs)
    if n < 3:
        return None
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = my - slope * mx
    syy = sum((y - my) ** 2 for y in ys)
    r2 = (sxy ** 2) / (sxx * syy) if syy > 0 else float("nan")
    return slope, intercept, r2


def comparable_sites(ref: str) -> int:
    """Sites where divergence can actually be observed.

    The reference is terminus-masked to match ncov's copy (first 100 and last
    50 bases are N), and `hamming_to_reference` skips any position that is not
    an unambiguous base on both sides. So the regression measures divergence
    over the unmasked sites, and the per-site rate has to divide by that count
    rather than by the raw genome length -- otherwise masking silently
    deflates the rate by the masked fraction.
    """
    return sum(1 for b in ref.upper() if b in "ACGT")


def per_site_per_year(slope_per_genome_per_day: float, n_sites: int) -> float:
    """Convert the regression slope to the units the literature uses.

    `n_sites` is the number of comparable sites, not the genome length. See
    comparable_sites().
    """
    return slope_per_genome_per_day * 365.0 / n_sites


def root_to_tip(fasta, reference, metadata, date_col="date", id_col="strain",
                max_records=None, quantity="mu_truth") -> dict:
    """Regress divergence-from-reference on sampling date.

    The estimator does not care whether the sequences are simulated or real:
    given an alignment and dates it returns the same quantity, which is what
    makes it the most directly comparable number available -- no TreeTime on
    either side. The caller names it (`mu_truth` for painted sequences,
    `mu_real_observed` for real ones).

    Real sequences must be ALIGNED to the reference first; unaligned records
    are skipped and counted, since silently comparing them would be wrong.

    Returns a dict of regression results, or one carrying `error` when there
    is not enough usable data to fit.
    """
    ref = read_reference(reference)
    dates = read_dates(metadata, date_col=date_col, id_col=id_col)
    if not dates:
        return {"error": f"no usable dates in {metadata}"}

    xs, ys = [], []
    n_seen = n_nodate = n_unaligned = 0
    for sid, seq in _iter_fasta(fasta):
        n_seen += 1
        d = dates.get(sid)
        if d is None:
            n_nodate += 1
            continue
        div = hamming_to_reference(seq, ref)
        if div is None:
            n_unaligned += 1
            continue
        xs.append(d.toordinal())
        ys.append(float(div))
        if max_records and len(xs) >= max_records:
            break

    if len(xs) < 3:
        msg = (f"only {len(xs)} usable records ({n_nodate} without a date, "
               f"{n_unaligned} not aligned to the reference, of {n_seen} sequences)")
        if n_unaligned > len(xs):
            msg += ("\n       Most records are not the reference's length. Real "
                    "sequences need aligning\n       first -- use ncov's "
                    "results/<build>/filtered.fasta, or `augur align`.")
        return {"error": msg}

    fit = linregress(xs, ys)
    if fit is None:
        return {"error": "all sampling dates identical; no temporal signal"}
    slope, intercept, r2 = fit
    span = (max(xs) - min(xs))
    n_sites = comparable_sites(ref)
    return {
        "quantity": quantity,
        "method": "root-to-tip vs date (Hamming from reference, tree-free)",
        "tips": len(xs),
        "slope_subs_per_genome_per_day": slope,
        "rate_subs_per_site_per_year": per_site_per_year(slope, n_sites),
        "r_squared": r2,
        "date_span_days": span,
        "mean_divergence": sum(ys) / len(ys),
        "genome_length": len(ref),
        "comparable_sites": n_sites,
        "skipped_no_date": n_nodate,
        "skipped_unaligned": n_unaligned,
    }


def operational_rate(branch_lengths_json) -> dict:
    """mu_operational: the rate ncov assumed, read back from its own output."""
    import json

    with _open_text(branch_lengths_json) as fh:
        data = json.load(fh)
    clock = data.get("clock", {}) or {}
    rate = clock.get("rate")
    if rate is None:
        return {"error": f"{branch_lengths_json} has no clock.rate "
                         f"(keys: {', '.join(list(data)[:8])})"}
    return {
        "quantity": "mu_operational",
        "method": "assumed by ncov (refine.clock_rate), not measured",
        "rate_subs_per_site_per_year": float(rate),
        "r_squared": clock.get("r_val") ** 2 if clock.get("r_val") is not None else None,
        "intercept": clock.get("intercept"),
    }


def format_row(arm: str, res: dict) -> dict:
    """Flatten one result into a CSV row."""
    row = {"arm": arm}
    row.update(res)
    return row
