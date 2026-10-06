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


# A, C, G, T as ASCII codes. Membership tested by four comparisons rather
# than np.isin, which is slower for a set this small.
_A, _C, _G, _T = (ord(c) for c in "ACGT")


def _codes(s: str):
    """Uppercased ASCII codes for a sequence, as a numpy array.

    bytes.upper() rather than str.upper() so the case fold happens in C, and
    frombuffer rather than fromiter so no per-character Python object is
    built. Bit tricks on 0x20 would be faster still but corrupt '-' (0x2D),
    which has to survive as an invalid character rather than become one.
    """
    import numpy as np

    return np.frombuffer(s.encode("ascii", "replace").upper(), dtype=np.uint8)


def _valid_mask(arr):
    return ((arr == _A) | (arr == _C) | (arr == _G) | (arr == _T))


def hamming_to_reference(seq: str, ref, ref_codes=None,
                         ref_valid=None) -> "int | None":
    """Observable differences from the reference.

    Ambiguity and gaps are skipped rather than counted: a real root-to-tip
    plot cannot see through an N either, and counting them would inflate
    divergence in proportion to sequence quality rather than to time.

    Returns None if the lengths differ, since that means the record is not
    aligned to the reference and silently comparing it would be wrong.

    Vectorised. The previous version summed a generator over ~30k characters
    with a set lookup each, which cost about 0.8 ms per record -- 70 minutes
    for a 5.35M-record painted set, and the reason `max_records` existed at
    all. `ref_codes`/`ref_valid` let a caller hoist the reference's own arrays
    out of the loop; without them they are recomputed per call.
    """
    import numpy as np

    if len(seq) != len(ref):
        return None
    a = _codes(seq)
    b = ref_codes if ref_codes is not None else _codes(ref)
    bv = ref_valid if ref_valid is not None else _valid_mask(b)
    return int(np.count_nonzero((a != b) & _valid_mask(a) & bv))


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
                max_records=None, quantity="mu_truth",
                date_min=None, date_max=None, exclude_ids=None,
                reservoir_seed=0) -> dict:
    """Regress divergence-from-reference on sampling date.

    The estimator does not care whether the sequences are simulated or real:
    given an alignment and dates it returns the same quantity, which is what
    makes it the most directly comparable number available -- no TreeTime on
    either side. The caller names it (`mu_truth` for painted sequences,
    `mu_real_observed` for real ones).

    Real sequences must be ALIGNED to the reference first; unaligned records
    are skipped and counted, since silently comparing them would be wrong.

    `date_min`/`date_max` bound the fitted window. They matter more than they
    look: a single record dated outside the simulated period -- the reference
    at Wuhan-Hu-1's real 2019-12-26 collection date, or a seed carrying its
    real-world date -- stretches `date_span_days` from the ~70 days actually
    studied to 537, which turns the one field you would glance at to judge
    the fit into a number describing data the fit barely contains. Excluded
    records are counted, not silently dropped.

    `exclude_ids` drops records by id, for the reference and any other
    context sequence injected into the alignment.

    Returns a dict of regression results, or one carrying `error` when there
    is not enough usable data to fit.
    """
    ref = read_reference(reference)
    dates = read_dates(metadata, date_col=date_col, id_col=id_col)
    if not dates:
        return {"error": f"no usable dates in {metadata}"}
    drop = {str(i) for i in (exclude_ids or ())}

    # Hoisted: the reference's codes and validity mask do not change per
    # record, and recomputing them was a third of the per-record work.
    _ref_codes = _codes(ref)
    _ref_valid = _valid_mask(_ref_codes)

    # max_records used to `break` after the first N usable records. The
    # painted FASTA is written in tick order, so that returned the earliest
    # infections and a compressed date range -- a 69-day span for a 299-day
    # study, which is how the same data gave 8.5e-04 capped and 2.8e-04
    # whole. Reservoir sampling instead: every usable record gets an equal
    # chance of being kept, so the retained sample spans the full window and
    # the slope is unbiased. Seeded, so a rerun reproduces.
    import random as _random

    reservoir = _random.Random(reservoir_seed)
    n_usable = 0

    xs, ys = [], []
    n_seen = n_nodate = n_unaligned = n_window = n_excluded = 0
    for sid, seq in _iter_fasta(fasta):
        n_seen += 1
        if sid in drop:
            n_excluded += 1
            continue
        d = dates.get(sid)
        if d is None:
            n_nodate += 1
            continue
        if (date_min is not None and d < date_min) or \
           (date_max is not None and d > date_max):
            n_window += 1
            continue
        div = hamming_to_reference(seq, ref, _ref_codes, _ref_valid)
        if div is None:
            n_unaligned += 1
            continue
        n_usable += 1
        if not max_records or len(xs) < max_records:
            xs.append(d.toordinal())
            ys.append(float(div))
        else:
            j = reservoir.randrange(n_usable)
            if j < max_records:
                xs[j] = d.toordinal()
                ys[j] = float(div)

    if len(xs) < 3:
        msg = (f"only {len(xs)} usable records ({n_nodate} without a date, "
               f"{n_unaligned} not aligned to the reference, "
               f"{n_window} outside {date_min or '-inf'}..{date_max or '+inf'}, "
               f"{n_excluded} excluded by id, of {n_seen} sequences)")
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
        "skipped_out_of_window": n_window,
        "skipped_excluded_id": n_excluded,
        "window": (f"{date_min or ''}..{date_max or ''}"
                   if (date_min or date_max) else "unbounded"),
        "first_date": _dt.date.fromordinal(min(xs)).isoformat(),
        "last_date": _dt.date.fromordinal(max(xs)).isoformat(),
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
