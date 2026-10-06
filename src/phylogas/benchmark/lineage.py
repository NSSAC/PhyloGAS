"""Decompose mu_truth into its within-chain and between-chain parts.

`clock.root_to_tip` regresses divergence-from-Wuhan on calendar date, pooled
across every tip. That is the right estimator for comparing against TreeTime,
which fits the same shape on the same tree root -- but it is not the painter's
per-lineage rate, for two reasons that this module separates.

The painter does not start from the reference. Every chain begins at a real
seed genome that already carries its own divergence from Wuhan (measured at
mean 41.6, sd 4.7 substitutions over 3,322 seeds for the Virginia Delta run),
and mutations accumulate along time since *that chain's* importation rather
than along calendar date. So the pooled slope mixes:

    mu_truth  =  w * mu_lineage  +  (1 - w) * mu_between_chains

where w is the share of date variance lying within chains. mu_lineage is the
fixed-effects within-chain slope -- one common rate, a free intercept per
chain -- and is the painter's own rate. mu_between_chains is the slope of
chain means on chain mean dates, tip-weighted, and is a property of the
importation regime rather than of evolution.

With importation continuing through a wave, every late date holds a mixture of
old and young chains, so the between-chain term is shallow and drags the
pooled slope down. On synthetic data with a known 8e-4 clock the pooled fit
came back at 2.5e-05, a 32x attenuation, while the within-chain fit recovered
the truth.

Everything streams. The painted set runs to millions of tips across hundreds
of thousands of chains, so each chain costs six running sums rather than a
list of points.
"""
from __future__ import annotations

import datetime as _dt

from .clock import (
    _codes,
    _iter_fasta,
    _open_text,
    _valid_mask,
    comparable_sites,
    hamming_to_reference,
    per_site_per_year,
    read_dates,
    read_reference,
)

# Values that mean "absent" in these tables. component_id uses -1 for an
# agent with no exposure event; real_strain is blank for non-index cases.
_EMPTY = {"", "none", "nan", "na", "null", "-1"}


class _Sums:
    """Running least-squares sums for one group."""
    __slots__ = ("n", "sx", "sy", "sxx", "sxy", "syy", "xmin", "ymin")

    def __init__(self):
        self.n = 0
        self.sx = self.sy = self.sxx = self.sxy = self.syy = 0.0
        self.xmin = None
        self.ymin = None

    def add(self, x, y):
        self.n += 1
        self.sx += x
        self.sy += y
        self.sxx += x * x
        self.sxy += x * y
        self.syy += y * y
        if self.xmin is None or x < self.xmin:
            self.xmin = x
        if self.ymin is None or y < self.ymin:
            self.ymin = y

    def centred(self):
        """(Sxx, Sxy, Syy) about this group's own means."""
        if self.n == 0:
            return 0.0, 0.0, 0.0
        mx, my = self.sx / self.n, self.sy / self.n
        return (self.sxx - self.n * mx * mx,
                self.sxy - self.n * mx * my,
                self.syy - self.n * my * my)


def _fit(sxx, sxy, syy):
    """Slope and R^2 from centred sums."""
    if sxx <= 0:
        return None, None
    slope = sxy / sxx
    r2 = (sxy * sxy) / (sxx * syy) if syy > 0 else float("nan")
    return slope, r2


def _weighted_fit(points):
    """[(x, y, w)] -> (slope, r2) by weighted least squares."""
    W = sum(w for _, _, w in points)
    if W <= 0 or len(points) < 3:
        return None, None
    mx = sum(x * w for x, _, w in points) / W
    my = sum(y * w for _, y, w in points) / W
    sxx = sum(w * (x - mx) ** 2 for x, _, w in points)
    sxy = sum(w * (x - mx) * (y - my) for x, y, w in points)
    syy = sum(w * (y - my) ** 2 for _, y, w in points)
    return _fit(sxx, sxy, syy)


def read_components(path, key_col="alias_pid", comp_col="component_id") -> dict:
    """Infection id -> chain id, from TwinSampler's all-events table.

    The painted metadata carries no chain column: the painter writes
    `real_strain` only for the index case, so a chain cannot be recovered from
    it. The all-events table has both `alias_pid` and `component_id`, which is
    the join.

    This is the memory high-water mark -- one entry per infection, so roughly
    a gigabyte for 5.3M. Repeated rows for one infection collapse (last wins),
    which is correct because component_id is assigned per transmission
    component, not per clinical state.
    """
    import csv

    out = {}
    with _open_text(path) as fh:
        first = fh.readline().rstrip("\r\n")
        delim = "\t" if first.count("\t") >= first.count(",") else ","
        cols = first.split(delim)
        for want in (key_col, comp_col):
            if want not in cols:
                raise ValueError(
                    f"{path} has no {want!r} column; has: "
                    f"{', '.join(cols[:12])}")
        ik, ic = cols.index(key_col), cols.index(comp_col)
        for row in csv.reader(fh, delimiter=delim):
            if len(row) > max(ik, ic):
                k, c = row[ik].strip(), row[ic].strip()
                if k and c.lower() not in _EMPTY:
                    out[k] = c
    return out


def read_meta_keys(metadata, id_col="strain", key_col="alias_pid",
                   seed_col="real_strain") -> tuple:
    """strain -> (infection id, seed id or None).

    `real_strain` names the seed this record was painted from, and is set only
    on index cases -- the painter passes a seed id to populate_sim_details
    only when it is seeding. That is what identifies the importations.
    """
    out = {}
    missing = []
    with _open_text(metadata) as fh:
        first = fh.readline().rstrip("\n").rstrip("\r")
        delim = "\t" if first.count("\t") >= first.count(",") else ","
        cols = first.split(delim)
        if id_col not in cols:
            return {}, [id_col]
        if key_col not in cols:
            missing.append(key_col)
        i_id = cols.index(id_col)
        i_key = cols.index(key_col) if key_col in cols else None
        i_seed = cols.index(seed_col) if seed_col in cols else None
        if i_seed is None:
            missing.append(seed_col)
        top = max(x for x in (i_id, i_key, i_seed) if x is not None)
        for line in fh:
            parts = line.rstrip("\n").rstrip("\r").split(delim)
            if len(parts) <= top:
                continue
            key = parts[i_key].strip() if i_key is not None else ""
            seed = parts[i_seed].strip() if i_seed is not None else ""
            out[parts[i_id].strip()] = (
                key or None,
                None if seed.lower() in _EMPTY else seed)
    return out, missing


def decompose(fasta, reference, metadata, components=None,
              date_col="date", id_col="strain", key_col="alias_pid",
              seed_col="real_strain", comp_col="component_id",
              date_min=None, date_max=None, exclude_ids=None,
              takeoff_min=10, max_records=None) -> list:
    """Split the pooled rate into within-chain and between-chain parts.

    Returns a list of result dicts shaped like `clock.root_to_tip`'s, one per
    quantity, so they go into clock_estimates.csv as ordinary rows.

    `components` is the mapping from `read_components`; without it only the
    seed-trend rows can be produced, since chain membership is unknown.
    """
    ref = read_reference(reference)
    ref_codes = _codes(ref)
    ref_valid = _valid_mask(ref_codes)
    n_sites = comparable_sites(ref)

    dates = read_dates(metadata, date_col=date_col, id_col=id_col)
    if not dates:
        return [{"quantity": "mu_lineage",
                 "error": f"no usable dates in {metadata}"}]
    keys, missing = read_meta_keys(metadata, id_col, key_col, seed_col)
    if not keys:
        return [{"quantity": "mu_lineage",
                 "error": f"{metadata} has no {id_col!r} column"}]

    drop = {str(i) for i in (exclude_ids or ())}
    chains: dict = {}
    seeds: list = []
    n_seen = n_used = n_nochain = n_unaligned = 0

    for sid, seq in _iter_fasta(fasta):
        n_seen += 1
        if max_records and n_used >= max_records:
            break
        if sid in drop:
            continue
        d = dates.get(sid)
        if d is None:
            continue
        if (date_min is not None and d < date_min) or \
           (date_max is not None and d > date_max):
            continue
        div = hamming_to_reference(seq, ref, ref_codes, ref_valid)
        if div is None:
            n_unaligned += 1
            continue
        x, y = float(d.toordinal()), float(div)
        n_used += 1

        infection, seed = keys.get(sid, (None, None))
        if seed is not None:
            # An index case: its divergence IS its seed's, and its date is
            # the importation date.
            seeds.append((x, y, components.get(infection)
                          if (components and infection) else None))
        if components is not None and infection is not None:
            chain = components.get(infection)
            if chain is None:
                n_nochain += 1
            else:
                chains.setdefault(chain, _Sums()).add(x, y)
        elif components is not None:
            n_nochain += 1

    diag = {
        "tips": n_used,
        "chains": len(chains),
        "skipped_unaligned": n_unaligned,
        "skipped_no_chain": n_nochain,
        "window": (f"{date_min or ''}..{date_max or ''}"
                   if (date_min or date_max) else "unbounded"),
        "genome_length": len(ref),
        "comparable_sites": n_sites,
    }
    rows = []
    rows.extend(_chain_rows(chains, n_sites, diag))
    rows.extend(_seed_rows(seeds, chains, n_sites, takeoff_min, diag))
    return rows


def _chain_rows(chains, n_sites, diag) -> list:
    """mu_lineage and mu_between_chains from the per-chain sums."""
    groups = [g for g in chains.values() if g.n > 0]
    N = sum(g.n for g in groups)
    if N < 3 or not groups:
        return [{"quantity": "mu_lineage", **diag,
                 "error": "no chain membership; pass the all-events table "
                          "so component_id can be joined on alias_pid"}]

    # Within: sum each chain's own centred sums. This is the fixed-effects
    # slope -- one shared rate, a free intercept per chain -- so each chain's
    # seed offset is absorbed rather than contributing to the fit.
    sxx_w = sxy_w = syy_w = 0.0
    for g in groups:
        a, b, c = g.centred()
        sxx_w += a
        sxy_w += b
        syy_w += c
    slope_w, r2_w = _fit(sxx_w, sxy_w, syy_w)

    # Between: chain means against chain mean dates, weighted by tips. The
    # weighting is what makes this the term that actually enters the pooled
    # slope -- a chain that never took off contributes almost nothing.
    means = [(g.sx / g.n, g.sy / g.n, float(g.n)) for g in groups]
    slope_b, r2_b = _weighted_fit(means)

    # Variance share, which is the weight in
    #   mu_truth = w * mu_lineage + (1 - w) * mu_between_chains
    X = sum(g.sx for g in groups) / N
    sxx_b = sum(g.n * (g.sx / g.n - X) ** 2 for g in groups)
    total = sxx_w + sxx_b
    w = (sxx_w / total) if total > 0 else float("nan")

    fitted = len(groups)
    sizes = sorted((g.n for g in groups), reverse=True)
    shared = {**diag, "chains_fitted": fitted,
              "largest_chain_tips": sizes[0] if sizes else 0,
              "median_chain_tips": sizes[len(sizes) // 2] if sizes else 0,
              "within_date_variance_share": w}

    rows = []
    if slope_w is not None:
        rows.append({
            "quantity": "mu_lineage",
            "method": ("within-chain slope, free intercept per chain "
                       "(the painter's per-lineage rate)"),
            "slope_subs_per_genome_per_day": slope_w,
            "rate_subs_per_site_per_year": per_site_per_year(slope_w, n_sites),
            "r_squared": r2_w,
            **shared,
        })
    if slope_b is not None:
        rows.append({
            "quantity": "mu_between_chains",
            "method": ("chain means vs chain mean dates, tip-weighted "
                       "(the importation regime, not evolution)"),
            "slope_subs_per_genome_per_day": slope_b,
            "rate_subs_per_site_per_year": per_site_per_year(slope_b, n_sites),
            "r_squared": r2_b,
            **shared,
        })
    return rows


def _seed_rows(seeds, chains, n_sites, takeoff_min, diag) -> list:
    """The seed trend: what was introduced, and what the founder lottery kept.

    Three fits over the same index cases, differing only in which seeds count
    and how much:

    mu_seed_trend             every importation, unweighted. The real-world
                              trend in what arrived, blind to whether it grew.
    mu_seed_trend_surviving   only chains reaching `takeoff_min` tips.
    mu_seed_trend_weighted    those, weighted by tip count -- the version that
                              actually enters the pooled slope.

    The gap between the first and the last is the founder lottery's effect on
    the trend, which is what separates this simulation's pooled slope from
    what real surveillance of the same wave would have seen. Real surveillance
    applies the same filter (you only sequence lineages that grew) but draws
    its own lottery, so the two need not agree even with identical seeds.
    """
    if len(seeds) < 3:
        return [{"quantity": "mu_seed_trend", **diag,
                 "error": (f"only {len(seeds)} index cases carried "
                           f"real_strain; the painted metadata sets it on "
                           f"seeded records only")}]

    shared = {**diag, "imports": len(seeds), "takeoff_min": takeoff_min}
    rows = []

    allpts = [(x, y, 1.0) for x, y, _ in seeds]
    slope, r2 = _weighted_fit(allpts)
    if slope is not None:
        rows.append({
            "quantity": "mu_seed_trend",
            "method": ("seed divergence vs import date, every importation, "
                       "unweighted (what arrived)"),
            "slope_subs_per_genome_per_day": slope,
            "rate_subs_per_site_per_year": per_site_per_year(slope, n_sites),
            "r_squared": r2,
            "mean_divergence": sum(y for _, y, _ in seeds) / len(seeds),
            **shared,
        })

    if not chains:
        return rows

    grew = [(x, y, c) for x, y, c in seeds
            if c is not None and c in chains and chains[c].n >= takeoff_min]
    if len(grew) >= 3:
        s2, r2b = _weighted_fit([(x, y, 1.0) for x, y, _ in grew])
        if s2 is not None:
            rows.append({
                "quantity": "mu_seed_trend_surviving",
                "method": (f"same, but only chains reaching {takeoff_min} "
                           f"tips (what took off)"),
                "slope_subs_per_genome_per_day": s2,
                "rate_subs_per_site_per_year": per_site_per_year(s2, n_sites),
                "r_squared": r2b,
                "surviving_imports": len(grew),
                **shared,
            })
        s3, r2c = _weighted_fit([(x, y, float(chains[c].n))
                                 for x, y, c in grew])
        if s3 is not None:
            rows.append({
                "quantity": "mu_seed_trend_weighted",
                "method": ("same, weighted by tips per chain (the term that "
                           "enters the pooled slope)"),
                "slope_subs_per_genome_per_day": s3,
                "rate_subs_per_site_per_year": per_site_per_year(s3, n_sites),
                "r_squared": r2c,
                "surviving_imports": len(grew),
                **shared,
            })
    return rows
