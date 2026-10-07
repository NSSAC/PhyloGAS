#!/usr/bin/env python3
"""Importation schedules and their seed sequences, from one cluster table.

The schedule (how many importations land on each day) and the seed genomes
(which real sequence founds each one) used to be produced by two different
code paths: `src/importation_analysis.py` wrote the schedules EpiHiper ran,
while `seq_prep.py --seed_mode` picked the seeds. Nothing tied them together,
and they disagreed in a way that mattered:

  * The schedule places a cluster's importation on its `earliest_date`, the
    collection date of its FIRST sample.
  * Seed picking took `samples.split(",")[0]`, but that list is ordered
    NEWEST FIRST. For 79% of Virginia's Delta seeds the chosen sequence was
    the cluster's LAST sample -- in the worst case 209 days after the
    importation it was supposed to found.

The painter hands out seed records in file order as importations come up, so
early importations received late genomes and vice versa (rank correlation
between importation date and genome collection date: -0.54). That inflates
root-to-tip divergence at the start of the study and is the likely cause of
the negative `mu_seed_trend` the clock benchmark reports.

This module derives both products from the same rows, so they cannot drift:
the seed of a cluster is its EARLIEST sample, and the FASTA is written in
schedule order.

Verified: `build_schedule` reproduces all five
`data/importations/schedules/<State>_schedule.csv` files exactly (rows and
per-day counts), and `absolute_ticks` reproduces all ten
`data/seed_schedule/<st>-<n>.csv` files EpiHiper actually consumed -- with
`outlier_method="none"`, which is what produced them.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

# A cluster's `earliest_date` can predate the variant's real emergence, from a
# misassigned or mis-dated sample. The schedule generator replaced such a date
# with the cluster's first sample after the threshold; reproducing the
# published schedules requires the same correction with the same dates.
VARIANT_START_DEFAULTS: Dict[str, str] = {
    "B.1.1.7": "2020-12-01",
    "B.1.617.2": "2021-03-01",
}

# EpiHiper's own tick 0, from cfg/exp1/config.json ("tickZero"). Shared across
# states deliberately, so every state's ticks are on one calendar.
EPIHIPER_TICK_ZERO = "2020-11-30"


def tick_zero_from_abm_config(path) -> pd.Timestamp:
    """Read tick 0 from an EpiHiper config JSON (its `tickZero` key).

    The ABM config stays in the format EpiHiper expects; PhyloGAS reads it
    rather than keeping a second copy of the date.
    """
    data = json.loads(Path(path).read_text())
    for key in ("tickZero", "tick_zero", "projectionBegDate"):
        if key in data:
            return pd.to_datetime(data[key])
    raise KeyError(f"no tickZero in {path}")


def parse_cluster_dates(values: pd.Series) -> pd.Series:
    """Parse `earliest_date` in either format the cluster tables use.

    UCSC's table writes `2021-May-12`; the filtered snapshot was re-saved as
    `2021-05-12`; some rows say `no-valid-date`. Left to infer, pandas picks
    one format from the first value and coerces everything else to NaT --
    which silently dropped 99% of Washington's Delta clusters while Virginia,
    whose first row happened to be unparseable, came through fine. So each
    format is tried explicitly.
    """
    v = values.astype(str).str.strip()
    out = pd.to_datetime(v, format="%Y-%b-%d", errors="coerce")
    miss = out.isna()
    if miss.any():
        out[miss] = pd.to_datetime(v[miss], format="%Y-%m-%d", errors="coerce")
    return out


def parse_samples(samples: str) -> List[Tuple[Optional[pd.Timestamp], str, str]]:
    """Split a cluster's `samples` field into (date, strain, accession), oldest first.

    Entries are `strain|accession|date`. Rows whose date will not parse sort
    last rather than being dropped: they are still usable as fallback seeds.
    """
    out = []
    for item in str(samples).split(","):
        item = item.strip()
        if not item:
            continue
        parts = item.split("|")
        strain = parts[0]
        acc = parts[1] if len(parts) > 1 else ""
        date = pd.to_datetime(parts[2], errors="coerce") if len(parts) > 2 else pd.NaT
        out.append((date, strain, acc))
    out.sort(key=lambda t: (pd.isna(t[0]), t[0] if not pd.isna(t[0]) else pd.Timestamp.min))
    return out


def regularize_lineages(df: pd.DataFrame, lineages: Sequence[str],
                        include_sublineages: bool, base_map_fn=None) -> pd.DataFrame:
    """Add `variant`: each row's requested ancestor lineage, or its own name."""
    df = df.copy()
    if include_sublineages and base_map_fn is not None:
        replace_map = base_map_fn(list(lineages))
        df["variant"] = df["annotation_2"].map(replace_map).fillna(df["annotation_2"])
    else:
        df["variant"] = df["annotation_2"]
    return df


def prepare_clusters(df: pd.DataFrame, state: str, lineages: Sequence[str],
                     include_sublineages: bool = True, base_map_fn=None,
                     thresholds: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """One row per importation: cluster, variant, intro date, ordered samples.

    Sorted by intro date, which is both the schedule's order and the order the
    seed FASTA is written in, so record i founds importation i.
    """
    if "region" not in df.columns or "annotation_2" not in df.columns:
        raise ValueError("cluster table needs 'region' and 'annotation_2' columns")
    thr = {k: pd.to_datetime(v) for k, v in
           (thresholds if thresholds is not None else VARIANT_START_DEFAULTS).items()}

    d = regularize_lineages(df, lineages, include_sublineages, base_map_fn)
    d = d[(d["region"] == state) & (d["variant"].isin(list(lineages)))].copy()
    if d.empty:
        return pd.DataFrame(columns=["cluster_id", "variant", "intro_date",
                                     "sample_count", "samples_ordered"])

    d["intro_date"] = parse_cluster_dates(d["earliest_date"])
    d["samples_ordered"] = d["samples"].map(parse_samples)

    # Threshold correction: an intro date before the variant existed is
    # replaced by the cluster's first sample after the threshold.
    def corrected(row):
        limit = thr.get(row["variant"])
        if limit is None or pd.isna(row["intro_date"]) or row["intro_date"] >= limit:
            return row["intro_date"]
        later = [dt for dt, _s, _a in row["samples_ordered"]
                 if not pd.isna(dt) and dt > limit]
        return later[0] if later else pd.NaT

    d["intro_date"] = d.apply(corrected, axis=1)
    d = d.dropna(subset=["intro_date"])

    # Keep only samples the corrected intro date admits. Where the correction
    # moved a cluster forward, its pre-threshold samples are the very ones we
    # distrusted, so they must not become its seed either -- otherwise the
    # genome predates the importation it founds. Undated samples are kept as
    # last-resort fallbacks.
    def admissible(row):
        return [(dt, st, ac) for dt, st, ac in row["samples_ordered"]
                if pd.isna(dt) or dt >= row["intro_date"]]

    d["samples_ordered"] = d.apply(admissible, axis=1)
    d = d[d["samples_ordered"].map(len) > 0]
    # The cluster's own Pango call (AY.44, AY.103, ...) before the roll-up to
    # the requested lineage: the default benchmark variant label.
    d["sublineage"] = d["annotation_2"].fillna(d["variant"]).astype(str)
    cols = ["cluster_id", "variant", "sublineage", "intro_date", "sample_count", "samples_ordered"]
    cols = [c for c in cols if c in d.columns]
    return d[cols].sort_values(["intro_date", "cluster_id"]).reset_index(drop=True)


def build_schedule(clusters: pd.DataFrame) -> pd.DataFrame:
    """Importations per day: tick (0 = first importation), date, variant, counts."""
    if clusters.empty:
        return pd.DataFrame(columns=["tick", "date", "variant", "importations", "seq_count"])
    c = clusters.copy()
    first = c["intro_date"].min()
    c["tick"] = (c["intro_date"] - first).dt.days.astype(int)
    agg = {"importations": ("cluster_id", "size")}
    if "sample_count" in c.columns:
        agg["seq_count"] = ("sample_count", "sum")
    out = (c.groupby(["tick", "intro_date", "variant"], as_index=False)
             .agg(**agg)
             .rename(columns={"intro_date": "date"})
             .sort_values(["tick", "variant"]))
    return out.reset_index(drop=True)


def variant_schedule(clusters: pd.DataFrame, tick_zero, label: str = "sublineage") -> pd.DataFrame:
    """Importation schedule for `phylogas assign-variants`, on absolute ticks.

    Columns tick, date, variant, clusters, sample_count -- the format of the
    hand-made overlay schedules (e.g. TwinSampler's Virginia_importation_
    schedule.csv), so either can be handed to assign-variants. `variant` is
    each cluster's sublineage by default, which gives `variant_benchmark`
    real co-circulating labels (AY.44, AY.103, ...) where the requested
    lineage alone would give every importation the same one.
    """
    if clusters.empty:
        return pd.DataFrame(columns=["tick", "date", "variant", "clusters", "sample_count"])
    c = clusters.copy()
    c["label"] = c[label] if label in c.columns else c["variant"]
    agg = {"clusters": ("label", "size")}
    agg["sample_count"] = (("sample_count", "sum") if "sample_count" in c.columns
                           else ("label", "size"))
    out = (c.groupby(["intro_date", "label"], as_index=False).agg(**agg)
             .rename(columns={"intro_date": "date", "label": "variant"}))
    out["tick"] = (pd.to_datetime(out["date"]) - pd.to_datetime(tick_zero)).dt.days.astype(int)
    out["date"] = pd.to_datetime(out["date"]).dt.strftime("%Y-%m-%d")
    out["sample_count"] = out["sample_count"].astype(int)
    return (out[["tick", "date", "variant", "clusters", "sample_count"]]
            .sort_values(["tick", "variant"]).reset_index(drop=True))


def absolute_ticks(schedule: pd.DataFrame, tick_zero) -> Dict[str, pd.DataFrame]:
    """Per-variant `tick,count` frames on the ABM's own tick axis.

    This is the shape EpiHiper's seeding template consumes (one file per
    variant, as data/seed_schedule/<st>-<n>.csv).
    """
    t0 = pd.to_datetime(tick_zero)
    out = {}
    for variant, g in schedule.groupby("variant"):
        f = pd.DataFrame({
            "tick": (pd.to_datetime(g["date"]) - t0).dt.days.astype(int),
            "count": g["importations"].astype(int),
        }).sort_values("tick").reset_index(drop=True)
        out[str(variant)] = f
    return out


def seed_plan(clusters: pd.DataFrame) -> pd.DataFrame:
    """The seed strain per importation: each cluster's EARLIEST sample.

    `fallbacks` holds the remaining samples, oldest first, for clusters whose
    preferred sequence cannot be retrieved.
    """
    rows = []
    for row in clusters.itertuples(index=False):
        ordered = row.samples_ordered
        if not ordered:
            continue
        rows.append({
            "cluster_id": getattr(row, "cluster_id", ""),
            "variant": row.variant,
            "sublineage": getattr(row, "sublineage", row.variant),
            "intro_date": row.intro_date,
            "strain": ordered[0][1],
            "seed_date": ordered[0][0],
            "fallbacks": [s for _d, s, _a in ordered[1:]],
        })
    return pd.DataFrame(rows)
