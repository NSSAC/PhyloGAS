#!/usr/bin/env python3
"""Does mu_truth measure the painter's clock, or the importation regime?

`mu_truth` regresses divergence-from-Wuhan on calendar date, pooled across
every tip. Two things make that not the painter's rate:

  1. Each chain starts at its seed's own divergence -- measured at mean 41.6,
     sd 4.7 substitutions for this run's 3,322 seeds -- so the y-axis carries
     a per-chain constant that is not elapsed time.
  2. The painter accumulates along time-since-importation, not calendar date.
     With importation continuing through the wave, every late date holds a
     mixture of old and young chains, so expected divergence stops growing
     with calendar time and the pooled slope attenuates.

This measures the same data three ways so the three can be compared:

  pooled        divergence vs calendar date            == what mu_truth does
  within        (divergence - chain's own baseline) vs days since chain start
  per_component the within fit done separately per large chain

If `within` lands near the literature ~8e-4 while `pooled` sits far below it,
the attenuation is the import mixture and `mu_truth` is a property of the
wave rather than of the painter. If `within` is also low, the painter's
emergent rate really is low and that is the finding.

Read-only. Nothing is written unless --out is given.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from phylogas.benchmark.clock import (            # noqa: E402
    _codes, _iter_fasta, _open_text, _valid_mask, hamming_to_reference,
)


def _header(fh):
    head = fh.readline().rstrip("\r\n")
    delim = "\t" if head.count("\t") >= head.count(",") else ","
    return head.split(delim), delim


def _read_components(path, key_col, comp_col):
    """join key -> component id, from TwinSampler's all-events table.

    The painted metadata carries no component: the painter writes
    `real_strain` only for the index case (populate_sim_details is passed a
    seed id only when seed_fasta is not None), so a chain cannot be recovered
    from it. The all-events table has one row per infection with both
    `alias_pid` and `component_id`, which is the join.

    This holds one dict entry per infection, so it is the memory high-water
    mark of the script. Only keys that also appear in the painted metadata
    matter, but filtering needs that set first, so the caller passes it in.
    """
    out = {}
    with _open_text(path) as fh:
        cols, delim = _header(fh)
        for want in (key_col, comp_col):
            if want not in cols:
                sys.exit(f"ERROR: {path} has no {want!r} column.\n"
                         f"       columns: {', '.join(cols[:14])}")
        i_k, i_c = cols.index(key_col), cols.index(comp_col)
        for row in csv.reader(fh, delimiter=delim):
            if len(row) > max(i_k, i_c):
                out[row[i_k]] = row[i_c]
    return out


def _read_meta(path, id_col, date_col, comp_col, join_col=None, comps=None):
    """strain -> (date, component). Streamed: the painted metadata is large.

    The component comes either from a column in this file, or -- for the
    painted metadata, which has none -- by joining `join_col` against the
    all-events mapping in `comps`.
    """
    import datetime

    out = {}
    n_nocomp = 0
    with _open_text(path) as fh:
        cols, delim = _header(fh)
        try:
            i_id, i_dt = cols.index(id_col), cols.index(date_col)
        except ValueError:
            sys.exit(f"ERROR: {path} needs {id_col!r} and {date_col!r}; "
                     f"has {', '.join(cols[:14])}")
        i_cp = cols.index(comp_col) if comp_col in cols else None
        i_jn = None
        if i_cp is None:
            if comps is None:
                sys.exit(
                    f"ERROR: {path} has no {comp_col!r} column, so chains "
                    f"cannot be identified.\n"
                    f"       The painted metadata never has one. Pass "
                    f"--components <allevents csv> to join on "
                    f"{join_col!r}.\n"
                    f"       columns: {', '.join(cols[:14])}")
            if join_col not in cols:
                sys.exit(f"ERROR: {path} has neither {comp_col!r} nor the "
                         f"join key {join_col!r}.\n"
                         f"       columns: {', '.join(cols[:14])}")
            i_jn = cols.index(join_col)
        for row in csv.reader(fh, delimiter=delim):
            if len(row) <= max(i_id, i_dt):
                continue
            try:
                d = datetime.date.fromisoformat(row[i_dt][:10])
            except ValueError:
                continue
            if i_cp is not None:
                comp = row[i_cp]
            else:
                comp = comps.get(row[i_jn]) if len(row) > i_jn else None
                if comp is None:
                    n_nocomp += 1
                    continue
            out[row[i_id]] = (d, comp)
    if n_nocomp:
        print(f"  note: {n_nocomp:,} metadata rows had no component after the "
              f"join and were dropped")
    return out


def _ols(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    icpt = my - slope * mx
    sst = sum((y - my) ** 2 for y in ys)
    sse = sum((y - (icpt + slope * x)) ** 2 for x, y in zip(xs, ys))
    return {"slope": slope, "intercept": icpt, "n": n,
            "r2": (1 - sse / sst) if sst else float("nan")}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-f", "--fasta", required=True, help="painted sequences")
    ap.add_argument("-m", "--metadata", required=True, help="painted metadata")
    ap.add_argument("-r", "--reference", required=True)
    ap.add_argument("--id-col", default="strain")
    ap.add_argument("--date-col", default="date")
    ap.add_argument("--component-col", default="component_id")
    ap.add_argument("--components", default=None,
                    help="all-events csv carrying component_id, joined when "
                         "the metadata has no component column (the painted "
                         "metadata never does)")
    ap.add_argument("--join-col", default="alias_pid",
                    help="key shared by the metadata and --components "
                         "(default: alias_pid)")
    ap.add_argument("--min-component", type=int, default=200,
                    help="chains smaller than this are pooled but not fitted "
                         "individually (default: 200)")
    ap.add_argument("--min-span", type=int, default=30,
                    help="skip per-chain fits spanning fewer days (default: 30)")
    ap.add_argument("--max-records", type=int, default=None)
    ap.add_argument("--out", default=None, help="write per-component fits here")
    args = ap.parse_args()

    ref = next(iter(_iter_fasta(args.reference)))[1]
    rc, rv = _codes(ref), _valid_mask(_codes(ref))
    sites = int(rv.sum())

    comps = None
    if args.components:
        comps = _read_components(args.components, args.join_col,
                                 args.component_col)
        print(f"components : {len(comps):,} {args.join_col} -> "
              f"{args.component_col} from {args.components}")
    meta = _read_meta(args.metadata, args.id_col, args.date_col,
                      args.component_col, args.join_col, comps)
    print(f"metadata   : {len(meta):,} rows with a usable date "
          f"and a component")

    by_comp = defaultdict(list)
    seen = kept = unaligned = nometa = 0
    for sid, seq in _iter_fasta(args.fasta):
        seen += 1
        hit = meta.get(sid)
        if hit is None:
            nometa += 1
            continue
        div = hamming_to_reference(seq, ref, rc, rv)
        if div is None:
            unaligned += 1
            continue
        d, comp = hit
        by_comp[comp].append((d.toordinal(), float(div)))
        kept += 1
        if args.max_records and kept >= args.max_records:
            break
    print(f"sequences  : {seen:,} seen, {kept:,} used "
          f"({nometa:,} no metadata, {unaligned:,} unaligned)")
    print(f"chains     : {len(by_comp):,}  "
          f"(>= {args.min_component} tips: "
          f"{sum(1 for v in by_comp.values() if len(v) >= args.min_component):,})")
    print(f"comparable sites: {sites:,}\n")

    def rate(slope_per_day):
        return slope_per_day / sites * 365.0

    # --- pooled: what mu_truth computes -------------------------------
    allpts = [p for v in by_comp.values() for p in v]
    pooled = _ols([x for x, _ in allpts], [y for _, y in allpts])

    # --- within: recentre each chain on its own start -----------------
    wx, wy = [], []
    for pts in by_comp.values():
        t0 = min(x for x, _ in pts)
        base = min(y for _, y in pts)
        for x, y in pts:
            wx.append(x - t0)
            wy.append(y - base)
    within = _ols(wx, wy)

    for label, fit, note in (
            ("pooled   (== mu_truth)", pooled, "divergence vs calendar date"),
            ("within   (recentred)", within, "excess divergence vs days since chain start")):
        if not fit:
            print(f"{label}: too few points"); continue
        print(f"{label}")
        print(f"    {note}")
        print(f"    rate      {rate(fit['slope']):.4e} subs/site/year")
        print(f"    slope     {fit['slope']:.5f} subs/genome/day")
        print(f"    intercept {fit['intercept']:.2f} subs")
        print(f"    n {fit['n']:,}   R2 {fit['r2']:.3f}\n")

    # --- per chain ----------------------------------------------------
    rows = []
    for comp, pts in by_comp.items():
        if len(pts) < args.min_component:
            continue
        span = max(x for x, _ in pts) - min(x for x, _ in pts)
        if span < args.min_span:
            continue
        fit = _ols([x for x, _ in pts], [y for _, y in pts])
        if fit:
            rows.append((comp, len(pts), span, rate(fit["slope"]), fit["r2"]))
    rows.sort(key=lambda r: -r[1])

    if rows:
        rs = sorted(r[3] for r in rows)
        mid = rs[len(rs) // 2]
        print(f"per-chain fits: {len(rows)}  "
              f"median {mid:.4e}  min {rs[0]:.4e}  max {rs[-1]:.4e}")
        print(f"\n{'component':>14}  {'tips':>7}  {'span':>5}  "
              f"{'rate':>11}  {'R2':>5}")
        for comp, n, span, r, r2 in rows[:15]:
            print(f"{comp:>14}  {n:>7,}  {span:>4}d  {r:>11.4e}  {r2:>5.3f}")
        if len(rows) > 15:
            print(f"{'...':>14}  {len(rows) - 15} more")
    else:
        print(f"no chain has >= {args.min_component} tips spanning "
              f">= {args.min_span} days; lower --min-component/--min-span")

    if args.out and rows:
        with open(args.out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["component", "tips", "span_days",
                        "rate_subs_per_site_per_year", "r_squared"])
            w.writerows(rows)
        print(f"\nWrote {args.out}  ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
