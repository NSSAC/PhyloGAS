#!/usr/bin/env python3
"""Does a sampling strategy escape the epidemic's demographic homophily?

The claim under test is that LASSO-Greedy uses demographic targeting as a
topological compass, while SURS stays trapped inside densely-ascertained,
demographically-uniform clusters. Graph coverage cannot distinguish these:
both can touch many nodes. What distinguishes them is the *character of the
transmission edges they retain*.

Method
------
1. Build the true transmission graph from the ABM all-events file: one
   undirected edge per infection event, joining `alias_contact` (the
   infector's infection id) to `alias_pid` (the infectee's). Root
   importations have no infector and contribute no edge.
2. Attach demographics to every infection through `sim_pid`.
3. For the FULL graph, compute Newman's categorical assortativity r for each
   attribute. This is the epidemic's natural homophily and is the only
   meaningful baseline -- an absolute r means nothing without it.
4. For each sampling strategy, induce the subgraph on sampled infections:
   keep an edge only when BOTH endpoints were sampled. Recompute r.

Reading the output
------------------
`r` is Newman's assortativity coefficient: +1 means every edge joins two
infections of the same category, 0 means category-blind mixing, negative
means edges preferentially cross categories.

Compare each strategy's r against the `__full_graph__` row, NOT against zero:

  r above baseline   the sample is MORE homophilic than the epidemic --
                     it settled inside demographic islands
  r below baseline   the sample spans demographic boundaries more than the
                     disease does -- it broke homophilic gravity

`n_edges_retained` and `edge_retention` are not optional context. A strategy
that scatters samples so thinly that few true edges have both endpoints
present will report a low or undefined r for a reason that is a failure, not
a success: it shattered the transmission tree rather than traversing it.
Always read r and edge_retention together. An r computed on a handful of
edges is noise; `r` is reported as empty below --min-edges.

`cross_county_edge_frac` and `n_counties_touched` are the spatial analogues,
comparable to the "boundary node" and "counties touched" figures from the
earlier cross-replicate sweep.

Usage
-----
    python assortativity.py \
        --infections results/va_delta_wave/05_benchmarks/allevents_variants.csv.xz \
        --population data/va/va_2_4_0_demographics.csv \
        --samples 'results/va_delta_wave/03_sampled_datasets/*_samples.csv.xz' \
        --state va \
        --out poster_project/data/assortativity_va.csv

Run it where the data is (the cluster); the output is a few dozen rows.
"""

from __future__ import annotations

import argparse
import glob as _glob
import os
import sys

import numpy as np
import pandas as pd

# Candidate demographic attributes, in preference order. Whichever are present
# in the population file get measured; BeyondBaseline's default stratifiers are
# age, race, county, sex.
# BeyondBaseline's STRAT_ALIAS maps age->age_group, race->smh_race, county->county,
# so the demographics file carries the right-hand spellings.
# The all-events file already carries the demographic columns, so the
# population file is usually unnecessary. rucc_code is the rural-urban
# continuum code: the urban/rural axis the homophily argument turns on.
CANDIDATE_ATTRS = ("age_group", "age", "smh_race", "race", "ethnicity",
                   "county", "county_fips", "fips", "rucc_code",
                   "ses_category", "sex", "gender")

ID_DTYPES = {"alias_pid": str, "alias_contact": str, "sim_pid": str, "pid": str}


def _pid_col(df: pd.DataFrame) -> str:
    for c in ("sim_pid", "pid"):
        if c in df.columns:
            return c
    raise SystemExit("no sim_pid or pid column found")


def newman_r(a: np.ndarray, b: np.ndarray) -> "float | None":
    """Newman's assortativity for a categorical attribute over edges (a_i, b_i).

    r = (Tr e - ||e^2||) / (1 - ||e^2||), with e the normalised mixing matrix
    symmetrised over both edge orientations. Returns None when the attribute
    takes a single value, where r is undefined rather than zero.
    """
    keep = (a >= 0) & (b >= 0)
    a, b = a[keep], b[keep]
    if a.size == 0:
        return None
    k = int(max(a.max(), b.max())) + 1
    if k < 2:
        return None
    # Count each edge in both directions so the mixing matrix is symmetric.
    e = np.bincount(a * k + b, minlength=k * k).astype(np.float64)
    e += np.bincount(b * k + a, minlength=k * k)
    total = e.sum()
    if total <= 0:
        return None
    e = (e / total).reshape(k, k)
    tr = float(np.trace(e))
    sq = float((e.sum(axis=0) * e.sum(axis=1)).sum())
    if abs(1.0 - sq) < 1e-12:
        return None
    return (tr - sq) / (1.0 - sq)


def label_from_path(path: str) -> str:
    """'1S-P__lasso_greedy_samples.csv.xz' -> '1S-P__lasso_greedy'."""
    base = os.path.basename(path)
    for suf in ("_samples.csv.xz", "_samples.csv.gz", "_samples.csv"):
        if base.endswith(suf):
            return base[: -len(suf)]
    return base.split(".")[0]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--infections", required=True,
                    help="ABM all-events file (after `phylogas assign-variants`)")
    ap.add_argument("--population", default=None,
                    help="demographics CSV; only needed for attributes the "
                         "all-events file does not already carry")
    ap.add_argument("--samples", required=True, nargs="+",
                    help="sample CSVs or globs, one per strategy")
    ap.add_argument("--state", default="", help="label for the output rows")
    ap.add_argument("--attributes", default="",
                    help="comma-separated attribute columns; default: auto-detect")
    ap.add_argument("--min-edges", type=int, default=200,
                    help="below this many retained edges, r is left empty (default 200)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    print(f"Infections : {args.infections}")
    inf = pd.read_csv(args.infections, dtype=ID_DTYPES, low_memory=False)
    inf.columns = [c.strip() for c in inf.columns]
    if "alias_pid" not in inf.columns or "alias_contact" not in inf.columns:
        raise SystemExit(f"need alias_pid and alias_contact; found {list(inf.columns)}")
    # One row per infection: alias_pid is 1:1 with a sequence, sim_pid is not.
    inf = inf.drop_duplicates(subset=["alias_pid"], keep="first")
    print(f"  {len(inf):,} infections")

    pop = None
    if args.population:
        print(f"Population : {args.population}")
        pop = pd.read_csv(args.population, dtype=ID_DTYPES, low_memory=False)
        pop.columns = [c.strip() for c in pop.columns]

    if args.attributes:
        attrs = [a.strip() for a in args.attributes.split(",") if a.strip()]
    else:
        have = set(inf.columns) | (set(pop.columns) if pop is not None else set())
        attrs = [a for a in CANDIDATE_ATTRS if a in have]
        # age and age_group are redundant; prefer the grouped form.
        if "age_group" in attrs and "age" in attrs:
            attrs.remove("age")
        if "smh_race" in attrs and "race" in attrs:
            attrs.remove("race")
        if "county" in attrs:
            for dup in ("county_fips", "fips"):
                if dup in attrs:
                    attrs.remove(dup)
    if not attrs:
        raise SystemExit(f"no demographic attributes found in {list(inf.columns)[:25]}...")
    print(f"  attributes: {', '.join(attrs)}")

    if pop is not None:
        from_pop = [a for a in attrs if a in pop.columns and a not in inf.columns]
        if from_pop:
            inf = inf.merge(pop[[_pid_col(pop)] + from_pop],
                            left_on=_pid_col(inf), right_on=_pid_col(pop),
                            how="left", suffixes=("", "_pop"))

    # ---- node index and attribute codes -----------------------------------
    inf = inf.reset_index(drop=True)
    node_of = pd.Series(np.arange(len(inf), dtype=np.int64), index=inf["alias_pid"])
    codes = {}
    for a in attrs:
        if a not in inf.columns:
            print(f"  note: {a} absent after join; skipping", file=sys.stderr)
            continue
        c, _ = pd.factorize(inf[a].astype("string"), use_na_sentinel=True)
        codes[a] = c.astype(np.int64)
    attrs = [a for a in attrs if a in codes]

    # ---- edges -------------------------------------------------------------
    src = node_of.reindex(inf["alias_contact"]).to_numpy()
    dst = np.arange(len(inf), dtype=np.float64)
    ok = ~np.isnan(src)
    src_i = src[ok].astype(np.int64)
    dst_i = dst[ok].astype(np.int64)
    # A root importation has no infector; a self-edge carries no mixing signal.
    keep = src_i != dst_i
    src_i, dst_i = src_i[keep], dst_i[keep]
    n_edges = src_i.size
    print(f"Graph      : {n_edges:,} transmission edges "
          f"({len(inf) - n_edges:,} roots/unmatched)")
    if n_edges == 0:
        raise SystemExit("no edges; check alias_contact values")

    county_code = codes.get("county", codes.get("county_fips", codes.get("fips")))

    def measure(name: str, mask: "np.ndarray | None", n_nodes: int) -> dict:
        s, d = (src_i, dst_i) if mask is None else (src_i[mask], dst_i[mask])
        row = {"state": args.state, "strategy": name,
               "n_nodes_sampled": n_nodes,
               "n_edges_retained": int(s.size),
               "edge_retention": round(s.size / n_edges, 6)}
        for a in attrs:
            c = codes[a]
            r = newman_r(c[s], c[d]) if s.size >= args.min_edges else None
            row[f"r_{a}"] = None if r is None else round(r, 6)
            # The mixing matrix diagonal, stated plainly: what share of the
            # retained edges joined two infections of the same category. Not
            # redundant with r -- r is normalised against the marginals and
            # saturates near 1 on a strongly clustered graph, where this
            # fraction still moves and is what a reader can interpret.
            cs_, cd_ = c[s], c[d]
            ok_ = (cs_ >= 0) & (cd_ >= 0)
            row[f"within_{a}_frac"] = (round(float((cs_[ok_] == cd_[ok_]).mean()), 6)
                                       if ok_.any() else None)
        if county_code is not None and s.size:
            cs, cd = county_code[s], county_code[d]
            valid = (cs >= 0) & (cd >= 0)
            row["cross_county_edge_frac"] = (round(float((cs[valid] != cd[valid]).mean()), 6)
                                             if valid.any() else None)
            touched = np.unique(np.concatenate([cs[valid], cd[valid]]))
            row["n_counties_touched"] = int(touched.size)
        else:
            row["cross_county_edge_frac"] = None
            row["n_counties_touched"] = None
        return row

    out = [measure("__full_graph__", None, len(inf))]
    print(f"  baseline: " + ", ".join(
        f"{a}={out[0][f'r_{a}']}" for a in attrs if out[0][f"r_{a}"] is not None))

    paths: list[str] = []
    for spec in args.samples:
        paths.extend(sorted(_glob.glob(spec)) or ([spec] if os.path.exists(spec) else []))
    if not paths:
        raise SystemExit(f"no sample files matched: {args.samples}")

    for p in paths:
        sdf = pd.read_csv(p, dtype=ID_DTYPES, low_memory=False)
        sdf.columns = [c.strip() for c in sdf.columns]
        if "alias_pid" not in sdf.columns:
            print(f"  WARNING: {os.path.basename(p)} has no alias_pid; skipped",
                  file=sys.stderr)
            continue
        idx = node_of.reindex(sdf["alias_pid"].dropna().unique()).dropna().astype(np.int64)
        sampled = np.zeros(len(inf), dtype=bool)
        sampled[idx.to_numpy()] = True
        mask = sampled[src_i] & sampled[dst_i]
        row = measure(label_from_path(p), mask, int(sampled.sum()))
        out.append(row)
        base = out[0]
        delta = ", ".join(
            f"{a}: {row[f'r_{a}']:+.3f} ({row[f'r_{a}'] - base[f'r_{a}']:+.3f} vs graph)"
            for a in attrs
            if row.get(f"r_{a}") is not None and base.get(f"r_{a}") is not None)
        print(f"  {row['strategy']:<28} n={row['n_nodes_sampled']:>7,} "
              f"edges={row['n_edges_retained']:>7,} "
              f"({row['edge_retention']*100:5.1f}%)  {delta}")

    df = pd.DataFrame(out)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"\n-> {args.out}  ({len(df)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
