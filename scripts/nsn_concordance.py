#!/usr/bin/env python3
"""Nearest-sampled-neighbour concordance: does a strategy walk bridges?

Replaces the sampled-subgraph assortativity measure, which fails at realistic
sampling fractions. Sampling p of the population retains only p^2 of the
transmission edges, so at p = 0.28% roughly 350 of Virginia's 5.3M edges have
both endpoints sampled -- statistical dust, and biased: a quota-filling
algorithm that samples heavily inside one cluster captures a few direct,
homophilic edges, while uniform random sampling's handful of retained edges
are chance pairings that cross boundaries by accident. The measure ranked the
algorithms backwards for that reason.

This measure needs no adjacency. For each sampled infection it finds the
nearest OTHER sampled infection along the true transmission tree and asks
whether the two sit on the same side of a demographic or geographic boundary.
Every sample contributes a pair, so N goes from ~350 edges to ~14,600 pairs.

It is also the right question. In real surveillance an infector-infectee pair
is almost never both sequenced; two sequenced tips are several transmission
hops apart and the phylogeny links them through an unobserved common ancestor.
The path from a sample to its nearest sampled neighbour is the topological
equivalent of that branch -- so asking whether it crosses a county or a
demographic boundary is asking what mugration has to infer.

Three families of output, from one graph traversal:

  nsn_hops_*            how far apart the samples sit on the tree.
                        "Pushed further apart" made measurable.
  cross_<attr>_frac     share of nearest-neighbour pairs whose endpoints
                        differ in that attribute. Demographic dissimilarity.
  purity_<attr>         Voronoi cell purity: assign every infection to its
                        nearest sample, then ask what share of them match the
                        demographic of the sample that owns them. A flow-style
                        partition of the graph into one basin per sample.

Direction of merit -- these do NOT all point the same way, and two of them
measure different things:

  nsn_hops_*            LARGER under the bridge-walking hypothesis. Samples
                        pushed apart on the tree rather than piled together.
  cross_<attr>_frac     LARGER under the hypothesis. The unobserved path
                        between two sequenced tips crosses a boundary, which
                        is the inference mugration has to make.
  purity_<attr>         LARGER is BETTER, and it is a representativeness
                        measure, not a homophily measure. Low purity means
                        infections are routinely nearest to a sample of some
                        other group -- i.e. whole groups have no sample near
                        them and get represented by a dissimilar one. A
                        strategy that piles its budget onto a few clusters
                        shows low purity because its basins must stretch
                        across the groups it ignored.

So the predicted signature of a bridge-walking strategy is high nsn_hops,
high cross_<attr>_frac AND high purity_<attr>: spread out, crossing
boundaries, and leaving no group unrepresented. A strategy can score high on
cross_<attr>_frac for the wrong reason -- scattering so thinly that nothing
is near anything -- and purity is what catches that.

All of these are relative comparisons between strategies on one shared graph,
so there is no baseline row to compute.

    python scripts/nsn_concordance.py \
        --infections results/va_delta_wave/05_benchmarks/allevents_variants.csv.xz \
        --samples 'results/va_delta_wave/03_sampled_datasets/*_samples.csv.xz' \
        --state va --out poster_project/data/nsn_va.csv
"""

from __future__ import annotations

import argparse
import glob as _glob
import os
import sys

import numpy as np
import pandas as pd

CANDIDATE_ATTRS = ("age_group", "age", "smh_race", "race", "ethnicity",
                   "county", "county_fips", "fips", "rucc_code",
                   "ses_category", "sex", "gender")
ID_DTYPES = {"alias_pid": str, "alias_contact": str, "sim_pid": str, "pid": str}


def build_csr(src: np.ndarray, dst: np.ndarray, n: int):
    """Undirected adjacency in CSR form, numpy only (no scipy dependency)."""
    a = np.concatenate([src, dst])
    b = np.concatenate([dst, src])
    order = np.argsort(a, kind="stable")
    indices = b[order].astype(np.int32)
    counts = np.bincount(a, minlength=n)
    indptr = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(counts, out=indptr[1:])
    return indptr, indices


def voronoi_bfs(indptr, indices, n: int, seeds: np.ndarray):
    """Multi-source BFS: distance to the nearest seed, and which seed that is.

    Level-synchronous and vectorised: each round expands the whole frontier at
    once, which is what keeps this tractable on a 5M-node graph.
    """
    dist = np.full(n, -1, dtype=np.int32)
    owner = np.full(n, -1, dtype=np.int32)
    dist[seeds] = 0
    owner[seeds] = seeds
    frontier = seeds
    d = 0
    while frontier.size:
        d += 1
        starts, ends = indptr[frontier], indptr[frontier + 1]
        deg = ends - starts
        if not deg.any():
            break
        # Flatten the frontier's neighbour lists without a Python loop.
        total = int(deg.sum())
        pos = np.repeat(starts, deg) + (np.arange(total) - np.repeat(np.cumsum(deg) - deg, deg))
        nbr = indices[pos]
        parent = np.repeat(frontier, deg)
        fresh = dist[nbr] < 0
        nbr, parent = nbr[fresh], parent[fresh]
        if nbr.size == 0:
            break
        # First writer wins; ties between seeds are arbitrary but consistent.
        uniq, first = np.unique(nbr, return_index=True)
        dist[uniq] = d
        owner[uniq] = owner[parent[first]]
        frontier = uniq
    return dist, owner


def nearest_other(indptr, indices, dist, owner, src, dst, seeds):
    """For each seed, the distance to the nearest OTHER seed, and which seed.

    A shortest path between two seeds must cross an edge whose endpoints
    belong to different Voronoi cells, and for that edge the path length is
    exactly dist[a] + 1 + dist[b]. So scanning the boundary edges once gives
    every seed's nearest other seed exactly -- no per-seed search.
    """
    best_d = np.full(owner.size, -1, dtype=np.int64)
    best_o = np.full(owner.size, -1, dtype=np.int64)
    oa, ob = owner[src], owner[dst]
    boundary = (oa != ob) & (oa >= 0) & (ob >= 0)
    if not boundary.any():
        # Every seed sits in its own component: no seed has a reachable other.
        return best_d, best_o
    a, b = src[boundary], dst[boundary]
    cand = (dist[a] + 1 + dist[b]).astype(np.int64)
    ow_a, ow_b = owner[a], owner[b]
    # Each boundary edge informs both of its owners.
    who = np.concatenate([ow_a, ow_b])
    other = np.concatenate([ow_b, ow_a])
    cc = np.concatenate([cand, cand])
    order = np.lexsort((cc, who))          # by owner, then by distance
    who, other, cc = who[order], other[order], cc[order]
    firsts = np.ones(who.size, dtype=bool)
    firsts[1:] = who[1:] != who[:-1]       # first row per owner is its minimum
    best_d[who[firsts]] = cc[firsts]
    best_o[who[firsts]] = other[firsts]
    return best_d, best_o


def label_from_path(path: str) -> str:
    base = os.path.basename(path)
    for suf in ("_samples.csv.xz", "_samples.csv.gz", "_samples.csv"):
        if base.endswith(suf):
            return base[: -len(suf)]
    return base.split(".")[0]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--infections", required=True)
    ap.add_argument("--samples", required=True, nargs="+")
    ap.add_argument("--state", default="")
    ap.add_argument("--attributes", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    print(f"Infections : {args.infections}")
    inf = pd.read_csv(args.infections, dtype=ID_DTYPES, low_memory=False)
    inf.columns = [c.strip() for c in inf.columns]
    inf = inf.drop_duplicates(subset=["alias_pid"], keep="first").reset_index(drop=True)
    n = len(inf)
    print(f"  {n:,} infections")

    if args.attributes:
        attrs = [a.strip() for a in args.attributes.split(",") if a.strip()]
    else:
        attrs = [a for a in CANDIDATE_ATTRS if a in inf.columns]
        if "age_group" in attrs and "age" in attrs:
            attrs.remove("age")
        if "smh_race" in attrs and "race" in attrs:
            attrs.remove("race")
        if "county" in attrs:
            for dup in ("county_fips", "fips"):
                if dup in attrs:
                    attrs.remove(dup)
    print(f"  attributes: {', '.join(attrs)}")
    codes = {a: pd.factorize(inf[a].astype('string'), use_na_sentinel=True)[0].astype(np.int32)
             for a in attrs}

    node_of = pd.Series(np.arange(n, dtype=np.int64), index=inf["alias_pid"])
    src = node_of.reindex(inf["alias_contact"]).to_numpy()
    ok = ~np.isnan(src)
    src_i = src[ok].astype(np.int64)
    dst_i = np.arange(n, dtype=np.int64)[ok]
    keep = src_i != dst_i
    src_i, dst_i = src_i[keep], dst_i[keep]
    print(f"Graph      : {src_i.size:,} transmission edges")

    indptr, indices = build_csr(src_i, dst_i, n)

    paths: list[str] = []
    for spec in args.samples:
        paths.extend(sorted(_glob.glob(spec)) or ([spec] if os.path.exists(spec) else []))
    if not paths:
        raise SystemExit(f"no sample files matched: {args.samples}")

    rows = []
    for p in paths:
        sdf = pd.read_csv(p, dtype=ID_DTYPES, low_memory=False)
        sdf.columns = [c.strip() for c in sdf.columns]
        idx = node_of.reindex(sdf["alias_pid"].dropna().unique()).dropna()
        seeds = np.unique(idx.to_numpy().astype(np.int32))
        dist, owner = voronoi_bfs(indptr, indices, n, seeds)
        best_d, best_o = nearest_other(indptr, indices, dist, owner, src_i, dst_i, seeds)

        row = {"state": args.state, "strategy": label_from_path(p),
               "n_samples": int(seeds.size)}
        have = best_d[seeds] > 0
        hops = best_d[seeds][have]
        partner = best_o[seeds][have].astype(np.int64)
        mine = seeds[have].astype(np.int64)
        row.update(n_pairs=int(hops.size),
                   pair_coverage=round(float(hops.size / seeds.size), 6),
                   nsn_hops_mean=round(float(hops.mean()), 4),
                   nsn_hops_median=float(np.median(hops)),
                   nsn_hops_p90=float(np.percentile(hops, 90)))
        for a in attrs:
            c = codes[a]
            va_, vb_ = c[mine], c[partner]
            valid = (va_ >= 0) & (vb_ >= 0)
            row[f"cross_{a}_frac"] = (round(float((va_[valid] != vb_[valid]).mean()), 6)
                                      if valid.any() else None)
            # Voronoi cell purity over every reached infection.
            reached = owner >= 0
            own_attr = c[owner[reached]]
            node_attr = c[reached]
            v2 = (own_attr >= 0) & (node_attr >= 0)
            row[f"purity_{a}"] = (round(float((own_attr[v2] == node_attr[v2]).mean()), 6)
                                  if v2.any() else None)
        rows.append(row)
        print(f"  {row['strategy']:<28} pairs={row['n_pairs']:>7,} "
              f"hops={row['nsn_hops_mean']:>6.2f} " +
              " ".join(f"{a[:8]}:{row[f'cross_{a}_frac']:.3f}" for a in attrs[:3]))

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"\n-> {args.out}  ({len(df)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
