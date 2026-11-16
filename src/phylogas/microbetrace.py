"""Export MicrobeTrace-ready sequence, node, and ground-truth link bundles.

MicrobeTrace (https://microbetrace.cdc.gov) computes pairwise genetic distances
and renders interactive 2D/3D networks, epi curves, and choropleth/pin maps in
the browser. Because browser-side alignment and O(N^2) distance calculation work
best on 200 to 1,500 sequences, this module extracts curated, self-contained
subsets from a PhyloGAS run and emits four files per bundle:

1. ``<prefix>_sequences.fasta``
   Uncompressed multi-FASTA whose record headers match the ``id`` column of the
   node and link tables 1:1 (no trailing header commentary that could trip
   browser parsers).

2. ``<prefix>_nodes.csv``
   Node metadata table keyed on ``id`` carrying epidemiological fields
   (``collection_date``, ``exposure_date``, ``county``, ``county_fips``,
   ``latitude``, ``longitude``, ``rucc_code``, ``urbanicity``, ``age_group``,
   ``sex``, ``race``, ``smh_race``, ``ses_category``) alongside ground-truth
   phylodynamic fields (``variant_benchmark``, ``true_component_id``,
   ``true_component_size``, ``subset_component_size``,
   ``sampled_in_surveillance``, ``true_infector_id``,
   ``nearest_subset_ancestor_id``, ``hops_to_subset_ancestor``,
   ``snp_distance_to_ancestor``).

3. ``<prefix>_truth_links.csv``
   Ground-truth transmission edge list with ``source`` and ``target`` columns
   ready to drop into MicrobeTrace as a link table. Includes both direct
   infector->infectee edges (``transmission_hops = 1``, ``edge_type = "direct"``)
   and ancestral bridge edges across unsampled intermediate hosts up to
   ``--max-hops`` (``transmission_hops = k``, ``edge_type = "ancestral"``),
   annotated with ``snp_distance``, ``tn93_distance``, ``substitutions``, and
   ``delta_exposure_days``.

4. ``<prefix>_summary.json``
   Bundle diagnostics, hop-stratified SNP divergence, and single-linkage SNP
   threshold clustering concordance (Precision, Recall, F1, Adjusted Rand Index
   at tau = 0, 1, 2, 3, 5 SNPs) comparing genetic clusters against true
   transmission chains.

Implemented using only the Python standard library so it runs in a bare Python
environment as well as inside the full Conda/Snakemake environment.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import gzip
import json
import lzma
import math
import random
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


_NO_CONTACT = {"", "-1", "nan", "None", "<NA>", "null"}

_DEFAULT_CONTEXT_IDS = {
    "Wuhan/Hu-1/2019",
    "MN908947.3",
    "NC_045512.2",
}


def _open_text(path: str | Path, mode: str = "rt"):
    """Open a plain, .xz, or .gz text file for reading or writing."""
    p = str(path)
    if p.endswith(".xz"):
        return lzma.open(p, mode, encoding="utf-8")
    if p.endswith(".gz"):
        return gzip.open(p, mode, encoding="utf-8")
    raw_mode = mode.replace("t", "")
    return open(p, raw_mode, encoding="utf-8", newline="" if "w" in raw_mode else None)


def _sniff_delimiter(header_line: str) -> str:
    return "\t" if header_line.count("\t") > header_line.count(",") else ","


def _parse_iso_date(val: str | None) -> _dt.date | None:
    if not val:
        return None
    s = str(val).strip()[:10]
    if not s or s.lower() in _NO_CONTACT:
        return None
    try:
        return _dt.date.fromisoformat(s)
    except ValueError:
        return None


def _alias_from_strain(strain: str) -> str | None:
    """Extract ``{pid}.{tick}`` from ``USA/VA-EHip-{pid}.{tick}/{year}``."""
    if not strain or "EHip-" not in strain:
        return None
    try:
        return strain.split("EHip-", 1)[1].rsplit("/", 1)[0]
    except IndexError:
        return None


def _urbanicity_label(rucc_raw: str | None) -> str:
    if not rucc_raw or str(rucc_raw).strip() in _NO_CONTACT:
        return "Unknown"
    try:
        code = int(float(str(rucc_raw).strip()))
    except ValueError:
        return "Unknown"
    if 1 <= code <= 3:
        return "Metro (RUCC 1-3)"
    if 4 <= code <= 9:
        return "Non-metro / Rural (RUCC 4-9)"
    return "Unknown"


def _snp_and_mutations(seq_source: str, seq_target: str, max_mut_list: int = 25) -> tuple[int, float, str]:
    """Compute exact ACGT SNP count, Tamura-Nei (TN93) distance, and mutation string.

    Ambiguous bases (N, gaps, IUPAC codes outside A/C/G/T) are ignored, matching
    MicrobeTrace's default pairwise SNP and TN93 comparison behavior on aligned
    SARS-CoV-2 sequences.
    """
    valid = {"A", "C", "G", "T"}
    purines = {"A", "G"}
    pyrimidines = {"C", "T"}

    snps = 0
    ts_pur = 0  # A <-> G transitions
    ts_pyr = 0  # C <-> T transitions
    tv = 0      # transversions
    compared = 0
    counts_a = 0
    counts_c = 0
    counts_g = 0
    counts_t = 0
    muts: list[str] = []

    s1 = seq_source.upper()
    s2 = seq_target.upper()
    for pos_0, (a, b) in enumerate(zip(s1, s2)):
        if a not in valid or b not in valid:
            continue
        compared += 1
        if a == "A":
            counts_a += 1
        elif a == "C":
            counts_c += 1
        elif a == "G":
            counts_g += 1
        else:
            counts_t += 1
        if b == "A":
            counts_a += 1
        elif b == "C":
            counts_c += 1
        elif b == "G":
            counts_g += 1
        else:
            counts_t += 1

        if a != b:
            snps += 1
            if len(muts) <= max_mut_list:
                muts.append(f"{a}{pos_0 + 1}{b}")
            if a in purines and b in purines:
                ts_pur += 1
            elif a in pyrimidines and b in pyrimidines:
                ts_pyr += 1
            else:
                tv += 1

    if compared == 0:
        return 0, 0.0, ""

    # TN93 distance formula (Tamura & Nei 1993), falling back to p-distance
    # if base frequencies are degenerate or logarithms saturate.
    total_bases = 2.0 * compared
    fa = counts_a / total_bases
    fc = counts_c / total_bases
    fg = counts_g / total_bases
    ft = counts_t / total_bases
    fr = fa + fg
    fy = fc + ft

    p1 = ts_pur / compared
    p2 = ts_pyr / compared
    q = tv / compared
    p_dist = snps / compared

    tn93 = p_dist
    if fa > 0 and fg > 0 and fc > 0 and ft > 0 and fr > 0 and fy > 0:
        k1 = 2.0 * fa * fg / fr
        k2 = 2.0 * fc * ft / fy
        k3 = 2.0 * fr * fy
        term1 = 1.0 - (p1 / k1) - (q / (2.0 * fr)) if k1 > 0 else -1.0
        term2 = 1.0 - (p2 / k2) - (q / (2.0 * fy)) if k2 > 0 else -1.0
        term3 = 1.0 - (q / k3) if k3 > 0 else -1.0
        if term1 > 0 and term2 > 0 and term3 > 0:
            tn93 = (
                -k1 * math.log(term1)
                - k2 * math.log(term2)
                - (k3 - k1 * fy - k2 * fr) * math.log(term3)
            )
            if tn93 < 0.0:
                tn93 = p_dist

    mut_str = ",".join(muts) if len(muts) <= max_mut_list else f"{','.join(muts[:max_mut_list])},...(+{snps - max_mut_list})"
    return snps, round(tn93, 8), mut_str


def _adjusted_rand_index(labels_true: list[str], labels_pred: list[int]) -> float:
    """Compute Adjusted Rand Index (Hubert & Arabie 1985) in pure Python."""
    n = len(labels_true)
    if n < 2:
        return 1.0

    contingency: dict[tuple[str, int], int] = defaultdict(int)
    a_counts: dict[str, int] = defaultdict(int)
    b_counts: dict[int, int] = defaultdict(int)
    for lt, lp in zip(labels_true, labels_pred):
        contingency[(lt, lp)] += 1
        a_counts[lt] += 1
        b_counts[lp] += 1

    def comb2(k: int) -> int:
        return k * (k - 1) // 2

    sum_nij = sum(comb2(v) for v in contingency.values())
    sum_ai = sum(comb2(v) for v in a_counts.values())
    sum_bj = sum(comb2(v) for v in b_counts.values())
    total_pairs = comb2(n)
    if total_pairs == 0:
        return 1.0

    expected = (sum_ai * sum_bj) / total_pairs
    max_idx = 0.5 * (sum_ai + sum_bj)
    denom = max_idx - expected
    if abs(denom) < 1e-12:
        return 1.0 if sum_nij == sum_ai == sum_bj else 0.0
    return float((sum_nij - expected) / denom)


def compute_snp_threshold_concordance(
    node_ids: list[str],
    seqs: dict[str, str],
    true_components: dict[str, str],
    thresholds: Iterable[int] = (0, 1, 2, 3, 5),
) -> list[dict]:
    """Evaluate MicrobeTrace single-linkage SNP clustering against true chains.

    For each SNP threshold ``tau``, connects any pair of sequenced nodes with
    ``snp_distance <= tau`` and compares the resulting connected components and
    same-cluster node pairs against ``true_component_id``.
    """
    present = [nid for nid in node_ids if nid in seqs and true_components.get(nid) not in (None, "", "nan", "-1")]
    n = len(present)
    if n < 2:
        return []

    # Precompute pairwise SNP distances for all pairs (O(N^2) on N <= 1500).
    pair_snps: list[tuple[int, int, int]] = []
    same_chain_pairs = 0
    total_pairs = n * (n - 1) // 2

    for i in range(n):
        si = seqs[present[i]]
        ci = true_components[present[i]]
        for j in range(i + 1, n):
            sj = seqs[present[j]]
            cj = true_components[present[j]]
            if ci == cj:
                same_chain_pairs += 1
            # Fast inline ACGT mismatch count
            d = sum(1 for a, b in zip(si, sj) if a != b and a in "ACGT" and b in "ACGT")
            pair_snps.append((i, j, d))

    true_labels = [str(true_components[nid]) for nid in present]
    results = []

    for tau in thresholds:
        # Union-Find for single-linkage clustering at distance <= tau
        parent = list(range(n))

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x: int, y: int) -> None:
            rx, ry = find(x), find(y)
            if rx != ry:
                if rx < ry:
                    parent[ry] = rx
                else:
                    parent[rx] = ry

        direct_links = 0
        true_direct_links = 0
        for i, j, d in pair_snps:
            if d <= tau:
                direct_links += 1
                if true_labels[i] == true_labels[j]:
                    true_direct_links += 1
                union(i, j)

        pred_labels = [find(i) for i in range(n)]
        cluster_sizes = Counter(pred_labels)
        n_clusters = len(cluster_sizes)
        n_multi_clusters = sum(1 for sz in cluster_sizes.values() if sz >= 2)

        # Pairwise cluster concordance (co-clustering in predicted vs true component)
        tp = 0
        fp = 0
        fn = 0
        for i in range(n):
            pi = pred_labels[i]
            ti = true_labels[i]
            for j in range(i + 1, n):
                same_pred = (pi == pred_labels[j])
                same_true = (ti == true_labels[j])
                if same_pred and same_true:
                    tp += 1
                elif same_pred and not same_true:
                    fp += 1
                elif not same_pred and same_true:
                    fn += 1

        precision = tp / (tp + fp) if (tp + fp) > 0 else 1.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = (2.0 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0
        link_precision = true_direct_links / direct_links if direct_links > 0 else 1.0
        ari = _adjusted_rand_index(true_labels, pred_labels)

        results.append({
            "snp_threshold": int(tau),
            "n_sequences": n,
            "n_genetic_links": direct_links,
            "link_same_chain_precision": round(link_precision, 4),
            "n_inferred_clusters": n_clusters,
            "n_multicase_clusters": n_multi_clusters,
            "largest_cluster_size": max(cluster_sizes.values()) if cluster_sizes else 0,
            "pairwise_precision": round(precision, 4),
            "pairwise_recall": round(recall, 4),
            "pairwise_f1": round(f1, 4),
            "adjusted_rand_index": round(ari, 4),
        })

    return results


def _load_candidate_samples(
    samples_path: Path | None,
) -> tuple[dict[str, dict[str, str]], set[str]]:
    """Read a sampled or ascertained line-list CSV into ``{alias_pid: row_dict}``."""
    if samples_path is None:
        return {}, set()

    by_alias: dict[str, dict[str, str]] = {}
    by_strain: set[str] = set()
    with _open_text(samples_path, "rt") as fh:
        first = fh.readline()
        if not first:
            return {}, set()
        delim = _sniff_delimiter(first)
        fh.seek(0)
        reader = csv.DictReader(fh, delimiter=delim)
        if not reader.fieldnames:
            return {}, set()
        reader.fieldnames = [c.strip() for c in reader.fieldnames]
        for row in reader:
            alias = (row.get("alias_pid") or "").strip()
            strain = (row.get("strain") or "").strip()
            if not alias and strain:
                alias = _alias_from_strain(strain) or ""
            if strain:
                by_strain.add(strain)
            if alias and alias not in by_alias:
                by_alias[alias] = row
    return by_alias, by_strain


def _scan_infections_graph(
    infections_path: Path | None,
    county_filter: str | None = None,
    date_min: _dt.date | None = None,
    date_max: _dt.date | None = None,
    variant_filter: str | None = None,
) -> tuple[
    dict[str, str],          # parent_of[alias_pid] = alias_contact
    dict[str, str],          # comp_of[alias_pid] = component_id
    Counter,                 # comp_total_size[component_id]
    dict[str, str],          # strain_of[alias_pid] = strain
    dict[str, int],          # tick_of[alias_pid] = exposure_tick
    dict[str, str],          # variant_of[alias_pid] = variant_benchmark
    set[str],                # filtered_aliases matching criteria
]:
    """Stream the ABM all-events table once to build the lightweight parent tree.

    Rather than holding all 46 string columns in memory for 2.8M-5.3M rows, this
    first pass extracts only the compact transmission topology and component IDs,
    then a second targeted pass loads full metadata rows only for the chosen
    subset (~600 rows).
    """
    parent_of: dict[str, str] = {}
    comp_of: dict[str, str] = {}
    comp_total_size: Counter = Counter()
    strain_of: dict[str, str] = {}
    tick_of: dict[str, int] = {}
    variant_of: dict[str, str] = {}
    filtered_aliases: set[str] = set()

    if infections_path is None or not infections_path.exists():
        return parent_of, comp_of, comp_total_size, strain_of, tick_of, variant_of, filtered_aliases

    county_q = county_filter.strip().lower() if county_filter else None
    variant_q = variant_filter.strip().lower() if variant_filter else None
    has_filters = bool(county_q or date_min or date_max or variant_q)

    with _open_text(infections_path, "rt") as fh:
        first = fh.readline()
        if not first:
            return parent_of, comp_of, comp_total_size, strain_of, tick_of, variant_of, filtered_aliases
        delim = _sniff_delimiter(first)
        cols = [c.strip() for c in first.rstrip("\r\n").split(delim)]
        col_idx = {c: i for i, c in enumerate(cols)}

        i_alias = col_idx.get("alias_pid")
        i_contact = col_idx.get("alias_contact")
        i_comp = col_idx.get("component_id")
        i_strain = col_idx.get("strain")
        i_tick = col_idx.get("exposure_tick", col_idx.get("sim_tick"))
        i_var = col_idx.get("variant_benchmark", col_idx.get("variant_label"))
        i_county = col_idx.get("county")
        i_fips = col_idx.get("county_fips", col_idx.get("FIPS"))
        i_date = col_idx.get("date", col_idx.get("collection_date", col_idx.get("exposure_date")))

        if i_alias is None:
            raise ValueError(f"Infections file {infections_path} has no 'alias_pid' column.")

        reader = csv.reader(fh, delimiter=delim)
        for row in reader:
            if len(row) <= i_alias:
                continue
            alias = row[i_alias].strip()
            if not alias or alias in parent_of:
                # De-duplicate to first row per alias_pid
                continue
            contact = row[i_contact].strip() if i_contact is not None and len(row) > i_contact else ""
            parent_of[alias] = contact if contact not in _NO_CONTACT and contact != alias else ""

            if i_comp is not None and len(row) > i_comp:
                comp = row[i_comp].strip()
                if comp and comp not in _NO_CONTACT:
                    comp_of[alias] = comp
                    comp_total_size[comp] += 1

            if i_strain is not None and len(row) > i_strain:
                st = row[i_strain].strip()
                if st:
                    strain_of[alias] = st

            if i_tick is not None and len(row) > i_tick:
                t_raw = row[i_tick].strip()
                if t_raw:
                    try:
                        tick_of[alias] = int(float(t_raw))
                    except ValueError:
                        pass

            v_raw = ""
            if i_var is not None and len(row) > i_var:
                v_raw = row[i_var].strip()
                if v_raw and v_raw not in _NO_CONTACT:
                    variant_of[alias] = v_raw

            if has_filters:
                if county_q:
                    c_val = row[i_county].strip().lower() if i_county is not None and len(row) > i_county else ""
                    f_val = row[i_fips].strip().lower() if i_fips is not None and len(row) > i_fips else ""
                    if county_q not in c_val and county_q != f_val:
                        continue
                if date_min or date_max:
                    d_str = row[i_date].strip() if i_date is not None and len(row) > i_date else ""
                    dv = _parse_iso_date(d_str)
                    if dv is None or (date_min and dv < date_min) or (date_max and dv > date_max):
                        continue
                if variant_q:
                    vv = v_raw.lower()
                    if vv != variant_q and not vv.startswith(variant_q):
                        continue
                filtered_aliases.add(alias)

    return parent_of, comp_of, comp_total_size, strain_of, tick_of, variant_of, filtered_aliases


def _find_nearest_subset_ancestor(
    alias: str,
    subset_aliases: set[str],
    parent_of: dict[str, str],
    max_hops: int,
) -> tuple[str | None, int | None]:
    """Walk parent pointers from ``alias`` up to ``max_hops`` steps to find the
    nearest ancestor present in ``subset_aliases``."""
    cur = alias
    visited = {alias}
    for hop in range(1, max_hops + 1):
        parent = parent_of.get(cur, "")
        if not parent or parent in visited:
            return None, None
        if parent in subset_aliases:
            return parent, hop
        visited.add(parent)
        cur = parent
    return None, None


def _build_k_hop_clusters(
    candidate_aliases: set[str],
    parent_of: dict[str, str],
    max_hops: int,
) -> list[list[str]]:
    """Partition ``candidate_aliases`` into connected components under <= max_hops
    ancestral links, sorted largest first, with members ordered in root-first BFS
    order so any prefix slice ``comp[:k]`` remains a connected subtree."""
    adj: dict[str, list[str]] = defaultdict(list)
    children_in_sub: dict[str, list[str]] = defaultdict(list)
    has_parent_in_sub: set[str] = set()
    for alias in sorted(candidate_aliases):
        anc, _ = _find_nearest_subset_ancestor(alias, candidate_aliases, parent_of, max_hops)
        if anc is not None:
            adj[alias].append(anc)
            adj[anc].append(alias)
            children_in_sub[anc].append(alias)
            has_parent_in_sub.add(alias)

    visited: set[str] = set()
    components: list[list[str]] = []
    for alias in sorted(candidate_aliases):
        if alias in visited:
            continue
        raw_comp = []
        queue = [alias]
        visited.add(alias)
        while queue:
            u = queue.pop()
            raw_comp.append(u)
            for v in adj.get(u, ()):
                if v not in visited:
                    visited.add(v)
                    queue.append(v)

        if len(raw_comp) <= 2:
            components.append(sorted(raw_comp))
            continue

        # Order raw_comp via BFS starting from its root(s) (nodes without an
        # ancestor in raw_comp) so prefix truncation preserves connectivity.
        comp_set = set(raw_comp)
        roots = sorted(u for u in raw_comp if u not in has_parent_in_sub)
        if not roots:
            roots = [min(raw_comp)]
        bfs_order: list[str] = []
        seen_bfs: set[str] = set()
        for r in roots:
            if r not in seen_bfs:
                seen_bfs.add(r)
                bfs_order.append(r)
        idx = 0
        while idx < len(bfs_order):
            curr = bfs_order[idx]
            idx += 1
            for ch in sorted(children_in_sub.get(curr, ())):
                if ch in comp_set and ch not in seen_bfs:
                    seen_bfs.add(ch)
                    bfs_order.append(ch)
        for u in sorted(raw_comp):
            if u not in seen_bfs:
                bfs_order.append(u)
        components.append(bfs_order)

    components.sort(key=lambda c: (-len(c), c[0]))
    return components


def _filter_rows_by_criteria(
    rows_by_alias: dict[str, dict[str, str]],
    comp_of: dict[str, str],
    variant_of: dict[str, str],
    county_filter: str | None = None,
    date_min: _dt.date | None = None,
    date_max: _dt.date | None = None,
    variant_filter: str | None = None,
) -> set[str]:
    """Filter candidate aliases by optional county, date range, and variant."""
    county_q = county_filter.strip().lower() if county_filter else None
    variant_q = variant_filter.strip().lower() if variant_filter else None

    kept: set[str] = set()
    for alias, row in rows_by_alias.items():
        if county_q:
            c_val = (row.get("county") or "").strip().lower()
            fips_val = (row.get("county_fips") or row.get("FIPS") or "").strip().lower()
            if county_q not in c_val and county_q != fips_val:
                continue
        if date_min or date_max:
            d_str = row.get("date") or row.get("collection_date") or row.get("exposure_date") or ""
            d_val = _parse_iso_date(d_str)
            if d_val is None:
                continue
            if date_min and d_val < date_min:
                continue
            if date_max and d_val > date_max:
                continue
        if variant_q:
            v_val = (
                row.get("variant_benchmark")
                or variant_of.get(alias)
                or row.get("variant_label")
                or ""
            ).strip().lower()
            if v_val != variant_q and not v_val.startswith(variant_q):
                continue
        kept.add(alias)
    return kept


def _preprune_large_candidate_pool(
    candidate_aliases: set[str],
    comp_of: dict[str, str],
    target_cap: int,
    rng: random.Random,
) -> set[str]:
    """When ``candidate_aliases`` is very large (>50,000, e.g. an entire multi-million
    infection table), pre-select top components plus a random background sample before
    running fine-grained subtree/cluster selection."""
    if len(candidate_aliases) <= 50_000:
        return candidate_aliases

    comp_counts: Counter = Counter()
    for al in candidate_aliases:
        cid = comp_of.get(al, "")
        if cid:
            comp_counts[cid] += 1

    # Pick medium/large components until we have ~10x target_cap (at least 15,000)
    keep_budget = max(target_cap * 10, 15_000)
    ranked_comps = [cid for cid, cnt in comp_counts.most_common() if cnt >= 2]
    chosen_comps: set[str] = set()
    accum = 0
    for cid in ranked_comps:
        chosen_comps.add(cid)
        accum += comp_counts[cid]
        if accum >= keep_budget:
            break

    pruned = {al for al in candidate_aliases if comp_of.get(al, "") in chosen_comps}
    if len(pruned) < keep_budget:
        rest = [al for al in candidate_aliases if al not in pruned]
        extra = rng.sample(rest, min(len(rest), keep_budget - len(pruned)))
        pruned.update(extra)
    return pruned


def _select_subset_aliases(
    mode: str,
    candidate_aliases: set[str],
    rows_by_alias: dict[str, dict[str, str]],
    parent_of: dict[str, str],
    comp_of: dict[str, str],
    max_sequences: int,
    max_hops: int,
    background_fraction: float,
    min_chain_size: int,
    max_chain_size: int,
    seed: int,
) -> list[str]:
    """Select the final list of ``alias_pid`` records according to ``mode``."""
    rng = random.Random(seed)

    if mode == "all" or max_sequences <= 0 or len(candidate_aliases) <= max_sequences:
        if mode != "complete-chains" or not parent_of:
            return sorted(candidate_aliases)

    if mode == "random":
        pool = sorted(candidate_aliases)
        return sorted(rng.sample(pool, min(max_sequences, len(pool))))

    if len(candidate_aliases) > 50_000 and mode in ("cluster-enriched", "focal-window"):
        candidate_aliases = _preprune_large_candidate_pool(
            candidate_aliases, comp_of, max_sequences, rng
        )

    if mode == "focal-window":
        # If multiple counties or a wide date window remain, focus on the county
        # and 56-day window with the highest density of multi-case transmission chains.
        clusters = _build_k_hop_clusters(candidate_aliases, parent_of, max_hops)
        multi_nodes = {u for c in clusters if len(c) >= 2 for u in c}
        county_counts: Counter = Counter()
        for u in (multi_nodes or candidate_aliases):
            cty = (rows_by_alias.get(u, {}).get("county") or "").strip()
            if cty:
                county_counts[cty] += 1
        focal_county = county_counts.most_common(1)[0][0] if county_counts else None
        if focal_county:
            focal_pool = [
                u for u in candidate_aliases
                if (rows_by_alias.get(u, {}).get("county") or "").strip() == focal_county
            ]
        else:
            focal_pool = sorted(candidate_aliases)

        if len(focal_pool) <= max_sequences:
            return sorted(focal_pool)

        # Slide an 8-week (56-day) window over focal_pool to find the densest peak
        dated: list[tuple[_dt.date, str]] = []
        for u in focal_pool:
            r = rows_by_alias.get(u, {})
            d = _parse_iso_date(r.get("date") or r.get("exposure_date"))
            if d is not None:
                dated.append((d, u))
        dated.sort()
        if dated:
            best_slice: list[str] = []
            j = 0
            for i in range(len(dated)):
                while j < len(dated) and (dated[j][0] - dated[i][0]).days <= 56:
                    j += 1
                window_aliases = [dated[k][1] for k in range(i, j)]
                if len(window_aliases) > len(best_slice):
                    best_slice = window_aliases
            if len(best_slice) > max_sequences:
                # Fall through to cluster-enriched inside best_slice
                candidate_aliases = set(best_slice)
            elif best_slice:
                return sorted(best_slice)

    if mode == "complete-chains":
        # Build directed children adjacency on the ground-truth graph restricted to
        # candidate_aliases (or all of parent_of if candidate_aliases is full truth)
        children_of: dict[str, list[str]] = defaultdict(list)
        pool_set = candidate_aliases if candidate_aliases else set(parent_of.keys())
        for u in pool_set:
            p = parent_of.get(u, "")
            if p and p in pool_set:
                children_of[p].append(u)

        # Extract contiguous BFS subtrees of size [min_chain_size, max_chain_size]
        # Prefer roots that have many descendants in pool_set
        roots = [u for u in pool_set if not parent_of.get(u) or parent_of.get(u) not in pool_set]
        roots.sort(key=lambda r: (-len(children_of.get(r, ())), r))

        selected: list[str] = []
        selected_set: set[str] = set()

        # First pass: collect medium-sized complete components or BFS subtrees
        by_comp: dict[str, list[str]] = defaultdict(list)
        for u in pool_set:
            cid = comp_of.get(u, "")
            if cid:
                by_comp[cid].append(u)

        eligible_comps = [
            (cid, members)
            for cid, members in by_comp.items()
            if min_chain_size <= len(members) <= max_chain_size
        ]
        # Sort by number of direct edges inside members
        def _internal_edges(members: list[str]) -> int:
            mset = set(members)
            return sum(1 for u in members if parent_of.get(u, "") in mset)

        eligible_comps.sort(key=lambda item: (-_internal_edges(item[1]), -len(item[1]), item[0]))

        for _, members in eligible_comps:
            if len(selected) + len(members) > max_sequences and len(selected) >= min_chain_size:
                continue
            for u in sorted(members):
                if u not in selected_set and len(selected) < max_sequences:
                    selected.append(u)
                    selected_set.add(u)
            if len(selected) >= max_sequences:
                break

        # If we still have room (e.g. components are huge > max_chain_size),
        # extract contiguous BFS subtrees of size up to max_chain_size
        if len(selected) < max_sequences:
            large_comps = [
                (cid, members)
                for cid, members in by_comp.items()
                if len(members) > max_chain_size
            ]
            large_comps.sort(key=lambda item: (-len(item[1]), item[0]))
            for _, members in large_comps:
                if len(selected) >= max_sequences:
                    break
                mset = set(members)
                comp_roots = [u for u in members if parent_of.get(u, "") not in mset]
                start_node = comp_roots[0] if comp_roots else sorted(members)[0]
                bfs_queue = [start_node]
                subtree: list[str] = []
                seen_sub = {start_node}
                idx_q = 0
                budget = min(max_chain_size, max_sequences - len(selected))
                while idx_q < len(bfs_queue) and len(subtree) < budget:
                    curr = bfs_queue[idx_q]
                    idx_q += 1
                    if curr not in selected_set:
                        subtree.append(curr)
                    for ch in sorted(children_of.get(curr, ())):
                        if ch in mset and ch not in seen_sub:
                            seen_sub.add(ch)
                            bfs_queue.append(ch)
                if len(subtree) >= min_chain_size or not selected:
                    for u in subtree:
                        if u not in selected_set and len(selected) < max_sequences:
                            selected.append(u)
                            selected_set.add(u)

        return sorted(selected)

    # Default: "cluster-enriched" (also used as final ranker for focal-window)
    # 1. Partition candidate_aliases into k-hop connected clusters in the true graph
    k_clusters = _build_k_hop_clusters(candidate_aliases, parent_of, max_hops)
    multi_clusters = [c for c in k_clusters if len(c) >= 2]
    singletons = [c[0] for c in k_clusters if len(c) == 1]

    # If k-hop clusters are sparse (e.g. parent_of was not provided or sampling is
    # very sparse), also group remaining singletons by true_component_id so same-chain
    # cases stay together.
    if multi_clusters:
        clustered_Flat = [u for c in multi_clusters for u in c]
    else:
        clustered_Flat = []

    clustered_set = set(clustered_Flat)
    remaining_by_comp: dict[str, list[str]] = defaultdict(list)
    true_singletons: list[str] = []
    for u in singletons:
        cid = comp_of.get(u) or (rows_by_alias.get(u, {}).get("component_id") or "").strip()
        if cid and cid not in _NO_CONTACT:
            remaining_by_comp[cid].append(u)
        else:
            true_singletons.append(u)

    comp_groups = sorted(
        [sorted(v) for v in remaining_by_comp.values() if len(v) >= 2],
        key=lambda g: (-len(g), g[0]),
    )
    comp_singletons = [v[0] for v in remaining_by_comp.values() if len(v) == 1] + true_singletons

    bg_Quota = int(round(max_sequences * max(0.0, min(1.0, background_fraction))))
    cluster_quota = max_sequences - bg_Quota

    selected: list[str] = []
    selected_set: set[str] = set()

    # Pass 1: Fill cluster quota across distinct k-hop connected subtrees and
    # multi-sample components, capping each cluster at `max_chain_size` (default 80)
    # as a root-first BFS subtree so a single giant lineage cannot crowd out
    # co-circulating outbreak chains.
    per_cluster_cap = max(8, max_chain_size) if max_chain_size > 0 else cluster_quota
    all_groups = multi_clusters + comp_groups
    for grp in all_groups:
        if len(selected) >= cluster_quota:
            break
        added_from_grp = 0
        for u in grp:
            if added_from_grp >= per_cluster_cap or len(selected) >= cluster_quota:
                break
            if u not in selected_set:
                selected.append(u)
                selected_set.add(u)
                added_from_grp += 1

    # Pass 2: If distinct clusters did not fill `cluster_quota`, top up from
    # remaining nodes of larger clusters before drawing background singletons.
    if len(selected) < cluster_quota:
        for grp in all_groups:
            if len(selected) >= cluster_quota:
                break
            for u in grp:
                if len(selected) >= cluster_quota:
                    break
                if u not in selected_set:
                    selected.append(u)
                    selected_set.add(u)

    # Fill remaining capacity from background singletons (or unused cluster nodes)
    bg_pool = [u for u in comp_singletons if u not in selected_set]
    rng.shuffle(bg_pool)
    for u in bg_pool:
        if len(selected) >= max_sequences:
            break
        selected.append(u)
        selected_set.add(u)

    # If still below max_sequences, top up with any remaining cluster nodes
    if len(selected) < max_sequences:
        leftover = [u for u in candidate_aliases if u not in selected_set]
        for u in sorted(leftover):
            if len(selected) >= max_sequences:
                break
            selected.append(u)
            selected_set.add(u)

    return sorted(selected)


def _fetch_full_metadata_for_subset(
    subset_aliases: set[str],
    samples_by_alias: dict[str, dict[str, str]],
    infections_path: Path | None,
) -> dict[str, dict[str, str]]:
    """Retrieve complete row dictionaries for ``subset_aliases`` from samples and/or
    the all-events infections file."""
    merged: dict[str, dict[str, str]] = {}
    for alias in subset_aliases:
        if alias in samples_by_alias:
            merged[alias] = dict(samples_by_alias[alias])

    # Enrich or backfill from infections_path if provided
    if infections_path is not None and infections_path.exists():
        with _open_text(infections_path, "rt") as fh:
            first = fh.readline()
            if first:
                delim = _sniff_delimiter(first)
                fh.seek(0)
                reader = csv.DictReader(fh, delimiter=delim)
                if reader.fieldnames:
                    reader.fieldnames = [c.strip() for c in reader.fieldnames]
                    for row in reader:
                        alias = (row.get("alias_pid") or "").strip()
                        if alias not in subset_aliases:
                            continue
                        if alias not in merged:
                            merged[alias] = dict(row)
                        else:
                            # Backfill any missing fields from allevents (e.g. exposure_date, variant_benchmark)
                            cur = merged[alias]
                            for k, v in row.items():
                                if k and (k not in cur or cur[k] in ("", None, "nan", "unassigned")) and v not in ("", None):
                                    cur[k] = v
    return merged


def _extract_sequences_for_subset(
    fasta_path: Path | None,
    subset_aliases: set[str],
    strain_by_alias: dict[str, str],
) -> dict[str, str]:
    """Stream a FASTA (.fasta, .fasta.xz, .fasta.gz) and return ``{alias_pid: seq}``."""
    if fasta_path is None or not fasta_path.exists():
        return {}

    alias_by_strain = {st: al for al, st in strain_by_alias.items() if st}
    seqs_by_alias: dict[str, str] = {}

    with _open_text(fasta_path, "rt") as fh:
        cur_alias: str | None = None
        chunks: list[str] = []
        for line in fh:
            if line.startswith(">"):
                if cur_alias is not None and chunks:
                    seqs_by_alias[cur_alias] = "".join(chunks)
                header_id = line[1:].strip().split()[0]
                if header_id in _DEFAULT_CONTEXT_IDS:
                    cur_alias = None
                    chunks = []
                    continue
                matched_alias = alias_by_strain.get(header_id)
                if not matched_alias:
                    if header_id in subset_aliases:
                        matched_alias = header_id
                    else:
                        extracted = _alias_from_strain(header_id)
                        if extracted and extracted in subset_aliases:
                            matched_alias = extracted
                if matched_alias and matched_alias in subset_aliases:
                    cur_alias = matched_alias
                    chunks = []
                else:
                    cur_alias = None
                    chunks = []
            elif cur_alias is not None:
                chunks.append(line.strip())
        if cur_alias is not None and chunks:
            seqs_by_alias[cur_alias] = "".join(chunks)

    return seqs_by_alias


def export_microbetrace_bundle(
    outdir: str | Path,
    prefix: str = "microbetrace",
    fasta_path: str | Path | None = None,
    samples_path: str | Path | None = None,
    infections_path: str | Path | None = None,
    selection_mode: str = "cluster-enriched",
    max_sequences: int = 600,
    max_hops: int = 5,
    background_fraction: float = 0.20,
    min_chain_size: int = 8,
    max_chain_size: int = 80,
    county_filter: str | None = None,
    date_min: str | None = None,
    date_max: str | None = None,
    variant_filter: str | None = None,
    id_column: str = "strain",
    seed: int = 43,
    allow_no_fasta: bool = False,
    compute_concordance: bool = True,
) -> dict:
    """Build and write a MicrobeTrace bundle (`_sequences.fasta`, `_nodes.csv`,
    `_truth_links.csv`, `_summary.json`)."""
    out_dir = Path(outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    f_path = Path(fasta_path) if fasta_path else None
    s_path = Path(samples_path) if samples_path else None
    i_path = Path(infections_path) if infections_path else None

    if s_path is None and i_path is None:
        raise ValueError("At least one of --samples/--metadata or --infections must be provided.")
    if f_path is not None and not f_path.exists() and not allow_no_fasta:
        raise FileNotFoundError(
            f"FASTA file not found: {f_path}\n"
            f"Pass --allow-no-fasta to export nodes and truth_links without sequences."
        )

    print(f"MicrobeTrace export ({selection_mode}, max_sequences={max_sequences}, max_hops={max_hops})")
    if s_path:
        print(f"  Samples / metadata : {s_path}")
    if i_path:
        print(f"  Ground-truth graph : {i_path}")
    if f_path:
        print(f"  Sequences FASTA    : {f_path}")

    # 1. Load candidate sample rows (if provided)
    samples_by_alias, _ = _load_candidate_samples(s_path) if s_path and s_path.exists() else ({}, set())
    sampled_surveillance_set = set(samples_by_alias.keys())
    if samples_by_alias:
        print(f"  Loaded {len(samples_by_alias):,} candidate sample records.")

    # 2. Scan ground-truth infections graph (if provided)
    d_min = _parse_iso_date(date_min)
    d_max = _parse_iso_date(date_max)
    has_filters = bool(county_filter or d_min or d_max or variant_filter)

    parent_of, comp_of, comp_total_size, strain_of, tick_of, variant_of, inf_filtered = (
        _scan_infections_graph(
            i_path,
            county_filter=county_filter if (not samples_by_alias or selection_mode == "complete-chains") else None,
            date_min=d_min if (not samples_by_alias or selection_mode == "complete-chains") else None,
            date_max=d_max if (not samples_by_alias or selection_mode == "complete-chains") else None,
            variant_filter=variant_filter if (not samples_by_alias or selection_mode == "complete-chains") else None,
        )
    )
    if parent_of:
        print(f"  Scanned {len(parent_of):,} ground-truth infections across {len(comp_total_size):,} components.")
    else:
        # Build fallback parent_of / comp_of / strain_of from samples_by_alias
        for alias, row in samples_by_alias.items():
            con = (row.get("alias_contact") or "").strip()
            parent_of[alias] = con if con not in _NO_CONTACT and con != alias else ""
            comp = (row.get("component_id") or "").strip()
            if comp and comp not in _NO_CONTACT:
                comp_of[alias] = comp
                comp_total_size[comp] += 1
            st = (row.get("strain") or "").strip()
            if st:
                strain_of[alias] = st
            vb = (row.get("variant_benchmark") or row.get("variant_label") or "").strip()
            if vb and vb not in _NO_CONTACT:
                variant_of[alias] = vb

    # 3. Determine initial candidate pool
    if selection_mode == "complete-chains" and i_path and i_path.exists():
        candidate_pool = inf_filtered if has_filters else set(parent_of.keys())
    elif samples_by_alias:
        candidate_pool = _filter_rows_by_criteria(
            samples_by_alias,
            comp_of=comp_of,
            variant_of=variant_of,
            county_filter=county_filter,
            date_min=d_min,
            date_max=d_max,
            variant_filter=variant_filter,
        )
    else:
        # Only infections_path was given (no separate samples CSV)
        candidate_pool = inf_filtered if has_filters else set(parent_of.keys())

    if not candidate_pool:
        raise ValueError("No records matched the requested filters.")

    # 4. Select the target subset of aliases
    selected_aliases = _select_subset_aliases(
        mode=selection_mode,
        candidate_aliases=candidate_pool,
        rows_by_alias=samples_by_alias,
        parent_of=parent_of,
        comp_of=comp_of,
        max_sequences=max_sequences,
        max_hops=max_hops,
        background_fraction=background_fraction,
        min_chain_size=min_chain_size,
        max_chain_size=max_chain_size,
        seed=seed,
    )
    selected_set = set(selected_aliases)
    print(f"  Selected {len(selected_aliases):,} records using mode '{selection_mode}'.")

    # 5. Fetch full metadata rows for the selected subset
    full_rows = _fetch_full_metadata_for_subset(selected_set, samples_by_alias, i_path)

    # Ensure strain_of is populated for every selected alias
    for al in selected_aliases:
        if al not in strain_of:
            st = (full_rows.get(al, {}).get("strain") or "").strip()
            if st:
                strain_of[al] = st

    # 6. Extract sequences from FASTA (if provided)
    seqs_by_alias: dict[str, str] = {}
    if f_path is not None and f_path.exists():
        seqs_by_alias = _extract_sequences_for_subset(f_path, selected_set, strain_of)
        print(f"  Extracted {len(seqs_by_alias):,} / {len(selected_aliases):,} sequences from {f_path.name}.")
        if seqs_by_alias and len(seqs_by_alias) < len(selected_aliases):
            # Restrict selected_aliases to those present in the FASTA so MicrobeTrace
            # receives a 1:1 matched node/sequence set.
            missing_n = len(selected_aliases) - len(seqs_by_alias)
            print(f"  Note: restricting node & link tables to the {len(seqs_by_alias):,} records present in FASTA "
                  f"({missing_n:,} records had no sequence in {f_path.name}).")
            selected_aliases = [al for al in selected_aliases if al in seqs_by_alias]
            selected_set = set(selected_aliases)

    # 7. Map each selected alias to its exported primary `id`
    id_of_alias: dict[str, str] = {}
    for al in selected_aliases:
        row = full_rows.get(al, {})
        if id_column == "alias_pid":
            id_of_alias[al] = al
        else:
            st = (row.get(id_column) or row.get("strain") or strain_of.get(al) or al).strip()
            id_of_alias[al] = st

    # Count subset component sizes
    # Count subset component sizes and true secondary cases (out-degree in full graph)
    subset_comp_counts: Counter = Counter()
    for al in selected_aliases:
        cid = comp_of.get(al) or (full_rows.get(al, {}).get("component_id") or "").strip()
        if cid and cid not in _NO_CONTACT:
            subset_comp_counts[cid] += 1

    true_out_degree: Counter = Counter()
    for _, par_al in parent_of.items():
        if par_al and par_al in selected_set:
            true_out_degree[par_al] += 1

    # Helper to classify shared person-to-place / household exposure setting
    def _classify_exposure_setting(src_row: dict, tgt_row: dict, same_county: bool) -> tuple[str, bool]:
        s_hid = (src_row.get("hid") or src_row.get("household_id") or "").strip()
        t_hid = (tgt_row.get("hid") or tgt_row.get("household_id") or "").strip()

        same_hh = bool(s_hid and s_hid not in _NO_CONTACT and s_hid == t_hid)
        if same_hh:
            return "Household", True

        s_rlid = (src_row.get("rlid") or src_row.get("residence_location_id") or "").strip()
        t_rlid = (tgt_row.get("rlid") or tgt_row.get("residence_location_id") or "").strip()
        if s_rlid and s_rlid not in _NO_CONTACT and s_rlid == t_rlid:
            return "Residence Location", False

        s_tract = (src_row.get("admin3") or src_row.get("census_tract") or "").strip()
        t_tract = (tgt_row.get("admin3") or tgt_row.get("census_tract") or "").strip()
        if same_county and s_tract and s_tract not in _NO_CONTACT and s_tract == t_tract:
            return "Same Census Tract", False

        if same_county:
            return "Community (Same County)", False
        return "Cross-County Import", False

    # 8. Build ground-truth links (direct + ancestral bridge edges within max_hops,
    #    plus coincidental household/location co-exposures from independent chains)
    truth_links: list[dict] = []
    ancestor_info_by_alias: dict[str, tuple[str, int, int | str]] = {}
    subset_out_degree: Counter = Counter()
    linked_pairs: set[tuple[str, str]] = set()

    for target_al in selected_aliases:
        anc_al, hops = _find_nearest_subset_ancestor(target_al, selected_set, parent_of, max_hops)
        if anc_al is None or hops is None:
            continue

        subset_out_degree[anc_al] += 1
        linked_pairs.add((anc_al, target_al))
        linked_pairs.add((target_al, anc_al))

        src_id = id_of_alias[anc_al]
        tgt_id = id_of_alias[target_al]
        src_row = full_rows.get(anc_al, {})
        tgt_row = full_rows.get(target_al, {})

        snp_dist: int | str = ""
        tn93_dist: float | str = ""
        mut_str = ""
        if anc_al in seqs_by_alias and target_al in seqs_by_alias:
            snp_dist, tn93_dist, mut_str = _snp_and_mutations(
                seqs_by_alias[anc_al], seqs_by_alias[target_al]
            )

        ancestor_info_by_alias[target_al] = (src_id, hops, snp_dist)

        # Exposure and collection date intervals
        src_exp = _parse_iso_date(src_row.get("exposure_date"))
        tgt_exp = _parse_iso_date(tgt_row.get("exposure_date"))
        if src_exp is not None and tgt_exp is not None:
            delta_exp: int | str = (tgt_exp - src_exp).days
        elif anc_al in tick_of and target_al in tick_of:
            delta_exp = tick_of[target_al] - tick_of[anc_al]
        else:
            delta_exp = ""

        src_col = _parse_iso_date(src_row.get("date") or src_row.get("collection_date"))
        tgt_col = _parse_iso_date(tgt_row.get("date") or tgt_row.get("collection_date"))
        delta_col: int | str = (tgt_col - src_col).days if (src_col and tgt_col) else ""

        src_county = (src_row.get("county") or "").strip()
        tgt_county = (tgt_row.get("county") or "").strip()
        same_county = bool(src_county and tgt_county and src_county == tgt_county)
        exp_setting, same_hh = _classify_exposure_setting(src_row, tgt_row, same_county)

        cid = comp_of.get(target_al) or (tgt_row.get("component_id") or "").strip()
        var = (
            tgt_row.get("variant_benchmark")
            or variant_of.get(target_al)
            or tgt_row.get("variant_label")
            or ""
        ).strip()

        truth_links.append({
            "source": src_id,
            "target": tgt_id,
            "edge_type": "direct" if hops == 1 else "ancestral",
            "exposure_setting": exp_setting,
            "same_household": same_hh,
            "same_true_chain": True,
            "transmission_hops": hops,
            "unsampled_intermediates": hops - 1,
            "snp_distance": snp_dist,
            "tn93_distance": tn93_dist,
            "substitutions": mut_str,
            "delta_exposure_days": delta_exp,
            "delta_collection_days": delta_col,
            "same_county": same_county,
            "source_county": src_county,
            "target_county": tgt_county,
            "true_component_id": cid,
            "variant_benchmark": var,
        })

    # Detect coincidental household / residence-location co-exposures: pairs of cases
    # sharing an explicit household_id (hid) or residence_location_id (rlid) within 21 days
    # that belong to DIFFERENT ground-truth imported chains (same_true_chain = False).
    # In MicrobeTrace, contact tracing would flag these as household/place contacts,
    # whereas overlaying the 0-2 SNP threshold rules them out as independent introductions.
    by_place: dict[str, list[str]] = defaultdict(list)
    for al in selected_aliases:
        r = full_rows.get(al, {})
        hid_v = (r.get("hid") or r.get("household_id") or "").strip()
        rlid_v = (r.get("rlid") or r.get("residence_location_id") or "").strip()
        if hid_v and hid_v not in _NO_CONTACT:
            by_place[f"hid:{hid_v}"].append(al)
        elif rlid_v and rlid_v not in _NO_CONTACT:
            by_place[f"rlid:{rlid_v}"].append(al)

    for _, members in sorted(by_place.items()):
        if len(members) < 2 or len(members) > 15:
            continue
        for idx_a in range(len(members)):
            for idx_b in range(idx_a + 1, len(members)):
                u_al = members[idx_a]
                v_al = members[idx_b]
                if (u_al, v_al) in linked_pairs:
                    continue
                u_row = full_rows.get(u_al, {})
                v_row = full_rows.get(v_al, {})
                u_cid = comp_of.get(u_al) or (u_row.get("component_id") or "").strip()
                v_cid = comp_of.get(v_al) or (v_row.get("component_id") or "").strip()
                if not u_cid or not v_cid or u_cid == v_cid:
                    continue
                u_exp = _parse_iso_date(u_row.get("exposure_date") or u_row.get("date"))
                v_exp = _parse_iso_date(v_row.get("exposure_date") or v_row.get("date"))
                if u_exp and v_exp and abs((v_exp - u_exp).days) > 21:
                    continue
                if u_exp and v_exp and v_exp < u_exp:
                    u_al, v_al = v_al, u_al
                    u_row, v_row = v_row, u_row
                    u_cid, v_cid = v_cid, u_cid
                    u_exp, v_exp = v_exp, u_exp

                snp_d: int | str = ""
                tn93_d: float | str = ""
                mut_s = ""
                if u_al in seqs_by_alias and v_al in seqs_by_alias:
                    snp_d, tn93_d, mut_s = _snp_and_mutations(
                        seqs_by_alias[u_al], seqs_by_alias[v_al]
                    )

                u_col = _parse_iso_date(u_row.get("date") or u_row.get("collection_date"))
                v_col = _parse_iso_date(v_row.get("date") or v_row.get("collection_date"))
                d_exp: int | str = (v_exp - u_exp).days if (u_exp and v_exp) else ""
                d_col: int | str = (v_col - u_col).days if (u_col and v_col) else ""
                u_cty = (u_row.get("county") or "").strip()
                v_cty = (v_row.get("county") or "").strip()
                same_cty = bool(u_cty and v_cty and u_cty == v_cty)
                exp_set, same_h = _classify_exposure_setting(u_row, v_row, same_cty)

                linked_pairs.add((u_al, v_al))
                linked_pairs.add((v_al, u_al))
                truth_links.append({
                    "source": id_of_alias[u_al],
                    "target": id_of_alias[v_al],
                    "edge_type": "coincidental_exposure",
                    "exposure_setting": exp_set,
                    "same_household": same_h,
                    "same_true_chain": False,
                    "transmission_hops": "",
                    "unsampled_intermediates": "",
                    "snp_distance": snp_d,
                    "tn93_distance": tn93_d,
                    "substitutions": mut_s,
                    "delta_exposure_days": d_exp,
                    "delta_collection_days": d_col,
                    "same_county": same_cty,
                    "source_county": u_cty,
                    "target_county": v_cty,
                    "true_component_id": f"{u_cid} vs {v_cid}",
                    "variant_benchmark": (u_row.get("variant_benchmark") or variant_of.get(u_al) or "").strip(),
                })

    # 9. Build node rows
    node_rows: list[dict] = []
    true_comp_by_node_id: dict[str, str] = {}

    for al in selected_aliases:
        row = full_rows.get(al, {})
        nid = id_of_alias[al]
        strain = (row.get("strain") or strain_of.get(al) or nid).strip()
        cid = comp_of.get(al) or (row.get("component_id") or "").strip()
        true_comp_by_node_id[nid] = cid

        infector_al = parent_of.get(al, "") or (row.get("alias_contact") or "").strip()
        if infector_al in _NO_CONTACT:
            infector_al = ""
        infector_id = (
            id_of_alias.get(infector_al)
            or strain_of.get(infector_al)
            or infector_al
        ) if infector_al else ""
        direct_in_subset = bool(infector_al and infector_al in selected_set)

        anc_id, anc_hops, anc_snps = ancestor_info_by_alias.get(al, ("", "", ""))
        rucc_raw = (row.get("rucc_code") or "").strip()
        lat = (row.get("latitude") or row.get("home_latitude") or "").strip()
        lon = (row.get("longitude") or row.get("home_longitude") or "").strip()
        fips = (row.get("county_fips") or row.get("FIPS") or "").strip()
        if fips and fips.isdigit() and len(fips) < 5:
            fips = fips.zfill(5)

        var_b = (
            row.get("variant_benchmark")
            or variant_of.get(al)
            or row.get("variant_label")
            or ""
        ).strip()

        in_surv = (al in sampled_surveillance_set) if sampled_surveillance_set else True
        true_sec = int(true_out_degree.get(al, 0))
        sub_sec = int(subset_out_degree.get(al, 0))

        node_rows.append({
            "id": nid,
            "strain": strain,
            "alias_pid": al,
            "sim_pid": (row.get("sim_pid") or row.get("pid") or al.split(".")[0]).strip(),
            "collection_date": (row.get("date") or row.get("collection_date") or "").strip(),
            "exposure_date": (row.get("exposure_date") or "").strip(),
            "exposure_tick": (row.get("exposure_tick") or str(tick_of.get(al, ""))).strip(),
            "division": (row.get("division") or "").strip(),
            "county": (row.get("county") or "").strip(),
            "county_fips": fips,
            "census_tract": (row.get("admin3") or row.get("census_tract") or "").strip(),
            "household_id": (row.get("hid") or row.get("household_id") or "").strip(),
            "residence_location_id": (row.get("rlid") or row.get("residence_location_id") or "").strip(),
            "occupation_socp": (row.get("occupation_socp") or "").strip(),
            "latitude": lat,
            "longitude": lon,
            "rucc_code": rucc_raw,
            "urbanicity": _urbanicity_label(rucc_raw),
            "age_group": (row.get("age_group") or "").strip(),
            "sex": (row.get("sex") or row.get("gender") or "").strip(),
            "race": (row.get("race") or "").strip(),
            "smh_race": (row.get("smh_race") or "").strip(),
            "latino": (row.get("latino") or "").strip(),
            "ses_category": (row.get("ses_category") or "").strip(),
            "hh_size": (row.get("hh_size") or "").strip(),
            "hh_income": (row.get("hh_income") or "").strip(),
            "asymptomatic": (row.get("asymptomatic") or "").strip(),
            "exit_state": (row.get("exit_state") or "").strip(),
            "variant_benchmark": var_b,
            "true_component_id": cid,
            "true_component_size": comp_total_size.get(cid, "") if cid else "",
            "subset_component_size": subset_comp_counts.get(cid, 1) if cid else 1,
            "true_secondary_cases": true_sec,
            "subset_secondary_cases": sub_sec,
            "is_superspreader": true_sec >= 5,
            "sampled_in_surveillance": in_surv,
            "true_infector_alias": infector_al,
            "true_infector_id": infector_id,
            "direct_infector_in_subset": direct_in_subset,
            "nearest_subset_ancestor_id": anc_id,
            "hops_to_subset_ancestor": anc_hops,
            "snp_distance_to_ancestor": anc_snps,
        })

    # 10. Write files
    fasta_out = out_dir / f"{prefix}_sequences.fasta"
    nodes_out = out_dir / f"{prefix}_nodes.csv"
    links_out = out_dir / f"{prefix}_truth_links.csv"
    summary_out = out_dir / f"{prefix}_summary.json"

    if seqs_by_alias:
        with open(fasta_out, "w", encoding="utf-8") as fh:
            for al in selected_aliases:
                if al in seqs_by_alias:
                    fh.write(f">{id_of_alias[al]}\n{seqs_by_alias[al]}\n")
        print(f"  Wrote FASTA        : {fasta_out} ({len(seqs_by_alias):,} sequences)")

    if node_rows:
        with open(nodes_out, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(node_rows[0].keys()))
            writer.writeheader()
            writer.writerows(node_rows)
        print(f"  Wrote Nodes CSV    : {nodes_out} ({len(node_rows):,} nodes)")

    link_fields = [
        "source", "target", "edge_type", "exposure_setting",
        "same_household", "same_true_chain", "transmission_hops",
        "unsampled_intermediates", "snp_distance", "tn93_distance",
        "substitutions", "delta_exposure_days", "delta_collection_days",
        "same_county", "source_county", "target_county",
        "true_component_id", "variant_benchmark",
    ]
    with open(links_out, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=link_fields)
        writer.writeheader()
        writer.writerows(truth_links)
    print(f"  Wrote Truth Links  : {links_out} ({len(truth_links):,} edges)")

    # 11. Compute summary & SNP threshold clustering concordance
    direct_edges = [e for e in truth_links if e["edge_type"] == "direct"]
    ancestral_edges = [e for e in truth_links if e["edge_type"] == "ancestral"]
    coincidental_edges = [e for e in truth_links if e["edge_type"] == "coincidental_exposure"]

    hop_stats: dict[str, dict] = {}
    for h in sorted({int(e["transmission_hops"]) for e in truth_links if e["transmission_hops"] != ""}):
        h_edges = [e for e in truth_links if e["transmission_hops"] != "" and int(e["transmission_hops"]) == h]
        h_snps = [int(e["snp_distance"]) for e in h_edges if e["snp_distance"] != ""]
        h_days = [float(e["delta_exposure_days"]) for e in h_edges if e["delta_exposure_days"] != ""]
        entry: dict = {"n_edges": len(h_edges)}
        if h_snps:
            entry["mean_snps"] = round(sum(h_snps) / len(h_snps), 3)
            entry["median_snps"] = round(sorted(h_snps)[len(h_snps) // 2], 3)
            entry["max_snps"] = max(h_snps)
        if h_days:
            entry["mean_delta_exposure_days"] = round(sum(h_days) / len(h_days), 2)
        hop_stats[str(h)] = entry

    seqs_by_node_id = {
        id_of_alias[al]: seqs_by_alias[al]
        for al in selected_aliases
        if al in seqs_by_alias
    }
    concordance = []
    if compute_concordance and len(seqs_by_node_id) >= 2 and len(seqs_by_node_id) <= 2000:
        concordance = compute_snp_threshold_concordance(
            node_ids=[id_of_alias[al] for al in selected_aliases],
            seqs=seqs_by_node_id,
            true_components=true_comp_by_node_id,
        )

    multi_comp_count = sum(1 for c, cnt in subset_comp_counts.items() if cnt >= 2)
    superspreader_count = sum(1 for r in node_rows if r["is_superspreader"])
    urban_counts = Counter(r["urbanicity"] for r in node_rows)
    age_counts = Counter(r["age_group"] for r in node_rows if r["age_group"])
    variant_counts = Counter(r["variant_benchmark"] for r in node_rows if r["variant_benchmark"])
    setting_counts = Counter(e["exposure_setting"] for e in truth_links if e["edge_type"] != "coincidental_exposure")

    summary = {
        "prefix": prefix,
        "selection_mode": selection_mode,
        "id_column": id_column,
        "max_sequences_requested": max_sequences,
        "max_hops": max_hops,
        "n_nodes": len(node_rows),
        "n_sequences": len(seqs_by_alias),
        "n_distinct_true_components": len(subset_comp_counts),
        "n_multicase_true_components": multi_comp_count,
        "largest_subset_component_size": max(subset_comp_counts.values()) if subset_comp_counts else 0,
        "n_superspreaders_in_subset": superspreader_count,
        "n_truth_links_total": len(truth_links),
        "n_truth_links_direct_1hop": len(direct_edges),
        "n_truth_links_ancestral_multihop": len(ancestral_edges),
        "n_coincidental_exposure_links": len(coincidental_edges),
        "transmission_exposure_settings": dict(setting_counts),
        "hop_stratified_divergence": hop_stats,
        "urbanicity_breakdown": dict(urban_counts),
        "age_group_breakdown": dict(age_counts),
        "top_variants": dict(variant_counts.most_common(10)),
        "snp_threshold_concordance": concordance,
        "files": {
            "sequences_fasta": str(fasta_out) if seqs_by_alias else None,
            "nodes_csv": str(nodes_out),
            "truth_links_csv": str(links_out),
            "summary_json": str(summary_out),
        },
    }

    with open(summary_out, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"  Wrote Summary JSON : {summary_out}")
    print(
        f"  Summary            : {len(node_rows):,} nodes across "
        f"{len(subset_comp_counts):,} true chains ({multi_comp_count:,} multi-case, "
        f"{superspreader_count:,} superspreaders); "
        f"{len(truth_links):,} links ({len(direct_edges):,} direct 1-hop, "
        f"{len(ancestral_edges):,} ancestral 2..{max_hops}-hop, "
        f"{len(coincidental_edges):,} coincidental household/place exposures)."
    )
    if concordance:
        print("  MicrobeTrace SNP Threshold Concordance (vs. true_component_id):")
        for row_c in concordance:
            print(
                f"    <= {row_c['snp_threshold']} SNPs : "
                f"{row_c['n_genetic_links']:5,d} links "
                f"(same-chain={100 * row_c['link_same_chain_precision']:5.1f}%), "
                f"F1={row_c['pairwise_f1']:.3f} "
                f"(P={row_c['pairwise_precision']:.3f}, R={row_c['pairwise_recall']:.3f}), "
                f"ARI={row_c['adjusted_rand_index']:.3f}"
            )

    return summary


def resolve_paths_from_config(cfg, args) -> tuple[str | None, str | None, str | None, str, str]:
    """Resolve (fasta, samples, infections, outdir, prefix) from a PhyloGAS Config
    and CLI arguments."""
    project = str(cfg.get("project_name", default="phylogas"))
    results_dir = Path(str(cfg.get("results_dir", default="results")))
    sample_dir = Path(str(cfg.get("sampling.outdir", default=str(results_dir / "03_sampled_datasets"))))
    bench_dir = Path(str(cfg.get("benchmark.outdir", default=str(results_dir / "05_benchmarks"))))
    comp = cfg.get("genetic_painter.compression", default="xz")
    ext = {"xz": ".xz", "bgzf": ".gz"}.get(comp, "")

    recipe = getattr(args, "recipe", None) or "1S__surs"
    # Resolve recipe alias if BeyondBaseline is installed
    try:
        from scenarios_simulation.recipes import aliases, all_recipes
        known = all_recipes()
        if recipe not in known and recipe in aliases():
            recipe = aliases()[recipe]
    except ImportError:
        pass

    # 1. Samples / metadata file
    samples = getattr(args, "samples", None) or getattr(args, "metadata", None)
    if not samples:
        if recipe == "all_infections":
            samples = None
        else:
            cand_s = sample_dir / f"{recipe}_samples.csv.xz"
            if cand_s.exists():
                samples = str(cand_s)
            else:
                # Fallback to ascertained linelist
                ll = cfg.get("ascertainment.output", default=None)
                if ll and Path(str(ll)).exists():
                    samples = str(ll)

    # 2. Ground-truth infections (prefer allevents_variants.csv.xz, then linelist_allevents.csv.xz)
    infections = getattr(args, "infections", None)
    if not infections:
        cand_var = bench_dir / "allevents_variants.csv.xz"
        if cand_var.exists():
            infections = str(cand_var)
        else:
            derived = cfg.ascertainment_outputs().get("allevents")
            if derived and Path(str(derived)).exists():
                infections = str(derived)

    # 3. FASTA file (prefer subset FASTA for the recipe unless complete-chains is requested)
    fasta = getattr(args, "fasta", None)
    mode = getattr(args, "selection_mode", "cluster-enriched")
    if not fasta:
        paint_prefix = cfg.get(
            "genetic_painter.output_prefix",
            default=str(results_dir / "01_synthetic_genomes" / project),
        )
        full_fasta = Path(f"{paint_prefix}.sequences.fasta{ext}")
        subset_fasta = sample_dir / f"{project}.{recipe}.fasta{ext}"
        if mode == "complete-chains" and full_fasta.exists():
            fasta = str(full_fasta)
        elif subset_fasta.exists():
            fasta = str(subset_fasta)
        elif full_fasta.exists():
            fasta = str(full_fasta)

    outdir = getattr(args, "outdir", None) or str(results_dir / "06_microbetrace")
    prefix = getattr(args, "prefix", None) or f"{project}_{recipe}_microbetrace"
    return fasta, samples, infections, outdir, prefix


def add_microbetrace_arguments(sp: argparse.ArgumentParser) -> None:
    """Attach CLI arguments for ``phylogas export-microbetrace``."""
    sp.add_argument(
        "--config", "-c", default=None, metavar="YAML",
        help="Pipeline configuration file (e.g. configs/va.yaml). Resolves "
             "--fasta, --samples, --infections, and --outdir automatically.",
    )
    sp.add_argument(
        "--recipe", "-r", default="1S__surs",
        help="Sampling recipe to export when using --config (default: 1S__surs; "
             "e.g. 1S-P__lasso_greedy, 4S-P__lasso_greedy, or all_infections).",
    )
    sp.add_argument(
        "--fasta", "-f", default=None,
        help="Input painted or subset FASTA (.fasta, .fasta.xz, .fasta.gz).",
    )
    sp.add_argument(
        "--samples", "--metadata", "-m", dest="samples", default=None,
        help="Sampled cases CSV (e.g. 1S__surs_samples.csv.xz) or ascertained line list.",
    )
    sp.add_argument(
        "--infections", "-i", default=None,
        help="Ground-truth ABM all-events table (allevents_variants.csv.xz or "
             "linelist_allevents.csv.xz) used to trace direct and multi-hop "
             "transmission links.",
    )
    sp.add_argument(
        "--outdir", "-o", default=None,
        help="Output directory (default: <results_dir>/06_microbetrace or ./microbetrace_out).",
    )
    sp.add_argument(
        "--prefix", "-p", default=None,
        help="Filename prefix for emitted files (default: <project>_<recipe>_microbetrace).",
    )
    sp.add_argument(
        "--selection-mode", dest="selection_mode", default="cluster-enriched",
        choices=["cluster-enriched", "focal-window", "complete-chains", "random", "all"],
        help="How to select a browser-friendly subset of cases (default: cluster-enriched).",
    )
    sp.add_argument(
        "--max-sequences", dest="max_sequences", type=int, default=600,
        help="Maximum number of nodes/sequences to export (default: 600; 0 for all).",
    )
    sp.add_argument(
        "--max-hops", dest="max_hops", type=int, default=5,
        help="Maximum transmission hops along the true tree to emit ancestral bridge "
             "links in <prefix>_truth_links.csv (default: 5).",
    )
    sp.add_argument(
        "--background-fraction", dest="background_fraction", type=float, default=0.20,
        help="Fraction of --max-sequences reserved for singleton background cases "
             "in cluster-enriched mode (default: 0.20).",
    )
    sp.add_argument(
        "--min-chain-size", dest="min_chain_size", type=int, default=8,
        help="Minimum chain/subtree size in complete-chains mode (default: 8).",
    )
    sp.add_argument(
        "--max-chain-size", dest="max_chain_size", type=int, default=80,
        help="Maximum chain/subtree size per component in complete-chains mode (default: 80).",
    )
    sp.add_argument(
        "--county", default=None,
        help="Optional county name or 5-digit FIPS filter (e.g. 'Fairfax VA' or '51059').",
    )
    sp.add_argument(
        "--date-min", dest="date_min", default=None,
        help="Optional earliest collection/exposure date (YYYY-MM-DD).",
    )
    sp.add_argument(
        "--date-max", dest="date_max", default=None,
        help="Optional latest collection/exposure date (YYYY-MM-DD).",
    )
    sp.add_argument(
        "--variant", default=None,
        help="Optional variant_benchmark filter (e.g. 'AY.44', 'AY.103').",
    )
    sp.add_argument(
        "--id-column", dest="id_column", default="strain",
        choices=["strain", "alias_pid"],
        help="Primary node identifier used in FASTA headers, nodes.csv, and "
             "truth_links.csv (default: strain).",
    )
    sp.add_argument(
        "--seed", type=int, default=43,
        help="Random seed for reproducible background/random subsampling (default: 43).",
    )
    sp.add_argument(
        "--allow-no-fasta", dest="allow_no_fasta", action="store_true",
        help="Export _nodes.csv, _truth_links.csv, and _summary.json even if the "
             "FASTA file is not present locally.",
    )
    sp.add_argument(
        "--no-concordance", dest="compute_concordance", action="store_false", default=True,
        help="Skip single-linkage SNP threshold clustering concordance calculation.",
    )


def run_from_args(args: argparse.Namespace) -> int:
    fasta = args.fasta
    samples = args.samples
    infections = args.infections
    outdir = args.outdir or "microbetrace_out"
    prefix = args.prefix or "microbetrace"

    if args.config:
        from .config import Config, ConfigError
        try:
            cfg = Config.load(args.config)
        except ConfigError as exc:
            sys.exit(f"ERROR: {exc}")
        fasta, samples, infections, outdir, prefix = resolve_paths_from_config(cfg, args)

    try:
        export_microbetrace_bundle(
            outdir=outdir,
            prefix=prefix,
            fasta_path=fasta,
            samples_path=samples,
            infections_path=infections,
            selection_mode=args.selection_mode,
            max_sequences=args.max_sequences,
            max_hops=args.max_hops,
            background_fraction=args.background_fraction,
            min_chain_size=args.min_chain_size,
            max_chain_size=args.max_chain_size,
            county_filter=args.county,
            date_min=args.date_min,
            date_max=args.date_max,
            variant_filter=args.variant,
            id_column=args.id_column,
            seed=args.seed,
            allow_no_fasta=args.allow_no_fasta,
            compute_concordance=args.compute_concordance,
        )
    except (FileNotFoundError, ValueError) as exc:
        sys.exit(f"ERROR: {exc}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="phylogas export-microbetrace",
        description="Export MicrobeTrace-ready FASTA, node metadata, and ground-truth transmission link bundles.",
    )
    add_microbetrace_arguments(parser)
    args = parser.parse_args(argv)
    return run_from_args(args)


if __name__ == "__main__":
    raise SystemExit(main())
