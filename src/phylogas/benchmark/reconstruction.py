"""Phylogenetic topology, divergence-dating, and transmission-tree benchmarks.

Evaluates sequence-inferred reconstructions across the three operational rungs
of genomic epidemiology against PhyloGAS's exact agent-based ground truth:

1. **Cluster / Outbreak Triage** (e.g. ``MicrobeTrace`` SNP / TN93 thresholds):
   - Adjusted Rand Index (``chain_ari``), pairwise co-clustering Precision /
     Recall / ``F1`` against both statewide ``component_id`` chains and local
     cluster introductions (``local_intro_ari``, ``local_intro_f1``).
   - Directed ``SNP + collection_date`` baseline donor heuristic.

2. **Rooted Viral Genealogy & Time-Tree** (e.g. ``Nextstrain`` ``IQ-TREE 2`` +
   ``TreeTime``, ``JUNIPER`` / ``TransPhylo`` / ``outbreaker2`` host trees):
   - Strict Normalized Robinson-Foulds (``strict_nrf``).
   - Polytomy-collapsed Clade Precision, Recall, ``F1``, and ``collapsed_nrf``
     (collapsing inferred branches shorter than ``min_mutations``, default 0.5
     SNPs, so forced binary bifurcations on 0-SNP SARS-CoV-2 polytomies are not
     falsely penalized).
   - Generalized Information-Theoretic Split Similarity (``info_rf_similarity``,
     clustering mutual information matching following Smith 2020 / BREATH).
   - Quartet Concordance (``quartet_concordance`` over sampled 4-tip quartets,
     polytomy-aware and insensitive to rogue tips).
   - Rooted Kendall-Colijn topological correlation & distance (``kc_spearman_r``,
     ``kc_norm_dist`` on pairwise MRCA depths).
   - Cross-chain spurious merge rate (``cross_chain_merge_pct``: % of internal
     nodes whose descendants mix >= 2 independent ``component_id`` chains).
   - Clade TMRCA MAE / bias / 95% CI coverage (``clade_tmrca_mae_days``) and
     within-chain pairwise divergence-duration MAE / bias
     (``pairwise_duration_mae_days``) against true ABM ``exposure_date``s.

3. **Who-Infected-Whom Transmission Reconstruction** (``JUNIPER``, ``TransPhylo``,
   ``outbreaker2``, ``SNP+Date``):
   - Evaluated separately for **Direct** (``k = 1``, collapsing any unsampled
     intermediate or outside-cluster source to ``"UNOBSERVED"``) and **Indirect /
     Ancestral** (``k >= 1``, nearest sampled ancestor in the true chain):
     * ``donor_prob_mean`` — Mean posterior probability assigned to true donor
       (Specht et al. 2025 "Donor-Recipient Accuracy").
     * ``top1_donor_acc`` — Fraction where highest-posterior donor matches truth
       (Specht et al. "Most-Likely Donor"; Campbell et al. 2018 MPC accuracy).
     * ``obs_donor_recall`` & ``obs_link_precision`` — Recall on recipients whose
       true donor is in the sample, and precision among predicted host->host links
       (preventing sparse-sampling inflation from predicting ``UNOBSERVED``).
     * ``unobs_root_specificity`` — Fraction of cluster introductions / unsampled-
       parent cases correctly attributed to ``UNOBSERVED``.
     * ``donor_95set_cov`` & ``donor_95set_med_size`` — 95% credible-set coverage
       and median credible-set size (Specht et al. "Donor Set Accuracy").
     * ``calibration_ece`` — Expected Calibration Error across posterior deciles
       (Colijn et al. 2026 BREATH posterior reliability check).

Implemented purely in the Python standard library for zero-dependency portability
and streaming low-memory execution on multi-million-row ABM tables.
"""

from __future__ import annotations

import csv
import datetime as _dt
import gzip
import json
import lzma
import math
import random
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

UNOBSERVED = "UNOBSERVED"

_NO_VAL = {"", "-1", "nan", "None", "<NA>", "null", "NA"}

_REFERENCE_IDS = {
    "Wuhan/Hu-1/2019",
    "hCoV-19/Wuhan/Hu-1/2019",
    "MN908947.3",
    "NC_045512.2",
    "ref_genome",
}


# ---------------------------------------------------------------------------
# Basic IO & ID normalization
# ---------------------------------------------------------------------------
def _open_text(path: str | Path, mode: str = "rt"):
    p = str(path)
    if p.endswith(".xz"):
        return lzma.open(p, mode, encoding="utf-8")
    if p.endswith(".gz"):
        return gzip.open(p, mode, encoding="utf-8")
    raw_mode = mode.replace("t", "")
    return open(p, raw_mode, encoding="utf-8", newline="" if "w" in raw_mode else None)


def alias_from_any_id(raw: str | None, lookup_map: dict[str, str] | None = None) -> str | None:
    """Normalize a tip/host identifier to canonical ``alias_pid`` (``{pid}.{tick}``).

    Handles:
      - Safe JUNIPER/outbreaker IDs (``S000001``) via ``lookup_map``
      - Auspice / Nextstrain strain IDs (``hCoV-19/USA/VA-EHip-123.45/2021`` or
        ``USA/VA-EHip-123.45/2021``)
      - Direct ``alias_pid`` strings (``123.45``)
      - Reference / outgroup IDs -> returns ``None``
    """
    if not raw:
        return None
    s = str(raw).strip().strip("'\"")
    if not s or s in _NO_VAL or s in _REFERENCE_IDS:
        return None
    if lookup_map and s in lookup_map:
        mapped = lookup_map[s]
        if mapped in _REFERENCE_IDS:
            return None
        s = mapped
    if s in _REFERENCE_IDS or "Wuhan/Hu-1" in s or s.startswith("HOST_"):
        return None
    if "EHip-" in s:
        tail = s.split("EHip-", 1)[1]
        return tail.split("/", 1)[0].strip()
    return s


def load_sample_lookup(lookup_csv: str | Path | None) -> tuple[dict[str, str], dict[str, _dt.date]]:
    """Load ``sample_lookup.csv`` mapping safe ID (``S000001``) -> ``alias_pid`` & collection date."""
    id_map: dict[str, str] = {}
    date_map: dict[str, _dt.date] = {}
    if not lookup_csv or not Path(lookup_csv).exists():
        return id_map, date_map
    with _open_text(lookup_csv, "rt") as fh:
        rd = csv.DictReader(fh)
        for row in rd:
            safe = (row.get("sample") or row.get("id") or "").strip()
            alias = (row.get("alias_pid") or "").strip()
            strain = (row.get("strain") or "").strip()
            target = alias if alias else (alias_from_any_id(strain) or strain)
            if safe and target:
                id_map[safe] = target
            if strain and target:
                id_map[strain] = target
            d_raw = (row.get("date") or row.get("collection_date") or "").strip()[:10]
            if target and d_raw:
                try:
                    date_map[target] = _dt.date.fromisoformat(d_raw)
                except ValueError:
                    pass
    return id_map, date_map


def date_to_decimal_year(d: _dt.date) -> float:
    y0 = _dt.date(d.year, 1, 1).toordinal()
    y1 = _dt.date(d.year + 1, 1, 1).toordinal()
    return d.year + (d.toordinal() - y0) / float(y1 - y0)


def decimal_year_to_ordinal(dec_yr: float) -> float:
    yr = int(math.floor(dec_yr))
    frac = dec_yr - yr
    y0 = _dt.date(yr, 1, 1).toordinal()
    y1 = _dt.date(yr + 1, 1, 1).toordinal()
    return y0 + frac * (y1 - y0)


# ---------------------------------------------------------------------------
# Lightweight Rooted Tree Representation
# ---------------------------------------------------------------------------
class PhyloTree:
    """Rooted multifurcating/bifurcating tree over canonical tip IDs.

    Nodes are integer indices ``0 .. n_nodes - 1`` with ``root = 0``.
    """

    def __init__(self) -> None:
        self.parent: list[int] = [-1]
        self.children: list[list[int]] = [[]]
        self.name: list[str] = ["ROOT"]
        self.tip_alias: list[str | None] = [None]
        self.branch_muts: list[float] = [0.0]
        self.branch_len: list[float] = [0.0]
        self.num_date: list[float | None] = [None]
        self.num_date_ci: list[tuple[float, float] | None] = [None]
        self.meta: dict[str, Any] = {}

    def add_node(
        self,
        parent_idx: int,
        name: str = "",
        tip_alias: str | None = None,
        branch_muts: float = 0.0,
        branch_len: float = 0.0,
        num_date: float | None = None,
        num_date_ci: tuple[float, float] | None = None,
    ) -> int:
        idx = len(self.parent)
        self.parent.append(parent_idx)
        self.children.append([])
        self.name.append(name)
        self.tip_alias.append(tip_alias)
        self.branch_muts.append(float(branch_muts))
        self.branch_len.append(float(branch_len))
        self.num_date.append(num_date)
        self.num_date_ci.append(num_date_ci)
        if parent_idx >= 0:
            self.children[parent_idx].append(idx)
        return idx

    @property
    def n_nodes(self) -> int:
        return len(self.parent)

    def postorder(self) -> list[int]:
        order: list[int] = []
        stack = [0]
        while stack:
            u = stack.pop()
            order.append(u)
            for c in self.children[u]:
                stack.append(c)
        order.reverse()
        return order

    def tips_map(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for i, al in enumerate(self.tip_alias):
            if al is not None:
                out[al] = i
        return out

    def prune_to_tips(self, keep_aliases: set[str]) -> "PhyloTree":
        """Return a cleanly pruned copy containing only ``keep_aliases`` leaves,
        collapsing unary internal nodes while summing branch lengths/mutations."""
        post = self.postorder()
        has_kept = [False] * self.n_nodes
        for u in post:
            al = self.tip_alias[u]
            if al is not None and al in keep_aliases:
                has_kept[u] = True
            for c in self.children[u]:
                if has_kept[c]:
                    has_kept[u] = True

        out = PhyloTree()
        out.meta = dict(self.meta)
        if not has_kept[0]:
            return out

        # Walk top-down, suppressing degree-1 non-tip intermediates
        def _effective_target(u: int, acc_muts: float, acc_len: float) -> tuple[int, float, float]:
            while True:
                kept_ch = [c for c in self.children[u] if has_kept[c]]
                is_kept_tip = self.tip_alias[u] is not None and self.tip_alias[u] in keep_aliases
                if is_kept_tip or len(kept_ch) != 1:
                    return u, acc_muts, acc_len
                c = kept_ch[0]
                acc_muts += self.branch_muts[c]
                acc_len += self.branch_len[c]
                u = c

        root_u, _, _ = _effective_target(0, 0.0, 0.0)
        out.name[0] = self.name[root_u]
        out.tip_alias[0] = self.tip_alias[root_u] if self.tip_alias[root_u] in keep_aliases else None
        out.num_date[0] = self.num_date[root_u]
        out.num_date_ci[0] = self.num_date_ci[root_u]

        queue: list[tuple[int, int]] = [(root_u, 0)]
        while queue:
            src_u, dst_u = queue.pop()
            for c in self.children[src_u]:
                if not has_kept[c]:
                    continue
                eff_c, m_sum, l_sum = _effective_target(c, self.branch_muts[c], self.branch_len[c])
                al = self.tip_alias[eff_c] if self.tip_alias[eff_c] in keep_aliases else None
                # If an internal node itself carries a kept tip AND still has kept children,
                # represent the tip as a zero-length pendant leaf so every sample is a leaf.
                if al is not None and any(has_kept[gc] for gc in self.children[eff_c]):
                    internal_idx = out.add_node(
                        dst_u,
                        name=f"INTERNAL_{self.name[eff_c]}",
                        tip_alias=None,
                        branch_muts=m_sum,
                        branch_len=l_sum,
                        num_date=self.num_date[eff_c],
                        num_date_ci=self.num_date_ci[eff_c],
                    )
                    out.add_node(
                        internal_idx,
                        name=self.name[eff_c],
                        tip_alias=al,
                        branch_muts=0.0,
                        branch_len=0.0,
                        num_date=self.num_date[eff_c],
                        num_date_ci=self.num_date_ci[eff_c],
                    )
                    queue.append((eff_c, internal_idx))
                else:
                    new_idx = out.add_node(
                        dst_u,
                        name=self.name[eff_c],
                        tip_alias=al,
                        branch_muts=m_sum,
                        branch_len=l_sum,
                        num_date=self.num_date[eff_c],
                        num_date_ci=self.num_date_ci[eff_c],
                    )
                    if al is None:
                        queue.append((eff_c, new_idx))
        return out


# ---------------------------------------------------------------------------
# Tree & Network Parsers
# ---------------------------------------------------------------------------
def parse_auspice_json(
    path: str | Path,
    lookup_map: dict[str, str] | None = None,
) -> PhyloTree:
    """Parse a Nextstrain Auspice v2 JSON into a ``PhyloTree``."""
    with _open_text(path, "rt") as fh:
        data = json.load(fh)
    root_obj = data.get("tree")
    if not root_obj:
        raise ValueError(f"No 'tree' key in Auspice JSON: {path}")

    tree = PhyloTree()

    def _extract_attrs(node_obj: dict, parent_div: float) -> tuple[float, float, float | None, tuple[float, float] | None]:
        na = node_obj.get("node_attrs") or {}
        ba = node_obj.get("branch_attrs") or {}
        div_val = na.get("div")
        cur_div = float(div_val) if isinstance(div_val, (int, float)) else parent_div
        nuc_muts = (ba.get("mutations") or {}).get("nuc")
        if isinstance(nuc_muts, list):
            b_muts = float(len(nuc_muts))
        else:
            b_muts = max(0.0, cur_div - parent_div)

        nd = na.get("num_date")
        num_date: float | None = None
        num_date_ci: tuple[float, float] | None = None
        if isinstance(nd, dict) and isinstance(nd.get("value"), (int, float)):
            num_date = float(nd["value"])
            ci = nd.get("confidence")
            if isinstance(ci, (list, tuple)) and len(ci) == 2:
                try:
                    num_date_ci = (float(ci[0]), float(ci[1]))
                except (TypeError, ValueError):
                    pass
        elif isinstance(nd, (int, float)):
            num_date = float(nd)
        return cur_div, b_muts, num_date, num_date_ci

    r_div, _, r_date, r_ci = _extract_attrs(root_obj, 0.0)
    r_name = str(root_obj.get("name") or "ROOT")
    tree.name[0] = r_name
    tree.num_date[0] = r_date
    tree.num_date_ci[0] = r_ci

    stack: list[tuple[dict, int, float]] = [(root_obj, 0, r_div)]
    while stack:
        obj, u_idx, u_div = stack.pop()
        chs = obj.get("children") or []
        for ch in reversed(chs):
            c_name = str(ch.get("name") or "")
            c_chs = ch.get("children") or []
            c_div, b_muts, c_date, c_ci = _extract_attrs(ch, u_div)
            p_date = tree.num_date[u_idx]
            b_len = max(0.0, c_date - p_date) if (c_date is not None and p_date is not None) else b_muts
            tip_al = alias_from_any_id(c_name, lookup_map) if not c_chs else None
            c_idx = tree.add_node(
                u_idx,
                name=c_name,
                tip_alias=tip_al,
                branch_muts=b_muts,
                branch_len=b_len,
                num_date=c_date,
                num_date_ci=c_ci,
            )
            if c_chs:
                stack.append((ch, c_idx, c_div))
    return tree


def parse_newick(
    newick_str: str,
    node_data: dict[str, Any] | None = None,
    lookup_map: dict[str, str] | None = None,
    genome_length: float = 29903.0,
) -> PhyloTree:
    """Iterative (non-recursive) Newick parser supporting optional Augur ``branch_lengths.json``."""
    s = newick_str.strip()
    if s.endswith(";"):
        s = s[:-1]

    nodes_meta = (node_data or {}).get("nodes") or {}
    clock_meta = (node_data or {}).get("clock") or {}

    tree = PhyloTree()
    if isinstance(clock_meta.get("rate"), (int, float)):
        tree.meta["clock_rate_yr"] = float(clock_meta["rate"])

    cur = 0
    n_chars = len(s)
    i = 0

    def _read_token(pos: int) -> tuple[str, int]:
        if pos < n_chars and s[pos] in ("'", '"'):
            quote = s[pos]
            pos += 1
            buf = []
            while pos < n_chars and s[pos] != quote:
                buf.append(s[pos])
                pos += 1
            return "".join(buf), min(n_chars, pos + 1)
        buf = []
        while pos < n_chars and s[pos] not in ("(", ")", ",", ":", ";"):
            buf.append(s[pos])
            pos += 1
        return "".join(buf).strip(), pos

    while i < n_chars:
        ch = s[i]
        if ch == "(":
            c_idx = tree.add_node(cur)
            cur = c_idx
            i += 1
        elif ch == ",":
            p_idx = tree.parent[cur]
            if p_idx < 0:
                p_idx = 0
            c_idx = tree.add_node(p_idx)
            cur = c_idx
            i += 1
        elif ch == ")":
            cur = tree.parent[cur]
            if cur < 0:
                cur = 0
            i += 1
        elif ch == ":":
            i += 1
            tok, i = _read_token(i)
            try:
                bl = float(tok)
            except ValueError:
                bl = 0.0
            tree.branch_len[cur] = bl
            # Heuristic: if branch length < 0.05 in a raw/iqtree Newick, it is
            # substitutions/site; convert to mutation count unless overridden by node_data.
            tree.branch_muts[cur] = bl * genome_length if bl < 0.1 else bl
        elif ch.isspace():
            i += 1
        else:
            tok, i = _read_token(i)
            if tok:
                tree.name[cur] = tok

    # Attach tip aliases and optional Augur branch_lengths.json attributes
    for idx in range(tree.n_nodes):
        nm = tree.name[idx]
        if not tree.children[idx]:
            tree.tip_alias[idx] = alias_from_any_id(nm, lookup_map)
        nd = nodes_meta.get(nm)
        if isinstance(nd, dict):
            if isinstance(nd.get("mutation_length"), (int, float)):
                tree.branch_muts[idx] = float(nd["mutation_length"])
            if isinstance(nd.get("branch_length"), (int, float)):
                tree.branch_len[idx] = float(nd["branch_length"])
            nd_val = nd.get("num_date") if nd.get("num_date") is not None else nd.get("numdate")
            if isinstance(nd_val, (int, float)):
                tree.num_date[idx] = float(nd_val)
            ci = nd.get("num_date_confidence") or nd.get("numdate_confidence")
            if isinstance(ci, (list, tuple)) and len(ci) == 2:
                try:
                    tree.num_date_ci[idx] = (float(ci[0]), float(ci[1]))
                except (TypeError, ValueError):
                    pass
    return tree


def parse_newick_file(
    newick_path: str | Path,
    branch_lengths_json: str | Path | None = None,
    lookup_map: dict[str, str] | None = None,
    genome_length: float = 29903.0,
) -> PhyloTree:
    with _open_text(newick_path, "rt") as fh:
        nwk = fh.read()
    nd = None
    if branch_lengths_json and Path(branch_lengths_json).exists():
        with _open_text(branch_lengths_json, "rt") as fh:
            nd = json.load(fh)
    return parse_newick(nwk, node_data=nd, lookup_map=lookup_map, genome_length=genome_length)


def export_dated_newick_for_transphylo(
    newick_path: str | Path,
    branch_lengths_json: str | Path,
    out_newick: str | Path,
    out_meta_json: str | Path,
    lookup_map: dict[str, str] | None = None,
    min_branch_years: float = 1e-5,
) -> dict[str, Any]:
    """Convert ``augur refine`` output (``tree.nwk`` + ``branch_lengths.json``) into
    a strictly bifurcating, outgroup-pruned Newick tree with branch lengths in
    **decimal years**, ready for ``TransPhylo::ptreeFromPhylo``."""
    raw_tree = parse_newick_file(newick_path, branch_lengths_json=branch_lengths_json, lookup_map=None)
    # Keep all non-reference leaves using their original FASTA sample names (e.g. S000001)
    keep_names: set[str] = set()
    tip_dates: list[float] = []
    for i in range(raw_tree.n_nodes):
        if not raw_tree.children[i]:
            nm = raw_tree.name[i]
            if nm and nm not in _REFERENCE_IDS and "Wuhan/Hu-1" not in nm:
                raw_tree.tip_alias[i] = nm
                keep_names.add(nm)
                if raw_tree.num_date[i] is not None:
                    tip_dates.append(float(raw_tree.num_date[i]))
            else:
                raw_tree.tip_alias[i] = None

    pruned = raw_tree.prune_to_tips(keep_names)

    # Recompute branch lengths in decimal years from parent/child num_date,
    # and binarize any multifurcations with epsilon = min_branch_years.
    children_bin: list[list[int]] = [list(ch) for ch in pruned.children]
    branch_yrs: list[float] = [0.0] * pruned.n_nodes
    for u in range(1, pruned.n_nodes):
        p = pruned.parent[u]
        du = pruned.num_date[u]
        dp = pruned.num_date[p]
        if du is not None and dp is not None:
            branch_yrs[u] = max(min_branch_years, float(du - dp))
        else:
            branch_yrs[u] = max(min_branch_years, pruned.branch_len[u])

    # Iterative serialization to Newick with inline binarization of >2 children
    def _serialize(u: int, bl: float | None) -> str:
        chs = children_bin[u]
        if not chs:
            lbl = pruned.name[u]
            return f"{lbl}:{bl:.8f}" if bl is not None else lbl
        sub_strs = [_serialize(c, branch_yrs[c]) for c in chs]
        expr = sub_strs[0]
        for nxt in sub_strs[1:]:
            expr = f"({expr},{nxt})"
            if nxt is not sub_strs[-1]:
                expr = f"{expr}:{min_branch_years:.8f}"
        return f"{expr}:{bl:.8f}" if bl is not None else expr

    nwk_out = _serialize(0, None) + ";\n"
    Path(out_newick).parent.mkdir(parents=True, exist_ok=True)
    Path(out_newick).write_text(nwk_out, encoding="utf-8")

    date_last = max(tip_dates) if tip_dates else 2021.75
    meta = {
        "n_tips": len(keep_names),
        "date_last_sample": round(date_last, 8),
        "clock_rate_yr": pruned.meta.get("clock_rate_yr"),
    }
    Path(out_meta_json).write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return meta


def parse_host_network_csv(
    nodes_csv: str | Path,
    lookup_csv: str | Path | None = None,
) -> tuple[PhyloTree, dict[str, dict[str, dict[str, float]]]]:
    """Parse a host-level transmission network CSV (``JUNIPER``, ``TransPhylo``,
    or ``outbreaker2`` ``nodes.csv``) into:
      1) A pendant-leaf ``PhyloTree`` (each observed host becomes an internal
         host node dated at ``infection_date`` with a pendant leaf dated at
         ``collection_date``, matching ``juniper_to_graph.py``).
      2) Point-estimate direct & indirect donor probability dictionaries.
    """
    id_map, coll_dates = load_sample_lookup(lookup_csv)
    rows: list[dict[str, str]] = []
    with _open_text(nodes_csv, "rt") as fh:
        for r in csv.DictReader(fh):
            rows.append(r)

    host_parent: dict[int, int | None] = {}
    host_inf_dec: dict[int, float | None] = {}
    host_alias: dict[int, str | None] = {}
    children_of: dict[int, list[int]] = defaultdict(list)

    for r in rows:
        idx = int(float(r["idx"]))
        p_raw = (r.get("infector_idx") or "").strip()
        p_idx: int | None = None
        if p_raw and p_raw not in _NO_VAL:
            cand = int(float(p_raw))
            if cand > 0 and cand != idx:
                p_idx = cand
        host_parent[idx] = p_idx
        if p_idx is not None:
            children_of[p_idx].append(idx)

        inf_d = (r.get("infection_date") or "").strip()[:10]
        dec_yr: float | None = None
        if inf_d and inf_d not in _NO_VAL:
            try:
                dec_yr = date_to_decimal_year(_dt.date.fromisoformat(inf_d))
            except ValueError:
                pass
        elif (r.get("infection_num_date") or "").strip() not in _NO_VAL:
            try:
                dec_yr = float(r["infection_num_date"])
            except ValueError:
                pass
        host_inf_dec[idx] = dec_yr

        obs_flag = str(r.get("observed") or "").strip().lower() in ("true", "1", "t", "yes")
        samp = (r.get("sample") or "").strip()
        host_alias[idx] = alias_from_any_id(samp, id_map) if (obs_flag and samp) else None

    # Build PhyloTree (supporting single root or multi-root forest under synthetic super-root)
    roots = [i for i, p in host_parent.items() if p is None or p not in host_parent]
    tree = PhyloTree()
    root_dates = [host_inf_dec[r] for r in roots if host_inf_dec.get(r) is not None]
    if root_dates:
        tree.num_date[0] = min(root_dates)

    host_node_idx: dict[int, int] = {}
    queue: list[tuple[int, int]] = []
    if len(roots) == 1 and host_alias.get(roots[0]) is None:
        r0 = roots[0]
        host_node_idx[r0] = 0
        tree.name[0] = f"HOST_{r0:06d}"
        tree.num_date[0] = host_inf_dec.get(r0)
        for ch in children_of.get(r0, ()):
            queue.append((ch, 0))
    else:
        for r0 in roots:
            queue.append((r0, 0))

    visited_hosts: set[int] = set(host_node_idx.keys())
    while queue:
        h_id, p_tree_idx = queue.pop(0)
        if h_id in visited_hosts:
            continue
        visited_hosts.add(h_id)
        h_date = host_inf_dec.get(h_id)
        p_date = tree.num_date[p_tree_idx]
        b_len = max(0.0, h_date - p_date) if (h_date is not None and p_date is not None) else 1.0
        h_tree_idx = tree.add_node(
            p_tree_idx,
            name=f"HOST_{h_id:06d}",
            tip_alias=None,
            branch_muts=1.0,  # host step treated as non-zero topological event
            branch_len=b_len,
            num_date=h_date,
        )
        host_node_idx[h_id] = h_tree_idx
        al = host_alias.get(h_id)
        if al is not None:
            c_date_obj = coll_dates.get(al)
            c_dec = date_to_decimal_year(c_date_obj) if c_date_obj else h_date
            pendant_len = max(0.0, c_dec - h_date) if (c_dec is not None and h_date is not None) else 0.0
            tree.add_node(
                h_tree_idx,
                name=al,
                tip_alias=al,
                branch_muts=0.0,
                branch_len=pendant_len,
                num_date=c_dec,
            )
        for ch in children_of.get(h_id, ()):
            queue.append((ch, h_tree_idx))

    all_tips = {al for al in host_alias.values() if al is not None}
    pruned_tree = tree.prune_to_tips(all_tips)

    # Point-estimate direct & indirect donors from host_parent
    direct_probs: dict[str, dict[str, float]] = {}
    indirect_probs: dict[str, dict[str, float]] = {}
    for h_id, rec_al in host_alias.items():
        if rec_al is None:
            continue
        p = host_parent.get(h_id)
        if p is not None and host_alias.get(p) is not None:
            direct_probs[rec_al] = {host_alias[p]: 1.0}  # type: ignore[dict-item]
        else:
            direct_probs[rec_al] = {UNOBSERVED: 1.0}

        # Walk ancestors until hitting an observed host
        cur = p
        seen = {h_id}
        anc_al: str | None = None
        while cur is not None and cur not in seen:
            seen.add(cur)
            if host_alias.get(cur) is not None:
                anc_al = host_alias[cur]
                break
            cur = host_parent.get(cur)
        indirect_probs[rec_al] = {anc_al if anc_al is not None else UNOBSERVED: 1.0}

    return pruned_tree, {"direct": direct_probs, "indirect": indirect_probs}


def parse_donor_probs_csv(
    probs_csv: str | Path,
    lookup_csv: str | Path | None = None,
) -> dict[str, dict[str, dict[str, float]]]:
    """Read ``<tool>_donor_probs.csv`` with columns:
    ``recipient_id, donor_id, direct_prob, indirect_prob``."""
    id_map, _ = load_sample_lookup(lookup_csv)
    direct: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    indirect: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))

    with _open_text(probs_csv, "rt") as fh:
        for r in csv.DictReader(fh):
            rec_raw = (r.get("recipient_id") or r.get("recipient") or "").strip()
            don_raw = (r.get("donor_id") or r.get("donor") or "").strip()
            rec = alias_from_any_id(rec_raw, id_map)
            if not rec:
                continue
            if don_raw in ("UNOBSERVED", "unobserved", "NA", "", "-1", "ref_genome"):
                don = UNOBSERVED
            else:
                don = alias_from_any_id(don_raw, id_map) or UNOBSERVED

            dp = float(r.get("direct_prob") or 0.0)
            ip = float(r.get("indirect_prob") or 0.0)
            if dp > 0.0:
                direct[rec][don] += dp
            if ip > 0.0:
                indirect[rec][don] += ip

    # Normalize each recipient's distribution to sum to 1.0 (attributing any residual to UNOBSERVED)
    out_dir: dict[str, dict[str, float]] = {}
    out_ind: dict[str, dict[str, float]] = {}
    all_recs = set(direct.keys()) | set(indirect.keys())
    for rec in all_recs:
        for src, dst in ((direct[rec], out_dir), (indirect[rec], out_ind)):
            tot = sum(src.values())
            if tot <= 0.0:
                dst[rec] = {UNOBSERVED: 1.0}
            elif tot < 0.999999:
                normed = dict(src)
                normed[UNOBSERVED] = normed.get(UNOBSERVED, 0.0) + (1.0 - tot)
                dst[rec] = normed
            else:
                dst[rec] = {k: v / tot for k, v in src.items()}
    return {"direct": out_dir, "indirect": out_ind}


# ---------------------------------------------------------------------------
# Ground-Truth ABM Forest Streaming & Exact Pruned Genealogy Construction
# ---------------------------------------------------------------------------
class TruthForest:
    """Compact ancestor DAG for a target tip universe."""

    def __init__(
        self,
        parent_of: dict[str, str],
        exp_ordinal: dict[str, int],
        comp_of: dict[str, str],
        variant_of: dict[str, str],
        fips_of: dict[str, str],
    ) -> None:
        self.parent_of = parent_of
        self.exp_ordinal = exp_ordinal
        self.comp_of = comp_of
        self.variant_of = variant_of
        self.fips_of = fips_of


def load_truth_forest(
    infections_path: str | Path,
    target_aliases: set[str] | None = None,
) -> TruthForest:
    """Stream ``allevents_variants.csv.xz`` once and retain only ancestors
    reachable from ``target_aliases`` (or all nodes if ``target_aliases`` is None)."""
    all_parent: dict[str, str] = {}
    all_ord: dict[str, int] = {}
    all_comp: dict[str, str] = {}
    all_var: dict[str, str] = {}
    all_fips: dict[str, str] = {}

    with _open_text(infections_path, "rt") as fh:
        hdr_line = fh.readline()
        delim = "\t" if hdr_line.count("\t") > hdr_line.count(",") else ","
        cols = [c.strip() for c in hdr_line.rstrip("\r\n").split(delim)]
        idx = {c: i for i, c in enumerate(cols)}

        ia = idx["alias_pid"]
        ip = idx.get("alias_contact")
        idt = idx.get("exposure_date", idx.get("date"))
        ic = idx.get("component_id")
        iv = idx.get("variant_benchmark", idx.get("variant_label"))
        ifips = idx.get("county_fips", idx.get("FIPS"))

        reader = csv.reader(fh, delimiter=delim)
        for row in reader:
            if len(row) <= ia:
                continue
            al = row[ia].strip()
            if not al or al in all_parent:
                continue
            p = row[ip].strip() if ip is not None and len(row) > ip else ""
            all_parent[al] = p if (p and p not in _NO_VAL and p != al) else ""
            if idt is not None and len(row) > idt:
                ds = row[idt].strip()[:10]
                if len(ds) == 10 and ds[4] == "-":
                    try:
                        all_ord[al] = _dt.date.fromisoformat(ds).toordinal()
                    except ValueError:
                        pass
            if ic is not None and len(row) > ic:
                c = row[ic].strip()
                if c and c not in _NO_VAL:
                    all_comp[al] = c
            if iv is not None and len(row) > iv:
                v = row[iv].strip()
                if v and v not in _NO_VAL:
                    all_var[al] = v
            if ifips is not None and len(row) > ifips:
                f = row[ifips].strip().split(".")[0]
                if f and f not in _NO_VAL:
                    all_fips[al] = f

    if target_aliases is None:
        return TruthForest(all_parent, all_ord, all_comp, all_var, all_fips)

    # Retain only target_aliases + all their ancestors up to root
    keep: set[str] = set()
    for t in target_aliases:
        cur = t
        while cur and cur not in keep:
            keep.add(cur)
            cur = all_parent.get(cur, "")

    return TruthForest(
        parent_of={k: all_parent.get(k, "") for k in keep},
        exp_ordinal={k: all_ord[k] for k in keep if k in all_ord},
        comp_of={k: all_comp[k] for k in keep if k in all_comp},
        variant_of={k: all_var[k] for k in keep if k in all_var},
        fips_of={k: all_fips[k] for k in keep if k in all_fips},
    )


def build_pruned_truth_genealogy(
    truth: TruthForest,
    tip_aliases: Iterable[str],
    tip_collection_dates: dict[str, _dt.date] | None = None,
) -> PhyloTree:
    """Construct the exact pruned ground-truth viral genealogy ``G*_S`` over
    ``tip_aliases``.

    Every sampled host becomes a leaf; if a sampled host ``u`` also has sampled
    descendants, ``u``'s infection event is an internal node dated at
    ``exposure_date(u)`` with a zero-mutation pendant leaf for ``u``. Multiple
    independent importation chains join at the synthetic root (index 0).
    """
    tips = sorted(set(tip_aliases))
    tip_set = set(tips)

    # Collect induced subgraph of all ancestor paths from tips
    active_nodes: set[str] = set()
    children_in_sub: dict[str, list[str]] = defaultdict(list)
    roots_in_sub: set[str] = set()

    for t in tips:
        cur = t
        visited_path: set[str] = set()
        while cur and cur not in visited_path:
            visited_path.add(cur)
            was_active = cur in active_nodes
            active_nodes.add(cur)
            p = truth.parent_of.get(cur, "")
            if p and p not in visited_path:
                if not was_active:
                    children_in_sub[p].append(cur)
                if p in active_nodes:
                    break
                cur = p
            else:
                roots_in_sub.add(cur)
                break

    for ch_list in children_in_sub.values():
        ch_list.sort()

    # Build full induced host tree, then prune_to_tips collapses unsampled degree-1 chains!
    raw = PhyloTree()
    root_ords = [truth.exp_ordinal[r] for r in roots_in_sub if r in truth.exp_ordinal]
    if root_ords:
        raw.num_date[0] = date_to_decimal_year(_dt.date.fromordinal(min(root_ords)))

    queue: list[tuple[str, int]] = [(r, 0) for r in sorted(roots_in_sub)]
    seen_host: set[str] = set()
    while queue:
        h, p_idx = queue.pop(0)
        if h in seen_host:
            continue
        seen_host.add(h)
        h_ord = truth.exp_ordinal.get(h)
        h_dec = date_to_decimal_year(_dt.date.fromordinal(h_ord)) if h_ord is not None else None
        p_dec = raw.num_date[p_idx]
        b_len = max(0.0, h_dec - p_dec) if (h_dec is not None and p_dec is not None) else 1.0
        chs = children_in_sub.get(h, [])

        if h in tip_set and chs:
            # Sampled host with sampled descendants -> internal node + pendant leaf
            int_idx = raw.add_node(
                p_idx,
                name=f"HOST_{h}",
                tip_alias=None,
                branch_muts=1.0,  # hop count stored in branch_muts on truth tree
                branch_len=b_len,
                num_date=h_dec,
            )
            c_dt = (tip_collection_dates or {}).get(h)
            c_dec = date_to_decimal_year(c_dt) if c_dt else h_dec
            raw.add_node(
                int_idx,
                name=h,
                tip_alias=h,
                branch_muts=0.0,
                branch_len=max(0.0, (c_dec - h_dec) if (c_dec and h_dec) else 0.0),
                num_date=c_dec,
            )
            for c in chs:
                queue.append((c, int_idx))
        elif h in tip_set:
            c_dt = (tip_collection_dates or {}).get(h)
            c_dec = date_to_decimal_year(c_dt) if c_dt else h_dec
            raw.add_node(
                p_idx,
                name=h,
                tip_alias=h,
                branch_muts=1.0,
                branch_len=b_len,
                num_date=c_dec,
            )
        else:
            int_idx = raw.add_node(
                p_idx,
                name=f"HOST_{h}",
                tip_alias=None,
                branch_muts=1.0,
                branch_len=b_len,
                num_date=h_dec,
            )
            for c in chs:
                queue.append((c, int_idx))

    return raw.prune_to_tips(tip_set)


def build_truth_donors(
    truth: TruthForest,
    tip_aliases: Iterable[str],
) -> tuple[dict[str, str], dict[str, str], dict[str, int], dict[str, str]]:
    """For each recipient ``j`` in ``tip_aliases``, compute:
      - ``true_direct[j]``: ``parent_of[j]`` if in ``tip_aliases`` else ``"UNOBSERVED"``
      - ``true_indirect[j]``: nearest ancestor of ``j`` in ``tip_aliases``, or ``"UNOBSERVED"``
      - ``hops_to_anc[j]``: transmission hops to ``true_indirect[j]`` (0 if ``UNOBSERVED``)
      - ``local_intro_of[j]``: ID of the local cluster founder (topmost sampled or
        unsampled ancestor inside the connected component induced by ``tip_aliases``
        or the local root in ``tip_aliases``).
    """
    tips = set(tip_aliases)
    true_direct: dict[str, str] = {}
    true_indirect: dict[str, str] = {}
    hops_to_anc: dict[str, int] = {}

    for j in tips:
        p = truth.parent_of.get(j, "")
        true_direct[j] = p if p in tips else UNOBSERVED

        cur = p
        hops = 1
        visited = {j}
        found_anc = UNOBSERVED
        found_hops = 0
        while cur and cur not in visited:
            visited.add(cur)
            if cur in tips:
                found_anc = cur
                found_hops = hops
                break
            cur = truth.parent_of.get(cur, "")
            hops += 1
        true_indirect[j] = found_anc
        hops_to_anc[j] = found_hops

    # Partition tips into connected components under `true_indirect` (each local root
    # whose `true_indirect == UNOBSERVED` defines a distinct local introduction subtree!)
    local_intro_of: dict[str, str] = {}
    for j in tips:
        cur = j
        seen = set()
        while true_indirect.get(cur, UNOBSERVED) != UNOBSERVED and cur not in seen:
            seen.add(cur)
            cur = true_indirect[cur]
        local_intro_of[j] = cur

    return true_direct, true_indirect, hops_to_anc, local_intro_of


# ---------------------------------------------------------------------------
# Phylogenetic Topology & Dating Metrics
# ---------------------------------------------------------------------------
def _extract_clades_fast(
    tree: PhyloTree,
    tip_zobrist: dict[str, int],
    tip_to_bit: dict[str, int] | None,
    min_branch_muts: float = 0.0,
) -> tuple[dict[int, int], dict[int, int] | None, list[int]]:
    """Extract non-trivial clades in strict ``O(N)`` time & memory via 128-bit
    Zobrist XOR hashing, plus optional exact bitmasks when ``N <= 2500`` (for
    ``InfoRF``)."""
    n_tips = len(tip_zobrist)
    z_hash = [0] * tree.n_nodes
    sub_sz = [0] * tree.n_nodes
    masks = [0] * tree.n_nodes if tip_to_bit is not None else None

    for u in tree.postorder():
        al = tree.tip_alias[u]
        if al is not None and al in tip_zobrist:
            z_hash[u] ^= tip_zobrist[al]
            sub_sz[u] += 1
            if masks is not None and tip_to_bit is not None:
                masks[u] |= 1 << tip_to_bit[al]
        for c in tree.children[u]:
            z_hash[u] ^= z_hash[c]
            sub_sz[u] += sub_sz[c]
            if masks is not None:
                masks[u] |= masks[c]

    clades_z: dict[int, int] = {}
    clades_bit: dict[int, int] | None = {} if masks is not None else None
    for u in range(1, tree.n_nodes):
        sz = sub_sz[u]
        if 2 <= sz < n_tips and tree.branch_muts[u] >= min_branch_muts:
            clades_z.setdefault(z_hash[u], u)
            if clades_bit is not None and masks is not None:
                clades_bit.setdefault(masks[u], u)
    return clades_z, clades_bit, sub_sz


def _split_entropy(k: int, n: int) -> float:
    if k <= 0 or k >= n:
        return 0.0
    p = k / float(n)
    q = 1.0 - p
    return -(p * math.log2(p) + q * math.log2(q))


_POPCOUNT = getattr(int, "bit_count", None)


def _bit_count(x: int) -> int:
    return _POPCOUNT(x) if _POPCOUNT is not None else bin(x).count("1")


def _mutual_clustering_info(m1: int, k1: int, m2: int, k2: int, n: int) -> float:
    """Clustering mutual information (in bits) between two splits on ``n`` leaves."""
    a = _bit_count(m1 & m2)    # in1 & in2
    b = k1 - a                 # in1 & out2
    c = k2 - a                 # out1 & in2
    d = n - k1 - k2 + a        # out1 & out2
    fn = float(n)
    mi = 0.0
    for cell, row_tot, col_tot in (
        (a, k1, k2),
        (b, k1, n - k2),
        (c, n - k1, k2),
        (d, n - k1, n - k2),
    ):
        if cell > 0 and row_tot > 0 and col_tot > 0:
            mi += (cell / fn) * math.log2((cell * fn) / (row_tot * col_tot))
    return max(0.0, mi)


def compute_info_rf_similarity(
    true_clades: dict[int, int],
    inf_clades: dict[int, int],
    n_tips: int,
    max_splits: int = 1200,
) -> float | None:
    """Information-theoretic Generalized Robinson-Foulds similarity in ``[0, 1]``
    (following Smith 2020 ``TreeDist::MutualClusteringInfo``).

    Returns ``None`` if ``max(|true_clades|, |inf_clades|) > max_splits`` to bound
    ``O(S^2)`` runtime on 15,000-tip statewide trees (where ``collapsed_clade_f1``
    and ``quartet_concordance`` remain ``O(N)``).
    """
    if not true_clades and not inf_clades:
        return 1.0
    if len(true_clades) > max_splits or len(inf_clades) > max_splits:
        return None

    # Canonicalize unrooted/rooted bipartitions by flipping masks where bit 0 is set
    full = (1 << n_tips) - 1

    def _canon(clade_keys: Iterable[int]) -> list[tuple[int, int]]:
        seen: set[int] = set()
        out: list[tuple[int, int]] = []
        for m in clade_keys:
            cm = (full ^ m) if (m & 1) else m
            k = _bit_count(cm)
            if 2 <= k <= n_tips - 2 and cm not in seen:
                seen.add(cm)
                out.append((cm, k))
        return out

    t_splits = _canon(true_clades.keys())
    i_splits = _canon(inf_clades.keys())
    if not t_splits and not i_splits:
        return 1.0

    h_true = sum(_split_entropy(k, n_tips) for _, k in t_splits)
    h_inf = sum(_split_entropy(k, n_tips) for _, k in i_splits)
    denom = 0.5 * (h_true + h_inf)
    if denom <= 1e-12:
        return 1.0

    # Exact match shortcut + greedy highest-MI bipartite matching on remaining splits
    i_dict = dict(i_splits)
    matched_t: set[int] = set()
    matched_i: set[int] = set()
    total_mi = 0.0

    for idx_t, (mt, kt) in enumerate(t_splits):
        if mt in i_dict:
            total_mi += _split_entropy(kt, n_tips)
            matched_t.add(idx_t)
            matched_i.add(mt)

    rem_t = [(idx_t, mt, kt) for idx_t, (mt, kt) in enumerate(t_splits) if idx_t not in matched_t]
    rem_i = [(mi_mask, ki) for mi_mask, ki in i_splits if mi_mask not in matched_i]

    cand_pairs: list[tuple[float, int, int]] = []
    for idx_t, mt, kt in rem_t:
        for idx_i, (mi_mask, ki) in enumerate(rem_i):
            val = _mutual_clustering_info(mt, kt, mi_mask, ki, n_tips)
            if val > 1e-9:
                cand_pairs.append((val, idx_t, idx_i))
    cand_pairs.sort(reverse=True)

    used_t: set[int] = set()
    used_i: set[int] = set()
    for val, idx_t, idx_i in cand_pairs:
        if idx_t not in used_t and idx_i not in used_i:
            total_mi += val
            used_t.add(idx_t)
            used_i.add(idx_i)

    return min(1.0, round(total_mi / denom, 6))


class _LCAHelper:
    """Fast ancestor-depth & LCA lookup on a ``PhyloTree`` via Euler tour + sparse table."""

    def __init__(self, tree: PhyloTree, collapse_below_muts: float = 0.0) -> None:
        self.tree = tree
        n = tree.n_nodes
        self.topo_depth = [0] * n
        self.eff_parent = list(tree.parent)

        # Optionally contract near-zero internal branches so unresolved polytomies
        # have equal LCA depth rather than arbitrary binary tie-breaking
        for u in range(1, n):
            p = self.eff_parent[u]
            if tree.tip_alias[u] is None and tree.branch_muts[u] < collapse_below_muts and p >= 0:
                self.eff_parent[u] = p

        eff_children: list[list[int]] = [[] for _ in range(n)]
        for u in range(1, n):
            p = tree.parent[u]
            while p > 0 and tree.tip_alias[p] is None and tree.branch_muts[p] < collapse_below_muts:
                p = tree.parent[p]
            self.eff_parent[u] = p
            eff_children[p].append(u)

        # Iterative Euler tour on contracted tree
        euler: list[int] = []
        depths: list[int] = []
        first_occ: dict[int, int] = {}

        stack: list[tuple[int, int, int]] = [(0, 0, 0)]  # (node, depth, next_child_pos)
        self.topo_depth[0] = 0
        while stack:
            u, d, pos = stack[-1]
            if pos == 0:
                first_occ[u] = len(euler)
                self.topo_depth[u] = d
            euler.append(u)
            depths.append(d)
            chs = eff_children[u]
            if pos < len(chs):
                stack[-1] = (u, d, pos + 1)
                stack.append((chs[pos], d + 1, 0))
            else:
                stack.pop()

        m = len(euler)
        k_max = max(1, m.bit_length())
        st_node = [list(euler)]
        st_depth = [list(depths)]
        for k in range(1, k_max):
            span = 1 << (k - 1)
            prev_n = st_node[k - 1]
            prev_d = st_depth[k - 1]
            limit = m - (1 << k) + 1
            row_n = [0] * limit
            row_d = [0] * limit
            for i in range(limit):
                if prev_d[i] <= prev_d[i + span]:
                    row_d[i] = prev_d[i]
                    row_n[i] = prev_n[i]
                else:
                    row_d[i] = prev_d[i + span]
                    row_n[i] = prev_n[i + span]
            st_node.append(row_n)
            st_depth.append(row_d)

        self._first = first_occ
        self._st_node = st_node
        self._st_depth = st_depth
        self.tip_idx = tree.tips_map()

    def lca(self, al_a: str, al_b: str) -> int:
        ia = self._first[self.tip_idx[al_a]]
        ib = self._first[self.tip_idx[al_b]]
        if ia > ib:
            ia, ib = ib, ia
        length = ib - ia + 1
        k = length.bit_length() - 1
        span = 1 << k
        if self._st_depth[k][ia] <= self._st_depth[k][ib - span + 1]:
            return self._st_node[k][ia]
        return self._st_node[k][ib - span + 1]

    def lca_depth(self, al_a: str, al_b: str) -> int:
        return self.topo_depth[self.lca(al_a, al_b)]


def _quartet_state(lca_helper: _LCAHelper, a: str, b: str, c: str, d: str) -> int:
    """Return quartet topology among 4 tips on a rooted/unrooted tree:
      0 = unresolved star/polytomy
      1 = ab | cd
      2 = ac | bd
      3 = ad | bc
    In any tree, among the 3 pairings (ab+cd, ac+bd, ad+bc), two tie for the
    minimum LCA-depth sum and the resolved split strictly exceeds that minimum!
    """
    dab = lca_helper.lca_depth(a, b) + lca_helper.lca_depth(c, d)
    dac = lca_helper.lca_depth(a, c) + lca_helper.lca_depth(b, d)
    dad = lca_helper.lca_depth(a, d) + lca_helper.lca_depth(b, c)
    mx = max(dab, dac, dad)
    winners = [i for i, v in ((1, dab), (2, dac), (3, dad)) if v == mx]
    return winners[0] if len(winners) == 1 else 0


def _spearman_r(xs: list[float], ys: list[float]) -> float | None:
    n = len(xs)
    if n < 3:
        return None

    def _ranks(vals: list[float]) -> list[float]:
        order = sorted(range(n), key=vals.__getitem__)
        rk = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and vals[order[j + 1]] == vals[order[i]]:
                j += 1
            avg = 0.5 * (i + j) + 1.0
            for k in range(i, j + 1):
                rk[order[k]] = avg
            i = j + 1
        return rk

    rx = _ranks(xs)
    ry = _ranks(ys)
    mx = sum(rx) / n
    my = sum(ry) / n
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx <= 0 or syy <= 0:
        return None
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    return float(sxy / math.sqrt(sxx * syy))


def _adjusted_rand_and_f1(true_labels: list[str], pred_labels: list[str]) -> dict[str, float]:
    n = len(true_labels)
    if n < 2:
        return {"ari": 1.0, "precision": 1.0, "recall": 1.0, "f1": 1.0}

    cont: Counter = Counter()
    ca: Counter = Counter()
    cb: Counter = Counter()
    for lt, lp in zip(true_labels, pred_labels):
        cont[(lt, lp)] += 1
        ca[lt] += 1
        cb[lp] += 1

    def c2(k: int) -> int:
        return k * (k - 1) // 2

    tp = sum(c2(v) for v in cont.values())
    sa = sum(c2(v) for v in ca.values())
    sb = sum(c2(v) for v in cb.values())
    tot = c2(n)

    prec = tp / float(sb) if sb > 0 else 1.0
    rec = tp / float(sa) if sa > 0 else 1.0
    f1 = (2.0 * prec * rec / (prec + rec)) if (prec + rec) > 0 else 0.0

    exp_val = (sa * sb) / float(tot) if tot > 0 else 0.0
    max_val = 0.5 * (sa + sb)
    denom = max_val - exp_val
    ari = 1.0 if abs(denom) < 1e-12 and tp == sa == sb else ((tp - exp_val) / denom if abs(denom) >= 1e-12 else 0.0)
    return {
        "ari": round(float(ari), 6),
        "precision": round(float(prec), 6),
        "recall": round(float(rec), 6),
        "f1": round(float(f1), 6),
    }


def score_phylogeny(
    inferred_tree: PhyloTree,
    truth: TruthForest,
    min_mutations_collapse: float = 0.5,
    n_quartets: int = 8000,
    max_kc_pairs: int = 20000,
    tip_collection_dates: dict[str, _dt.date] | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    """Score an inferred ``PhyloTree`` against the pruned ABM truth genealogy."""
    inf_tips = set(inferred_tree.tips_map().keys()) & set(truth.parent_of.keys())
    n = len(inf_tips)
    if n < 4:
        return {"n_tips_scored": n, "error": "fewer than 4 shared tips with ground truth"}

    inf_pruned = inferred_tree.prune_to_tips(inf_tips)
    true_tree = build_pruned_truth_genealogy(truth, inf_tips, tip_collection_dates=tip_collection_dates)

    sorted_tips = sorted(inf_tips)
    z_rng = random.Random(0x5F3759DF)
    tip_zobrist = {al: z_rng.getrandbits(128) for al in sorted_tips}
    tip_to_bit = {al: i for i, al in enumerate(sorted_tips)} if n <= 2500 else None

    # 1. Split / Clade recovery (Strict & Polytomy-Collapsed) in O(N)
    true_clades, true_bits, _ = _extract_clades_fast(true_tree, tip_zobrist, tip_to_bit, min_branch_muts=0.0)
    inf_strict, _, _ = _extract_clades_fast(inf_pruned, tip_zobrist, None, min_branch_muts=0.0)
    inf_coll, inf_bits, _ = _extract_clades_fast(inf_pruned, tip_zobrist, tip_to_bit, min_branch_muts=min_mutations_collapse)

    t_keys = set(true_clades.keys())
    is_keys = set(inf_strict.keys())
    ic_keys = set(inf_coll.keys())

    strict_denom = len(t_keys) + len(is_keys)
    strict_nrf = (len(t_keys ^ is_keys) / float(strict_denom)) if strict_denom > 0 else 0.0

    tp_coll = len(t_keys & ic_keys)
    coll_prec = tp_coll / float(len(ic_keys)) if ic_keys else (1.0 if not t_keys else 0.0)
    coll_rec = tp_coll / float(len(t_keys)) if t_keys else 1.0
    coll_f1 = (2.0 * coll_prec * coll_rec / (coll_prec + coll_rec)) if (coll_prec + coll_rec) > 0 else 0.0
    coll_denom = len(t_keys) + len(ic_keys)
    coll_nrf = (len(t_keys ^ ic_keys) / float(coll_denom)) if coll_denom > 0 else 0.0

    info_rf = (
        compute_info_rf_similarity(true_bits, inf_bits, n)
        if (true_bits is not None and inf_bits is not None)
        else None
    )

    # 2. Quartet concordance & Kendall-Colijn rooted MRCA depth similarity
    lca_true = _LCAHelper(true_tree, collapse_below_muts=0.0)
    lca_inf = _LCAHelper(inf_pruned, collapse_below_muts=min_mutations_collapse)

    rng = random.Random(seed)
    q_match = 0
    q_resolved_true = 0
    q_resolved_match = 0
    for _ in range(n_quartets):
        a, b, c, d = rng.sample(sorted_tips, 4)
        qt = _quartet_state(lca_true, a, b, c, d)
        qi = _quartet_state(lca_inf, a, b, c, d)
        if qt == qi:
            q_match += 1
        if qt != 0:
            q_resolved_true += 1
            if qi == qt:
                q_resolved_match += 1

    quartet_all = q_match / float(n_quartets)
    quartet_res = (q_resolved_match / float(q_resolved_true)) if q_resolved_true > 0 else quartet_all

    total_pairs = n * (n - 1) // 2
    if total_pairs <= max_kc_pairs:
        pair_list = [(sorted_tips[i], sorted_tips[j]) for i in range(n) for j in range(i + 1, n)]
    else:
        pair_list = []
        for _ in range(max_kc_pairs):
            u, v = rng.sample(sorted_tips, 2)
            pair_list.append((u, v))

    d_true_vec = [float(lca_true.lca_depth(u, v)) for u, v in pair_list]
    d_inf_vec = [float(lca_inf.lca_depth(u, v)) for u, v in pair_list]
    kc_r = _spearman_r(d_true_vec, d_inf_vec)
    kc_l2 = math.sqrt(sum((a - b) ** 2 for a, b in zip(d_true_vec, d_inf_vec)) / len(pair_list))

    # 3. Cross-chain spurious merge rate in O(N) & Clade-cut Chain ARI
    first_comp: list[str | None] = [None] * inf_pruned.n_nodes
    is_mixed: list[bool] = [False] * inf_pruned.n_nodes
    for u in inf_pruned.postorder():
        al = inf_pruned.tip_alias[u]
        if al is not None and al in inf_tips:
            first_comp[u] = truth.comp_of.get(al, al)
        for c in inf_pruned.children[u]:
            if is_mixed[c]:
                is_mixed[u] = True
            elif first_comp[c] is not None:
                if first_comp[u] is None:
                    first_comp[u] = first_comp[c]
                elif first_comp[u] != first_comp[c]:
                    is_mixed[u] = True

    n_internal = 0
    n_merge_nodes = 0
    for u in range(1, inf_pruned.n_nodes):
        if inf_pruned.children[u]:
            n_internal += 1
            if is_mixed[u]:
                n_merge_nodes += 1
    cross_merge_pct = (100.0 * n_merge_nodes / n_internal) if n_internal > 0 else 0.0

    # Partition tree into clusters by cutting deep multi-mutation branches near root
    # via top-down propagation in O(N)
    node_cluster = ["C0"] * inf_pruned.n_nodes
    cid_counter = 0
    stack_cut = [0]
    while stack_cut:
        u = stack_cut.pop()
        for c in inf_pruned.children[u]:
            if u == 0 or (lca_inf.topo_depth[u] <= 3 and inf_pruned.branch_muts[c] >= 4.0):
                cid_counter += 1
                node_cluster[c] = f"C{cid_counter}"
            else:
                node_cluster[c] = node_cluster[u]
            stack_cut.append(c)

    inf_tip_idx_map = inf_pruned.tips_map()
    true_chain_vec = [truth.comp_of.get(al, al) for al in sorted_tips]
    pred_chain_vec = [node_cluster[inf_tip_idx_map[al]] for al in sorted_tips]
    chain_part = _adjusted_rand_and_f1(true_chain_vec, pred_chain_vec)

    # 4. Divergence Dating: Matched Clade TMRCA MAE & Pairwise Within-Chain Duration MAE
    tmrca_abs_errs: list[float] = []
    tmrca_signed_errs: list[float] = []
    ci_hits = 0
    ci_total = 0

    for mask in (t_keys & ic_keys):
        u_t = true_clades[mask]
        u_i = inf_coll[mask]
        dt_t = true_tree.num_date[u_t]
        dt_i = inf_pruned.num_date[u_i]
        if dt_t is not None and dt_i is not None:
            err_days = decimal_year_to_ordinal(dt_i) - decimal_year_to_ordinal(dt_t)
            tmrca_signed_errs.append(err_days)
            tmrca_abs_errs.append(abs(err_days))
            ci = inf_pruned.num_date_ci[u_i]
            if ci is not None:
                ci_total += 1
                lo_d = decimal_year_to_ordinal(ci[0])
                hi_d = decimal_year_to_ordinal(ci[1])
                true_d = decimal_year_to_ordinal(dt_t)
                if lo_d - 0.5 <= true_d <= hi_d + 0.5:
                    ci_hits += 1

    # Within-chain pairwise evolutionary duration error (Bryan's Analysis 2 equivalent)
    dur_abs_errs: list[float] = []
    dur_signed_errs: list[float] = []
    for u, v in pair_list:
        if truth.comp_of.get(u, u) != truth.comp_of.get(v, "__none__"):
            continue
        anc_t = lca_true.lca(u, v)
        anc_i = lca_inf.lca(u, v)
        t_mrca_true = true_tree.num_date[anc_t]
        t_mrca_inf = inf_pruned.num_date[anc_i]
        tu_true = truth.exp_ordinal.get(u)
        tv_true = truth.exp_ordinal.get(v)
        tu_inf = inf_pruned.num_date[inf_pruned.tips_map()[u]]
        tv_inf = inf_pruned.num_date[inf_pruned.tips_map()[v]]
        if None not in (t_mrca_true, t_mrca_inf, tu_true, tv_true, tu_inf, tv_inf):
            mrca_t_ord = decimal_year_to_ordinal(t_mrca_true)  # type: ignore[arg-type]
            mrca_i_ord = decimal_year_to_ordinal(t_mrca_inf)   # type: ignore[arg-type]
            ui_ord = decimal_year_to_ordinal(tu_inf)           # type: ignore[arg-type]
            vi_ord = decimal_year_to_ordinal(tv_inf)           # type: ignore[arg-type]
            true_dur = (tu_true - mrca_t_ord) + (tv_true - mrca_t_ord)  # type: ignore[operator]
            inf_dur = (ui_ord - mrca_i_ord) + (vi_ord - mrca_i_ord)
            diff = inf_dur - true_dur
            dur_signed_errs.append(diff)
            dur_abs_errs.append(abs(diff))

    return {
        "n_tips_scored": n,
        "n_true_clades": len(t_keys),
        "n_inferred_clades_strict": len(is_keys),
        "n_inferred_clades_collapsed": len(ic_keys),
        "strict_nrf": round(strict_nrf, 6),
        "collapsed_nrf": round(coll_nrf, 6),
        "collapsed_clade_precision": round(coll_prec, 6),
        "collapsed_clade_recall": round(coll_rec, 6),
        "collapsed_clade_f1": round(coll_f1, 6),
        "info_rf_similarity": info_rf,
        "quartet_concordance": round(quartet_all, 6),
        "quartet_resolved_concordance": round(quartet_res, 6),
        "kc_spearman_r": round(kc_r, 6) if kc_r is not None else None,
        "kc_norm_dist": round(kc_l2, 4),
        "cross_chain_merge_pct": round(cross_merge_pct, 4),
        "tree_chain_ari": chain_part["ari"],
        "tree_chain_f1": chain_part["f1"],
        "clade_tmrca_n": len(tmrca_abs_errs),
        "clade_tmrca_mae_days": round(sum(tmrca_abs_errs) / len(tmrca_abs_errs), 3) if tmrca_abs_errs else None,
        "clade_tmrca_bias_days": round(sum(tmrca_signed_errs) / len(tmrca_signed_errs), 3) if tmrca_signed_errs else None,
        "clade_tmrca_95ci_cov": round(ci_hits / float(ci_total), 4) if ci_total > 0 else None,
        "pairwise_duration_n": len(dur_abs_errs),
        "pairwise_duration_mae_days": round(sum(dur_abs_errs) / len(dur_abs_errs), 3) if dur_abs_errs else None,
        "pairwise_duration_bias_days": round(sum(dur_signed_errs) / len(dur_signed_errs), 3) if dur_signed_errs else None,
    }


# ---------------------------------------------------------------------------
# Transmission Who-Infected-Whom Scoring (JUNIPER / BREATH / outbreaker2)
# ---------------------------------------------------------------------------
def _score_single_donor_mode(
    true_donor_of: dict[str, str],
    post_probs_of: dict[str, dict[str, float]],
    prefix: str,
) -> dict[str, Any]:
    recipients = sorted(true_donor_of.keys())
    n = len(recipients)
    if n == 0:
        return {}

    prob_true_sum = 0.0
    top1_hits = 0
    obs_true_n = 0
    obs_true_hits = 0
    obs_pred_n = 0
    obs_pred_tp = 0
    unobs_true_n = 0
    unobs_true_hits = 0
    set95_hits = 0
    set95_sizes: list[int] = []

    # Calibration decile bins [0..9] over candidate predictions > 0
    bin_pred_sum = [0.0] * 10
    bin_true_sum = [0.0] * 10
    bin_counts = [0] * 10

    for j in recipients:
        t_don = true_donor_of[j]
        dist = post_probs_of.get(j) or {UNOBSERVED: 1.0}

        # Normalize defensively
        tot = sum(dist.values())
        if tot <= 0.0:
            dist = {UNOBSERVED: 1.0}
        elif abs(tot - 1.0) > 1e-6:
            dist = {k: v / tot for k, v in dist.items()}

        p_true = dist.get(t_don, 0.0)
        prob_true_sum += p_true

        ranked = sorted(dist.items(), key=lambda kv: (-kv[1], kv[0] != UNOBSERVED, kv[0]))
        top_don = ranked[0][0]
        is_hit = top_don == t_don
        if is_hit:
            top1_hits += 1

        if t_don != UNOBSERVED:
            obs_true_n += 1
            if is_hit:
                obs_true_hits += 1
        else:
            unobs_true_n += 1
            if is_hit:
                unobs_true_hits += 1

        if top_don != UNOBSERVED:
            obs_pred_n += 1
            if is_hit:
                obs_pred_tp += 1

        # Smallest 95% credible donor set (Specht et al. 2025)
        cum = 0.0
        cset: set[str] = set()
        for cand_d, cand_p in ranked:
            if cand_p <= 0.0:
                break
            cset.add(cand_d)
            cum += cand_p
            if cum >= 0.95 - 1e-9:
                break
        set95_sizes.append(len(cset))
        if t_don in cset:
            set95_hits += 1

        # Posterior calibration across observed candidate edges with p >= 0.01 (BREATH Fig. 4)
        for cand_d, cand_p in dist.items():
            if cand_p < 0.01:
                continue
            b = min(9, int(cand_p * 10.0))
            bin_pred_sum[b] += cand_p
            bin_true_sum[b] += 1.0 if cand_d == t_don else 0.0
            bin_counts[b] += 1

    set95_sizes.sort()
    med_95 = set95_sizes[len(set95_sizes) // 2] if set95_sizes else 0

    tot_cal = sum(bin_counts)
    ece = 0.0
    if tot_cal > 0:
        for b in range(10):
            if bin_counts[b] > 0:
                avg_p = bin_pred_sum[b] / bin_counts[b]
                emp_p = bin_true_sum[b] / bin_counts[b]
                ece += (bin_counts[b] / float(tot_cal)) * abs(avg_p - emp_p)

    return {
        f"{prefix}_n_recipients": n,
        f"{prefix}_n_true_observed_donors": obs_true_n,
        f"{prefix}_donor_prob_mean": round(prob_true_sum / n, 6),
        f"{prefix}_top1_donor_acc": round(top1_hits / n, 6),
        f"{prefix}_obs_donor_recall": round(obs_true_hits / obs_true_n, 6) if obs_true_n > 0 else None,
        f"{prefix}_obs_link_precision": round(obs_pred_tp / obs_pred_n, 6) if obs_pred_n > 0 else None,
        f"{prefix}_unobs_root_specificity": round(unobs_true_hits / unobs_true_n, 6) if unobs_true_n > 0 else None,
        f"{prefix}_donor_95set_cov": round(set95_hits / n, 6),
        f"{prefix}_donor_95set_med_size": med_95,
        f"{prefix}_calibration_ece": round(ece, 6),
    }


def score_transmission(
    donor_probs: dict[str, dict[str, dict[str, float]]],
    truth: TruthForest,
    tip_aliases: Iterable[str] | None = None,
) -> dict[str, Any]:
    """Score posterior or point who-infected-whom distributions against ground truth."""
    all_recs = set(tip_aliases) if tip_aliases is not None else (
        set(donor_probs.get("direct", {}).keys()) | set(donor_probs.get("indirect", {}).keys())
    )
    valid_recs = sorted(all_recs & set(truth.parent_of.keys()))
    if not valid_recs:
        return {"error": "no recipients matched ground truth"}

    true_direct, true_indirect, _, _ = build_truth_donors(truth, valid_recs)
    out: dict[str, Any] = {}
    out.update(_score_single_donor_mode(true_direct, donor_probs.get("direct", {}), "direct"))
    out.update(_score_single_donor_mode(true_indirect, donor_probs.get("indirect", {}), "indirect"))
    return out


# ---------------------------------------------------------------------------
# MicrobeTrace SNP-Network & SNP+Date Heuristic Baseline
# ---------------------------------------------------------------------------
def score_snp_network_and_heuristic(
    fasta_path: str | Path,
    truth: TruthForest,
    collection_dates: dict[str, _dt.date],
    lookup_map: dict[str, str] | None = None,
    snp_thresholds: Iterable[int] = (0, 1, 2, 3),
    direct_max_snps: int = 1,
    indirect_max_snps: int = 3,
) -> tuple[list[dict[str, Any]], dict[str, dict[str, dict[str, float]]]]:
    """Score MicrobeTrace single-linkage SNP clustering against both statewide
    ``component_id`` and local cluster introductions, and build a naive field-epi
    ``SNP + collection_date`` donor probability baseline."""
    seqs: dict[str, str] = {}
    with _open_text(fasta_path, "rt") as fh:
        cur_id: str | None = None
        buf: list[str] = []
        for line in fh:
            if line.startswith(">"):
                if cur_id is not None:
                    al = alias_from_any_id(cur_id, lookup_map)
                    if al and al in truth.parent_of:
                        seqs[al] = "".join(buf).upper()
                cur_id = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur_id is not None:
            al = alias_from_any_id(cur_id, lookup_map)
            if al and al in truth.parent_of:
                seqs[al] = "".join(buf).upper()

    tips = sorted(seqs.keys())
    n = len(tips)
    if n < 2:
        return [], {"direct": {}, "indirect": {}}

    _, _, _, local_intro_of = build_truth_donors(truth, tips)

    # Precompute ACGT pairwise SNP matrix
    dist_mat: dict[tuple[int, int], int] = {}
    for i in range(n):
        si = seqs[tips[i]]
        for j in range(i + 1, n):
            sj = seqs[tips[j]]
            d = sum(1 for a, b in zip(si, sj) if a != b and a in "ACGT" and b in "ACGT")
            dist_mat[(i, j)] = d

    true_comp = [truth.comp_of.get(u, u) for u in tips]
    true_local = [local_intro_of.get(u, u) for u in tips]

    cluster_rows: list[dict[str, Any]] = []
    for tau in snp_thresholds:
        parent_uf = list(range(n))

        def _find(x: int) -> int:
            while parent_uf[x] != x:
                parent_uf[x] = parent_uf[parent_uf[x]]
                x = parent_uf[x]
            return x

        n_edges = 0
        same_comp_edges = 0
        for (i, j), d in dist_mat.items():
            if d <= tau:
                n_edges += 1
                if true_comp[i] == true_comp[j]:
                    same_comp_edges += 1
                ri, rj = _find(i), _find(j)
                if ri != rj:
                    parent_uf[max(ri, rj)] = min(ri, rj)

        pred = [f"G{_find(i)}" for i in range(n)]
        comp_stats = _adjusted_rand_and_f1(true_comp, pred)
        loc_stats = _adjusted_rand_and_f1(true_local, pred)
        cluster_rows.append({
            "snp_threshold": int(tau),
            "n_tips": n,
            "n_edges": n_edges,
            "edge_same_chain_precision": round(same_comp_edges / float(n_edges), 6) if n_edges > 0 else 1.0,
            "n_clusters": len(set(pred)),
            "chain_ari": comp_stats["ari"],
            "chain_f1": comp_stats["f1"],
            "chain_precision": comp_stats["precision"],
            "chain_recall": comp_stats["recall"],
            "local_intro_ari": loc_stats["ari"],
            "local_intro_f1": loc_stats["f1"],
        })

    # Field-epi SNP + collection_date donor heuristic:
    # Candidate donors i for recipient j must have collection_date(i) <= collection_date(j)
    # (with strict preference for earlier dates + minimal SNPs, weighted by Gamma GI ~4.7d)
    direct_probs: dict[str, dict[str, float]] = {}
    indirect_probs: dict[str, dict[str, float]] = {}

    for j_idx, rec in enumerate(tips):
        dj = collection_dates.get(rec)
        dir_weights: dict[str, float] = {}
        ind_weights: dict[str, float] = {}
        for i_idx, don in enumerate(tips):
            if i_idx == j_idx:
                continue
            d_snp = dist_mat[(min(i_idx, j_idx), max(i_idx, j_idx))]
            if d_snp > indirect_max_snps:
                continue
            di = collection_dates.get(don)
            if dj is not None and di is not None:
                dt_days = (dj - di).days
            else:
                dt_days = 4
            # Plausible serial interval window: -1 to +21 days for direct, -1 to +45 for indirect
            if -1 <= dt_days <= 21 and d_snp <= direct_max_snps:
                eff_dt = max(0.5, float(dt_days + 1.5))
                w = math.exp(-1.2 * d_snp) * (eff_dt ** 2.0) * math.exp(-eff_dt / 2.2)
                dir_weights[don] = w
            if -1 <= dt_days <= 45 and d_snp <= indirect_max_snps:
                eff_dt = max(0.5, float(dt_days + 1.5))
                w_ind = math.exp(-0.9 * d_snp) * math.exp(-abs(eff_dt - 5.0 * max(1, d_snp)) / 8.0)
                ind_weights[don] = w_ind

        if not dir_weights:
            direct_probs[rec] = {UNOBSERVED: 1.0}
        else:
            # Allow some mass on UNOBSERVED when closest candidate is >=1 SNP or >10d away
            w_tot = sum(dir_weights.values())
            unobs_prior = 0.25 * w_tot
            denom = w_tot + unobs_prior
            direct_probs[rec] = {k: v / denom for k, v in dir_weights.items()}
            direct_probs[rec][UNOBSERVED] = unobs_prior / denom

        if not ind_weights:
            indirect_probs[rec] = {UNOBSERVED: 1.0}
        else:
            w_tot = sum(ind_weights.values())
            unobs_prior = 0.15 * w_tot
            denom = w_tot + unobs_prior
            indirect_probs[rec] = {k: v / denom for k, v in ind_weights.items()}
            indirect_probs[rec][UNOBSERVED] = unobs_prior / denom

    return cluster_rows, {"direct": direct_probs, "indirect": indirect_probs}
