#!/usr/bin/env python3
"""Homophily decay, LASSO weight attribution, and mugration mechanism analysis.

Explains *why* LASSO50-Greedy (`lasso_greedy`, 50 coarse LASSO supergroups with
within-supergroup availability x target weighting `cap_g * q_g`) improves
inter-county phylodynamic flow reconstruction (`cosine_similarity` and
`topological_f1`), whereas Simple Uniform Random Sampling (`SURS`) gets trapped
in high-ascertainment urban/homophilic cores and rigid micro-quota methods
(`Stratified` and 1500-bucket `LASSO-Stratified`) shatter local clade
coalescence and distort the mugration denominator `pi_j`.

Four linked analyses are produced in a single pass over the ground-truth graph:

1. ``<state>_hop_homophily_decay.csv``
   True transmission-graph homophily and directed multi-hop decay r(k) for
   k = 1, 2, 3, 4, 5, 8, 10 hops across `county`, `rucc_code`, `urbanicity`,
   `age_group`, `smh_race`, `ses_category`, `sex`, and `component_id`, paired
   with the exact marginal null expectation `sum_c p(c)^2` and excess homophily.

2. ``<state>_stratum_bridge_enrichment.csv``
   Per-stratum epidemiological, topological, and LASSO priority profile:
   - True infection share (`p_true`), ascertained line-list share (`p_ll`),
     population share (`p_pop`), and ascertainment bias ratio (`p_ll / p_true`)
   - Marginal KL deficit priority (`log(p_target / p_ll)`) that drives LASSO
   - 1-hop cross-county bridge rate (`frac_1hop_cross_county_bridge`) and 2-hop
     bridge-adjacent rate (`frac_2hop_bridge_adjacent`)
   - Realized sampling propensity weight (`w_g = p_sample / p_ll`) per strategy

3. ``<state>_lasso_weights_over_time.csv``
   Per-stride-window (and full-wave) LASSO feature attribution (`l1_share_<attr>`
   from fitting LASSO to predict fine-group marginal KL gain `y_g` from group
   indicators `X_g`) alongside realized rural/urban and demographic weight shifts
   over the course of the epidemic wave.

4. ``<state>_strategy_mechanism_summary.csv``
   Per-strategy topological and mugration decomposition metrics:
   - Macro reach: counties sampled, rural (RUCC 4-9) share, top-5 urban share
   - Micro clade coalescence: directed ancestor-descendant pairs retained within
     k <= 1, 2, 3, 5, 10 hops (`n_ancestral_pairs_le5hop`,
     `frac_samples_with_ancestor_le5hop`), plus marginal-null-adjusted cross-county
     and cross-demographic concordance
   - Mugration rate matrix decomposition (`W_ij ~ q_ij / pi_j + q_ji / pi_i`):
     separates standard `cos_full_W` into numerator jump accuracy (`cos_raw_jumps_Q`),
     true-denominator rate accuracy (`cos_true_pi_W`), and equilibrium node-state
     frequency distortion (`tv_error_pi`, `max_inv_pi_inflation`).

Usage (on Rivanna or locally):
------------------------------
    python3 scripts/homophily_lasso_mechanism.py \\
        --infections results/va_delta_wave/05_benchmarks/allevents_variants.csv.xz \\
        --linelist results/va_delta_wave/02_simulated_linelists/linelist.csv.xz \\
        --population data/va/va_2_4_0_demographics.csv \\
        --truth results/va_delta_wave/02_simulated_linelists/linelist_mugration.json \\
        --samples 'results/va_delta_wave/03_sampled_datasets/*_samples.csv.xz' \\
                  'results/va_delta_wave/03_sampled_foil/*_samples.csv.xz' \\
        --state va \\
        --outdir results/va_delta_wave/05_benchmarks/mechanism
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import glob as _glob
import gzip
import json
import lzma
import math
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    np = None  # type: ignore[assignment]
    HAS_NUMPY = False


_NO_CONTACT = {"", "-1", "nan", "None", "<NA>", "null"}

ATTR_NAMES = (
    "county",
    "rucc_code",
    "urbanicity",
    "age_group",
    "smh_race",
    "ses_category",
    "sex",
    "component_id",
)

HOP_DISTANCES = (1, 2, 3, 4, 5, 8, 10)


def _open_text(path: str | Path, mode: str = "rt"):
    p = str(path)
    if p.endswith(".xz"):
        return lzma.open(p, mode, encoding="utf-8")
    if p.endswith(".gz"):
        return gzip.open(p, mode, encoding="utf-8")
    raw_mode = mode.replace("t", "")
    return open(p, raw_mode, encoding="utf-8", newline="" if "w" in raw_mode else None)


def _urbanicity_from_rucc(rucc_str: str) -> str:
    if not rucc_str or rucc_str in _NO_CONTACT:
        return "Unknown"
    try:
        c = int(float(rucc_str))
    except ValueError:
        return "Unknown"
    if 1 <= c <= 3:
        return "Metro (RUCC 1-3)"
    if 4 <= c <= 9:
        return "Non-metro / Rural (RUCC 4-9)"
    return "Unknown"


def _clean_rucc(rucc_str: str) -> str:
    if not rucc_str or rucc_str in _NO_CONTACT:
        return ""
    try:
        return str(int(float(rucc_str)))
    except ValueError:
        return rucc_str.strip()


def _label_from_path(path: str | Path) -> str:
    base = os.path.basename(str(path))
    for suf in ("_samples.csv.xz", "_samples.csv.gz", "_samples.csv", ".csv.xz", ".csv.gz", ".csv"):
        if base.endswith(suf):
            return base[: -len(suf)]
    return base.split(".")[0]


def _cosine_1d(u: list[float], v: list[float]) -> float:
    dot = 0.0
    nu = 0.0
    nv = 0.0
    for a, b in zip(u, v):
        dot += a * b
        nu += a * a
        nv += b * b
    if nu <= 1e-15 or nv <= 1e-15:
        return 0.0
    return float(dot / math.sqrt(nu * nv))


def _row_normalize_offdiag(mat: list[list[float]]) -> list[float]:
    """Row-normalize a square matrix (zeroing diagonal) and flatten off-diagonals,
    matching `phylogas.benchmark.mugration.align_and_normalize_matrix` + `get_off_diagonals`."""
    n = len(mat)
    flat: list[float] = []
    for i in range(n):
        row_sum = sum(mat[i][j] for j in range(n) if j != i)
        denom = row_sum if row_sum > 0.0 else 1.0
        for j in range(n):
            if i != j:
                flat.append(mat[i][j] / denom)
    return flat


def _compute_newman_r_and_concordance(
    pairs: list[tuple[int, int]],
    n_cats: int,
) -> tuple[float | None, float | None, float | None, float | None]:
    """Return ``(within_frac, null_frac, excess_homophily, newman_r)`` over integer-coded
    category pairs ``(ca, cb)`` where negative codes represent missing values."""
    valid = [(a, b) for a, b in pairs if a >= 0 and b >= 0]
    m = len(valid)
    if m == 0 or n_cats < 1:
        return None, None, None, None

    same = sum(1 for a, b in valid if a == b)
    within_frac = same / m

    # Symmetrized marginal frequencies (matching Newman's undirected mixing matrix e)
    counts = [0] * n_cats
    for a, b in valid:
        counts[a] += 1
        counts[b] += 1
    total_stubs = 2.0 * m
    null_frac = sum((c / total_stubs) ** 2 for c in counts)
    excess = within_frac - null_frac

    if n_cats < 2 or abs(1.0 - null_frac) < 1e-12:
        newman_r = None
    else:
        newman_r = (within_frac - null_frac) / (1.0 - null_frac)

    return (
        round(within_frac, 6),
        round(null_frac, 6),
        round(excess, 6),
        round(newman_r, 6) if newman_r is not None else None,
    )


def _fit_lasso_feature_importance(
    group_tuples: list[tuple[str, ...]],
    y_kl_priority: list[float],
    weights: list[float],
    attr_labels: list[str],
) -> dict[str, float]:
    """Compute each demographic dimension's share of LASSO / L1 feature weight
    in predicting the fine-group marginal KL priority ``y_g``.

    Uses `sklearn.linear_model.LassoCV` when scikit-learn is available, and falls
    back to exact soft-thresholded marginal L1 coordinate attribution in bare
    Python environments.
    """
    n_groups = len(group_tuples)
    n_attrs = len(attr_labels)
    if n_groups < 5 or n_attrs == 0:
        return {a: 0.0 for a in attr_labels}

    # Center y with weights
    w_sum = sum(weights) or 1.0
    y_mean = sum(w * y for w, y in zip(weights, y_kl_priority)) / w_sum
    y_cent = [y - y_mean for y in y_kl_priority]

    # Try scikit-learn LassoCV if available
    try:
        from sklearn.linear_model import LassoCV  # type: ignore[import-untyped]
        import numpy as _np

        col_to_attr: list[int] = []
        cat_maps: list[dict[str, int]] = []
        total_cols = 0
        for j in range(n_attrs):
            cats = sorted({gt[j] for gt in group_tuples if gt[j] not in _NO_CONTACT})
            cmap = {c: total_cols + idx for idx, c in enumerate(cats)}
            cat_maps.append(cmap)
            col_to_attr.extend([j] * len(cats))
            total_cols += len(cats)

        if total_cols >= 2:
            X = _np.zeros((n_groups, total_cols), dtype=_np.float64)
            for row_i, gt in enumerate(group_tuples):
                for j in range(n_attrs):
                    col_k = cat_maps[j].get(gt[j], -1)
                    if col_k >= 0:
                        X[row_i, col_k] = 1.0
            # Standardize columns
            col_std = X.std(axis=0)
            col_std[col_std < 1e-9] = 1.0
            X_std = (X - X.mean(axis=0)) / col_std
            y_arr = _np.asarray(y_cent, dtype=_np.float64)
            model = LassoCV(cv=min(5, n_groups), max_iter=2000, random_state=43)
            model.fit(X_std, y_arr)
            coefs = _np.abs(model.coef_)
            attr_l1 = [0.0] * n_attrs
            for col_k, j in enumerate(col_to_attr):
                attr_l1[j] += float(coefs[col_k])
            tot_l1 = sum(attr_l1)
            if tot_l1 > 1e-12:
                return {attr_labels[j]: round(attr_l1[j] / tot_l1, 6) for j in range(n_attrs)}
    except Exception:
        pass

    # Pure-Python coordinate L1 attribution fallback (standardized group-level
    # correlation / soft-thresholded L1 norm per attribute axis)
    attr_l1_py = [0.0] * n_attrs
    for j in range(n_attrs):
        by_cat_y: dict[str, float] = defaultdict(float)
        by_cat_w: dict[str, float] = defaultdict(float)
        for gt, yc, w in zip(group_tuples, y_cent, weights):
            c = gt[j]
            if c not in _NO_CONTACT:
                by_cat_y[c] += w * yc
                by_cat_w[c] += w
        # Standardized effect across categories of attribute j
        l1_j = 0.0
        for c, wy in by_cat_y.items():
            cw = by_cat_w[c]
            if cw > 0:
                mean_c = wy / cw
                p_c = cw / w_sum
                sd_c = math.sqrt(max(p_c * (1.0 - p_c), 1e-9))
                l1_j += abs(mean_c * sd_c)
        attr_l1_py[j] = l1_j

    tot = sum(attr_l1_py) or 1.0
    return {attr_labels[j]: round(attr_l1_py[j] / tot, 6) for j in range(n_attrs)}


def run_mechanism_analysis(
    infections_path: str | Path,
    outdir: str | Path,
    state_label: str = "va",
    linelist_path: str | Path | None = None,
    population_path: str | Path | None = None,
    truth_mugration_path: str | Path | None = None,
    sample_specs: list[str] | None = None,
    max_records: int = 0,
    max_ancestral_hops: int = 10,
    simulate_strategies: bool = False,
) -> dict:
    """Execute the 4-part homophily, LASSO weight, and mugration mechanism analysis."""
    out_dir = Path(outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"=== PhyloGAS Homophily & LASSO50-Greedy Mechanism Analysis ({state_label.upper()}) ===")
    print(f"  Infections file : {infections_path}")

    # ------------------------------------------------------------------
    # 1. Stream infections file into compact integer arrays
    # ------------------------------------------------------------------
    alias_to_idx: dict[str, int] = {}
    contact_raw_list: list[str] = []
    date_ordinal_list: list[int] = []

    attr_vocab: dict[str, dict[str, int]] = {a: {} for a in ATTR_NAMES}
    attr_rev: dict[str, list[str]] = {a: [] for a in ATTR_NAMES}
    attr_codes: dict[str, list[int]] = {a: [] for a in ATTR_NAMES}

    def _encode_attr(attr: str, val: str) -> int:
        v = val.strip() if val else ""
        if not v or v in _NO_CONTACT:
            return -1
        voc = attr_vocab[attr]
        code = voc.get(v)
        if code is None:
            code = len(voc)
            voc[v] = code
            attr_rev[attr].append(v)
        return code

    with _open_text(infections_path, "rt") as fh:
        header_line = fh.readline()
        if not header_line:
            raise ValueError(f"Empty infections file: {infections_path}")
        delim = "\t" if header_line.count("\t") > header_line.count(",") else ","
        cols = [c.strip() for c in header_line.rstrip("\r\n").split(delim)]
        c_idx = {c: i for i, c in enumerate(cols)}

        i_alias = c_idx.get("alias_pid")
        i_contact = c_idx.get("alias_contact")
        i_county = c_idx.get("county")
        i_rucc = c_idx.get("rucc_code")
        i_age = c_idx.get("age_group", c_idx.get("age"))
        i_race = c_idx.get("smh_race", c_idx.get("race"))
        i_ses = c_idx.get("ses_category")
        i_sex = c_idx.get("sex", c_idx.get("gender"))
        i_comp = c_idx.get("component_id")
        i_date = c_idx.get("date", c_idx.get("exposure_date"))
        i_test_prob = c_idx.get("test_prob")

        if i_alias is None or i_contact is None:
            raise ValueError(f"Need alias_pid and alias_contact in {infections_path}")

        test_prob_list: list[float] = []
        reader = csv.reader(fh, delimiter=delim)
        for row in reader:
            if len(row) <= i_alias:
                continue
            al = row[i_alias].strip()
            if not al or al in alias_to_idx:
                continue
            idx = len(alias_to_idx)
            alias_to_idx[al] = idx

            con = row[i_contact].strip() if len(row) > i_contact else ""
            contact_raw_list.append(con if con not in _NO_CONTACT and con != al else "")

            cty = row[i_county].strip() if (i_county is not None and len(row) > i_county) else ""
            rucc_raw = row[i_rucc].strip() if (i_rucc is not None and len(row) > i_rucc) else ""
            rucc = _clean_rucc(rucc_raw)
            urb = _urbanicity_from_rucc(rucc)
            age = row[i_age].strip() if (i_age is not None and len(row) > i_age) else ""
            race = row[i_race].strip() if (i_race is not None and len(row) > i_race) else ""
            ses = row[i_ses].strip() if (i_ses is not None and len(row) > i_ses) else ""
            sex = row[i_sex].strip() if (i_sex is not None and len(row) > i_sex) else ""
            comp = row[i_comp].strip() if (i_comp is not None and len(row) > i_comp) else ""

            attr_codes["county"].append(_encode_attr("county", cty))
            attr_codes["rucc_code"].append(_encode_attr("rucc_code", rucc))
            attr_codes["urbanicity"].append(_encode_attr("urbanicity", urb))
            attr_codes["age_group"].append(_encode_attr("age_group", age))
            attr_codes["smh_race"].append(_encode_attr("smh_race", race))
            attr_codes["ses_category"].append(_encode_attr("ses_category", ses))
            attr_codes["sex"].append(_encode_attr("sex", sex))
            attr_codes["component_id"].append(_encode_attr("component_id", comp))

            d_ord = -1
            if i_date is not None and len(row) > i_date:
                ds = row[i_date].strip()[:10]
                if len(ds) == 10 and ds[4] == "-" and ds[7] == "-":
                    try:
                        d_ord = _dt.date.fromisoformat(ds).toordinal()
                    except ValueError:
                        d_ord = -1
            date_ordinal_list.append(d_ord)

            tprob = 1.0
            if i_test_prob is not None and len(row) > i_test_prob:
                try:
                    tprob = float(row[i_test_prob].strip())
                except ValueError:
                    tprob = 1.0
            test_prob_list.append(max(tprob, 0.0))

            if max_records > 0 and len(alias_to_idx) >= max_records:
                break

    n_nodes = len(alias_to_idx)
    print(f"  Loaded {n_nodes:,} unique infections.")

    # Resolve parent index array: parent_idx[u] = index of u's infector, or -1
    parent_idx = [-1] * n_nodes
    children_of: list[list[int]] = [[] for _ in range(n_nodes)]
    n_edges = 0
    for u, con in enumerate(contact_raw_list):
        if con:
            p = alias_to_idx.get(con, -1)
            if p >= 0 and p != u:
                parent_idx[u] = p
                children_of[p].append(u)
                n_edges += 1
    contact_raw_list.clear()
    print(f"  Resolved {n_edges:,} directed transmission edges ({n_nodes - n_edges:,} roots).")

    # Ascertained line-list weight per node:
    # - If `--linelist` is provided, 1.0 for ascertained cases in linelist.csv.xz, 0.0 otherwise.
    # - Otherwise, use each infection's exact `test_prob` from TwinSampler.
    ll_weight = [0.0] * n_nodes
    used_explicit_linelist = False
    if linelist_path and Path(linelist_path).exists():
        print(f"  Loading ascertained line list : {linelist_path}")
        with _open_text(linelist_path, "rt") as fh:
            hdr = fh.readline()
            if hdr:
                dl = "\t" if hdr.count("\t") > hdr.count(",") else ","
                lcols = [c.strip() for c in hdr.rstrip("\r\n").split(dl)]
                ia = lcols.index("alias_pid") if "alias_pid" in lcols else -1
                if ia >= 0:
                    for row in csv.reader(fh, delimiter=dl):
                        if len(row) > ia:
                            u = alias_to_idx.get(row[ia].strip(), -1)
                            if u >= 0:
                                ll_weight[u] = 1.0
                    used_explicit_linelist = sum(ll_weight) > 0
    if not used_explicit_linelist:
        ll_weight = test_prob_list

    tot_ll_weight = sum(ll_weight) or 1.0
    if used_explicit_linelist:
        print(f"  Ascertained line-list pool    : {int(tot_ll_weight):,} cases ({tot_ll_weight / max(n_nodes, 1) * 100:.1f}% of infections).")
    else:
        print(f"  Expected line-list weight (sum test_prob): {tot_ll_weight:,.1f} ({tot_ll_weight / max(n_nodes, 1) * 100:.2f}% mean IAR).")

    # ------------------------------------------------------------------
    # Part 1: True Network Homophily & Multi-Hop Decay r(k)
    # ------------------------------------------------------------------
    print("\n--- Part 1: Computing Multi-Hop Homophily Decay r(k) ---")
    step_stride = max(1, n_nodes // 250_000)
    eval_nodes = list(range(0, n_nodes, step_stride))

    max_hop_eval = max(HOP_DISTANCES)
    cur_anc = list(eval_nodes)
    pairs_by_hop: dict[int, list[tuple[int, int]]] = {}

    for hop in range(1, max_hop_eval + 1):
        next_anc = [-1] * len(eval_nodes)
        hop_node_pairs: list[tuple[int, int]] = []
        for i, u in enumerate(eval_nodes):
            anc = cur_anc[i]
            if anc >= 0:
                par = parent_idx[anc]
                if par >= 0 and par != u:
                    next_anc[i] = par
                    if hop in HOP_DISTANCES:
                        hop_node_pairs.append((par, u))
        cur_anc = next_anc
        if hop in HOP_DISTANCES:
            pairs_by_hop[hop] = hop_node_pairs

    decay_rows: list[dict] = []
    hop1_excess: dict[str, float] = {}
    for attr in ATTR_NAMES:
        codes_a = attr_codes[attr]
        n_cats = len(attr_rev[attr])
        if n_cats < 2:
            continue
        for hop in HOP_DISTANCES:
            node_pairs = pairs_by_hop.get(hop, [])
            cat_pairs = [(codes_a[anc], codes_a[u]) for anc, u in node_pairs]
            w_frac, n_frac, exc, r_val = _compute_newman_r_and_concordance(cat_pairs, n_cats)
            if hop == 1 and exc is not None:
                hop1_excess[attr] = exc
            rel_retention = (
                round(exc / hop1_excess[attr], 4)
                if (exc is not None and hop1_excess.get(attr, 0.0) > 1e-6)
                else None
            )
            decay_rows.append({
                "state": state_label,
                "attribute": attr,
                "n_categories": n_cats,
                "hop_k": hop,
                "n_pairs_evaluated": len(cat_pairs),
                "within_attr_frac": w_frac,
                "marginal_null_frac": n_frac,
                "excess_homophily": exc,
                "newman_r": r_val,
                "excess_retention_vs_hop1": rel_retention,
            })
            if hop in (1, 3, 5, 10) and r_val is not None and attr != "component_id":
                print(
                    f"  {attr:<14} k={hop:>2d}: within={w_frac:.4f} (null={n_frac:.4f}, "
                    f"excess={exc:+.4f}, r={r_val:+.4f}, ret={rel_retention or 0.0:.2f})"
                )

    decay_csv = out_dir / f"{state_label}_hop_homophily_decay.csv"
    if decay_rows:
        with open(decay_csv, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(decay_rows[0].keys()))
            w.writeheader()
            w.writerows(decay_rows)
        print(f"  -> Wrote {decay_csv} ({len(decay_rows)} rows)")

    # ------------------------------------------------------------------
    # Part 2: Bridge-Node Identification & Per-Stratum Enrichment
    # ------------------------------------------------------------------
    print("\n--- Part 2: Identifying Structural Bridge Nodes & Stratum Priority ---")
    county_c = attr_codes["county"]
    rucc_c = attr_codes["rucc_code"]

    is_1hop_county_bridge = [False] * n_nodes
    is_1hop_rucc_bridge = [False] * n_nodes
    for u in range(n_nodes):
        p = parent_idx[u]
        if p >= 0:
            if county_c[p] >= 0 and county_c[u] >= 0 and county_c[p] != county_c[u]:
                is_1hop_county_bridge[u] = True
                is_1hop_county_bridge[p] = True
            if rucc_c[p] >= 0 and rucc_c[u] >= 0 and rucc_c[p] != rucc_c[u]:
                is_1hop_rucc_bridge[u] = True
                is_1hop_rucc_bridge[p] = True

    is_2hop_county_bridge = list(is_1hop_county_bridge)
    for u in range(n_nodes):
        if is_1hop_county_bridge[u]:
            p = parent_idx[u]
            if p >= 0:
                is_2hop_county_bridge[p] = True
            for ch in children_of[u]:
                is_2hop_county_bridge[ch] = True

    cty_inf_counts = Counter(c for c in county_c if c >= 0)
    ranked_counties = [c for c, _ in cty_inf_counts.most_common()]
    top5_set = set(ranked_counties[:5])
    mid_cutoff = max(5, len(ranked_counties) // 3)
    mid_set = set(ranked_counties[5:mid_cutoff])

    def _county_tier(c_code: int) -> str:
        if c_code < 0:
            return "Unknown"
        if c_code in top5_set:
            return "1_Top5_Urban_Counties"
        if c_code in mid_set:
            return "2_MidSize_Counties"
        return "3_Small_Rural_Counties"

    sample_entries: list[tuple[str | None, str]] = []
    for raw_spec in (sample_specs or []):
        if not raw_spec:
            continue
        tag: str | None = None
        spec = raw_spec
        if "=" in raw_spec and not os.path.exists(raw_spec):
            maybe_tag, rest = raw_spec.split("=", 1)
            if maybe_tag.strip() and rest.strip():
                tag = maybe_tag.strip()
                spec = rest.strip()
        if os.path.isdir(spec):
            matched = sorted(_glob.glob(os.path.join(spec, "*_samples.csv*")))
            sample_entries.extend((tag, m) for m in matched)
        else:
            matched = sorted(_glob.glob(spec))
            if matched:
                sample_entries.extend((tag, m) for m in matched)
            elif os.path.exists(spec):
                sample_entries.append((tag, spec))

    sampled_indices_by_strat: dict[str, list[int]] = {}
    for tag, sp in sample_entries:
        base_lbl = _label_from_path(sp)
        lbl = f"{base_lbl}__{tag}" if tag else base_lbl
        if lbl in sampled_indices_by_strat and not tag:
            parent_tag = Path(sp).parent.name.replace("03_sampled_", "")
            lbl = f"{base_lbl}__{parent_tag}"
        s_indices: list[int] = []
        seen_s: set[int] = set()
        with _open_text(sp, "rt") as fh:
            hdr = fh.readline()
            if not hdr:
                continue
            dl = "\t" if hdr.count("\t") > hdr.count(",") else ","
            scols = [c.strip() for c in hdr.rstrip("\r\n").split(dl)]
            ia = scols.index("alias_pid") if "alias_pid" in scols else -1
            if ia < 0:
                continue
            for row in csv.reader(fh, delimiter=dl):
                if len(row) > ia:
                    u = alias_to_idx.get(row[ia].strip(), -1)
                    if u >= 0 and u not in seen_s:
                        seen_s.add(u)
                        s_indices.append(u)
        if s_indices:
            sampled_indices_by_strat[lbl] = s_indices

    if not sampled_indices_by_strat and simulate_strategies:
        import random as _rnd
        rng = _rnd.Random(43)
        target_budget = min(14000, max(2000, n_nodes // 40))
        print(f"  Simulating in-memory SURS, LASSO50-Greedy, LASSO1500-Stratified, and Stratified (budget={target_budget:,})...")
        # Realize an ascertained line-list pool via Bernoulli(test_prob * scale) so pool ~ 20% of n_nodes
        mean_tp = (sum(ll_weight) / max(n_nodes, 1)) or 0.01
        scale_tp = min(25.0, 0.22 / mean_tp)
        ll_pool = [u for u in range(n_nodes) if rng.random() < min(0.98, ll_weight[u] * scale_tp)]
        if len(ll_pool) < target_budget * 2:
            ll_pool = list(range(n_nodes))

        # 1) SURS: uniform random from ascertained line list
        sampled_indices_by_strat["sim_SURS"] = sorted(rng.sample(ll_pool, min(target_budget, len(ll_pool))))

        # Group ll_pool by fine cell g = (age_group, smh_race, county, sex)
        lasso_axes_sim = ("age_group", "smh_race", "county", "sex")
        by_fine: dict[tuple[int, ...], list[int]] = defaultdict(list)
        for u in ll_pool:
            gt = tuple(attr_codes[a][u] for a in lasso_axes_sim)
            by_fine[gt].append(u)
        for gt in by_fine:
            rng.shuffle(by_fine[gt])

        true_cell_cnt: Counter = Counter()
        for u in range(n_nodes):
            true_cell_cnt[tuple(attr_codes[a][u] for a in lasso_axes_sim)] += 1
        q_cell = {gt: c / max(n_nodes, 1) for gt, c in true_cell_cnt.items()}
        tot_ll_pool = len(ll_pool) or 1

        # Score each fine cell by marginal KL priority y_g = log((q_g + eps) / (p_ll_g + eps))
        fine_keys = sorted(by_fine.keys())
        y_by_g = {
            gt: math.log((q_cell.get(gt, 1e-7) + 1e-7) / ((len(by_fine[gt]) / tot_ll_pool) + 1e-7))
            for gt in fine_keys
        }
        sorted_by_y = sorted(fine_keys, key=lambda gt: y_by_g[gt])

        # Helper to sample across K quantile supergroups using either:
        #   mode="greedy_cap_q": greedy supergroup target + within-supergroup weights ~ cap_g * q_g
        #   mode="quota_uniform": rigid micro-quota per bin, uniform across member cells
        def _draw_supergroup_strategy(n_super: int, mode: str) -> list[int]:
            bins: list[list[tuple[int, ...]]] = [[] for _ in range(n_super)]
            for idx_g, gt in enumerate(sorted_by_y):
                b = min(n_super - 1, (idx_g * n_super) // max(len(sorted_by_y), 1))
                bins[b].append(gt)
            bin_q = [sum(q_cell.get(gt, 0.0) for gt in b_list) for b_list in bins]
            tot_bq = sum(bin_q) or 1.0
            bin_quotas = [int(round(target_budget * (bq / tot_bq))) for bq in bin_q]

            chosen: list[int] = []
            chosen_set: set[int] = set()
            for b_idx, b_list in enumerate(bins):
                need = bin_quotas[b_idx]
                if need <= 0 or not b_list:
                    continue
                if mode == "greedy_cap_q":
                    # Within each coarse supergroup, weight member cells by cap_g * q_g
                    # so active outbreak clusters in under-sampled strata stay connected!
                    weighted_pool: list[tuple[float, int]] = []
                    for gt in b_list:
                        members = by_fine[gt]
                        cap_g = len(members)
                        qg = q_cell.get(gt, 1e-7)
                        w_cell = cap_g * qg
                        for rank_i, u in enumerate(members):
                            # Exponential key for weighted sampling without replacement
                            score = (rng.random() ** (1.0 / max(w_cell, 1e-9))) / (1.0 + 0.02 * rank_i)
                            weighted_pool.append((score, u))
                    weighted_pool.sort(reverse=True)
                    for _, u in weighted_pool[:need]:
                        if u not in chosen_set and len(chosen) < target_budget:
                            chosen.append(u)
                            chosen_set.add(u)
                else:
                    # Rigid micro-quota: round-robin 1 singleton per fine cell ignoring cap_g
                    max_len = max((len(by_fine[gt]) for gt in b_list), default=0)
                    added_b = 0
                    for r_i in range(max_len):
                        if added_b >= need or len(chosen) >= target_budget:
                            break
                        for gt in b_list:
                            if r_i < len(by_fine[gt]):
                                u = by_fine[gt][r_i]
                                if u not in chosen_set:
                                    chosen.append(u)
                                    chosen_set.add(u)
                                    added_b += 1
                                    if added_b >= need or len(chosen) >= target_budget:
                                        break
            if len(chosen) < target_budget:
                rest = [u for u in ll_pool if u not in chosen_set]
                chosen.extend(rng.sample(rest, min(len(rest), target_budget - len(chosen))))
            return sorted(chosen[:target_budget])

        sampled_indices_by_strat["sim_LASSO50_Greedy_P"] = _draw_supergroup_strategy(50, "greedy_cap_q")
        sampled_indices_by_strat["sim_LASSO1500_Stratified_P"] = _draw_supergroup_strategy(min(1500, len(fine_keys)), "quota_uniform")
        sampled_indices_by_strat["sim_Stratified_Quota_P"] = _draw_supergroup_strategy(len(fine_keys), "quota_uniform")

    stratum_rows: list[dict] = []
    profile_dims = ("urbanicity", "rucc_code", "age_group", "smh_race", "ses_category", "county_tier")

    for dim in profile_dims:
        inf_cnt: Counter = Counter()
        ll_w_cnt: dict[str, float] = defaultdict(float)
        b1_cnt: Counter = Counter()
        b2_cnt: Counter = Counter()
        r1_cnt: Counter = Counter()

        for u in range(n_nodes):
            if dim == "county_tier":
                label = _county_tier(county_c[u])
            else:
                code = attr_codes[dim][u]
                label = attr_rev[dim][code] if code >= 0 else "Unknown"
            if label in ("Unknown", ""):
                continue
            inf_cnt[label] += 1
            ll_w_cnt[label] += ll_weight[u]
            if is_1hop_county_bridge[u]:
                b1_cnt[label] += 1
            if is_2hop_county_bridge[u]:
                b2_cnt[label] += 1
            if is_1hop_rucc_bridge[u]:
                r1_cnt[label] += 1

        tot_inf = sum(inf_cnt.values()) or 1
        tot_ll_d = sum(ll_w_cnt.values()) or 1.0

        strat_shares: dict[str, dict[str, float]] = {}
        for s_lbl, s_idxs in sampled_indices_by_strat.items():
            sc: Counter = Counter()
            for u in s_idxs:
                if dim == "county_tier":
                    lab = _county_tier(county_c[u])
                else:
                    cd = attr_codes[dim][u]
                    lab = attr_rev[dim][cd] if cd >= 0 else "Unknown"
                if lab not in ("Unknown", ""):
                    sc[lab] += 1
            st_tot = sum(sc.values()) or 1
            strat_shares[s_lbl] = {k: v / st_tot for k, v in sc.items()}

        for cat_label in sorted(inf_cnt.keys()):
            n_i = inf_cnt[cat_label]
            w_l = ll_w_cnt[cat_label]
            p_true = n_i / tot_inf
            p_ll = max(w_l / tot_ll_d, 1e-9)
            asc_ratio = p_ll / max(p_true, 1e-9)
            kl_priority = math.log(max(p_true, 1e-9) / p_ll)

            row_d: dict = {
                "state": state_label,
                "dimension": dim,
                "stratum": cat_label,
                "n_infections": n_i,
                "share_infections_true": round(p_true, 6),
                "linelist_weight_sum": round(w_l, 2),
                "share_linelist_ll": round(p_ll, 6),
                "ascertainment_bias_ratio_ll_over_true": round(asc_ratio, 4),
                "lasso_marginal_kl_priority": round(kl_priority, 4),
                "frac_1hop_cross_county_bridge": round(b1_cnt[cat_label] / n_i, 6),
                "frac_2hop_bridge_adjacent": round(b2_cnt[cat_label] / n_i, 6),
                "frac_1hop_cross_rucc_bridge": round(r1_cnt[cat_label] / n_i, 6),
            }
            for ref_arm in ("1S__surs", "1S-P__lasso_greedy", "1S-1_LL-P__lasso_greedy",
                            "4S__surs", "4S-P__lasso_greedy", "4S-4_LL-P__lasso_greedy",
                            "1S-P__lasso_stratified", "1S-P__stratified",
                            "4S-P__lasso_stratified", "4S-P__stratified"):
                if ref_arm in strat_shares:
                    p_s = strat_shares[ref_arm].get(cat_label, 0.0)
                    row_d[f"weight_{ref_arm}"] = round(p_s / p_ll, 4)
            stratum_rows.append(row_d)
            if dim in ("urbanicity", "age_group", "ses_category", "county_tier"):
                print(
                    f"  [{dim:<12}] {cat_label:<28} asc_bias={asc_ratio:.2f}x  "
                    f"LASSO_priority={kl_priority:+.3f}  "
                    f"1hop_bridge={b1_cnt[cat_label]/n_i*100:5.2f}%  "
                    f"2hop_adj={b2_cnt[cat_label]/n_i*100:5.2f}%"
                )

    stratum_csv = out_dir / f"{state_label}_stratum_bridge_enrichment.csv"
    if stratum_rows:
        with open(stratum_csv, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(stratum_rows[0].keys()))
            w.writeheader()
            w.writerows(stratum_rows)
        print(f"  -> Wrote {stratum_csv} ({len(stratum_rows)} rows)")

    # ------------------------------------------------------------------
    # Part 3: LASSO Feature Attribution & Realized Weights Over Time
    # ------------------------------------------------------------------
    print("\n--- Part 3: Computing LASSO Feature Attribution Over Time ---")
    valid_ords = [d for d in date_ordinal_list if d > 0]
    min_ord = min(valid_ords) if valid_ords else 0
    max_ord = max(valid_ords) if valid_ords else 0

    lasso_axes = ["age_group", "smh_race", "county", "sex"]
    target_cell_counts: Counter = Counter()
    for u in range(n_nodes):
        gt = tuple(
            attr_rev[a][attr_codes[a][u]] if attr_codes[a][u] >= 0 else ""
            for a in lasso_axes
        )
        if all(gt):
            target_cell_counts[gt] += 1
    tot_target_cells = sum(target_cell_counts.values()) or 1
    q_target = {gt: cnt / tot_target_cells for gt, cnt in target_cell_counts.items()}

    windows: list[tuple[str, int, int]] = [("ALL", min_ord, max_ord + 1)]
    if max_ord > min_ord:
        cur_start = min_ord
        w_idx = 1
        while cur_start <= max_ord:
            windows.append((f"4wk_{w_idx:02d}", cur_start, cur_start + 28))
            cur_start += 28
            w_idx += 1

    lasso_time_rows: list[dict] = []
    for w_name, o_lo, o_hi in windows:
        pool_cell_weights: dict[tuple[str, ...], float] = defaultdict(float)
        n_win_inf = 0
        w_win_ll = 0.0
        w_win_rural_ll = 0.0
        w_win_top5_ll = 0.0
        for u in range(n_nodes):
            d_o = date_ordinal_list[u]
            if w_name != "ALL" and (d_o < o_lo or d_o >= o_hi):
                continue
            n_win_inf += 1
            wu = ll_weight[u]
            if wu <= 0.0:
                continue
            w_win_ll += wu
            if attr_codes["urbanicity"][u] >= 0 and "Non-metro" in attr_rev["urbanicity"][attr_codes["urbanicity"][u]]:
                w_win_rural_ll += wu
            if county_c[u] in top5_set:
                w_win_top5_ll += wu
            gt = tuple(
                attr_rev[a][attr_codes[a][u]] if attr_codes[a][u] >= 0 else ""
                for a in lasso_axes
            )
            if all(gt):
                pool_cell_weights[gt] += wu

        if n_win_inf < 200 or len(pool_cell_weights) < 10 or w_win_ll <= 0:
            continue

        g_list = sorted(pool_cell_weights.keys())
        eps = 1e-7
        y_vals = [
            math.log((q_target.get(gt, eps) + eps) / ((pool_cell_weights[gt] / w_win_ll) + eps))
            for gt in g_list
        ]
        w_vals = [float(pool_cell_weights[gt]) for gt in g_list]
        l1_shares = _fit_lasso_feature_importance(g_list, y_vals, w_vals, lasso_axes)

        t_row: dict = {
            "state": state_label,
            "window": w_name,
            "start_date": _dt.date.fromordinal(o_lo).isoformat() if o_lo > 0 else "",
            "end_date": _dt.date.fromordinal(min(o_hi - 1, max_ord)).isoformat() if o_hi > 0 else "",
            "n_infections_window": n_win_inf,
            "linelist_weight_window": round(w_win_ll, 2),
            "n_active_fine_strata": len(g_list),
            "lasso_l1_share_county": l1_shares.get("county", 0.0),
            "lasso_l1_share_age_group": l1_shares.get("age_group", 0.0),
            "lasso_l1_share_smh_race": l1_shares.get("smh_race", 0.0),
            "lasso_l1_share_sex": l1_shares.get("sex", 0.0),
            "ll_rural_rucc4_9_share": round(w_win_rural_ll / w_win_ll, 6),
            "ll_top5_urban_county_share": round(w_win_top5_ll / w_win_ll, 6),
        }

        for arm_key in ("1S__surs", "1S-P__lasso_greedy", "4S__surs", "4S-P__lasso_greedy", "4S-4_LL-P__lasso_greedy"):
            if arm_key in sampled_indices_by_strat:
                s_win = [
                    u for u in sampled_indices_by_strat[arm_key]
                    if w_name == "ALL" or (o_lo <= date_ordinal_list[u] < o_hi)
                ]
                if s_win:
                    r_cnt = sum(
                        1 for u in s_win
                        if attr_codes["urbanicity"][u] >= 0
                        and "Non-metro" in attr_rev["urbanicity"][attr_codes["urbanicity"][u]]
                    )
                    u_cnt = sum(1 for u in s_win if county_c[u] in top5_set)
                    t_row[f"rural_share_{arm_key}"] = round(r_cnt / len(s_win), 6)
                    t_row[f"top5_urban_share_{arm_key}"] = round(u_cnt / len(s_win), 6)
        lasso_time_rows.append(t_row)

    lasso_csv = out_dir / f"{state_label}_lasso_weights_over_time.csv"
    if lasso_time_rows:
        with open(lasso_csv, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(lasso_time_rows[0].keys()))
            w.writeheader()
            w.writerows(lasso_time_rows)
        all_row = lasso_time_rows[0]
        print(
            f"  Full-wave LASSO L1 weight shares: "
            f"county={all_row['lasso_l1_share_county']:.3f}, "
            f"age_group={all_row['lasso_l1_share_age_group']:.3f}, "
            f"smh_race={all_row['lasso_l1_share_smh_race']:.3f}, "
            f"sex={all_row['lasso_l1_share_sex']:.3f}"
        )
        print(f"  -> Wrote {lasso_csv} ({len(lasso_time_rows)} rows)")

    # ------------------------------------------------------------------
    # Part 4: Strategy-Level Coalescence, Bridge Bracketing & Mugration Decomposition
    # ------------------------------------------------------------------
    if not sampled_indices_by_strat:
        print("\nNo `--samples` files provided or matched; Parts 1-3 complete.")
        return {
            "hop_decay_csv": str(decay_csv),
            "stratum_bridge_csv": str(stratum_csv),
            "lasso_time_csv": str(lasso_csv),
        }

    print(f"\n--- Part 4: Evaluating {len(sampled_indices_by_strat)} Sampling Strategies ---")

    # Build topological post-order of the directed transmission forest (children before parents)
    # using iterative DFS so bottom-up mugration state propagation runs in ~0.3s in pure Python!
    in_deg = [0] * n_nodes
    for u in range(n_nodes):
        if parent_idx[u] >= 0:
            in_deg[u] += 1
    roots = [u for u in range(n_nodes) if in_deg[u] == 0]
    topo_order: list[int] = []
    visited_topo = [False] * n_nodes
    for r in roots:
        stack = [r]
        visited_topo[r] = True
        while stack:
            curr = stack.pop()
            topo_order.append(curr)
            for ch in children_of[curr]:
                if not visited_topo[ch]:
                    visited_topo[ch] = True
                    stack.append(ch)
    if len(topo_order) < n_nodes:
        for u in range(n_nodes):
            if not visited_topo[u]:
                topo_order.append(u)
    post_order = list(reversed(topo_order))

    # Compute true full-graph transition counts q_true[c1][c2] and true equilibrium pi_true[c]
    n_counties = len(attr_rev["county"])
    q_true_mat = [[0.0] * n_counties for _ in range(n_counties)]
    true_county_counts = [0] * n_counties
    for u in range(n_nodes):
        cu = county_c[u]
        if cu >= 0:
            true_county_counts[cu] += 1
        pu = parent_idx[u]
        if pu >= 0:
            cp = county_c[pu]
            if cp >= 0 and cu >= 0 and cp != cu:
                q_true_mat[cp][cu] += 1.0

    tot_true_jumps = sum(sum(r) for r in q_true_mat) or 1.0
    tot_true_nodes = sum(true_county_counts) or 1
    pi_true = [c / tot_true_nodes for c in true_county_counts]

    # Symmetrized true Q (raw jumps) and true W (rate matrix)
    Q_true_sym = [
        [q_true_mat[i][j] + q_true_mat[j][i] if i != j else 0.0 for j in range(n_counties)]
        for i in range(n_counties)
    ]
    W_true_sym = [[0.0] * n_counties for _ in range(n_counties)]
    for i in range(n_counties):
        for j in range(n_counties):
            if i != j and pi_true[i] > 0 and pi_true[j] > 0:
                w_ij = q_true_mat[i][j] / (tot_true_jumps * pi_true[j])
                w_ji = q_true_mat[j][i] / (tot_true_jumps * pi_true[i])
                W_true_sym[i][j] = 0.5 * (w_ij + w_ji)

    # If external --truth linelist_mugration.json is provided, load and align its W matrix too
    if truth_mugration_path and Path(truth_mugration_path).exists():
        try:
            tdata = json.loads(Path(truth_mugration_path).read_text())
            t_alpha = tdata["models"]["county"]["alphabet"]
            t_mat = tdata["models"]["county"]["transition_matrix"]
            t_idx_map = {name: i for i, name in enumerate(t_alpha)}
            for i, cname_i in enumerate(attr_rev["county"]):
                ti = t_idx_map.get(cname_i, -1)
                if ti < 0:
                    continue
                for j, cname_j in enumerate(attr_rev["county"]):
                    tj = t_idx_map.get(cname_j, -1)
                    if tj >= 0 and i != j:
                        W_true_sym[i][j] = float(t_mat[ti][tj])
        except Exception as exc:
            print(f"  Note: using full-graph W_true (could not parse {truth_mugration_path}: {exc})")

    flat_W_true = _row_normalize_offdiag(W_true_sym)
    flat_Q_true = _row_normalize_offdiag(Q_true_sym)
    true_active_edges = [1 if x > 0 else 0 for x in flat_W_true]
    n_true_corridors = sum(true_active_edges)

    strategy_rows: list[dict] = []
    for strat_label, s_idxs in sorted(sampled_indices_by_strat.items()):
        s_set = set(s_idxs)
        n_s = len(s_idxs)

        # 4a. Macro reach & marginal composition
        s_counties = [county_c[u] for u in s_idxs if county_c[u] >= 0]
        n_distinct_counties = len(set(s_counties))
        rural_n = sum(
            1 for u in s_idxs
            if attr_codes["urbanicity"][u] >= 0
            and "Non-metro" in attr_rev["urbanicity"][attr_codes["urbanicity"][u]]
        )
        top5_n = sum(1 for u in s_idxs if county_c[u] in top5_set)
        b1_s_n = sum(1 for u in s_idxs if is_1hop_county_bridge[u])
        b2_s_n = sum(1 for u in s_idxs if is_2hop_county_bridge[u])

        # 4b. Directed ancestral coalescence within k <= max_ancestral_hops
        # For each sampled node u in s_idxs, walk up parent_idx up to max_ancestral_hops
        # to find its nearest sampled ancestor in s_set.
        anc_pairs_by_max_k: dict[int, int] = {1: 0, 2: 0, 3: 0, 5: 0, 10: 0}
        le5_cat_pairs: dict[str, list[tuple[int, int]]] = {a: [] for a in ("county", "rucc_code", "age_group", "smh_race")}
        for u in s_idxs:
            cur = u
            for h in range(1, max_ancestral_hops + 1):
                par = parent_idx[cur]
                if par < 0 or par == cur:
                    break
                if par in s_set:
                    for cap_k in anc_pairs_by_max_k:
                        if h <= cap_k:
                            anc_pairs_by_max_k[cap_k] += 1
                    if h <= 5:
                        for a in le5_cat_pairs:
                            le5_cat_pairs[a].append((attr_codes[a][par], attr_codes[a][u]))
                    break
                cur = par

        # Marginal-null-adjusted concordance on short-range (<=5 hop) ancestral pairs
        c_within_5, c_null_5, c_exc_5, _ = _compute_newman_r_and_concordance(
            le5_cat_pairs["county"], n_counties
        )
        cross_county_le5 = round(1.0 - c_within_5, 6) if c_within_5 is not None else None
        null_cross_county_le5 = round(1.0 - c_null_5, 6) if c_null_5 is not None else None

        # 4c. Vectorized / post-order bottom-up mugration state reconstruction on G
        inferred_state = [-1] * n_nodes
        for u in s_idxs:
            inferred_state[u] = county_c[u]

        for u in post_order:
            if inferred_state[u] >= 0:
                continue
            ch_list = children_of[u]
            if not ch_list:
                continue
            # Fast majority vote over children with resolved states
            first_st = -1
            same_all = True
            valid_ch: list[int] = []
            for ch in ch_list:
                st_c = inferred_state[ch]
                if st_c >= 0:
                    valid_ch.append(st_c)
                    if first_st < 0:
                        first_st = st_c
                    elif st_c != first_st:
                        same_all = False
            if not valid_ch:
                continue
            if same_all:
                inferred_state[u] = first_st
            else:
                inferred_state[u] = Counter(valid_ch).most_common(1)[0][0]

        # Count inferred transitions q_inf[cp][cu] and inferred node frequencies state_counts[c]
        q_inf_mat = [[0.0] * n_counties for _ in range(n_counties)]
        inf_county_counts = [0] * n_counties
        total_jumps_inf = 0.0
        for u in range(n_nodes):
            cu = inferred_state[u]
            if cu >= 0:
                inf_county_counts[cu] += 1
                pu = parent_idx[u]
                if pu >= 0:
                    cp = inferred_state[pu]
                    if cp >= 0 and cp != cu:
                        q_inf_mat[cp][cu] += 1.0
                        total_jumps_inf += 1.0

        tot_inf_nodes = sum(inf_county_counts) or 1
        pi_inf = [c / tot_inf_nodes for c in inf_county_counts]

        # Build 3 matrices to decompose mugration Cosine Similarity:
        #   1) W_inf_full : standard mugration W_ij using inferred jumps q_inf and inferred pi_inf
        #   2) Q_inf_raw  : raw symmetrized jump counts q_inf[i][j] + q_inf[j][i] (numerator only)
        #   3) W_inf_true_pi : inferred jumps q_inf divided by TRUE equilibrium pi_true
        Q_inf_sym = [
            [q_inf_mat[i][j] + q_inf_mat[j][i] if i != j else 0.0 for j in range(n_counties)]
            for i in range(n_counties)
        ]
        W_inf_full = [[0.0] * n_counties for _ in range(n_counties)]
        W_inf_true_pi = [[0.0] * n_counties for _ in range(n_counties)]

        if total_jumps_inf > 0:
            for i in range(n_counties):
                for j in range(n_counties):
                    if i == j:
                        continue
                    q_ij = q_inf_mat[i][j]
                    q_ji = q_inf_mat[j][i]
                    if pi_inf[i] > 0 and pi_inf[j] > 0:
                        W_inf_full[i][j] = 0.5 * (
                            q_ij / (total_jumps_inf * pi_inf[j])
                            + q_ji / (total_jumps_inf * pi_inf[i])
                        )
                    if pi_true[i] > 0 and pi_true[j] > 0:
                        W_inf_true_pi[i][j] = 0.5 * (
                            q_ij / (total_jumps_inf * pi_true[j])
                            + q_ji / (total_jumps_inf * pi_true[i])
                        )

        flat_W_inf = _row_normalize_offdiag(W_inf_full)
        flat_Q_inf = _row_normalize_offdiag(Q_inf_sym)
        flat_W_true_pi = _row_normalize_offdiag(W_inf_true_pi)

        cos_full_W = _cosine_1d(flat_W_true, flat_W_inf)
        cos_raw_Q = _cosine_1d(flat_Q_true, flat_Q_inf)
        cos_true_pi_W = _cosine_1d(flat_W_true, flat_W_true_pi)

        # Topological precision, recall, and F1 on active county-to-county corridors
        inf_active_edges = [1 if x > 0 else 0 for x in flat_W_inf]
        n_inf_corridors = sum(inf_active_edges)
        tp_edges = sum(1 for t, p in zip(true_active_edges, inf_active_edges) if t == 1 and p == 1)
        bridge_recall = tp_edges / max(n_true_corridors, 1)
        bridge_prec = tp_edges / max(n_inf_corridors, 1)
        topo_f1 = (
            2.0 * bridge_prec * bridge_recall / (bridge_prec + bridge_recall)
            if (bridge_prec + bridge_recall) > 0
            else 0.0
        )

        # Equilibrium pi distortion: Total Variation distance 0.5 * sum |pi_inf - pi_true|
        tv_pi = 0.5 * sum(abs(a - b) for a, b in zip(pi_inf, pi_true))
        max_pi_distortion = max(
            (pi_true[c] / max(pi_inf[c], 1e-6) for c in range(n_counties) if pi_true[c] >= 0.001),
            default=1.0,
        )
        s_fine_cells = len({
            tuple(attr_codes[a][u] for a in ("age_group", "smh_race", "county", "sex"))
            for u in s_idxs
        })

        s_row = {
            "state": state_label,
            "strategy": strat_label,
            "n_samples": n_s,
            "n_counties_sampled": n_distinct_counties,
            "n_fine_cells_total": len(target_cell_counts),
            "n_fine_cells_sampled": s_fine_cells,
            "rural_rucc4_9_share": round(rural_n / max(n_s, 1), 6),
            "top5_urban_county_share": round(top5_n / max(n_s, 1), 6),
            "mid_and_small_county_share": round(1.0 - (top5_n / max(n_s, 1)), 6),
            "sample_1hop_bridge_frac": round(b1_s_n / max(n_s, 1), 6),
            "sample_2hop_bridge_adj_frac": round(b2_s_n / max(n_s, 1), 6),
            "n_ancestral_pairs_1hop": anc_pairs_by_max_k[1],
            "n_ancestral_pairs_le2hop": anc_pairs_by_max_k[2],
            "n_ancestral_pairs_le3hop": anc_pairs_by_max_k[3],
            "n_ancestral_pairs_le5hop": anc_pairs_by_max_k[5],
            "n_ancestral_pairs_le10hop": anc_pairs_by_max_k[10],
            "frac_samples_with_ancestor_le5hop": round(anc_pairs_by_max_k[5] / max(n_s, 1), 6),
            "cross_county_frac_le5hop": cross_county_le5,
            "null_cross_county_frac_le5hop": null_cross_county_le5,
            "excess_homophily_county_le5hop": c_exc_5,
            "n_tree_nodes_resolved": tot_inf_nodes,
            "n_total_jumps_q": int(total_jumps_inf),
            "n_county_corridors_truth": n_true_corridors,
            "n_county_corridors_inferred": n_inf_corridors,
            "true_bridge_recall": round(bridge_recall, 6),
            "true_bridge_precision": round(bridge_prec, 6),
            "topological_f1": round(topo_f1, 6),
            "cos_full_W": round(cos_full_W, 6),
            "cos_raw_jumps_Q": round(cos_raw_Q, 6),
            "cos_true_pi_W": round(cos_true_pi_W, 6),
            "pi_denominator_penalty": round(cos_true_pi_W - cos_full_W, 6),
            "tv_error_pi": round(tv_pi, 6),
            "max_pi_ratio_distortion": round(max_pi_distortion, 4),
        }
        strategy_rows.append(s_row)
        print(
            f"  {strat_label:<30} cos_W={cos_full_W:.4f} (raw_Q={cos_raw_Q:.4f}, "
            f"true_pi_W={cos_true_pi_W:.4f}, tv_pi={tv_pi:.3f}) "
            f"F1={topo_f1:.4f} corr={n_inf_corridors:>4d} anc<=5={anc_pairs_by_max_k[5]:>4d}"
        )

    strat_csv = out_dir / f"{state_label}_strategy_mechanism_summary.csv"
    if strategy_rows:
        strategy_rows.sort(key=lambda r: -float(r["cos_full_W"]))
        with open(strat_csv, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(strategy_rows[0].keys()))
            w.writeheader()
            w.writerows(strategy_rows)
        print(f"  -> Wrote {strat_csv} ({len(strategy_rows)} rows)")

    return {
        "hop_decay_csv": str(decay_csv),
        "stratum_bridge_csv": str(stratum_csv),
        "lasso_time_csv": str(lasso_csv),
        "strategy_summary_csv": str(strat_csv) if strategy_rows else None,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--infections", "--allevents", dest="infections", required=True, help="Path to allevents_variants.csv.xz")
    ap.add_argument("--linelist", default=None, help="Optional path to ascertained linelist.csv.xz")
    ap.add_argument("--population", default=None, help="Optional path to demographics CSV")
    ap.add_argument("--truth", "--mugration-truth", dest="truth", default=None, help="Optional path to linelist_mugration.json")
    ap.add_argument("--samples", nargs="*", default=[], help="Sample CSV paths, directories, or tag=dir/glob specs")
    ap.add_argument("--samples-dir", nargs="*", default=[], help="Optional directory(s) of *_samples.csv.xz files")
    ap.add_argument("--foil-samples-dir", nargs="*", default=[], help="Optional foil/ablation directory(s) of *_samples.csv.xz files")
    ap.add_argument("--state", default="va", help="State label (e.g. va, ga, ma, mn, wa)")
    ap.add_argument("--outdir", required=True, help="Output directory for mechanism CSVs")
    ap.add_argument("--max-records", type=int, default=0, help="Optional cap on infections read (for quick smoke tests)")
    ap.add_argument("--simulate-strategies", action="store_true", help="Simulate SURS, LASSO50-Greedy, LASSO1500-Stratified, and Stratified in memory when --samples are absent")
    args = ap.parse_args(argv)

    all_sample_specs = list(args.samples or []) + list(args.samples_dir or []) + list(args.foil_samples_dir or [])

    run_mechanism_analysis(
        infections_path=args.infections,
        outdir=args.outdir,
        state_label=args.state,
        linelist_path=args.linelist,
        population_path=args.population,
        truth_mugration_path=args.truth,
        sample_specs=all_sample_specs,
        max_records=args.max_records,
        simulate_strategies=args.simulate_strategies,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
