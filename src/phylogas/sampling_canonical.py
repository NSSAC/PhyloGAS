"""Canonical public-health specimen-submission pools and line-list-only subsamplers.

Implements two complementary layers for benchmarking real-world public-health
genomic surveillance designs (e.g., Oltean et al. 2023 Washington State DOH
sentinel evaluation; Brito et al. 2022 ``Subsampler``; Nextstrain ``augur filter
--group-by``) **without requiring an external Census population model** and
using only the Python standard library so it runs in both bare Python and
full Conda HPC environments:

1. **Specimen-Submission Pool Regimes (Type-A Laboratory Catchment Filters)**
   Separates the *Epidemiologic Reported Case Line List* (all ascertained cases
   in ``linelist.csv.xz``, known as county case counts to state epidemiologists)
   from the *Sequencable Swab Pool* physically routed to sequencing laboratories:

   * ``ship_convenience`` (*Pre-Sentinel Convenience Status Quo*):
     Concentrated in ``K = 2`` flagship metropolitan academic hubs with steep
     distance decay, enriched ``2.2x`` for hospitalized/critical inpatients
     (``hM*``/``dM*``) and ``2.0x`` for large outbreak components, with higher
     ``Ct > 30`` dropout on asymptomatic/presymptomatic cases (``0.40x``).
   * ``ship_sentinel`` (*Regional Clinical Sentinel Network*, Oltean et al. 2023):
     Selects ``K = 6`` geographically dispersed regional clinical lab hubs via
     population-weighted minimax ``k``-center placement (``cases^alpha * min_dist``),
     models spatial catchment decay plus a ``2.5%`` statewide commercial-lab
     baseline (CDC NS3), and applies state-dependent ``Ct < 30`` / hospital-lab
     routing pass multipliers (outpatient ``I* = 1.00x``, asymptomatic
     ``A*/P* = 0.50x``, hospitalized inpatient ``hM*/dM* = 0.45x`` -> ``O/E ~ 0.57``).

2. **Line-List-Only Allocation Strategies (Fixed Weekly / 4-Weekly Budget B)**
   * ``convenience``        : Uniform random draw (SURS) from ``ship_convenience == 1``.
   * ``sentinel``           : Uniform random draw (SURS) from ``ship_sentinel == 1``.
   * ``sentinel_caseprop``  : Case-proportional allocation (Brito ``Subsampler``)
                              drawing swabs from ``ship_sentinel == 1`` to match
                              **full statewide reported county case shares** in
                              ``linelist``, with iterative water-filling when
                              under-covered counties exhaust shipped swabs.
   * ``sentinel_geocap``    : Nextstrain ``augur filter --group-by county`` equal
                              county quota on ``ship_sentinel == 1`` with
                              iterative water-filling redistribution.
   * ``statewide_geocap``   : Equal county quota on ``100%`` mandatory statewide
                              line list (isolates pure equi-county sampling).
   * ``statewide_caseprop`` : Case-proportional county quota on ``100%``
                              mandatory statewide line list.
"""

from __future__ import annotations

import csv
import datetime as _dt
import gzip
import json
import lzma
import math
import random
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

EARTH_RADIUS_KM = 6371.0088

CANONICAL_STRATEGIES = (
    "convenience",
    "sentinel",
    "sentinel_caseprop",
    "sentinel_geocap",
    "statewide_surs",
    "statewide_caseprop",
    "statewide_geocap",
)


@dataclass(frozen=True)
class PoolConfig:
    """Parameters governing Type-A laboratory specimen-submission regimes."""

    # Sentinel network (post-March 2021 WA DOH style: ~12-18% statewide pool)
    sentinel_k_hubs: int = 6
    sentinel_hub_pop_alpha: float = 0.25
    sentinel_min_hub_share: float = 0.008
    sentinel_p_hub: float = 0.42
    sentinel_decay_km: float = 45.0
    sentinel_p_floor: float = 0.020
    sentinel_mult_outpatient: float = 1.00
    sentinel_mult_asymp: float = 0.50
    sentinel_mult_inpatient: float = 0.45

    # Pre-sentinel convenience (academic metro hubs + inpatient/outbreak skew)
    convenience_k_hubs: int = 2
    convenience_p_hub: float = 0.22
    convenience_decay_km: float = 28.0
    convenience_p_floor: float = 0.003
    convenience_mult_outpatient: float = 1.00
    convenience_mult_asymp: float = 0.40
    convenience_mult_inpatient: float = 2.40
    convenience_mult_outbreak: float = 1.60
    convenience_outbreak_min_cases: int = 50


def open_text_auto(path: Path | str, mode: str = "rt"):
    p = str(path)
    if p.endswith(".xz"):
        return lzma.open(p, mode, encoding="utf-8", newline="")
    if p.endswith(".gz"):
        return gzip.open(p, mode, encoding="utf-8", newline="")
    return open(p, mode, encoding="utf-8", newline="")


def haversine_km_scalar(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Great-circle distance in kilometers between two WGS84 points."""
    r_lat1 = math.radians(lat1)
    r_lon1 = math.radians(lon1)
    r_lat2 = math.radians(lat2)
    r_lon2 = math.radians(lon2)
    dlat = r_lat2 - r_lat1
    dlon = r_lon2 - r_lon1
    a = math.sin(dlat / 2.0) ** 2 + math.cos(r_lat1) * math.cos(r_lat2) * math.sin(dlon / 2.0) ** 2
    return 2.0 * EARTH_RADIUS_KM * math.asin(min(1.0, math.sqrt(max(0.0, a))))


def _parse_iso_ordinal(ds: str) -> int | None:
    s = ds.strip()[:10]
    if len(s) == 10 and s[4] == "-" and s[7] == "-":
        try:
            return _dt.date.fromisoformat(s).toordinal()
        except ValueError:
            return None
    return None


def select_regional_hubs(
    county_cases: dict[str, int],
    county_coords: dict[str, tuple[float, float]],
    k_hubs: int = 6,
    pop_alpha: float = 0.35,
    min_hub_share: float = 0.005,
) -> list[str]:
    """Deterministically select ``k_hubs`` geographically dispersed regional lab hubs."""
    if not county_cases or k_hubs <= 0:
        return []
    ordered = sorted(county_cases.items(), key=lambda kv: (-kv[1], kv[0]))
    k = min(k_hubs, len(ordered))

    if len(county_coords) < len(ordered):
        return [c for c, _ in ordered[:k]]

    tot = float(sum(county_cases.values()))
    min_req = tot * min_hub_share

    elig = [c for c, cnt in ordered if cnt >= min_req and c in county_coords]
    if len(elig) < k:
        elig = [c for c, _ in ordered if c in county_coords]
    if not elig:
        return [c for c, _ in ordered[:k]]

    chosen = [elig[0]]
    min_d = {c: float("inf") for c in elig}
    while len(chosen) < k:
        last = chosen[-1]
        la0, lo0 = county_coords[last]
        for c in elig:
            if c in chosen:
                min_d[c] = -1.0
            else:
                la1, lo1 = county_coords[c]
                d = haversine_km_scalar(la0, lo0, la1, lo1)
                if d < min_d[c]:
                    min_d[c] = d

        best_c = None
        best_score = -1.0
        for c in elig:
            if min_d[c] <= 0:
                continue
            sc = (float(county_cases[c]) ** pop_alpha) * (min_d[c] ** 1.20)
            if sc > best_score:
                best_score = sc
                best_c = c
        if best_c is None:
            break
        chosen.append(best_c)
    return chosen


def allocate_waterfill_quotas(
    capacities: dict[str, int],
    weights: dict[str, float] | None,
    budget: int,
    rng: random.Random,
) -> dict[str, int]:
    """Allocate integer quotas summing to ``min(budget, sum(capacities))`` subject
    to ``0 <= quota[g] <= capacities[g]``.

    * ``weights is None`` -> **GeoCap** (equal quota per active group with water-filling).
    * ``weights`` given   -> **CaseProp** (proportional to ``weights[g]``, iteratively
                             redistributing excess quota whenever an under-shipped
                             county hits its available swab cap).
    """
    avail = {g: int(c) for g, c in capacities.items() if int(c) > 0}
    if not avail or budget <= 0:
        return {g: 0 for g in capacities}

    tot_cap = sum(avail.values())
    if tot_cap <= budget:
        return {g: int(capacities.get(g, 0)) for g in capacities}

    if weights is None:
        w = {g: 1.0 for g in avail}
    else:
        w = {g: max(float(weights.get(g, 0.0)), 0.0) for g in avail}
        if sum(w.values()) <= 0.0:
            w = {g: 1.0 for g in avail}

    locked: dict[str, int] = {}
    rem_budget = float(budget)
    active = set(avail.keys())

    while active:
        w_sum = sum(w[g] for g in active)
        if w_sum <= 0.0:
            for g in active:
                w[g] = 1.0
            w_sum = float(len(active))

        over_capped = [g for g in sorted(active) if (rem_budget * (w[g] / w_sum)) >= avail[g]]
        if not over_capped:
            break
        for g in over_capped:
            locked[g] = avail[g]
            rem_budget -= avail[g]
            active.remove(g)
        if rem_budget <= 1e-9:
            break

    quotas = {g: 0 for g in capacities}
    quotas.update(locked)

    int_rem = budget - sum(locked.values())
    if int_rem > 0 and active:
        active_sorted = sorted(active)
        w_sum = sum(w[g] for g in active_sorted)
        ideals = {g: int_rem * (w[g] / w_sum) for g in active_sorted}
        floors = {g: min(avail[g], int(math.floor(ideals[g]))) for g in active_sorted}
        for g, fl in floors.items():
            quotas[g] = fl
        leftover = int_rem - sum(floors.values())
        if leftover > 0:
            jitter = {g: rng.random() * 1e-9 for g in active_sorted}
            rem_order = sorted(
                [g for g in active_sorted if quotas[g] < avail[g]],
                key=lambda g: (ideals[g] - floors[g]) + jitter[g],
                reverse=True,
            )
            for g in rem_order[:leftover]:
                quotas[g] += 1
    return quotas


def _draw_window_indices(
    win_indices: list[int],
    row_county: list[str],
    ship_sent: list[int],
    ship_conv: list[int],
    strategy: str,
    budget: int,
    rng: random.Random,
) -> list[int]:
    """Select up to ``budget`` row indices without replacement for one stride window."""
    if not win_indices or budget <= 0:
        return []

    if strategy == "convenience":
        pool = [i for i in win_indices if ship_conv[i] == 1]
        return sorted(pool) if len(pool) <= budget else sorted(rng.sample(pool, budget))

    if strategy == "sentinel":
        pool = [i for i in win_indices if ship_sent[i] == 1]
        return sorted(pool) if len(pool) <= budget else sorted(rng.sample(pool, budget))

    if strategy == "statewide_surs":
        pool = list(win_indices)
        return sorted(pool) if len(pool) <= budget else sorted(rng.sample(pool, budget))

    if strategy in ("sentinel_caseprop", "sentinel_geocap"):
        pool = [i for i in win_indices if ship_sent[i] == 1]
    elif strategy in ("statewide_geocap", "statewide_caseprop"):
        pool = list(win_indices)
    else:
        raise ValueError(f"Unknown canonical strategy: {strategy}")

    if len(pool) <= budget:
        return sorted(pool)

    by_cty: dict[str, list[int]] = defaultdict(list)
    for i in pool:
        by_cty[row_county[i]].append(i)
    caps = {c: len(idxs) for c, idxs in by_cty.items()}

    if strategy in ("sentinel_caseprop", "statewide_caseprop"):
        # Epidemiologist's target weights come from ALL reported cases in `win_indices`
        all_counts = Counter(row_county[i] for i in win_indices)
        weights: dict[str, float] | None = {c: float(cnt) for c, cnt in all_counts.items()}
    else:
        weights = None

    quotas = allocate_waterfill_quotas(caps, weights, budget, rng)
    picked: list[int] = []
    for c in sorted(quotas):
        q = quotas[c]
        if q <= 0:
            continue
        c_idxs = by_cty[c]
        if len(c_idxs) <= q:
            picked.extend(c_idxs)
        else:
            picked.extend(rng.sample(c_idxs, q))
    return sorted(picked)


def run_canonical_from_csv(
    linelist_path: Path | str,
    outdir: Path | str,
    start_date_iso: str,
    min_pool: int = 50,
    weekly_budget: int = 400,
    strides: Iterable[int] = (1, 4),
    strategies: Iterable[str] = CANONICAL_STRATEGIES,
    county_col: str = "county",
    date_field: str = "date",
    pool_cfg: PoolConfig = PoolConfig(),
    seed: int = 43,
    diagnostics_dir: Path | str | None = None,
    state_label: str = "state",
    synthesize_ascertainment_from_test_prob: bool = False,
) -> tuple[list[Path], dict]:
    """Stream ``linelist_path`` once, assign ``ship_sentinel`` & ``ship_convenience``
    pools, execute all canonical subsamplers across strides, write compressed
    ``*_samples.csv.xz`` outputs, and emit catchment/O-E diagnostics.

    If ``synthesize_ascertainment_from_test_prob`` is True (used when testing directly
    against ``allevents_variants.csv.xz`` on machines without ``02_simulated_linelists``),
    rows are first filtered via reproducible ``Bernoulli(test_prob)`` draws.
    """
    out_path = Path(outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    start_ord = _parse_iso_ordinal(start_date_iso)
    if start_ord is None:
        raise ValueError(f"Invalid ISO start_date: {start_date_iso}")
    week0_ord = start_ord - 7

    rows_raw: list[list[str]] = []
    row_county: list[str] = []
    row_wk: list[int] = []
    row_rucc_ge4: list[bool] = []
    row_comp: list[str] = []
    row_exit_prefix: list[str] = []
    row_age: list[str] = []

    cty_lat_sum: dict[str, float] = defaultdict(float)
    cty_lon_sum: dict[str, float] = defaultdict(float)
    cty_coord_n: dict[str, int] = defaultdict(int)
    seen_alias: set[str] = set()
    rng_asc = random.Random(seed + 17)

    with open_text_auto(linelist_path, "rt") as fh:
        first_line = fh.readline()
        if not first_line:
            raise ValueError(f"Empty linelist file: {linelist_path}")
        delim = "\t" if first_line.count("\t") > first_line.count(",") else ","
        header = [c.strip() for c in next(csv.reader([first_line], delimiter=delim))]
        c_idx = {c: i for i, c in enumerate(header)}

        i_alias = c_idx.get("alias_pid")
        i_cty = c_idx.get(county_col)
        i_date = c_idx.get(date_field)
        if i_cty is None or i_date is None:
            raise ValueError(f"Missing required '{county_col}' or '{date_field}' in {linelist_path}")

        i_lat, i_lon = None, None
        for la, lo in (("latitude", "longitude"), ("home_latitude", "home_longitude")):
            if la in c_idx and lo in c_idx:
                i_lat, i_lon = c_idx[la], c_idx[lo]
                break

        i_rucc = c_idx.get("rucc_code")
        i_comp = c_idx.get("component_id")
        i_exit = c_idx.get("exit_state")
        i_age = c_idx.get("age_group", c_idx.get("age"))
        i_tprob = c_idx.get("test_prob")

        for row in csv.reader(fh, delimiter=delim):
            if len(row) <= max(i_cty, i_date):
                continue
            if i_alias is not None and len(row) > i_alias:
                al = row[i_alias].strip()
                if al in seen_alias:
                    continue
                seen_alias.add(al)

            if synthesize_ascertainment_from_test_prob and i_tprob is not None and len(row) > i_tprob:
                try:
                    tp = float(row[i_tprob])
                except ValueError:
                    tp = 0.20
                if rng_asc.random() >= tp:
                    continue

            d_ord = _parse_iso_ordinal(row[i_date])
            if d_ord is None:
                continue
            wk = (d_ord - week0_ord) // 7
            cty = row[i_cty].strip()
            if not cty:
                continue

            rows_raw.append(row)
            row_county.append(cty)
            row_wk.append(wk)

            if i_lat is not None and i_lon is not None and len(row) > max(i_lat, i_lon):
                try:
                    la_v = float(row[i_lat])
                    lo_v = float(row[i_lon])
                    cty_lat_sum[cty] += la_v
                    cty_lon_sum[cty] += lo_v
                    cty_coord_n[cty] += 1
                except ValueError:
                    pass

            is_r = False
            if i_rucc is not None and len(row) > i_rucc and row[i_rucc].strip():
                try:
                    is_r = float(row[i_rucc]) >= 4.0
                except ValueError:
                    is_r = False
            row_rucc_ge4.append(is_r)

            row_comp.append(row[i_comp].strip() if (i_comp is not None and len(row) > i_comp) else "")
            ep = row[i_exit].strip().split("_", 1)[0] if (i_exit is not None and len(row) > i_exit) else ""
            row_exit_prefix.append(ep)
            row_age.append(row[i_age].strip().lower() if (i_age is not None and len(row) > i_age) else "")

    n_total = len(rows_raw)
    wk_counts = Counter(w for w in row_wk if w >= 0)
    viable_wks = [w for w, cnt in wk_counts.items() if cnt >= min_pool]
    n_weeks = (max(viable_wks) + 1) if viable_wks else 0
    if n_weeks <= 0:
        raise RuntimeError(f"No week >= 0 holds >= {min_pool} cases in {linelist_path}")

    county_cases = Counter(row_county)
    county_coords = {
        c: (cty_lat_sum[c] / cty_coord_n[c], cty_lon_sum[c] / cty_coord_n[c])
        for c in county_cases
        if cty_coord_n[c] > 0
    }

    sentinel_hubs = select_regional_hubs(
        dict(county_cases),
        county_coords,
        k_hubs=pool_cfg.sentinel_k_hubs,
        pop_alpha=pool_cfg.sentinel_hub_pop_alpha,
        min_hub_share=pool_cfg.sentinel_min_hub_share,
    )
    convenience_hubs = select_regional_hubs(
        dict(county_cases),
        county_coords,
        k_hubs=pool_cfg.convenience_k_hubs,
        pop_alpha=0.60,
        min_hub_share=0.02,
    )

    cty_dist_sent: dict[str, float] = {}
    cty_p_sent: dict[str, float] = {}
    cty_p_conv: dict[str, float] = {}
    ordered_ctys = [c for c, _ in county_cases.most_common()]

    for rank_i, cty in enumerate(ordered_ctys):
        if cty in county_coords and sentinel_hubs and all(h in county_coords for h in sentinel_hubs):
            la, lo = county_coords[cty]
            d_s = min(haversine_km_scalar(la, lo, *county_coords[h]) for h in sentinel_hubs)
            d_c = min(haversine_km_scalar(la, lo, *county_coords[h]) for h in convenience_hubs)
        else:
            d_s = 0.0 if cty in sentinel_hubs else 30.0 * (rank_i + 1)
            d_c = 0.0 if cty in convenience_hubs else 45.0 * (rank_i + 1)
        cty_dist_sent[cty] = round(d_s, 2)
        cty_p_sent[cty] = pool_cfg.sentinel_p_floor + (
            pool_cfg.sentinel_p_hub - pool_cfg.sentinel_p_floor
        ) * math.exp(-d_s / max(pool_cfg.sentinel_decay_km, 1e-6))
        cty_p_conv[cty] = pool_cfg.convenience_p_floor + (
            pool_cfg.convenience_p_hub - pool_cfg.convenience_p_floor
        ) * math.exp(-d_c / max(pool_cfg.convenience_decay_km, 1e-6))

    # Clinical tier classification (`asymp`, `outpatient`, `inpatient`)
    has_explicit_ih = any(ep.startswith(("I", "hM", "dM")) for ep in row_exit_prefix)
    rng_tier = random.Random(seed + 7919)
    row_tier: list[str] = []
    for ep, ag in zip(row_exit_prefix, row_age):
        if has_explicit_ih:
            if ep.startswith(("hM", "dM")):
                row_tier.append("inpatient")
            elif ep.startswith(("A", "P")):
                row_tier.append("asymp")
            else:
                row_tier.append("outpatient")
        else:
            if ep.startswith("A"):
                row_tier.append("asymp")
            else:
                hp = 0.24 if ("65" in ag or "senior" in ag) else (
                     0.11 if ("50-64" in ag or "older" in ag) else (
                     0.02 if ("0-4" in ag or "5-17" in ag or "school" in ag or "preschool" in ag) else 0.055))
                row_tier.append("inpatient" if rng_tier.random() < hp else "outpatient")

    comp_counts = Counter(c for c in row_comp if c)
    rng_s = random.Random(seed + 101)
    rng_c = random.Random(seed + 202)
    ship_sent: list[int] = [0] * n_total
    ship_conv: list[int] = [0] * n_total

    for i in range(n_total):
        cty = row_county[i]
        tier = row_tier[i]
        ms = (
            pool_cfg.sentinel_mult_inpatient if tier == "inpatient"
            else (pool_cfg.sentinel_mult_asymp if tier == "asymp" else pool_cfg.sentinel_mult_outpatient)
        )
        mc = (
            pool_cfg.convenience_mult_inpatient if tier == "inpatient"
            else (pool_cfg.convenience_mult_asymp if tier == "asymp" else pool_cfg.convenience_mult_outpatient)
        )
        if row_comp[i] and comp_counts[row_comp[i]] >= pool_cfg.convenience_outbreak_min_cases:
            mc *= pool_cfg.convenience_mult_outbreak

        ps = min(0.98, cty_p_sent[cty] * ms)
        pc = min(0.98, cty_p_conv[cty] * mc)
        if rng_s.random() < ps:
            ship_sent[i] = 1
        if rng_c.random() < pc:
            ship_conv[i] = 1

    # Compute clinical-tier O/E diagnostics
    tot_sent = max(sum(ship_sent), 1)
    tot_conv = max(sum(ship_conv), 1)
    tier_oe: dict[str, dict[str, float]] = {}
    for tname in ("outpatient", "inpatient", "asymp"):
        n_t = sum(1 for t in row_tier if t == tname)
        p_all = n_t / max(n_total, 1)
        p_s = sum(1 for i, t in enumerate(row_tier) if t == tname and ship_sent[i] == 1) / tot_sent
        p_c = sum(1 for i, t in enumerate(row_tier) if t == tname and ship_conv[i] == 1) / tot_conv
        tier_oe[tname] = {
            "linelist_share": round(p_all, 4),
            "sentinel_share": round(p_s, 4),
            "sentinel_OE": round(p_s / p_all, 3) if p_all > 0 else 0.0,
            "convenience_share": round(p_c, 4),
            "convenience_OE": round(p_c / p_all, 3) if p_all > 0 else 0.0,
        }

    # Index rows by week 0 .. n_weeks - 1
    by_wk: list[list[int]] = [[] for _ in range(n_weeks)]
    for i, w in enumerate(row_wk):
        if 0 <= w < n_weeks:
            by_wk[w].append(i)

    out_header = header + ["ship_sentinel", "ship_convenience"]
    written_files: list[Path] = []
    strategy_summaries: list[dict] = []
    all_inpat_share = tier_oe["inpatient"]["linelist_share"]

    print(f"  State={state_label} | {n_total:,} cases across {len(county_cases)} counties, {n_weeks} active weeks")
    print(f"  Sentinel K={len(sentinel_hubs)} hubs    : {', '.join(sentinel_hubs)} ({sum(ship_sent):,} swabs, {sum(ship_sent)/max(n_total,1)*100:.1f}%)")
    print(f"  Convenience K={len(convenience_hubs)} hubs : {', '.join(convenience_hubs)} ({sum(ship_conv):,} swabs, {sum(ship_conv)/max(n_total,1)*100:.1f}%)")

    for stride in strides:
        stride_i = int(stride)
        for strat_idx, strat in enumerate(strategies):
            rng = random.Random(seed + 1000 * stride_i + 37 * (strat_idx + 1))
            picked_all: list[int] = []
            for w_lo in range(0, n_weeks, stride_i):
                w_hi = min(w_lo + stride_i - 1, n_weeks - 1)
                span_wks = w_hi - w_lo + 1
                win_budget = int(weekly_budget) * span_wks
                win_idxs: list[int] = []
                for w in range(w_lo, w_hi + 1):
                    win_idxs.extend(by_wk[w])
                picked_all.extend(
                    _draw_window_indices(win_idxs, row_county, ship_sent, ship_conv, strat, win_budget, rng)
                )

            recipe_id = f"{stride_i}S__{strat}"
            fpath = out_path / f"{recipe_id}_samples.csv.xz"
            with lzma.open(fpath, "wt", encoding="utf-8", newline="", preset=1) as out_fh:
                writer = csv.writer(out_fh)
                writer.writerow(out_header)
                for idx in picked_all:
                    writer.writerow(rows_raw[idx] + [str(ship_sent[idx]), str(ship_conv[idx])])
            written_files.append(fpath)

            n_s = max(len(picked_all), 1)
            ctys_hit = len({row_county[i] for i in picked_all})
            rur_sh = round(sum(1 for i in picked_all if row_rucc_ge4[i]) / n_s, 4) if picked_all else 0.0
            inp_sh = sum(1 for i in picked_all if row_tier[i] == "inpatient") / n_s if picked_all else 0.0
            inp_oe = round(inp_sh / all_inpat_share, 3) if all_inpat_share > 0 else 0.0

            strategy_summaries.append({
                "state": state_label,
                "recipe": recipe_id,
                "stride_weeks": stride_i,
                "strategy": strat,
                "n_samples": len(picked_all),
                "counties_sampled": ctys_hit,
                "counties_total": len(county_cases),
                "rural_rucc4_9_share": rur_sh,
                "inpatient_OE_vs_linelist": inp_oe,
                "samples_file": str(fpath),
            })
            print(
                f"  [{recipe_id:<24}] {len(picked_all):>6,} samples | "
                f"{ctys_hit:>3d}/{len(county_cases)} counties | "
                f"rural={rur_sh*100:5.1f}% | inpatient O/E={inp_oe:.2f}"
            )

    diag = {
        "state": state_label,
        "n_linelist_cases": n_total,
        "n_weeks": n_weeks,
        "n_counties_total": len(county_cases),
        "sentinel_hubs": sentinel_hubs,
        "convenience_hubs": convenience_hubs,
        "n_ship_sentinel": sum(ship_sent),
        "sentinel_pool_fraction": round(sum(ship_sent) / max(n_total, 1), 4),
        "n_ship_convenience": sum(ship_conv),
        "convenience_pool_fraction": round(sum(ship_conv) / max(n_total, 1), 4),
        "clinical_tier_OE": tier_oe,
        "county_nearest_sentinel_km": cty_dist_sent,
        "strategy_summaries": strategy_summaries,
    }

    if diagnostics_dir is not None:
        ddir = Path(diagnostics_dir)
        ddir.mkdir(parents=True, exist_ok=True)
        (ddir / f"{state_label}_sentinel_hubs.json").write_text(json.dumps(diag, indent=2), encoding="utf-8")
        if strategy_summaries:
            with open(ddir / f"{state_label}_canonical_pool_summary.csv", "w", encoding="utf-8", newline="") as sfh:
                w = csv.DictWriter(sfh, fieldnames=list(strategy_summaries[0].keys()))
                w.writeheader()
                w.writerows(strategy_summaries)

    return written_files, diag
