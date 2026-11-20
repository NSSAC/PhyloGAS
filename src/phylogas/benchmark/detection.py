"""Detection lead-time benchmark (T95, UT95, ΔT95, and true cases at T95).

Adapts the T95 / UT95 variant-detection framework (Wohl et al. 2022; Moon et al.
WSC 2024) to PhyloGAS's weekly/stride genomic surveillance architecture:

* **Target unit (default: established importation chains):** Any column in
  ``allevents_variants.csv.xz`` (default ``component_id``, optionally
  ``variant_benchmark``) filtered by an establishment threshold
  (``--min-chain-size``, default 100 true infections) and restricted to
  introductions originating within the active surveillance window
  (``origin_week >= 0``).
* **UT95 (homogeneous baseline):** Earliest week offset ``W - w0`` at which
  cumulative homogeneous detection probability reaches ``alpha`` (default 0.95)
  under true statewide weekly prevalence ``pi_v(w) = I_v(w) / I_all(w)`` and
  the strategy's realized weekly sequencing volume ``B_w``:

      P_UT(W; v) = 1 - prod_{w=w0..W} (1 - pi_v(w)) ** B_w

* **T95 (stratified hypergeometric):** Exact zero-noise detection time for a
  strategy's realized weekly stratum allocations ``n_{w,g}`` drawn without
  replacement from ascertained line-list stratum pools ``N_{w,g}`` containing
  ``k_{w,g}(v)`` cases of target ``v``:

      q_w(v) = prod_g C(N_{w,g} - k_{w,g}(v), n_{w,g}) / C(N_{w,g}, n_{w,g})
      P_strat(W; v) = 1 - prod_{w=w0..W} q_w(v)

* **Paired ΔT95 = T95 - UT95:** Computed per target first, then summarized by
  median and IQR overall and stratified by seeding geography (``metro`` vs.
  ``rural`` via USDA RUCC), both among detected targets and right-censored at
  the surveillance horizon across all established targets (preventing
  survivorship bias when convenience/hub designs miss rural chains), alongside
  cumulative true infections at T95 and empirical first-detection lag.
"""

from __future__ import annotations

import csv
import datetime as _dt
import glob as _glob
import gzip
import lzma
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


_EMPTY_TOKENS = {"", "none", "nan", "na", "null", "-1", "background", "unassigned"}


def _open_text(path: Path | str):
    p = str(path)
    if p.endswith(".xz"):
        return lzma.open(p, "rt", encoding="utf-8", newline="")
    if p.endswith(".gz"):
        return gzip.open(p, "rt", encoding="utf-8", newline="")
    return open(p, "rt", encoding="utf-8", newline="")


def _parse_iso_ordinal(raw: str | None) -> int | None:
    if not raw:
        return None
    s = raw.strip()[:10]
    if len(s) == 10 and s[4] == "-" and s[7] == "-":
        try:
            return _dt.date(int(s[0:4]), int(s[5:7]), int(s[8:10])).toordinal()
        except ValueError:
            return None
    return None


def _week_index_from_ordinal(day_ord: int, start_ord: int) -> int:
    """Match BeyondBaseline / truth_runner calendar: week 0 = [start-7d, start-1d]."""
    week0_ord = start_ord - 7
    return (day_ord - week0_ord) // 7


def _strategy_label_from_path(path: Path) -> str:
    """Recover clean recipe/algorithm label matching truth_runner._label_from_path."""
    name = path.name
    for suffix in (".csv.xz", ".csv.gz", ".csv"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    name = re.sub(r"_samples$", "", name)
    m = re.search(r"scenario\d+_(.+)$", name)
    if m:
        return m.group(1)
    return name


def _parse_rucc(raw: str | None) -> int | None:
    if not raw:
        return None
    s = raw.strip()
    if not s:
        return None
    try:
        return int(float(s))
    except ValueError:
        return None


def _hypergeom_miss_prob(N: int, k: int, n: int) -> float:
    """Exact probability of drawing 0 target items in n draws without replacement
    from a pool of size N containing k target items."""
    if k <= 0 or n <= 0 or N <= 0:
        return 1.0
    n_eff = min(n, N)
    if k >= N or n_eff > (N - k):
        return 0.0
    log_p = (
        math.lgamma(N - k + 1)
        - math.lgamma(N - k - n_eff + 1)
        - math.lgamma(N + 1)
        + math.lgamma(N - n_eff + 1)
    )
    return min(1.0, max(0.0, math.exp(log_p)))


def _quantile(values: list[float], q: float) -> float:
    """Linear-interpolation quantile matching numpy/pandas default."""
    if not values:
        return float("nan")
    xs = sorted(values)
    n = len(xs)
    if n == 1:
        return float(xs[0])
    pos = q * (n - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return float(xs[lo])
    frac = pos - lo
    return float(xs[lo] * (1.0 - frac) + xs[hi] * frac)


def _stratum_key_positional(row: list[str], col_indices: tuple[int, ...]) -> str:
    if not col_indices:
        return "__all__"
    if len(col_indices) == 1:
        idx = col_indices[0]
        return row[idx].strip() if idx < len(row) else ""
    n = len(row)
    return "|".join(row[i].strip() if i < n else "" for i in col_indices)


@dataclass
class TargetMeta:
    target_id: str
    origin_week: int
    total_true_cases: int
    seed_county: str
    seed_rucc: int | None
    geo_group: str  # "metro", "rural", or "unknown"


@dataclass
class TruthIndex:
    start_date_iso: str
    start_ord: int
    n_weeks: int
    pid_to_target: dict[str, str]
    weekly_total_true: dict[int, int]
    weekly_target_true: dict[str, dict[int, int]]
    targets: dict[str, TargetMeta]


def load_truth_index(
    infections_path: Path,
    *,
    target_col: str = "component_id",
    date_field: str = "date",
    start_date: str | None = None,
    min_chain_size: int = 100,
    metro_rucc_max: int = 3,
    n_weeks_override: int | None = None,
) -> TruthIndex:
    """Stream ``allevents`` in a single positional CSV pass, deduplicate on
    ``alias_pid``, and index true weekly counts and established target metadata."""
    explicit_start_ord: int | None = None
    explicit_start_iso: str | None = None
    if start_date and str(start_date).strip().lower() not in ("", "auto"):
        s = str(start_date).strip()[:10]
        explicit_start_ord = _parse_iso_ordinal(s)
        if explicit_start_ord is None:
            raise ValueError(f"Invalid --start-date '{start_date}'; expected YYYY-MM-DD.")
        explicit_start_iso = s

    # Deduplicate by alias_pid in one pass: value = (sort_key, d_ord, target_id, county, rucc)
    seen_pid: dict[str, tuple[int, int, str, str, int | None]] = {}
    min_ord: int | None = None

    with _open_text(infections_path) as fh:
        reader = csv.reader(fh)
        try:
            raw_hdr = next(reader)
        except StopIteration:
            raise ValueError(f"Empty infections file: {infections_path}")
        fnames = [c.strip() for c in raw_hdr]
        col_map = {c: i for i, c in enumerate(fnames)}

        if date_field not in col_map:
            raise ValueError(f"Date column '{date_field}' not found in {infections_path}")
        if target_col not in col_map:
            raise ValueError(
                f"Target column '{target_col}' not found in {infections_path}. "
                f"Available columns: {fnames}"
            )

        date_idx = col_map[date_field]
        tgt_idx = col_map[target_col]
        pid_idx = col_map.get("alias_pid", col_map.get("sim_pid", col_map.get("pid", -1)))
        tick_idx = col_map.get("exposure_tick", col_map.get("sim_tick", -1))
        cty_idx = col_map.get("county", -1)
        rucc_idx = col_map.get("rucc_code", -1)

        for row_idx, row in enumerate(reader):
            if date_idx >= len(row):
                continue
            d_ord = _parse_iso_ordinal(row[date_idx])
            if d_ord is None:
                continue
            if min_ord is None or d_ord < min_ord:
                min_ord = d_ord

            pid = (row[pid_idx].strip() if 0 <= pid_idx < len(row) else "") or f"__row_{row_idx}"
            if 0 <= tick_idx < len(row):
                ts = row[tick_idx].strip()
                if ts:
                    try:
                        sort_key = int(ts) if ts.isdigit() else int(float(ts))
                    except ValueError:
                        sort_key = d_ord
                else:
                    sort_key = d_ord
            else:
                sort_key = d_ord

            tgt = row[tgt_idx].strip() if tgt_idx < len(row) else ""
            cty = row[cty_idx].strip() if 0 <= cty_idx < len(row) else ""
            rucc = _parse_rucc(row[rucc_idx]) if 0 <= rucc_idx < len(row) else None

            prev = seen_pid.get(pid)
            if prev is None or sort_key < prev[0]:
                seen_pid[pid] = (sort_key, d_ord, tgt, cty, rucc)

    if explicit_start_ord is not None and explicit_start_iso is not None:
        start_ord = explicit_start_ord
        start_iso = explicit_start_iso
    else:
        if min_ord is None:
            raise ValueError(f"No parseable '{date_field}' values in {infections_path}")
        start_ord = min_ord + 7
        start_iso = _dt.date.fromordinal(start_ord).isoformat()

    week0_ord = start_ord - 7
    max_wk = -1
    for _, d_ord, _, _, _ in seen_pid.values():
        wk = (d_ord - week0_ord) // 7
        if wk > max_wk:
            max_wk = wk
    n_weeks = n_weeks_override if n_weeks_override is not None else (max_wk + 1)

    pid_to_target: dict[str, str] = {}
    weekly_total_true: dict[int, int] = Counter()
    weekly_target_true: dict[str, dict[int, int]] = defaultdict(Counter)
    target_earliest: dict[str, tuple[int, int, str, int | None]] = {}
    target_window_total: dict[str, int] = Counter()

    for pid, (sort_key, d_ord, tgt, cty, rucc) in seen_pid.items():
        wk = (d_ord - week0_ord) // 7
        if 0 <= wk < n_weeks:
            weekly_total_true[wk] += 1
        if tgt.lower() in _EMPTY_TOKENS:
            continue
        pid_to_target[pid] = tgt
        if 0 <= wk < n_weeks:
            weekly_target_true[tgt][wk] += 1
            target_window_total[tgt] += 1
        earl = target_earliest.get(tgt)
        if earl is None or sort_key < earl[0]:
            target_earliest[tgt] = (sort_key, wk, cty, rucc)

    targets: dict[str, TargetMeta] = {}
    for tgt, tot in target_window_total.items():
        if tot < min_chain_size:
            continue
        _, origin_wk, seed_cty, seed_rucc = target_earliest[tgt]
        # Exclude chains already circulating before surveillance week 0
        if origin_wk < 0 or origin_wk >= n_weeks:
            continue
        if seed_rucc is None:
            geo = "unknown"
        elif seed_rucc <= metro_rucc_max:
            geo = "metro"
        else:
            geo = "rural"
        targets[tgt] = TargetMeta(
            target_id=tgt,
            origin_week=origin_wk,
            total_true_cases=tot,
            seed_county=seed_cty,
            seed_rucc=seed_rucc,
            geo_group=geo,
        )

    return TruthIndex(
        start_date_iso=start_iso,
        start_ord=start_ord,
        n_weeks=n_weeks,
        pid_to_target=pid_to_target,
        weekly_total_true=dict(weekly_total_true),
        weekly_target_true={t: dict(weekly_target_true[t]) for t in targets},
        targets=targets,
    )


@dataclass
class LinelistPoolIndex:
    n_weeks: int
    # weekly_stratum_N[wk][stratum] = total ascertained cases in stratum
    weekly_stratum_N: dict[int, dict[str, int]]
    # weekly_target_stratum_k[target][wk][stratum] = ascertained target cases in stratum
    weekly_target_stratum_k: dict[str, dict[int, dict[str, int]]]


def load_linelist_pool(
    linelist_path: Path,
    truth: TruthIndex,
    *,
    stratum_cols: tuple[str, ...] = ("county",),
    date_field: str = "date",
    target_col: str = "component_id",
    min_pool: int = 50,
) -> LinelistPoolIndex:
    """Index ascertained weekly stratum pools ``N_{w,g}`` and target counts
    ``k_{w,g}(v)`` from ``linelist.csv.xz`` in a single positional pass."""
    weekly_stratum_N: dict[int, dict[str, int]] = defaultdict(Counter)
    weekly_target_stratum_k: dict[str, dict[int, dict[str, int]]] = defaultdict(
        lambda: defaultdict(Counter)
    )
    weekly_counts: dict[int, int] = Counter()
    week0_ord = truth.start_ord - 7
    valid_targets = set(truth.targets.keys())
    pid_to_target = truth.pid_to_target

    with _open_text(linelist_path) as fh:
        reader = csv.reader(fh)
        raw_hdr = next(reader, [])
        fnames = [c.strip() for c in raw_hdr]
        col_map = {c: i for i, c in enumerate(fnames)}

        date_idx = col_map.get(date_field, -1)
        pid_idx = col_map.get("alias_pid", col_map.get("sim_pid", col_map.get("pid", -1)))
        tgt_idx = col_map.get(target_col, -1)
        if not stratum_cols or stratum_cols == ("__all__",):
            strat_indices: tuple[int, ...] = ()
        else:
            strat_indices = tuple(col_map.get(c, -1) for c in stratum_cols)

        for row in reader:
            if date_idx < 0 or date_idx >= len(row):
                continue
            d_ord = _parse_iso_ordinal(row[date_idx])
            if d_ord is None:
                continue
            wk = (d_ord - week0_ord) // 7
            if wk < 0:
                continue
            weekly_counts[wk] += 1
            strat = _stratum_key_positional(row, strat_indices) if strat_indices else "__all__"
            weekly_stratum_N[wk][strat] += 1

            pid = row[pid_idx].strip() if 0 <= pid_idx < len(row) else ""
            tgt = pid_to_target.get(pid)
            if tgt is None and 0 <= tgt_idx < len(row):
                cand = row[tgt_idx].strip()
                if cand in valid_targets:
                    tgt = cand
            if tgt is not None and tgt in valid_targets:
                weekly_target_stratum_k[tgt][wk][strat] += 1

    viable = [w for w, c in weekly_counts.items() if c >= min_pool]
    n_weeks = (max(viable) + 1) if viable else truth.n_weeks

    return LinelistPoolIndex(
        n_weeks=n_weeks,
        weekly_stratum_N={w: dict(c) for w, c in weekly_stratum_N.items()},
        weekly_target_stratum_k={
            t: {w: dict(sc) for w, sc in wmap.items()}
            for t, wmap in weekly_target_stratum_k.items()
        },
    )


def compute_ut95_weeks(
    weekly_target_true: dict[int, int],
    weekly_total_true: dict[int, int],
    weekly_budget: dict[int, int],
    origin_week: int,
    n_weeks: int,
    alpha: float = 0.95,
) -> int | None:
    """Homogeneous-mixing detection lag (weeks after ``origin_week``) to reach
    cumulative detection probability ``alpha`` under true prevalence."""
    surv = 1.0
    for wk in range(origin_week, n_weeks):
        bw = weekly_budget.get(wk, 0)
        tot = weekly_total_true.get(wk, 0)
        kv = weekly_target_true.get(wk, 0)
        if bw > 0 and tot > 0 and kv > 0:
            p = kv / tot
            if p >= 1.0:
                surv = 0.0
            else:
                surv *= (1.0 - p) ** bw
        if (1.0 - surv) >= alpha - 1e-12:
            return wk - origin_week
    return None


def compute_stratified_t95_weeks(
    target_stratum_k: dict[int, dict[str, int]],
    weekly_stratum_N: dict[int, dict[str, int]],
    weekly_stratum_n: dict[int, dict[str, int]],
    origin_week: int,
    n_weeks: int,
    alpha: float = 0.95,
) -> int | None:
    """Exact stratified hypergeometric detection lag (weeks after ``origin_week``)
    to reach cumulative detection probability ``alpha``."""
    surv = 1.0
    for wk in range(origin_week, n_weeks):
        k_map = target_stratum_k.get(wk)
        n_map = weekly_stratum_n.get(wk)
        if k_map and n_map:
            N_map = weekly_stratum_N.get(wk, {})
            for strat, kv in k_map.items():
                nv = n_map.get(strat, 0)
                if kv > 0 and nv > 0:
                    Nv = N_map.get(strat, 0)
                    if Nv > 0:
                        surv *= _hypergeom_miss_prob(Nv, kv, nv)
                        if surv <= 0.0:
                            break
        if (1.0 - surv) >= alpha - 1e-12:
            return wk - origin_week
    return None


def _score_one_sample_file(
    sample_path: Path,
    truth: TruthIndex,
    pool: LinelistPoolIndex | None,
    *,
    stratum_cols: tuple[str, ...],
    date_field: str,
    target_col: str,
    alpha: float,
) -> tuple[str, list[dict[str, object]]]:
    strategy = _strategy_label_from_path(sample_path)
    n_weeks = pool.n_weeks if pool is not None else truth.n_weeks
    week0_ord = truth.start_ord - 7
    valid_targets = set(truth.targets.keys())
    pid_to_target = truth.pid_to_target

    weekly_budget: dict[int, int] = Counter()
    weekly_stratum_n: dict[int, dict[str, int]] = defaultdict(Counter)
    emp_first_wk: dict[str, int] = {}

    with _open_text(sample_path) as fh:
        reader = csv.reader(fh)
        raw_hdr = next(reader, [])
        fnames = [c.strip() for c in raw_hdr]
        col_map = {c: i for i, c in enumerate(fnames)}

        date_idx = col_map.get(date_field, -1)
        pid_idx = col_map.get("alias_pid", col_map.get("sim_pid", col_map.get("pid", -1)))
        tgt_idx = col_map.get(target_col, -1)
        if not stratum_cols or stratum_cols == ("__all__",):
            strat_indices: tuple[int, ...] = ()
        else:
            strat_indices = tuple(col_map.get(c, -1) for c in stratum_cols)

        for row in reader:
            if date_idx < 0 or date_idx >= len(row):
                continue
            d_ord = _parse_iso_ordinal(row[date_idx])
            if d_ord is None:
                continue
            wk = (d_ord - week0_ord) // 7
            if wk < 0 or wk >= n_weeks:
                continue
            weekly_budget[wk] += 1
            strat = _stratum_key_positional(row, strat_indices) if strat_indices else "__all__"
            weekly_stratum_n[wk][strat] += 1

            pid = row[pid_idx].strip() if 0 <= pid_idx < len(row) else ""
            tgt = pid_to_target.get(pid)
            if tgt is None and 0 <= tgt_idx < len(row):
                cand = row[tgt_idx].strip()
                if cand in valid_targets:
                    tgt = cand
            if tgt is not None and tgt in valid_targets:
                prev_w = emp_first_wk.get(tgt)
                if prev_w is None or wk < prev_w:
                    emp_first_wk[tgt] = wk

    events: list[dict[str, object]] = []
    for tgt_id, meta in sorted(truth.targets.items()):
        w0 = meta.origin_week
        if w0 >= n_weeks:
            continue
        horizon_wks = n_weeks - w0
        true_wmap = truth.weekly_target_true[tgt_id]
        ut95 = compute_ut95_weeks(
            true_wmap, truth.weekly_total_true, weekly_budget, w0, n_weeks, alpha=alpha
        )

        if pool is not None:
            k_wmap = pool.weekly_target_stratum_k.get(tgt_id, {})
            t95 = compute_stratified_t95_weeks(
                k_wmap, pool.weekly_stratum_N, weekly_stratum_n, w0, n_weeks, alpha=alpha
            )
        else:
            t95 = None

        delta_t95 = (t95 - ut95) if (t95 is not None and ut95 is not None) else None
        # Right-censored lag at surveillance horizon (avoids survivorship bias on missed chains)
        t95_cens = t95 if t95 is not None else horizon_wks
        ut95_cens = ut95 if ut95 is not None else horizon_wks
        delta_t95_cens = (t95_cens - ut95_cens) if pool is not None else None

        emp_wk = emp_first_wk.get(tgt_id)
        emp_lag = (emp_wk - w0) if (emp_wk is not None and emp_wk >= w0) else None

        # True cumulative infections and variant share at T95 (or empirical first det if no pool)
        ref_lag = t95 if t95 is not None else emp_lag
        if ref_lag is not None:
            ref_wk = w0 + ref_lag
            cum_true = sum(cnt for w, cnt in true_wmap.items() if w0 <= w <= ref_wk)
            tot_ref = truth.weekly_total_true.get(ref_wk, 0)
            p95_share = (true_wmap.get(ref_wk, 0) / tot_ref) if tot_ref > 0 else float("nan")
        else:
            cum_true = None
            p95_share = None

        events.append({
            "strategy": strategy,
            "target_id": tgt_id,
            "geo_group": meta.geo_group,
            "seed_county": meta.seed_county,
            "seed_rucc": meta.seed_rucc if meta.seed_rucc is not None else "",
            "origin_week": w0,
            "horizon_weeks": horizon_wks,
            "total_true_cases": meta.total_true_cases,
            "ut95_weeks": ut95 if ut95 is not None else "",
            "t95_weeks": t95 if t95 is not None else "",
            "delta_t95_weeks": delta_t95 if delta_t95 is not None else "",
            "t95_censored_weeks": t95_cens if pool is not None else "",
            "delta_t95_censored_weeks": delta_t95_cens if delta_t95_cens is not None else "",
            "emp_first_det_weeks": emp_lag if emp_lag is not None else "",
            "true_cases_at_t95": cum_true if cum_true is not None else "",
            "p95_true_share": f"{p95_share:.6f}" if isinstance(p95_share, float) and not math.isnan(p95_share) else "",
        })

    return strategy, events


def _summarize_events(strategy: str, events: list[dict[str, object]]) -> list[dict[str, object]]:
    """Aggregate per-target event rows into ``all``, ``metro``, and ``rural`` rows."""
    out: list[dict[str, object]] = []
    groups = [("all", events)]
    for g in ("metro", "rural"):
        sub = [e for e in events if e["geo_group"] == g]
        if sub:
            groups.append((g, sub))

    def _fmt(v: float | None, digits: int = 2) -> str:
        if v is None or (isinstance(v, float) and math.isnan(v)):
            return ""
        return f"{v:.{digits}f}"

    for geo_label, sub in groups:
        n_tgt = len(sub)
        ut95_vals = [float(e["ut95_weeks"]) for e in sub if e["ut95_weeks"] != ""]
        t95_vals = [float(e["t95_weeks"]) for e in sub if e["t95_weeks"] != ""]
        delta_vals = [float(e["delta_t95_weeks"]) for e in sub if e["delta_t95_weeks"] != ""]
        t95_cens_vals = [float(e["t95_censored_weeks"]) for e in sub if e["t95_censored_weeks"] != ""]
        delta_cens_vals = [float(e["delta_t95_censored_weeks"]) for e in sub if e["delta_t95_censored_weeks"] != ""]
        emp_vals = [float(e["emp_first_det_weeks"]) for e in sub if e["emp_first_det_weeks"] != ""]
        cases_vals = [float(e["true_cases_at_t95"]) for e in sub if e["true_cases_at_t95"] != ""]
        share_vals = [float(e["p95_true_share"]) for e in sub if e["p95_true_share"] != ""]

        iqr_delta = (
            _quantile(delta_vals, 0.75) - _quantile(delta_vals, 0.25)
            if delta_vals
            else float("nan")
        )
        out.append({
            "strategy": strategy,
            "geo_group": geo_label,
            "n_targets": n_tgt,
            "detection_rate_95": _fmt(len(t95_vals) / n_tgt if n_tgt else float("nan"), 4),
            "detection_rate_emp": _fmt(len(emp_vals) / n_tgt if n_tgt else float("nan"), 4),
            "median_ut95_wk": _fmt(_quantile(ut95_vals, 0.5), 2),
            "median_t95_wk": _fmt(_quantile(t95_vals, 0.5), 2),
            "median_delta_t95_wk": _fmt(_quantile(delta_vals, 0.5), 2),
            "iqr_delta_t95_wk": _fmt(iqr_delta, 2),
            "median_t95_censored_wk": _fmt(_quantile(t95_cens_vals, 0.5), 2),
            "median_delta_t95_censored_wk": _fmt(_quantile(delta_cens_vals, 0.5), 2),
            "median_emp_lag_wk": _fmt(_quantile(emp_vals, 0.5), 2),
            "median_true_cases_at_t95": _fmt(_quantile(cases_vals, 0.5), 1),
            "median_p95_share": _fmt(_quantile(share_vals, 0.5), 5),
        })
    return out


def _expand_sample_paths(samples_globs: Iterable[str]) -> list[Path]:
    paths: list[Path] = []
    seen: set[Path] = set()
    for pat in samples_globs:
        matches = sorted(_glob.glob(str(pat)))
        if not matches and Path(pat).is_file():
            matches = [str(pat)]
        for m in matches:
            p = Path(m).resolve()
            if p not in seen and p.is_file():
                seen.add(p)
                paths.append(p)
    return paths


def benchmark_detection(
    samples_globs: list[str],
    infections: str | Path,
    *,
    linelist: str | Path | None = None,
    target_col: str = "component_id",
    date_field: str = "date",
    start_date: str | None = None,
    min_pool: int = 50,
    min_chain_size: int = 100,
    stratum_cols: list[str] | tuple[str, ...] = ("county",),
    metro_rucc_max: int = 3,
    alpha: float = 0.95,
    out_csv: str | Path | None = None,
    out_events_csv: str | Path | None = None,
    truth_cache: TruthIndex | None = None,
    pool_cache: LinelistPoolIndex | None = None,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Run the lean T95 / UT95 / ΔT95 detection benchmark across sample CSVs."""
    sample_paths = _expand_sample_paths(samples_globs)
    if not sample_paths:
        raise FileNotFoundError(f"No sample files matched: {samples_globs}")

    inf_path = Path(infections)
    strat_tuple = tuple(stratum_cols) if stratum_cols else ("county",)

    print(f"Detection benchmark : target='{target_col}' (min_size={min_chain_size}, alpha={alpha})")
    if truth_cache is not None:
        truth = truth_cache
    else:
        print(f"Infections          : {inf_path}")
        truth = load_truth_index(
            inf_path,
            target_col=target_col,
            date_field=date_field,
            start_date=start_date,
            min_chain_size=min_chain_size,
            metro_rucc_max=metro_rucc_max,
        )
        n_metro = sum(1 for m in truth.targets.values() if m.geo_group == "metro")
        n_rural = sum(1 for m in truth.targets.values() if m.geo_group == "rural")
        print(
            f"  Calendar start={truth.start_date_iso}, {truth.n_weeks} weeks; "
            f"{len(truth.targets)} established targets (metro={n_metro}, rural={n_rural})"
        )

    pool: LinelistPoolIndex | None = pool_cache
    if pool is None and linelist:
        ll_path = Path(linelist)
        print(f"Linelist pool       : {ll_path} (strata={list(strat_tuple)})")
        pool = load_linelist_pool(
            ll_path,
            truth,
            stratum_cols=strat_tuple,
            date_field=date_field,
            target_col=target_col,
            min_pool=min_pool,
        )
        print(f"  Active sampling window: {pool.n_weeks} weeks (min_pool={min_pool})")

    summary_rows: list[dict[str, object]] = []
    all_events: list[dict[str, object]] = []
    for sp in sample_paths:
        strategy, events = _score_one_sample_file(
            sp,
            truth,
            pool,
            stratum_cols=strat_tuple,
            date_field=date_field,
            target_col=target_col,
            alpha=alpha,
        )
        all_events.extend(events)
        summary_rows.extend(_summarize_events(strategy, events))

    if out_csv:
        out_p = Path(out_csv)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        if summary_rows:
            with open(out_p, "w", encoding="utf-8", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()))
                w.writeheader()
                w.writerows(summary_rows)
        print(f"Wrote summary ({len(summary_rows)} rows) -> {out_p}")

    if out_events_csv:
        ev_p = Path(out_events_csv)
        ev_p.parent.mkdir(parents=True, exist_ok=True)
        if all_events:
            with open(ev_p, "w", encoding="utf-8", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(all_events[0].keys()))
                w.writeheader()
                w.writerows(all_events)
        print(f"Wrote events  ({len(all_events)} rows) -> {ev_p}")

    # Compact console preview (all targets)
    print("\nSummary (geo_group = all):")
    hdr = f"{'strategy':30s} {'n':>5s} {'det95':>7s} {'UT95':>6s} {'T95':>6s} {'dT95':>6s} {'dT95_all':>8s} {'cases@T95':>10s}"
    print(hdr)
    print("-" * len(hdr))
    for r in summary_rows:
        if r["geo_group"] == "all":
            print(
                f"{str(r['strategy']):30s} {int(r['n_targets']):5d} "
                f"{str(r['detection_rate_95']):>7s} {str(r['median_ut95_wk']):>6s} "
                f"{str(r['median_t95_wk']):>6s} {str(r['median_delta_t95_wk']):>6s} "
                f"{str(r['median_delta_t95_censored_wk']):>8s} "
                f"{str(r['median_true_cases_at_t95']):>10s}"
            )

    return summary_rows, all_events
