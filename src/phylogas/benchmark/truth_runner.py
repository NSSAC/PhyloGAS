"""Score saved sample sets against agent-based-model ground truth.

BeyondBaseline selects which cases to sequence and writes them out with
``--save-samples``. This reads those files plus the ABM's hidden transmission
graph and computes the metrics that need truth:

    stride_variant_prevalence_error  prevalence estimation error, per stride
    stride_component_coverage        distinct sampled / distinct true
                                     transmission components, per stride
    coverage_size_{0,10,100,1000}    Mean Reciprocal Distance, cumulative,
                                     over cases in components larger than N
    8_week_rolling_tree_coverage     Mean Reciprocal Distance, 8-week window
    equity_<age group>               Mean Reciprocal Distance, cumulative,
                                     one series per age group

Mean Reciprocal Distance is (1/|Pt|) * sum_u 1/(d(u,S) + 1): for each
infection, the reciprocal of one plus its distance along the true
transmission tree to the nearest sampled case, averaged. It is 1 when every
infection is sampled and falls toward 0 as infections sit further from any
sample.
Definitions follow BeyondBaseline's run_all_scenarios.py (figures F, I-L, M
and N) as of bb8f2fd^, with two deliberate changes. Distances and tree sizes
use the true transmission graph from the all-events file, where
BeyondBaseline built its graph from the line list's rows and so measured
fragments broken at every unreported person. And the coverage metrics are
averaged over every infection in the simulated window, where BeyondBaseline
averaged over reported cases only -- which folded ascertainment bias into the
denominator instead of measuring it. The metric names follow its
eval_names.py registry.

    cumulative_infections            KL(sample || true infections) over the
                                     demographic groups, cumulative
    stride_window_infections         the same over the stride window

Output matches BeyondBaseline's historical ``AUC_rankings.csv`` schema, so
existing aggregation and plotting keep working:

    eval_type, algorithm, scenario_id, scenario_label, weeks, auc

Every metric uses the sampler's calendar, so a scoring window covers exactly
the weeks the sampler drew from: week w is the seven days ending the day
before start_date + w weeks (week 0 is the week *before* start_date), and the
run ends at the first week whose line list holds fewer than min_pool cases.
Each recipe is scored at its own stride. start_date and min_pool default to
BeyondBaseline's own defaults, which is what the sampler uses when the
pipeline does not pass them.

The transmission graph is built once and reused across every sample set, so
scoring N strategies costs barely more than scoring one.
"""

from __future__ import annotations

import glob as _glob
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from . import truth_metrics as T
from .variants import VARIANT_COLUMN, LEGACY_VARIANT_COLUMN


def _label_from_path(p: Path) -> tuple[str, str]:
    """Recover (scenario, algorithm) from a BeyondBaseline samples filename."""
    name = p.name
    for suffix in (".csv.xz", ".csv.gz", ".csv"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    name = re.sub(r"_samples$", "", name)
    m = re.search(r"scenario(\d+)_(.+)$", name)
    if m:
        return m.group(1), m.group(2)
    return "1", name


def _read_samples(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype={"alias_pid": str, "alias_contact": str,
                                  "sim_pid": str, "pid": str})
    df.columns = [c.strip() for c in df.columns]
    return df


# Tree-size thresholds for coverage_size_<N>: BeyondBaseline's figures I-L.
COVERAGE_SIZE_THRESHOLDS = (0, 10, 100, 1000)
ROLLING_TREE_WEEKS = 8

# Metrics where a larger value is better. Everything else here is an error.
_HIGHER_IS_BETTER_PREFIXES = ("stride_component_coverage", "coverage_size_",
                              "8_week_rolling_tree_coverage", "equity_")


def _higher_is_better(eval_type: str) -> bool:
    return str(eval_type).startswith(_HIGHER_IS_BETTER_PREFIXES)


def _recipe_stride(recipe: str, default: int) -> int:
    """Sampling stride in weeks from a recipe id: '4S__surs' -> 4, '1S-P__...' -> 1.

    BeyondBaseline evaluated each scenario at its own stride, so a weekly
    sampler is scored every week and a four-weekly one every fourth.
    """
    m = re.match(r"(\d+)S", str(recipe))
    return int(m.group(1)) if m else default


def _equity_name(age_group: str) -> str:
    """BeyondBaseline's cleaning of an age-group label into an eval_type."""
    ag = str(age_group)
    for ch in (" ", "/"):
        ag = ag.replace(ch, "_")
    for ch in ("(", ")"):
        ag = ag.replace(ch, "")
    return f"equity_{ag}"


def _series_auc(xs, ys) -> float:
    """Trapezoidal AUC, normalised by span, ignoring NaNs."""
    pairs = [(x, y) for x, y in zip(xs, ys) if y is not None and not (isinstance(y, float) and np.isnan(y))]
    if len(pairs) < 2:
        return float("nan")
    xs2 = [p[0] for p in pairs]
    ys2 = [p[1] for p in pairs]
    span = xs2[-1] - xs2[0]
    trapz_fn = getattr(np, "trapezoid", None) or np.trapz   # trapezoid is NumPy >= 2.0
    return float(trapz_fn(ys2, xs2) / span) if span else float("nan")


DEFAULT_STRATIFIERS = ("age", "race", "county", "sex")   # BeyondBaseline's default


def _bb_calendar_defaults():
    """BeyondBaseline's start date and minimum weekly pool, or (None, None)."""
    try:
        from scenarios_simulation.scenarios_config import (
            START_DATE_DEFAULT, MINIMUM_POOL_SIZE_DEFAULT)
        return pd.Timestamp(START_DATE_DEFAULT), int(MINIMUM_POOL_SIZE_DEFAULT)
    except ImportError:
        return None, None


def _week_index(dates: pd.Series, start_date: pd.Timestamp) -> pd.Series:
    """BeyondBaseline's calendar week: 0 is [start_date - 7d, start_date - 1d]."""
    week0 = start_date - pd.Timedelta(days=7)
    return (pd.to_datetime(dates, errors="coerce") - week0).dt.days // 7


def _n_weeks(ll_wk: pd.Series, min_pool: int) -> int:
    """Weeks the sampler ran: week 0 through the last week with >= min_pool cases.

    Must match BeyondBaseline's sampling_week_count: min_pool marks the end of
    the wave, and sparser weeks inside the window are sampled too.
    """
    counts = ll_wk[ll_wk >= 0].value_counts()
    viable = counts[counts >= min_pool]
    return int(viable.index.max()) + 1 if len(viable) else 0


def resolve_sampling_start(allevents: str, linelist: str, start_date: str = "auto",
                           min_pool: int | None = None, date_field: str = "date") -> pd.Timestamp:
    """The sampler's --start-date, shared by sampling and `benchmark truth`.

    "auto" places week 0 on the first week of the all-events file -- TwinSampler
    writes it for the simulated window only, so its first date is where the
    wave being benchmarked begins. An explicit date is used as given. Sparse
    opening weeks are sampled, not skipped: BeyondBaseline runs through the
    last week with min_pool cases and takes the whole pool whenever the budget
    covers it.
    """
    if min_pool is None:
        min_pool = _bb_calendar_defaults()[1] or 50
    ll_dates = pd.to_datetime(pd.read_csv(linelist, usecols=[date_field])[date_field],
                              errors="coerce").dropna()

    if str(start_date).strip().lower() in ("", "auto"):
        first = pd.to_datetime(pd.read_csv(allevents, usecols=[date_field])[date_field],
                               errors="coerce").min()
        if pd.isna(first):
            sys.exit(f"ERROR: no parseable '{date_field}' in {allevents}")
        start = first.normalize() + pd.Timedelta(days=7)
        how = f"all-events begins {first.date()}"
    else:
        start = pd.Timestamp(start_date)
        how = "configured"

    wk = _week_index(ll_dates, start)
    wk = wk[wk.notna()].astype(int)
    n = _n_weeks(wk, min_pool)
    counts = wk.value_counts()
    print(f"Sampling start: {start.date()} ({how}); week 0 = "
          f"{(start - pd.Timedelta(days=7)).date()} .. {(start - pd.Timedelta(days=1)).date()}")
    if n == 0:
        sys.exit(f"ERROR: no week from {start.date()} on holds {min_pool} or more reported "
                 f"cases, so the sampler has no window. Lower sampling.min_pool.")
    sparse = [w for w in range(n) if counts.get(w, 0) < min_pool]
    print(f"  {n} weeks, through the last week with {min_pool}+ reported cases"
          + (f"; {len(sparse)} sparser week(s) inside the window are sampled too "
             f"(first weeks: {', '.join(str(int(counts.get(w, 0))) for w in range(min(4, n)))} cases)"
             if sparse else ""))
    return start


def _weekly_counts(df: pd.DataFrame, value_col: str, n: int) -> list:
    """value_counts of `value_col` for calendar weeks 0..n-1 (empty weeks kept)."""
    grouped = {int(w): g[value_col].value_counts()
               for w, g in df[df["_wk"].between(0, n - 1)].groupby("_wk")}
    return [grouped.get(w, pd.Series(dtype=float)) for w in range(n)]


def _sum_window(hist: list, lo: int, hi: int) -> pd.Series:
    out = pd.Series(dtype=float)
    for j in range(max(0, lo), min(hi, len(hist) - 1) + 1):
        out = out.add(hist[j], fill_value=0)
    return out


def score_samples(
    samples_globs: list[str],
    infections: str,
    linelist: str | None = None,
    date_field: str = "date",
    out_csv: str | None = None,
    stride_weeks: int = 4,
    population: str | None = None,
    start_date: str | None = None,
    min_pool: int | None = None,
    stratifiers: list[str] | None = None,
) -> pd.DataFrame:
    """Score every matching sample set against ABM truth."""
    bb_start, bb_pool = _bb_calendar_defaults()
    start = pd.Timestamp(start_date) if start_date else bb_start
    pool = int(min_pool) if min_pool is not None else (bb_pool if bb_pool is not None else 50)
    if start is None:
        sys.exit("ERROR: no --start-date given and BeyondBaseline is not importable for its "
                 "default.\n       Pass the start date the sampler used.")
    print(f"Calendar   : week 0 = {(start - pd.Timedelta(days=7)).date()} .. "
          f"{(start - pd.Timedelta(days=1)).date()}, min pool {pool}"
          + ("" if start_date else "  (BeyondBaseline defaults)"))

    print(f"Infections : {infections}")
    inf = pd.read_csv(infections, dtype={"alias_pid": str, "alias_contact": str,
                                         "sim_pid": str, "pid": str})
    inf.columns = [c.strip() for c in inf.columns]
    print(f"  {len(inf):,} rows")
    inf = T._one_row_per_infection(inf, "infections")
    if date_field in inf.columns:
        inf["_wk"] = _week_index(inf[date_field], start)

    vcol = T.variant_column(inf)
    if vcol is None:
        print(f"  WARNING: no '{VARIANT_COLUMN}' (or legacy '{LEGACY_VARIANT_COLUMN}') "
              f"column; prevalence error will be skipped.\n"
              f"           Run `phylogas assign-variants` first.", file=sys.stderr)
    else:
        print(f"  variant column: '{vcol}'")

    # Line list: the sampler's pool. It fixes the number of weeks and is the
    # population the coverage metrics are scored over.
    ll = None
    if linelist:
        ll = pd.read_csv(linelist, dtype={"alias_pid": str, "alias_contact": str})
        ll.columns = [c.strip() for c in ll.columns]
        print(f"Linelist   : {linelist} ({len(ll):,} rows)")

    # The true transmission graph, from the all-events file: every infection
    # and its infector. Distances and tree sizes come from it. BeyondBaseline
    # built its graph from the line list's rows, which carry an edge only for
    # reported cases, so every chain broke at its first unreported person and
    # "tree size" measured those fragments.
    graph = None
    if {"alias_pid", "alias_contact"} <= set(inf.columns):
        graph = T.TreeGraph(inf["alias_pid"], inf["alias_contact"])
        print(f"Graph      : {len(graph):,} infections in {graph.n_trees:,} "
              f"transmission trees (largest {int(graph.size.max()):,})")

    # The line list is the sampler's pool: it fixes how many weeks it ran.
    if ll is not None and date_field in ll.columns:
        ll_wk = _week_index(ll[date_field], start)
    elif date_field in inf.columns:
        ll_wk = inf["_wk"]
    else:
        ll_wk = pd.Series(dtype=float)
    n_weeks = _n_weeks(ll_wk[ll_wk.notna()].astype(int), pool)

    # Who the coverage metrics are scored over: every infection in the
    # simulated window (the all-events file), not only the reported ones. A
    # denominator of reported cases would measure how well sequencing covers
    # what surveillance happened to see, ascertainment bias and all; against
    # all infections the bias shows up in the score, which is the point of
    # having ground truth.
    targets, age_groups = None, []
    if "_wk" in inf.columns:
        tg = inf.loc[inf["_wk"].notna()].copy()
        tg["_wk"] = tg["_wk"].astype(int)
        tg["alias_pid"] = tg["alias_pid"].astype(str)
        tg = T.normalize_age_group_col(tg, "age_group")
        if "age_group" in tg.columns:
            age_groups = sorted(ag for ag in tg["age_group"].dropna().unique() if str(ag) != "nan")
        keep = ["alias_pid", "_wk"] + [c for c in ("age_group", "component_id") if c in tg.columns]
        targets = tg[keep].drop_duplicates("alias_pid")
        if graph is not None:
            targets["_node"] = graph.codes(targets["alias_pid"])
            targets["_size"] = np.where(targets["_node"] >= 0,
                                        graph.size[targets["_node"].clip(lower=0)], 1)
        print(f"Coverage   : {len(targets):,} infections, {n_weeks} sampled weeks"
              + (f", {len(age_groups)} age groups" if age_groups else ""))
    if n_weeks == 0:
        sys.exit(f"ERROR: no week from {start.date()} holds {pool} or more cases; check "
                 f"--start-date / --min-pool against the line list.")

    # True weekly variant counts, on the same calendar.
    true_variants = None
    if vcol and "_wk" in inf.columns:
        true_variants = _weekly_counts(inf, vcol, n_weeks)

    # Demographic groups for the KL metrics, built by BeyondBaseline's own
    # loader so the population recoding and the group keys are identical to
    # what the sampler targeted.
    true_groups, ll_groups, kl_dist = None, None, None
    if population and linelist:
        try:
            from scenarios_simulation.run_all_scenarios import (
                load_linelist_and_population, _normalize_stratifiers)
            from scenarios_simulation.sampling_algorithms import kl_dist
            feats = _normalize_stratifiers(list(stratifiers or DEFAULT_STRATIFIERS))
            ll_g, pop_df, _, _ = load_linelist_and_population(
                linelist, population, date_field, start, pool, feats)
            keys = [k for k in ("alias_pid", "sim_tick") if k in ll_g.columns]
            ll_groups = ll_g[keys + ["group"]].astype({k: str for k in keys}).drop_duplicates(keys)
            pid_col = next((c for c in ("sim_pid", "pid") if c in inf.columns), None)
            pop_pid = "sim_pid" if "sim_pid" in pop_df.columns else "pid"
            if pid_col and "_wk" in inf.columns:
                g = inf[[pid_col, "_wk"]].merge(
                    pop_df[[pop_pid, "group"]].dropna().astype({pop_pid: str}),
                    left_on=pid_col, right_on=pop_pid, how="inner")
                true_groups = _weekly_counts(g, "group", n_weeks)
                print(f"Groups     : {', '.join(feats)}; {len(g):,} infections grouped")
        except ImportError as exc:
            print(f"  WARNING: BeyondBaseline not importable ({exc}); "
                  f"cumulative_infections and stride_window_infections skipped.",
                  file=sys.stderr)
    elif not population:
        print("  note: no --population; cumulative_infections and "
              "stride_window_infections skipped.", file=sys.stderr)

    paths: list[Path] = []
    for spec in samples_globs:
        paths.extend(sorted(Path(p) for p in _glob.glob(spec)))
    if not paths:
        sys.exit(f"ERROR: no sample files matched: {', '.join(samples_globs)}")

    rows = []
    for p in paths:
        scen, algo = _label_from_path(p)
        sdf = _read_samples(p)
        if date_field not in sdf.columns:
            print(f"  WARNING: {p.name} has no '{date_field}' column; skipping",
                  file=sys.stderr)
            continue
        # Only the weeks the sampler drew from, as BeyondBaseline's
        # split_samples_by_calendar_week kept them.
        sdf["_wk"] = _week_index(sdf[date_field], start)
        sdf = sdf.loc[sdf["_wk"].between(0, n_weeks - 1)].copy()
        sdf["_wk"] = sdf["_wk"].astype(int)

        def _add(eval_type, xs, ys):
            if xs:
                rows.append(dict(eval_type=eval_type, algorithm=algo,
                                 scenario_id=scen, scenario_label=scen,
                                 weeks=len(xs), auc=_series_auc(xs, ys),
                                 samples_file=str(p)))

        stride = _recipe_stride(algo, stride_weeks)
        eval_idx = list(range(stride - 1, n_weeks, stride))
        s_col = next((c for c in ("alias_pid", "pid") if c in sdf.columns), None)
        if s_col:
            sdf[s_col] = sdf[s_col].astype(str)

        # --- cumulative_infections / stride_window_infections (figures B, C)
        if true_groups is not None:
            keys = [k for k in ("alias_pid", "sim_tick") if k in ll_groups.columns and k in sdf.columns]
            sg = sdf[keys + ["_wk"]].astype({k: str for k in keys}).merge(ll_groups, on=keys, how="inner")
            sample_groups = _weekly_counts(sg, "group", n_weeks)
            for name, window in (("cumulative_infections", None), ("stride_window_infections", stride)):
                xs, ys = [], []
                for end in eval_idx:
                    lo = 0 if window is None else end - window + 1
                    s_cnt = _sum_window(sample_groups, lo, end)
                    t_cnt = _sum_window(true_groups, lo, end)
                    xs.append(end + 1)
                    ys.append(kl_dist(s_cnt / s_cnt.sum(), t_cnt / t_cnt.sum())
                              if s_cnt.sum() and t_cnt.sum() else float("nan"))
                _add(name, xs, ys)

        # --- stride_variant_prevalence_error (figure E) ---------------------
        svcol = T.variant_column(sdf)
        if true_variants is not None and svcol:
            xs, ys = [], []
            for end in eval_idx:
                lo = max(0, end - stride + 1)
                p_hat = sdf.loc[sdf["_wk"].between(lo, end), svcol].value_counts() \
                           .drop("background", errors="ignore")
                p_true = _sum_window(true_variants, lo, end).drop("background", errors="ignore")
                p_hat = p_hat / p_hat.sum() if p_hat.sum() else p_hat.astype(float)
                p_true = p_true / p_true.sum() if p_true.sum() else p_true.astype(float)
                idx = sorted(set(p_hat.index) | set(p_true.index))
                xs.append(end + 1)
                ys.append(float((p_hat.reindex(idx, fill_value=0.0)
                                 - p_true.reindex(idx, fill_value=0.0)).abs().sum()))
            _add("stride_variant_prevalence_error", xs, ys)

        # --- stride_component_coverage (figure F) ---------------------------
        # Distinct transmission components among the stride window's samples
        # over distinct components among all its infections. A ratio, not a
        # distance.
        if targets is not None and "component_id" in targets.columns and "component_id" in sdf.columns:
            xs, ys = [], []
            for end in eval_idx:
                lo = max(0, end - stride + 1)
                num = sdf.loc[sdf["_wk"].between(lo, end), "component_id"].dropna().nunique()
                den = targets.loc[targets["_wk"].between(lo, end), "component_id"].dropna().nunique()
                xs.append(end + 1)
                ys.append(num / den if den else float("nan"))
            _add("stride_component_coverage", xs, ys)

        # --- Mean Reciprocal Distance family (figures I-L, M, N) -----------
        if targets is not None and graph is not None and s_col:
            def _mrd(rows, sampled, dist):
                return T.coverage_from_codes(rows["_node"].to_numpy(),
                                             rows["alias_pid"].to_numpy(), sampled, dist)

            by_size = {t: ([], []) for t in COVERAGE_SIZE_THRESHOLDS}
            by_age = {ag: ([], []) for ag in age_groups}
            roll_x, roll_y = [], []
            for end in eval_idx:
                # Cumulative: every sample, and every infection up to this
                # week's end (including any before week 0, as BeyondBaseline
                # counted cases).
                # One search serves every size threshold and every age group.
                s_cum = set(sdf.loc[sdf["_wk"] <= end, s_col])
                t_cum = targets[targets["_wk"] <= end]
                dist = graph.distances(s_cum)
                for t in COVERAGE_SIZE_THRESHOLDS:
                    by_size[t][0].append(end + 1)
                    by_size[t][1].append(_mrd(t_cum[t_cum["_size"] > t], s_cum, dist))
                for ag in age_groups:
                    by_age[ag][0].append(end + 1)
                    by_age[ag][1].append(_mrd(t_cum[t_cum["age_group"] == ag], s_cum, dist))

                # Rolling: samples and infections both restricted to the last 8 weeks.
                lo = max(0, end - ROLLING_TREE_WEEKS + 1)
                s_win = set(sdf.loc[sdf["_wk"].between(lo, end), s_col])
                roll_x.append(end + 1)
                roll_y.append(_mrd(targets[targets["_wk"].between(lo, end)],
                                   s_win, graph.distances(s_win)))

            for t in COVERAGE_SIZE_THRESHOLDS:
                _add(f"coverage_size_{t}", *by_size[t])
            _add("8_week_rolling_tree_coverage", roll_x, roll_y)
            for ag in age_groups:
                _add(_equity_name(ag), *by_age[ag])

        print(f"  scored {p.name}  (scenario {scen}, {algo}, stride {stride}, {len(sdf):,} samples)")

    if not rows:
        sys.exit("ERROR: nothing could be scored. Check that the sample files carry\n"
                 "       a date column and that variants have been assigned.")

    df = pd.DataFrame(rows)
    # Rank 1 is best. Coverage metrics are higher-is-better and errors are
    # lower-is-better, so the direction is per eval_type; ranking every metric
    # ascending put the worst coverage first.
    asc = ~df["eval_type"].map(_higher_is_better)
    df["rank_overall"] = df["auc"].where(asc, -df["auc"]).groupby(df["eval_type"]).rank(method="dense")
    df = df.sort_values(["eval_type", "rank_overall", "algorithm"])
    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)
        print(f"\nWrote {len(df)} rows -> {out_csv}")
    return df
