"""Score saved sample sets against agent-based-model ground truth.

BeyondBaseline selects which cases to sequence and writes them out with
``--save-samples``. This reads those files plus the ABM's hidden transmission
graph and computes the metrics that need truth:

    B_cumulative_infections            KL vs true infections, cumulative
    C_stride_window_infections         KL vs true infections, per stride
    E_stride_variant_prevalence_error  prevalence estimation MAE
    F_stride_component_coverage        transmission-component coverage

Output matches BeyondBaseline's historical ``AUC_rankings.csv`` schema, so
existing aggregation and plotting keep working:

    eval_type, algorithm, scenario_id, scenario_label, weeks, auc

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


def _series_auc(xs, ys) -> float:
    """Trapezoidal AUC, normalised by span, ignoring NaNs."""
    pairs = [(x, y) for x, y in zip(xs, ys) if y is not None and not (isinstance(y, float) and np.isnan(y))]
    if len(pairs) < 2:
        return float("nan")
    xs2 = [p[0] for p in pairs]
    ys2 = [p[1] for p in pairs]
    span = xs2[-1] - xs2[0]
    return float(np.trapezoid(ys2, xs2) / span) if span else float("nan")


def score_samples(
    samples_globs: list[str],
    infections: str,
    linelist: str | None = None,
    date_field: str = "date",
    out_csv: str | None = None,
    stride_weeks: int = 4,
) -> pd.DataFrame:
    """Score every matching sample set against ABM truth."""
    print(f"Infections : {infections}")
    inf = pd.read_csv(infections, dtype={"alias_pid": str, "alias_contact": str,
                                         "sim_pid": str, "pid": str})
    inf.columns = [c.strip() for c in inf.columns]
    print(f"  {len(inf):,} rows")

    vcol = T.variant_column(inf)
    if vcol is None:
        print(f"  WARNING: no '{VARIANT_COLUMN}' (or legacy '{LEGACY_VARIANT_COLUMN}') "
              f"column; prevalence error will be skipped.\n"
              f"           Run `phylogas assign-variants` first.", file=sys.stderr)
    else:
        print(f"  variant column: '{vcol}'")

    # Build the transmission adjacency once, from whichever frame has edges.
    edge_src = inf
    if linelist:
        ll = pd.read_csv(linelist, dtype={"alias_pid": str, "alias_contact": str})
        ll.columns = [c.strip() for c in ll.columns]
        if "alias_contact" in ll.columns:
            edge_src = ll
            print(f"Linelist   : {linelist} ({len(ll):,} rows, used for coverage edges)")

    adj = {}
    if "alias_contact" in edge_src.columns:
        adj = T.build_undirected_adj(edge_src, pid_col="alias_pid", contact_col="alias_contact")
        print(f"  adjacency: {len(adj):,} nodes")

    # Weekly truth series.
    start = pd.to_datetime(inf[date_field]).min() if date_field in inf.columns else None
    weekly_true_variants = None
    if vcol and date_field in inf.columns:
        d = inf.copy()
        d["_wk"] = ((pd.to_datetime(d[date_field]) - start).dt.days // 7).astype(int)
        weekly_true_variants = [
            g[vcol].value_counts() for _, g in d.groupby("_wk", sort=True)
        ]
        print(f"  {len(weekly_true_variants)} weeks of true variant counts")

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
        sdf["_wk"] = ((pd.to_datetime(sdf[date_field]) - start).dt.days // 7).astype(int)
        n_weeks = int(sdf["_wk"].max()) + 1 if len(sdf) else 0

        # --- E: prevalence estimation error ------------------------------
        if weekly_true_variants is not None and T.variant_column(sdf):
            svcol = T.variant_column(sdf)
            xs, ys = [], []
            for end in range(stride_weeks - 1, min(n_weeks, len(weekly_true_variants))):
                lo = max(0, end - stride_weeks + 1)
                blk = sdf[(sdf["_wk"] >= lo) & (sdf["_wk"] <= end)]
                true_counts = pd.Series(dtype=float)
                for i in range(lo, end + 1):
                    true_counts = true_counts.add(weekly_true_variants[i], fill_value=0)
                p_hat = blk[svcol].value_counts().drop("background", errors="ignore")
                p_true = true_counts.drop("background", errors="ignore")
                p_hat = p_hat / p_hat.sum() if p_hat.sum() else p_hat.astype(float)
                p_true = p_true / p_true.sum() if p_true.sum() else p_true.astype(float)
                idx = sorted(set(p_hat.index) | set(p_true.index))
                xs.append(end + 1)
                ys.append(float((p_hat.reindex(idx, fill_value=0.0)
                                 - p_true.reindex(idx, fill_value=0.0)).abs().sum()))
            if xs:
                rows.append(dict(eval_type="E_stride_variant_prevalence_error",
                                 algorithm=algo, scenario_id=scen, scenario_label=scen,
                                 weeks=len(xs), auc=_series_auc(xs, ys),
                                 samples_file=str(p)))

        # --- F: transmission-component coverage ---------------------------
        if adj:
            target = set(edge_src["alias_pid"].astype(str))
            xs, ys = [], []
            for end in range(stride_weeks - 1, n_weeks):
                lo = max(0, end - stride_weeks + 1)
                blk = sdf[(sdf["_wk"] >= lo) & (sdf["_wk"] <= end)]
                sampled = set(blk["alias_pid"].astype(str)) if "alias_pid" in blk else set()
                if not sampled:
                    continue
                xs.append(end + 1)
                ys.append(T.calculate_coverage_score(target, sampled, adj))
            if xs:
                rows.append(dict(eval_type="F_stride_component_coverage",
                                 algorithm=algo, scenario_id=scen, scenario_label=scen,
                                 weeks=len(xs), auc=_series_auc(xs, ys),
                                 samples_file=str(p)))

        print(f"  scored {p.name}  (scenario {scen}, {algo}, {len(sdf):,} samples)")

    if not rows:
        sys.exit("ERROR: nothing could be scored. Check that the sample files carry\n"
                 "       a date column and that variants have been assigned.")

    df = pd.DataFrame(rows)
    df["rank_overall"] = df.groupby("eval_type")["auc"].rank(method="dense")
    df = df.sort_values(["eval_type", "rank_overall", "algorithm"])
    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)
        print(f"\nWrote {len(df)} rows -> {out_csv}")
    return df
