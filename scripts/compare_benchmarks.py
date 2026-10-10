#!/usr/bin/env python3
"""Compare benchmark results between runs: does the story still hold?

Reads only the raw per-run outputs the tools write -- never a notebook's
aggregated tables -- so every number here can be traced to a pipeline run.

Two input layouts are recognised, per --run:

  BeyondBaseline replicate directory   replicate_*/AUC_rankings.csv
  (e.g. the CSTE lasso_greedy50 run)   replicate_*/Mugration_Metrics.csv

  PhyloGAS project directory           05_benchmarks/AUC_truth_rankings.csv
  (e.g. results/va_delta_wave)         05_benchmarks/Mugration_Metrics.csv
                                       03_sampled_datasets/AUC_rankings.csv

Old panel codes (I_coverage_size_0), algorithm names and scenario ids are
translated to current metric names and recipe ids with BeyondBaseline's own
registry (eval_names.canonicalize, recipes.recipe_id), so both sides speak
the same vocabulary.

Recipes are ranked WITHIN a sampling stride (1-weekly against 1-weekly,
4-weekly against 4-weekly). BeyondBaseline's AUCs are raw trapezoid areas
while PhyloGAS divides by the span, and the CSVs do not record the x values
needed to convert one to the other; within a stride every recipe shares the
same evaluation weeks, so ranks there need no conversion. It is also the fair
comparison: a weekly and a four-weekly sampler spend the same budget on a
different cadence.

Outputs (--out):
  long.csv            every value, one row per run/state/replicate/recipe/metric
  head_to_head.csv    the --focus recipes side by side, per run and metric
  rank_agreement.csv  per metric and stride: Spearman between the first run's
                      ordering of recipes and each other run's (common states)
  paired.csv          the --focus pair replicate by replicate: how often the
                      second recipe wins, the median difference, and an exact
                      two-sided sign test (needs several replicates)
  tradeoff.csv        per run, state and stride: Spearman across recipes
                      between each coverage metric and mugration cosine.
                      Negative = recipes with more coverage reconstruct flow
                      worse (the abstract's claim).

Example:
  python scripts/compare_benchmarks.py \\
      --run cste=../BeyondBaseline/scripts/scenarios_simulation/analysis/lasso_greedy50:va \\
      --run rerun='results/*_delta_wave' \\
      --out results/compare_cste_vs_rerun
"""

from __future__ import annotations

import argparse
import glob
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Vocabulary: BeyondBaseline's registry when installed, a faithful fallback
# otherwise (the fallback cannot map scenario ids to recipe ids, so old-format
# runs need BeyondBaseline).
# --------------------------------------------------------------------------
try:
    from scenarios_simulation.eval_names import canonicalize as _canon, higher_is_better as _hib
except ImportError:                                    # pragma: no cover
    def _canon(e):
        n = re.sub(r"^[A-N]_", "", str(e))
        return "kl_targets" if n == "targets" else n

    def _hib(e):
        return str(_canon(e)).startswith(("stride_component_coverage", "coverage_size_",
                                          "8_week_rolling_tree_coverage", "equity_"))
try:
    from scenarios_simulation.recipes import recipe_id as _recipe_id
except ImportError:                                    # pragma: no cover
    _recipe_id = None

MUGRATION = {"cosine_similarity": True, "topological_f1": True,
             "pearson_r": True, "masked_mae": False}
COVERAGE_PREFIXES = ("coverage_size_", "8_week_rolling_tree_coverage", "equity_",
                     "stride_component_coverage")


def higher_is_better(metric: str) -> bool:
    return MUGRATION[metric] if metric in MUGRATION else bool(_hib(metric))


def stride_of(recipe: str) -> str:
    m = re.match(r"(\d+)S", str(recipe))
    return f"{m.group(1)}S" if m else "?"


def recipe_from(algorithm, scenario_id) -> str:
    if _recipe_id is None:
        sys.exit("ERROR: BeyondBaseline (scenarios_simulation) is needed to map old "
                 "scenario ids to recipe ids. Activate the PhyloGAS environment.")
    return _recipe_id(int(scenario_id), str(algorithm))


# --------------------------------------------------------------------------
# Loaders -> long rows: run, state, replicate, recipe, metric, value
# --------------------------------------------------------------------------
def _bb_auc_rows(csv: Path, run, state, rep):
    df = pd.read_csv(csv)
    out = []
    for r in df.itertuples(index=False):
        out.append(dict(run=run, state=state, replicate=rep,
                        recipe=recipe_from(r.algorithm, r.scenario_id),
                        metric=_canon(r.eval_type), value=float(r.auc)))
    return out


def _bb_mugration_rows(csv: Path, run, state, rep):
    df = pd.read_csv(csv)
    out = []
    for r in df.itertuples(index=False):
        rec = recipe_from(r.algorithm, r.scenario_id)
        for m in MUGRATION:
            if hasattr(r, m):
                out.append(dict(run=run, state=state, replicate=rep, recipe=rec,
                                metric=m, value=float(getattr(r, m))))
    return out


def load_bb_replicates(path: Path, run: str, state: str) -> list:
    rows = []
    for rep_dir in sorted(path.glob("replicate_*")):
        rep = rep_dir.name.split("_", 1)[1]
        if (rep_dir / "AUC_rankings.csv").exists():
            rows += _bb_auc_rows(rep_dir / "AUC_rankings.csv", run, state, rep)
        if (rep_dir / "Mugration_Metrics.csv").exists():
            rows += _bb_mugration_rows(rep_dir / "Mugration_Metrics.csv", run, state, rep)
    return rows


def load_phylogas_project(path: Path, run: str, state: str, rep: str) -> list:
    rows = []
    bench = path / "05_benchmarks"
    f = bench / "AUC_truth_rankings.csv"
    if f.exists():
        df = pd.read_csv(f)
        for r in df.itertuples(index=False):
            rows.append(dict(run=run, state=state, replicate=rep, recipe=str(r.algorithm),
                             metric=_canon(r.eval_type), value=float(r.auc)))
    f = bench / "Mugration_Metrics.csv"
    if f.exists():
        df = pd.read_csv(f)
        for r in df.itertuples(index=False):
            for m in MUGRATION:
                if hasattr(r, m):
                    rows.append(dict(run=run, state=state, replicate=rep, recipe=str(r.label),
                                     metric=m, value=float(getattr(r, m))))
    # BeyondBaseline's own rankings (kl_targets), written next to the samples
    for f in (path / "03_sampled_datasets").glob("AUC_rankings.csv"):
        rows += _bb_auc_rows(f, run, state, rep)
    return rows


def load_run(spec: str) -> list:
    """LABEL=PATH[:STATE]; PATH may be a glob of several project directories."""
    if "=" not in spec:
        sys.exit(f"ERROR: --run takes LABEL=PATH[:STATE], got {spec!r}")
    label, rest = spec.split("=", 1)
    m = re.match(r"^(.*?)(?::([a-z]{2}))?$", rest)
    pattern, forced_state = m.group(1), m.group(2)
    paths = sorted(Path(p) for p in glob.glob(pattern)) or [Path(pattern)]
    rows = []
    for p in paths:
        if not p.is_dir():
            print(f"  WARNING: {p} is not a directory; skipped", file=sys.stderr)
            continue
        # project names look like va_delta_wave or va_delta_wave_ci0.25_r07
        pm = re.match(r"^([a-z]{2})_.*?(?:_r(\d+))?$", p.name)
        state = forced_state or (pm.group(1) if pm else None)
        if (p / "05_benchmarks").is_dir():
            rep = (pm.group(2) if pm and pm.group(2) else "0")
            got = load_phylogas_project(p, label, state or "?", rep)
            kind = "PhyloGAS project"
        elif any(p.glob("replicate_*/AUC_rankings.csv")) or any(p.glob("replicate_*/Mugration_Metrics.csv")):
            if not state:
                sys.exit(f"ERROR: {p} is a BeyondBaseline replicate directory; give its state "
                         f"as {label}={p}:va")
            got = load_bb_replicates(p, label, state)
            kind = "BeyondBaseline replicates"
        else:
            print(f"  WARNING: {p}: no recognised benchmark files; skipped", file=sys.stderr)
            continue
        reps = len({r["replicate"] for r in got})
        print(f"  {label:10} {state or '?':3} {kind:26} {len(got):6,} values, "
              f"{reps} replicate(s)  {p}")
        rows += got
    return rows


# --------------------------------------------------------------------------
# Analyses
# --------------------------------------------------------------------------
def summarise(long: pd.DataFrame) -> pd.DataFrame:
    """Median over replicates (and over duplicate recipe rows, e.g. old SURS
    scored once per scenario sharing its stride)."""
    agg = (long.groupby(["run", "state", "recipe", "metric"], as_index=False)["value"]
               .agg(value="median", n="size"))
    agg["stride"] = agg["recipe"].map(stride_of)
    agg["higher_is_better"] = agg["metric"].map(higher_is_better)
    # rank 1 = best, within run/state/metric/stride
    signed = np.where(agg["higher_is_better"], -agg["value"], agg["value"])
    agg["rank"] = (pd.Series(signed, index=agg.index)
                     .groupby([agg["run"], agg["state"], agg["metric"], agg["stride"]])
                     .rank(method="average"))
    return agg


def _spearman(a: pd.Series, b: pd.Series) -> float:
    ok = a.notna() & b.notna()
    if ok.sum() < 3:
        return float("nan")
    ra, rb = a[ok].rank(), b[ok].rank()
    if ra.std() == 0 or rb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def head_to_head(agg: pd.DataFrame, focus: list) -> pd.DataFrame:
    sub = agg[agg["recipe"].isin(focus)]
    wide = sub.pivot_table(index=["run", "state", "metric"], columns="recipe",
                           values="value").reset_index()
    if all(f in wide.columns for f in focus[:2]):
        a, b = focus[0], focus[1]
        hib = wide["metric"].map(higher_is_better)
        better_b = np.where(hib, wide[b] > wide[a], wide[b] < wide[a])
        wide["winner"] = np.where(wide[a].isna() | wide[b].isna(), "",
                                  np.where(better_b, b, a))
    return wide


def _sign_test_p(k: int, n: int) -> float:
    """Exact two-sided sign test: k successes out of n non-tied pairs."""
    from math import comb
    if n == 0:
        return float("nan")
    tail = sum(comb(n, i) for i in range(0, min(k, n - k) + 1)) / 2 ** n
    return min(1.0, 2 * tail)


def paired(long: pd.DataFrame, focus: list) -> pd.DataFrame:
    """Second focus recipe against the first, one comparison per replicate."""
    a, b = focus[0], focus[1]
    rows = []
    for (run, state, metric), g in long.groupby(["run", "state", "metric"]):
        w = g.groupby(["replicate", "recipe"])["value"].mean().unstack()
        if a not in w.columns or b not in w.columns:
            continue
        d = (w[b] - w[a]).dropna()
        if len(d) < 2:
            continue
        better = d > 0 if higher_is_better(metric) else d < 0
        n = int((d != 0).sum())
        k = int((better & (d != 0)).sum())
        rows.append(dict(run=run, state=state, metric=metric, baseline=a, challenger=b,
                         replicates=len(d), challenger_wins=k, non_tied=n,
                         median_diff=float(d.median()), sign_test_p=_sign_test_p(k, n)))
    return pd.DataFrame(rows)


def rank_agreement(agg: pd.DataFrame) -> pd.DataFrame:
    runs = list(dict.fromkeys(agg["run"]))
    if len(runs) < 2:
        return pd.DataFrame()
    base, rows = runs[0], []
    for other in runs[1:]:
        for (state, metric, stride), g in agg.groupby(["state", "metric", "stride"]):
            x = g[g["run"] == base].set_index("recipe")["rank"]
            y = g[g["run"] == other].set_index("recipe")["rank"]
            common = x.index.intersection(y.index)
            if len(common) < 3:
                continue
            rows.append(dict(baseline=base, run=other, state=state, metric=metric,
                             stride=stride, n_recipes=len(common),
                             spearman=_spearman(x[common], y[common]),
                             best_baseline=x[common].idxmin(), best_run=y[common].idxmin()))
    return pd.DataFrame(rows)


def tradeoff(agg: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (run, state, stride), g in agg.groupby(["run", "state", "stride"]):
        cos = g[g["metric"] == "cosine_similarity"].set_index("recipe")["value"]
        if len(cos) < 3:
            continue
        for metric, gm in g[g["metric"].str.startswith(COVERAGE_PREFIXES)].groupby("metric"):
            cov = gm.set_index("recipe")["value"]
            common = cos.index.intersection(cov.index)
            rows.append(dict(run=run, state=state, stride=stride, coverage_metric=metric,
                             n_recipes=len(common),
                             spearman_vs_cosine=_spearman(cov[common], cos[common])))
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split("\n\n", 1)[1])
    ap.add_argument("--run", action="append", required=True,
                    help="LABEL=PATH[:STATE]; repeat. The first is the baseline.")
    ap.add_argument("--focus", nargs=2, default=["1S__surs", "1S-P__lasso_greedy"],
                    help="Two recipes for the head-to-head (default: SURS vs "
                         "LASSO-Greedy-50 with a population target, both weekly)")
    ap.add_argument("--out", required=True, help="Output directory")
    args = ap.parse_args(argv)

    print("Runs:")
    rows = []
    for spec in args.run:
        rows += load_run(spec)
    if not rows:
        sys.exit("ERROR: nothing loaded.")
    long = pd.DataFrame(rows)
    agg = summarise(long)
    h2h = head_to_head(agg, args.focus)
    agree = rank_agreement(agg)
    trade = tradeoff(agg)
    pair = paired(long, args.focus)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    long.to_csv(out / "long.csv", index=False)
    agg.to_csv(out / "summary.csv", index=False)
    h2h.to_csv(out / "head_to_head.csv", index=False)
    agree.to_csv(out / "rank_agreement.csv", index=False)
    trade.to_csv(out / "tradeoff.csv", index=False)
    pair.to_csv(out / "paired.csv", index=False)

    pd.set_option("display.width", 200)
    a, b = args.focus
    print(f"\nHead to head, {a} vs {b} (median over replicates):")
    show = ["cosine_similarity", "topological_f1", "coverage_size_0", "coverage_size_100",
            "8_week_rolling_tree_coverage", "stride_component_coverage", "kl_targets",
            "cumulative_infections"]
    if {a, b} <= set(h2h.columns):
        print(h2h[h2h["metric"].isin(show)].round(4).to_string(index=False))
    if not pair.empty:
        print(f"\n{b} vs {a}, replicate by replicate (exact sign test):")
        print(pair[pair["metric"].isin(show)][["run", "state", "metric", "challenger_wins", "non_tied",
              "median_diff", "sign_test_p"]].round(4).to_string(index=False))
    if not agree.empty:
        print("\nRank agreement with the baseline (Spearman across recipes, within stride):")
        print(agree.pivot_table(index=["metric"], columns=["run", "state", "stride"],
                                values="spearman").round(2).to_string())
    if not trade.empty:
        print("\nCoverage vs cosine across recipes (negative = the abstract's trade-off):")
        print(trade.pivot_table(index="coverage_metric", columns=["run", "state", "stride"],
                                values="spearman_vs_cosine").round(2).to_string())
    print("\nNote: coverage metrics changed definition between the CSTE runs and now "
          "(true tree, all infections, wave-start calendar), so only orderings and "
          "signs are comparable, not values. Mugration metrics did not change.")
    print(f"\nWrote {out}/{{long,summary,head_to_head,paired,rank_agreement,tradeoff}}.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
