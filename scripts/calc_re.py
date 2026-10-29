#!/usr/bin/env python3
"""Calculate overall, early-wave R_0, and time-varying cohort R_e(t) from EpiHiper output.

Ported from TwinSampler's ``scripts/linelist_generation/calc_re.py`` and
``label_components.py`` (``create_graph_classic`` / ``calculate_Re``), with:
  1. Tuple-unpacking fix and streaming out-degree counting (avoids building a
     10M-node NetworkX DiGraph in memory; works even without pandas/matplotlib).
  2. Early-wave R_0 summary metrics (first 1k / 10k / 50k cases and peak 7-tick
     rolling cohort R_e) alongside the full-wave R_e(t) CSV and PNG plot
     (note: overall R_e over a completed wave is tautologically ~1.0 because
     every non-seed infection has exactly 1 parent).
  3. Optional ``--compare`` mode to compare exp1 vs. exp2 across all 5 states.

Examples
--------
Single file (drop-in replacement for TwinSampler's calc_re.py):
    python scripts/calc_re.py --epihiper data/example_data/0.25_crossinfection/output.csv.gz \
        --variant 1 --out results/ga_e1_re.csv

Side-by-side exp1 vs. exp2 comparison on Rivanna:
    python scripts/calc_re.py --compare --variant 1 --cell 0.25 --out-dir results/re_exp1_vs_exp2
"""
from __future__ import annotations

import argparse
import collections
import csv
import gzip
from pathlib import Path

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


DEFAULT_EXP1_ROOT = (
    "/project/bii_nssac/epihiper-simulations/pipeline-jc/run/"
    "20250120_1/output_root/proj/20250120_1/batch_1"
)
DEFAULT_EXP2_ROOT = (
    "/project/bii_nssac/epihiper-simulations/pipeline-jc/run/"
    "20261022_2/output_root/proj/20261022_2/batch_1"
)
DEFAULT_EXP3_ROOT = (
    "/project/bii_nssac/epihiper-simulations/pipeline-jc/run/"
    "20261026_1/output_root/proj/20261026_1/batch_1"
)
STATES = ("ga", "ma", "mn", "va", "wa")
FIRST_V2_SEED = {"ga": 130, "ma": 103, "mn": 124, "va": 128, "wa": 124}


def compute_cohort_re(
    epihiper_path: str | Path,
    variant: str | None = None,
    start_tick: int | None = None,
    stop_tick: int | None = None,
) -> tuple[float, list[dict], dict]:
    """Trace contact_pid -> pid across Exposure ('E*') events and compute cohort R_e(t).

    Matches ``create_graph_classic`` + ``calculate_Re`` in TwinSampler's
    ``label_components.py``: each exposure event gets ``alias_pid = f"{pid}.{tick}"``,
    resolves ``contact_pid`` to the infector's most recent ``alias_pid``, and
    counts the out-degree (secondary infections) of every ``alias_pid``.
    """
    prefix = f"E{variant}_" if variant is not None else "E"
    open_fn = gzip.open if str(epihiper_path).endswith(".gz") else open

    events_tick: list[int] = []
    out_degree: collections.Counter[int] = collections.Counter()
    alias_map: dict[str, int] = {}
    n_seeds = 0

    with open_fn(epihiper_path, "rt") as f:
        reader = csv.DictReader(f)
        for row in reader:
            st = row["exit_state"]
            if not st.startswith(prefix):
                continue
            t = int(row["tick"])
            if start_tick is not None and t < start_tick:
                continue
            if stop_tick is not None and t > stop_tick:
                continue

            pid = row["pid"].strip()
            cpid = row["contact_pid"].strip()
            idx = len(events_tick)
            events_tick.append(t)

            if cpid not in ("", "-1", "0") and cpid in alias_map:
                parent_idx = alias_map[cpid]
                out_degree[parent_idx] += 1
            else:
                n_seeds += 1

            alias_map[pid] = idx

    n_total = len(events_tick)
    if n_total == 0:
        return 0.0, [], {}

    ks = [out_degree[i] for i in range(n_total)]
    overall_re = sum(ks) / n_total

    tick_sum: collections.Counter[int] = collections.Counter()
    tick_cnt: collections.Counter[int] = collections.Counter()
    for t, k in zip(events_tick, ks):
        tick_sum[t] += k
        tick_cnt[t] += 1

    ticks_sorted = sorted(tick_cnt)
    re_rows: list[dict] = []
    window_vals: collections.deque[float] = collections.deque(maxlen=7)
    for t in ticks_sorted:
        re_val = tick_sum[t] / tick_cnt[t]
        window_vals.append(re_val)
        roll7 = sum(window_vals) / len(window_vals)
        re_rows.append(
            {
                "tick": t,
                "R_e": re_val,
                "R_e_rolling7": roll7,
                "cases": tick_cnt[t],
            }
        )

    # Early-wave R_0 metrics (skipping the first 50 sparse importations)
    r0_1k = sum(ks[50:1050]) / 1000 if n_total >= 1050 else float("nan")
    r0_10k = sum(ks[50:10050]) / 10000 if n_total >= 10050 else float("nan")
    re_50k = sum(ks[50:50050]) / 50000 if n_total >= 50050 else float("nan")

    substantial = [r["R_e_rolling7"] for r in re_rows if r["cases"] >= 50]
    peak_roll7 = max(substantial) if substantial else float("nan")
    peak_row = max(re_rows, key=lambda r: r["cases"])

    summary = {
        "n_total": n_total,
        "n_seeds": n_seeds,
        "overall_Re": overall_re,
        "R0_1k": r0_1k,
        "R0_10k": r0_10k,
        "Re_50k": re_50k,
        "peak_rolling7_Re": peak_roll7,
        "peak_cases_tick": peak_row["tick"],
        "peak_daily_cases": peak_row["cases"],
    }
    return overall_re, re_rows, summary


def write_re_csv(re_rows: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["tick", "R_e", "R_e_rolling7", "cases"])
        writer.writeheader()
        writer.writerows(re_rows)


def plot_re(re_rows: list[dict], plot_file: Path, title: str = "Time-Varying Cohort Reproduction Number ($R_c$) from Simulation") -> bool:
    if plt is None or not re_rows:
        return False
    ticks = [r["tick"] for r in re_rows]
    re_vals = [r["R_e"] for r in re_rows]
    roll7 = [r["R_e_rolling7"] for r in re_rows]

    plt.figure(figsize=(12, 6))
    plt.plot(ticks, re_vals, marker="o", linestyle="-", color="dodgerblue", alpha=0.4, markersize=4, label="Daily Cohort $R_e$")
    if len(re_rows) > 7:
        plt.plot(ticks, roll7, color="firebrick", linewidth=2.5, label="7-Tick Rolling Average")
    plt.axhline(y=1.0, color="black", linestyle="--", linewidth=1.5, label="Epidemic Threshold ($R_e=1$)")
    plt.title(title, fontsize=14)
    plt.xlabel("Simulation Tick (Time of Exposure)", fontsize=12)
    plt.ylabel("Effective Reproduction Number ($R_e$)", fontsize=12)
    plt.grid(True, linestyle=":", alpha=0.7)
    plt.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)
    plt.close()
    return True


def run_check_gi(args: argparse.Namespace) -> None:
    """Measure realized GI (E1 & E2), reinfection rates, and inter-infection gaps."""
    exp_map = {"exp1": args.exp1_root, "exp2": args.exp2_root, "exp3": args.exp3_root}
    root = exp_map.get(args.check_gi, args.check_gi)

    for st in args.states:
        path = Path(root) / args.cell / st / f"replicate_{args.replicate}" / "output.csv.gz"
        if not path.exists():
            print(f"{st.upper()}: MISSING {path}")
            continue
        t0, t1 = FIRST_V2_SEED.get(st, 125), FIRST_V2_SEED.get(st, 125) + 300
        last_exp: dict[str, tuple[int, str]] = {}
        gi_sum = {"E1": 0, "E2": 0, "E2_win": 0}
        gi_cnt = {"E1": 0, "E2": 0, "E2_win": 0}
        e2_ticks = collections.defaultdict(list)
        e2_win_cnt: collections.Counter[str] = collections.Counter()
        max_tick = 0
        e1_total = e2_total = 0

        with gzip.open(path, "rt") as f:
            reader = csv.DictReader(f)
            for row in reader:
                t = int(row["tick"])
                if t > max_tick:
                    max_tick = t
                es = row["exit_state"]
                if not (es.startswith("E1_") or es.startswith("E2_")):
                    continue
                var = "E1" if es.startswith("E1_") else "E2"
                pid = row["pid"].strip()
                cpid = row["contact_pid"].strip()

                if var == "E1":
                    e1_total += 1
                else:
                    e2_total += 1
                    e2_ticks[pid].append(t)
                    if t0 <= t <= t1:
                        e2_win_cnt[pid] += 1

                if cpid in last_exp:
                    ctick, cvar = last_exp[cpid]
                    dt = t - ctick
                    if cvar == var and dt >= 0:
                        gi_sum[var] += dt
                        gi_cnt[var] += 1
                        if var == "E2" and t0 <= t <= t1:
                            gi_sum["E2_win"] += dt
                            gi_cnt["E2_win"] += 1
                last_exp[pid] = (t, var)

        gaps = [ts[i + 1] - ts[i] for ts in e2_ticks.values() if len(ts) > 1 for i in range(len(ts) - 1)]
        gaps.sort()
        n_u = len(e2_ticks)
        n_r = sum(1 for ts in e2_ticks.values() if len(ts) > 1)
        n_wu = len(e2_win_cnt)
        n_wr = sum(1 for c in e2_win_cnt.values() if c > 1)
        max_inf = max((len(ts) for ts in e2_ticks.values()), default=0)
        g30 = sum(1 for g in gaps if g <= 30) / len(gaps) * 100 if gaps else 0.0
        g60 = sum(1 for g in gaps if g <= 60) / len(gaps) * 100 if gaps else 0.0

        print(
            f"{st.upper()} (cell {args.cell}) | max_tick={max_tick} | "
            f"E1={e1_total:,} (GI={gi_sum['E1']/max(1,gi_cnt['E1']):.2f}d) | "
            f"E2={e2_total:,} (GI={gi_sum['E2']/max(1,gi_cnt['E2']):.2f}d, win_edges={gi_cnt['E2_win']:,}) | "
            f"E2 unique={n_u:,}, reinf_all={n_r/max(1,n_u)*100:.1f}%, reinf_win={n_wr/max(1,n_wu)*100:.1f}%, max_inf={max_inf} | "
            f"gaps<=30d={g30:.2f}%, <=60d={g60:.2f}%"
        )


def run_comparison(args: argparse.Namespace) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    var_label = f"E{args.variant}" if args.variant else "All"

    print(
        f"{'State':<5} {'Exp':<4} | {'Total ' + var_label:>11} | "
        f"{'R0(1k)':>7} {'R0(10k)':>8} {'Re(50k)':>8} {'Peak7dRe':>9} | "
        f"{'PeakTick':>8} {'Peak/day':>9}"
    )
    print("-" * 86)

    exp_roots = [("exp1", args.exp1_root), ("exp2", args.exp2_root), ("exp3", args.exp3_root)]
    exp_colors = {"exp1": "#E66101", "exp2": "#008891", "exp3": "#5E3C99"}

    rows = []
    for st in args.states:
        series_map = {}
        for exp_name, root in exp_roots:
            path = Path(root) / args.cell / st / f"replicate_{args.replicate}" / "output.csv.gz"
            if not path.exists():
                continue
            _, re_rows, summ = compute_cohort_re(path, variant=args.variant, start_tick=args.start_tick, stop_tick=args.stop_tick)
            re_csv = out_dir / f"{st}_{exp_name}_v{args.variant or 'all'}_Re.csv"
            write_re_csv(re_rows, re_csv)
            series_map[exp_name] = re_rows

            print(
                f"{st.upper():<5} {exp_name:<4} | {summ['n_total']:>11,d} | "
                f"{summ['R0_1k']:>7.2f} {summ['R0_10k']:>8.2f} {summ['Re_50k']:>8.2f} {summ['peak_rolling7_Re']:>9.2f} | "
                f"tick {summ['peak_cases_tick']:>3d} {summ['peak_daily_cases']:>9,d}"
            )
            rows.append({"state": st, "experiment": exp_name, **summ})
        print("-" * 86)

        if len(series_map) >= 2 and plt is not None:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
            for exp_name, rlist in series_map.items():
                color = exp_colors.get(exp_name, "black")
                sub = [r for r in rlist if r["cases"] >= 20]
                ax1.plot([r["tick"] for r in sub], [r["R_e_rolling7"] for r in sub], color=color, linewidth=2.2, label=f"{exp_name} (7-tick rolling $R_e$, cases>=20)")
                ax2.plot([r["tick"] for r in rlist], [r["cases"] for r in rlist], color=color, linewidth=2.0, label=f"{exp_name} daily exposures")
            ax1.axhline(1.0, color="black", linestyle="--", linewidth=1.2, label="$R_e = 1.0$")
            ax1.axhline(4.0, color="gray", linestyle=":", linewidth=1.2, label="Target $R_0 = 4.0$")
            ax1.set_ylabel(f"Cohort $R_e(t)$ ({var_label})")
            ax1.set_title(f"{st.upper()} (cell {args.cell}) — {var_label} Cohort $R_e(t)$ & Daily Incidence: {' vs. '.join(series_map.keys())}")
            ax1.legend()
            ax1.grid(True, linestyle=":", alpha=0.7)
            ax2.set_ylabel("Daily Exposures")
            ax2.set_xlabel("Simulation Tick")
            ax2.legend()
            ax2.grid(True, linestyle=":", alpha=0.7)
            fig.tight_layout()
            fig.savefig(out_dir / f"{st}_v{args.variant or 'all'}_Re_comparison.png", dpi=250)
            plt.close(fig)

    if rows:
        sum_csv = out_dir / f"summary_v{args.variant or 'all'}_Re.csv"
        with open(sum_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"Saved summary table (and plots if matplotlib available) to: {out_dir}")


def run_plot_grid(args: argparse.Namespace) -> None:
    """Plot 5-state x 4-cell daily incidence curves (E1 dashed, E2 solid) across experiments."""
    if plt is None:
        print("Error: matplotlib is required for --plot-grid")
        return
    cells = ("0.25", "0.5", "0.75", "1.0")
    exp_roots = [("exp1", args.exp1_root, "#E66101"), ("exp2", args.exp2_root, "#008891"), ("exp3", args.exp3_root, "#5E3C99")]
    fig, axes = plt.subplots(len(args.states), len(cells), figsize=(20, 3.2 * len(args.states)), sharex=True)
    if len(args.states) == 1:
        axes = [axes]

    for i, st in enumerate(args.states):
        for j, cell in enumerate(cells):
            ax = axes[i][j]
            for exp_name, root, color in exp_roots:
                path = Path(root) / cell / st / f"replicate_{args.replicate}" / "output.csv.gz"
                if not path.exists():
                    continue
                e1_cnt: collections.Counter[int] = collections.Counter()
                e2_cnt: collections.Counter[int] = collections.Counter()
                with gzip.open(path, "rt") as f:
                    for row in csv.DictReader(f):
                        es = row["exit_state"]
                        if es.startswith("E1_"):
                            e1_cnt[int(row["tick"])] += 1
                        elif es.startswith("E2_"):
                            e2_cnt[int(row["tick"])] += 1
                t1 = sorted(e1_cnt)
                t2 = sorted(e2_cnt)
                if t1:
                    ax.plot(t1, [e1_cnt[t] for t in t1], linestyle="--", color=color, alpha=0.85, linewidth=1.4, label=f"{exp_name} E1 (Alpha)")
                if t2:
                    ax.plot(t2, [e2_cnt[t] for t in t2], linestyle="-", color=color, alpha=0.9, linewidth=1.6, label=f"{exp_name} E2 (Delta)")
            ax.set_xlim(0, 500)
            ax.set_title(f"{st.upper()} | cell {cell}", fontsize=10)
            ax.grid(True, linestyle=":", alpha=0.5)
            if j == 0:
                ax.set_ylabel("Exposures / day")
            if i == len(args.states) - 1:
                ax.set_xlabel("Tick")
            if i == 0 and j == len(cells) - 1:
                ax.legend(fontsize=7, loc="upper right")

    fig.suptitle("EpiHiper Daily Exposures: exp1 vs. exp2 vs. exp3 (replicate_0, ticks 0-500)", fontsize=14, y=0.995)
    fig.tight_layout()
    out_file = Path(args.plot_grid)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=250)
    plt.close(fig)
    print(f"Saved 5-state x 4-cell daily incidence grid to: {out_file}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate overall, early-wave R_0, and time-varying cohort R_e(t) from EpiHiper output.")
    parser.add_argument("--epihiper", help="Path to raw EpiHiper output.csv.gz (single-file mode)")
    parser.add_argument("--out", help="Output path for the Re_over_time.csv file (single-file mode)")
    parser.add_argument("--variant", type=str, default=None, help="Specific variant number to filter for (e.g., '1' for E1, '2' for E2)")
    parser.add_argument("--start_tick", type=int, default=None, help="Only consider infections occurring on or after this tick")
    parser.add_argument("--stop_tick", type=int, default=None, help="Only consider infections occurring on or before this tick")

    # Side-by-side comparison & GI diagnostic flags
    parser.add_argument("--compare", action="store_true", help="Compare exp1 vs. exp2 vs. exp3 across states")
    parser.add_argument("--check-gi", metavar="EXP", help="Run GI & reinfection-gap diagnostic for 'exp1', 'exp2', 'exp3', or a batch_1 root path")
    parser.add_argument("--plot-grid", metavar="OUT_PNG", help="Plot 5-state x 4-cell daily incidence curves (exp1 vs. exp2 vs. exp3)")
    parser.add_argument("--exp1-root", default=DEFAULT_EXP1_ROOT, help="Root path for exp1 batch_1 directory")
    parser.add_argument("--exp2-root", default=DEFAULT_EXP2_ROOT, help="Root path for exp2 batch_1 directory")
    parser.add_argument("--exp3-root", default=DEFAULT_EXP3_ROOT, help="Root path for exp3 batch_1 directory")
    parser.add_argument("--cell", default="0.25", help="Cross-infection cell folder (default: 0.25)")
    parser.add_argument("--replicate", type=int, default=0, help="Replicate number (default: 0)")
    parser.add_argument("--states", nargs="+", default=list(STATES), help="States to compare in --compare / --check-gi mode")
    parser.add_argument("--out-dir", default="results/re_comparison", help="Output directory for --compare mode")

    args = parser.parse_args()

    if args.check_gi:
        run_check_gi(args)
        return

    if args.plot_grid:
        run_plot_grid(args)
        return

    if args.compare:
        run_comparison(args)
        return

    if not args.epihiper or not args.out:
        parser.error("Provide --epihiper and --out for single-file mode, or use --compare.")

    print(f"Loading EpiHiper data from {args.epihiper}...")
    overall_re, re_rows, summ = compute_cohort_re(
        args.epihiper,
        variant=args.variant,
        start_tick=args.start_tick,
        stop_tick=args.stop_tick,
    )

    print("\n====================================")
    print(f"Total Exposures        : {summ.get('n_total', 0):,}")
    print(f"Importation Seeds      : {summ.get('n_seeds', 0):,}")
    print(f"Early R_0 (first 1k)   : {summ.get('R0_1k', float('nan')):.3f}")
    print(f"Early R_0 (first 10k)  : {summ.get('R0_10k', float('nan')):.3f}")
    print(f"Early R_e (first 50k)  : {summ.get('Re_50k', float('nan')):.3f}")
    print(f"Peak 7-Tick Rolling R_e: {summ.get('peak_rolling7_Re', float('nan')):.3f} (cases>=50/d)")
    print(f"OVERALL EPIDEMIC R_e   : {overall_re:.3f}")
    print("====================================\n")

    out_path = Path(args.out)
    write_re_csv(re_rows, out_path)
    print(f"Time-varying R_e successfully saved to: {out_path}")

    plot_file = out_path.with_suffix(".png")
    if plot_re(re_rows, plot_file):
        print(f"Companion plot saved to: {plot_file}")


if __name__ == "__main__":
    main()
