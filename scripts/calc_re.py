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
STATES = ("ga", "ma", "mn", "va", "wa")


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

    rows = []
    for st in args.states:
        series_pair = {}
        for exp_name, root in (("exp1", args.exp1_root), ("exp2", args.exp2_root)):
            path = Path(root) / args.cell / st / f"replicate_{args.replicate}" / "output.csv.gz"
            if not path.exists():
                print(f"{st.upper():<5} {exp_name:<4} | MISSING: {path}")
                continue
            _, re_rows, summ = compute_cohort_re(path, variant=args.variant, start_tick=args.start_tick, stop_tick=args.stop_tick)
            re_csv = out_dir / f"{st}_{exp_name}_v{args.variant or 'all'}_Re.csv"
            write_re_csv(re_rows, re_csv)
            series_pair[exp_name] = re_rows

            print(
                f"{st.upper():<5} {exp_name:<4} | {summ['n_total']:>11,d} | "
                f"{summ['R0_1k']:>7.2f} {summ['R0_10k']:>8.2f} {summ['Re_50k']:>8.2f} {summ['peak_rolling7_Re']:>9.2f} | "
                f"tick {summ['peak_cases_tick']:>3d} {summ['peak_daily_cases']:>9,d}"
            )
            rows.append({"state": st, "experiment": exp_name, **summ})
        print("-" * 86)

        if len(series_pair) == 2 and plt is not None:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
            for exp_name, color in (("exp1", "#E66101"), ("exp2", "#008891")):
                rlist = series_pair[exp_name]
                sub = [r for r in rlist if r["cases"] >= 20]
                ax1.plot([r["tick"] for r in sub], [r["R_e_rolling7"] for r in sub], color=color, linewidth=2.2, label=f"{exp_name} (7-tick rolling $R_e$, cases>=20)")
                ax2.plot([r["tick"] for r in rlist], [r["cases"] for r in rlist], color=color, linewidth=2.0, label=f"{exp_name} daily exposures")
            ax1.axhline(1.0, color="black", linestyle="--", linewidth=1.2, label="$R_e = 1.0$")
            ax1.axhline(4.0, color="gray", linestyle=":", linewidth=1.2, label="Target $R_0 = 4.0$")
            ax1.set_ylabel(f"Cohort $R_e(t)$ ({var_label})")
            ax1.set_title(f"{st.upper()} (cell {args.cell}) — {var_label} Cohort $R_e(t)$ & Daily Incidence: exp1 vs. exp2")
            ax1.legend()
            ax1.grid(True, linestyle=":", alpha=0.7)
            ax2.set_ylabel("Daily Exposures")
            ax2.set_xlabel("Simulation Tick")
            ax2.legend()
            ax2.grid(True, linestyle=":", alpha=0.7)
            fig.tight_layout()
            fig.savefig(out_dir / f"{st}_v{args.variant or 'all'}_exp1_vs_exp2_Re.png", dpi=250)
            plt.close(fig)

    if rows:
        sum_csv = out_dir / f"summary_v{args.variant or 'all'}_Re.csv"
        with open(sum_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"Saved summary table (and plots if matplotlib available) to: {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate overall, early-wave R_0, and time-varying cohort R_e(t) from EpiHiper output.")
    parser.add_argument("--epihiper", help="Path to raw EpiHiper output.csv.gz (single-file mode)")
    parser.add_argument("--out", help="Output path for the Re_over_time.csv file (single-file mode)")
    parser.add_argument("--variant", type=str, default=None, help="Specific variant number to filter for (e.g., '1' for E1, '2' for E2)")
    parser.add_argument("--start_tick", type=int, default=None, help="Only consider infections occurring on or after this tick")
    parser.add_argument("--stop_tick", type=int, default=None, help="Only consider infections occurring on or before this tick")

    # Side-by-side comparison mode flags
    parser.add_argument("--compare", action="store_true", help="Compare exp1 vs. exp2 across states")
    parser.add_argument("--exp1-root", default=DEFAULT_EXP1_ROOT, help="Root path for exp1 batch_1 directory")
    parser.add_argument("--exp2-root", default=DEFAULT_EXP2_ROOT, help="Root path for exp2 batch_1 directory")
    parser.add_argument("--cell", default="0.25", help="Cross-infection cell folder (default: 0.25)")
    parser.add_argument("--replicate", type=int, default=0, help="Replicate number (default: 0)")
    parser.add_argument("--states", nargs="+", default=list(STATES), help="States to compare in --compare mode")
    parser.add_argument("--out-dir", default="results/re_exp1_vs_exp2", help="Output directory for --compare mode")

    args = parser.parse_args()

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
