#!/usr/bin/env python3
"""Generate exp3 epidemiological, molecular-clock, and phylodynamic characterization tables & figures.

Produces publication-ready CSV tables, GitHub-flavored Markdown tables, native
vector SVGs, and (when Google Chrome is available) high-DPI PNG rasters in
``docs/characterization/``. Uses only the Python standard library so it runs
identically on local macOS workstations and on Rivanna.

Tables generated:
  1. exp3_epi_characterization.{csv,md}
     5-state epidemiological characterization (E1 Alpha & E2 Delta total cases,
     generation intervals, Gamma(4, 45.75) waning immunity repeat-infection
     gaps, cohort R0(1k)/R0(10k)/Re(50k)/Peak7dRe, peak dates, pre-peak seeds).
  2. exp3_molecular_clock_closure.{csv,md}
     Molecular-clock decomposition across all 5 states (exp1 baseline vs.
     pre-run exp3 predictions vs. measured exp3 results from
     sequence_divergence.csv and clock_estimates.csv; baseline TreeTime clock
     reported on 1S__surs separately from 1S-P__lasso_greedy).
  3. exp3_mugration_benchmarks.{csv,md}
     Inter-county phylogeographic reconstruction metrics (cosine similarity,
     topological F1, Pearson r, masked MAE) across all 8 sampling strategies.
  4. exp3_auc_truth_benchmarks.{csv,md}
     Tree-free and phylogenetic truth-coverage AUC/MRD metrics across all 8
     sampling strategies.

Figures generated:
  1. fig1_exp3_epi_dynamics.{svg,png}
     3-panel epidemiological characterization: (A) Generation Interval & waning
     reinfection gap calibration (exp1 vs. exp2 vs. exp3 vs. Hart et al. 2022),
     (B) Cohort R0/Re progression across the wave, (C) Cumulative empirical E2
     importation schedules vs. wave peak timing (showing elimination of exp2's
     June founder bottleneck).
  2. fig2_exp3_clock_closure.{svg,png}
     3-panel molecular-clock closure: (A) Two-way orthogonality in (GI, lambda)
     space across exp1(poisson) -> exp1(lognormal) -> exp3(lognormal),
     (B) State-by-state clock rate hierarchy (lambda/GI, mu_lineage, mu_truth,
     and 1S__surs TreeTime mu_sim vs. real SARS-CoV-2 8.0e-4),
     (C) Root-to-tip R^2 temporal signal jump from exp1 to exp3.
  3. fig3_exp3_surveillance_benchmarks.{svg,png}
     2-panel surveillance benchmark summary across 5 states and 8 strategies:
     (A) Inter-county flow Cosine Similarity & Topological F1,
     (B) 8-week rolling tree coverage MRD & Senior 65+ demographic equity MRD.
  4. fig4_exp3_state_timelines.{svg,png}
     5-state timeline panels (GA, MA, MN, VA, WA; cell 0.25 only, no empty
     subplot columns) showing empirical E1 & E2 importation schedules,
     vaccination rollout, 300-day painted window, wave peak alignment (exp1 vs.
     exp2 vs. exp3), and (when results/re_comparison/incidence_cache/ is
     present) daily E1 & E2 exposure curves.

Usage:
    python3 scripts/utils/characterize_exp3.py \
        --results-dir results \
        --exp1-archive results/_archive/2026-10-16_exp1_lognormal \
        --out-dir docs/characterization
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
import os
import shutil
import statistics as stats
import subprocess
from pathlib import Path

STATES = ["ga", "ma", "mn", "va", "wa"]
STATE_NAMES = {
    "ga": "Georgia (GA)",
    "ma": "Massachusetts (MA)",
    "mn": "Minnesota (MN)",
    "va": "Virginia (VA)",
    "wa": "Washington (WA)",
}
TICK_ZERO = dt.date(2020, 11, 30)
GENOME_LEN = 29750.0
DAYS_PER_YEAR = 365.25
REAL_SC2_CLOCK = 8.0e-4  # subs/site/yr

# Accessible palette matching scripts/poster_panels.py
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
YELLOW = "#eda100"
PURPLE = "#5e3c99"
RED = "#c93b2b"
GRAY = "#7f7e7a"
INK = "#0b0b0b"
INK2 = "#4a4946"
INK3 = "#7d7c78"
GRID = "#e5e4e0"
BG_PANEL = "#faf9f6"
FONT = "Arial, Helvetica, sans-serif"

# Validated full-ABM epidemiological metrics from scripts/calc_re.py --compare
# on Rivanna (cell 0.25, replicate 0) recorded in design/exp3.md Section 4.
ABM_EPI_METRICS = {
    "ga": {
        "e1_total": 4_408_224, "e1_gi": 5.30, "e1_r0_1k": 3.41, "e1_r0_10k": 2.01,
        "e1_re_50k": 2.22, "e1_peak_re": 5.23, "e1_peak_tick": 94, "e1_peak_cases": 109_739,
        "e2_total": 5_342_451, "e2_gi": 4.74, "e2_r0_1k": 2.18, "e2_r0_10k": 1.71,
        "e2_re_50k": 1.49, "e2_peak_re": 2.31, "e2_peak_tick": 239, "e2_peak_cases": 44_795,
        "e2_win_edges": 3_917_729, "e2_unique_hosts": 4_156_497,
        "reinf_all_pct": 27.6, "reinf_win_pct": 16.0, "max_inf": 4,
        "gaps_le30_pct": 0.03, "gaps_le60_pct": 0.44,
        "exp1_e2_peak_tick": 254, "exp2_e2_peak_tick": 207,
        "exp1_e1_peak_tick": 125, "exp2_e1_peak_tick": 59,
        "exp2_e1_gi": 4.96, "exp2_e2_gi": 4.56,
    },
    "ma": {
        "e1_total": 3_254_688, "e1_gi": 5.23, "e1_r0_1k": 2.94, "e1_r0_10k": 2.39,
        "e1_re_50k": 2.57, "e1_peak_re": 3.87, "e1_peak_tick": 104, "e1_peak_cases": 116_827,
        "e2_total": 3_615_712, "e2_gi": 4.71, "e2_r0_1k": 1.92, "e2_r0_10k": 1.52,
        "e2_re_50k": 1.50, "e2_peak_re": 2.19, "e2_peak_tick": 243, "e2_peak_cases": 29_200,
        "e2_win_edges": 2_103_538, "e2_unique_hosts": 2_868_302,
        "reinf_all_pct": 24.9, "reinf_win_pct": 11.4, "max_inf": 4,
        "gaps_le30_pct": 0.03, "gaps_le60_pct": 0.63,
        "exp1_e2_peak_tick": 230, "exp2_e2_peak_tick": 190,
        "exp1_e1_peak_tick": 109, "exp2_e1_peak_tick": 88,
        "exp2_e1_gi": 4.94, "exp2_e2_gi": 4.58,
    },
    "mn": {
        "e1_total": 3_513_248, "e1_gi": 5.09, "e1_r0_1k": 3.06, "e1_r0_10k": 2.67,
        "e1_re_50k": 2.72, "e1_peak_re": 4.26, "e1_peak_tick": 97, "e1_peak_cases": 116_934,
        "e2_total": 4_108_415, "e2_gi": 4.62, "e2_r0_1k": 2.12, "e2_r0_10k": 1.55,
        "e2_re_50k": 1.46, "e2_peak_re": 2.98, "e2_peak_tick": 243, "e2_peak_cases": 46_196,
        "e2_win_edges": 3_031_622, "e2_unique_hosts": 3_194_733,
        "reinf_all_pct": 27.4, "reinf_win_pct": 16.9, "max_inf": 5,
        "gaps_le30_pct": 0.03, "gaps_le60_pct": 0.62,
        "exp1_e2_peak_tick": 262, "exp2_e2_peak_tick": 216,
        "exp1_e1_peak_tick": 115, "exp2_e1_peak_tick": 94,
        "exp2_e1_gi": 4.93, "exp2_e2_gi": 4.57,
    },
    "va": {
        "e1_total": 3_808_552, "e1_gi": 5.20, "e1_r0_1k": 3.45, "e1_r0_10k": 2.34,
        "e1_re_50k": 2.22, "e1_peak_re": 4.84, "e1_peak_tick": 116, "e1_peak_cases": 110_774,
        "e2_total": 4_617_163, "e2_gi": 4.70, "e2_r0_1k": 1.70, "e2_r0_10k": 1.65,
        "e2_re_50k": 1.73, "e2_peak_re": 1.82, "e2_peak_tick": 266, "e2_peak_cases": 36_253,
        "e2_win_edges": 3_213_682, "e2_unique_hosts": 3_639_635,
        "reinf_all_pct": 26.0, "reinf_win_pct": 15.1, "max_inf": 4,
        "gaps_le30_pct": 0.03, "gaps_le60_pct": 0.56,
        "exp1_e2_peak_tick": 247, "exp2_e2_peak_tick": 217,
        "exp1_e1_peak_tick": 122, "exp2_e1_peak_tick": 102,
        "exp2_e1_gi": 4.93, "exp2_e2_gi": 4.59,
    },
    "wa": {
        "e1_total": 3_100_280, "e1_gi": 5.29, "e1_r0_1k": 3.88, "e1_r0_10k": 1.90,
        "e1_re_50k": 2.18, "e1_peak_re": 5.45, "e1_peak_tick": 105, "e1_peak_cases": 78_092,
        "e2_total": 3_548_756, "e2_gi": 4.74, "e2_r0_1k": 2.12, "e2_r0_10k": 1.59,
        "e2_re_50k": 1.36, "e2_peak_re": 2.25, "e2_peak_tick": 248, "e2_peak_cases": 24_119,
        "e2_win_edges": 2_362_479, "e2_unique_hosts": 2_841_641,
        "reinf_all_pct": 24.1, "reinf_win_pct": 12.5, "max_inf": 4,
        "gaps_le30_pct": 0.03, "gaps_le60_pct": 0.55,
        "exp1_e2_peak_tick": 239, "exp2_e2_peak_tick": 200,
        "exp1_e1_peak_tick": 122, "exp2_e1_peak_tick": 95,
        "exp2_e1_gi": 4.92, "exp2_e2_gi": 4.58,
    },
}

# Pre-QC Poisson painter baseline on exp1 (from docs/abm_contribution.md)
EXP1_POISSON = {
    "ga": {"lambda": 0.5800, "gi": 6.93, "mu_truth": 9.01e-4},
    "ma": {"lambda": 0.5819, "gi": 6.88, "mu_truth": 9.57e-4},
    "mn": {"lambda": 0.5830, "gi": 6.76, "mu_truth": 8.85e-4},
    "va": {"lambda": 0.5818, "gi": 6.88, "mu_truth": 9.39e-4},
    "wa": {"lambda": 0.5831, "gi": 6.91, "mu_truth": 9.21e-4},
}


# ---------------------------------------------------------------------------
# SVG & File Helpers
# ---------------------------------------------------------------------------
def esc(s: object) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def svg_open(w: int, h: int) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
        f'viewBox="0 0 {w} {h}">\n'
        f'<rect width="{w}" height="{h}" fill="#ffffff"/>\n'
    )


def txt(
    x: float,
    y: float,
    s: object,
    size: int = 14,
    fill: str = INK,
    anchor: str = "start",
    weight: str = "normal",
    style: str = "",
    rotate: float | None = None,
) -> str:
    extra = f' font-style="{style}"' if style else ""
    if rotate is not None:
        extra += f' transform="rotate({rotate},{x:.1f},{y:.1f})"'
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
        f'fill="{fill}" text-anchor="{anchor}" font-weight="{weight}"{extra}>'
        f"{esc(s)}</text>"
    )


def line(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    stroke: str = INK,
    width: float = 1.5,
    dash: str = "",
    opacity: float = 1.0,
) -> str:
    d = f' stroke-dasharray="{dash}"' if dash else ""
    op = f' stroke-opacity="{opacity:.2f}"' if opacity < 1.0 else ""
    return (
        f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
        f'stroke="{stroke}" stroke-width="{width}"{d}{op}/>'
    )


def rect(
    x: float,
    y: float,
    w: float,
    h: float,
    fill: str = "none",
    stroke: str = "none",
    width: float = 1.0,
    rx: float = 0.0,
    opacity: float = 1.0,
) -> str:
    r_attr = f' rx="{rx:.1f}"' if rx > 0 else ""
    op_attr = f' fill-opacity="{opacity:.2f}"' if opacity < 1.0 else ""
    return (
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{width}"{r_attr}{op_attr}/>'
    )


def circle(
    cx: float,
    cy: float,
    r: float,
    fill: str = BLUE,
    stroke: str = "#ffffff",
    width: float = 1.5,
    opacity: float = 1.0,
) -> str:
    op = f' fill-opacity="{opacity:.2f}"' if opacity < 1.0 else ""
    return (
        f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{width}"{op}/>'
    )


def tick_to_date_str(tick: float) -> str:
    d = TICK_ZERO + dt.timedelta(days=int(round(tick)))
    return d.strftime("%Y-%m-%d")


def find_chrome() -> str | None:
    candidates = [
        "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
        "google-chrome",
        "chromium",
        "chromium-browser",
    ]
    for c in candidates:
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return c
        found = shutil.which(c)
        if found:
            return found
    return None


def convert_svg_to_png(svg_path: Path, width: int, height: int) -> Path | None:
    png_path = svg_path.with_suffix(".png")
    chrome = find_chrome()
    if chrome:
        cmd = [
            chrome,
            "--headless=new",
            "--disable-gpu",
            "--hide-scrollbars",
            f"--window-size={width},{height}",
            f"--screenshot={png_path}",
            f"file://{svg_path.resolve()}",
        ]
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if png_path.exists():
                return png_path
        except Exception:
            pass
    sips = shutil.which("sips") or ("/usr/bin/sips" if os.path.exists("/usr/bin/sips") else None)
    if sips:
        try:
            subprocess.run(
                [sips, "-s", "format", "png", str(svg_path), "--out", str(png_path)],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            if png_path.exists():
                return png_path
        except Exception:
            pass
    return None


# ---------------------------------------------------------------------------
# Data Loading
# ---------------------------------------------------------------------------
def load_seed_schedule(seed_dir: Path, state: str, variant: int) -> list[tuple[int, int]]:
    path = seed_dir / f"{state}-{variant}.csv"
    if not path.exists():
        return []
    rows = []
    with open(path, "r", newline="") as f:
        for r in csv.DictReader(f):
            rows.append((int(r["tick"]), int(r["count"])))
    return sorted(rows)


def load_vax_schedule(vax_dir: Path, state: str) -> list[tuple[int, int]]:
    path = vax_dir / f"{state}.csv"
    if not path.exists():
        return []
    by_tick: dict[int, int] = {}
    with open(path, "r", newline="") as f:
        for r in csv.DictReader(f):
            t = int(float(r["tick"]))
            cnt = int(float(r.get("count", r.get(" vaccinations", 0))))
            by_tick[t] = by_tick.get(t, 0) + cnt
    return sorted(by_tick.items())


def read_single_row_csv(path: Path) -> dict[str, str]:
    with open(path, "r", newline="") as f:
        rows = list(csv.DictReader(f))
        return rows[0] if rows else {}


def load_benchmarks(results_dir: Path, exp1_archive: Path) -> dict:
    data: dict[str, dict] = {"exp3": {}, "exp1": {}}
    for st in STATES:
        b3 = results_dir / f"{st}_delta_wave" / "05_benchmarks"
        b1 = exp1_archive / st
        data["exp3"][st] = {
            "seq_div": read_single_row_csv(b3 / "sequence_divergence.csv"),
            "clock": list(csv.DictReader(open(b3 / "clock_estimates.csv"))),
            "mugration": list(csv.DictReader(open(b3 / "Mugration_Metrics.csv"))),
            "auc": list(csv.DictReader(open(b3 / "AUC_truth_rankings.csv"))),
        }
        data["exp1"][st] = {
            "seq_div": read_single_row_csv(b1 / "sequence_divergence.csv"),
            "clock": list(csv.DictReader(open(b1 / "clock_estimates.csv"))),
            "mugration": list(csv.DictReader(open(b1 / "Mugration_Metrics.csv")))
            if (b1 / "Mugration_Metrics.csv").exists()
            else [],
            "auc": list(csv.DictReader(open(b1 / "AUC_truth_rankings.csv")))
            if (b1 / "AUC_truth_rankings.csv").exists()
            else [],
        }
    return data


# ---------------------------------------------------------------------------
# Table 1: Epidemiological Characterization (CSV + Markdown)
# ---------------------------------------------------------------------------
def build_table1_epi(seed_dir: Path, out_dir: Path) -> list[dict]:
    rows = []
    for st in STATES:
        m = ABM_EPI_METRICS[st]
        s2 = load_seed_schedule(seed_dir, st, 2)
        total_s2 = sum(c for _, c in s2)
        pre_exp1 = sum(c for t, c in s2 if t <= m["exp1_e2_peak_tick"])
        pre_exp2 = sum(c for t, c in s2 if t <= m["exp2_e2_peak_tick"])
        pre_exp3 = sum(c for t, c in s2 if t <= m["e2_peak_tick"])

        rows.append({
            "state": st.upper(),
            "e1_total": m["e1_total"],
            "e1_gi_days": f"{m['e1_gi']:.2f}",
            "e1_r0_1k": f"{m['e1_r0_1k']:.2f}",
            "e1_r0_10k": f"{m['e1_r0_10k']:.2f}",
            "e1_re_50k": f"{m['e1_re_50k']:.2f}",
            "e1_peak_7d_re": f"{m['e1_peak_re']:.2f}",
            "e1_peak_tick": m["e1_peak_tick"],
            "e1_peak_date": tick_to_date_str(m["e1_peak_tick"]),
            "e1_peak_daily": m["e1_peak_cases"],
            "e2_total": m["e2_total"],
            "e2_win_edges": m["e2_win_edges"],
            "e2_unique_hosts": m["e2_unique_hosts"],
            "e2_gi_days": f"{m['e2_gi']:.2f}",
            "e2_r0_1k": f"{m['e2_r0_1k']:.2f}",
            "e2_r0_10k": f"{m['e2_r0_10k']:.2f}",
            "e2_re_50k": f"{m['e2_re_50k']:.2f}",
            "e2_peak_7d_re": f"{m['e2_peak_re']:.2f}",
            "e2_peak_tick": m["e2_peak_tick"],
            "e2_peak_date": tick_to_date_str(m["e2_peak_tick"]),
            "e2_peak_daily": m["e2_peak_cases"],
            "e2_seeds_total": total_s2,
            "e2_seeds_le_peak_exp3": pre_exp3,
            "e2_seeds_le_peak_exp3_pct": f"{pre_exp3 / total_s2 * 100:.1f}",
            "e2_seeds_le_peak_exp2": pre_exp2,
            "e2_seeds_le_peak_exp2_pct": f"{pre_exp2 / total_s2 * 100:.1f}",
            "e2_seeds_le_peak_exp1": pre_exp1,
            "e2_seeds_le_peak_exp1_pct": f"{pre_exp1 / total_s2 * 100:.1f}",
            "reinf_all_pct": f"{m['reinf_all_pct']:.1f}",
            "reinf_win_pct": f"{m['reinf_win_pct']:.1f}",
            "max_inf_per_host": m["max_inf"],
            "gaps_le30d_pct": f"{m['gaps_le30_pct']:.2f}",
            "gaps_le60d_pct": f"{m['gaps_le60_pct']:.2f}",
        })

    # Mean row
    mean_e1_tick = stats.mean(ABM_EPI_METRICS[s]["e1_peak_tick"] for s in STATES)
    mean_e2_tick = stats.mean(ABM_EPI_METRICS[s]["e2_peak_tick"] for s in STATES)
    rows.append({
        "state": "MEAN",
        "e1_total": int(round(stats.mean(ABM_EPI_METRICS[s]["e1_total"] for s in STATES))),
        "e1_gi_days": f"{stats.mean(ABM_EPI_METRICS[s]['e1_gi'] for s in STATES):.2f}",
        "e1_r0_1k": f"{stats.mean(ABM_EPI_METRICS[s]['e1_r0_1k'] for s in STATES):.2f}",
        "e1_r0_10k": f"{stats.mean(ABM_EPI_METRICS[s]['e1_r0_10k'] for s in STATES):.2f}",
        "e1_re_50k": f"{stats.mean(ABM_EPI_METRICS[s]['e1_re_50k'] for s in STATES):.2f}",
        "e1_peak_7d_re": f"{stats.mean(ABM_EPI_METRICS[s]['e1_peak_re'] for s in STATES):.2f}",
        "e1_peak_tick": f"{mean_e1_tick:.1f}",
        "e1_peak_date": tick_to_date_str(mean_e1_tick),
        "e1_peak_daily": int(round(stats.mean(ABM_EPI_METRICS[s]["e1_peak_cases"] for s in STATES))),
        "e2_total": int(round(stats.mean(ABM_EPI_METRICS[s]["e2_total"] for s in STATES))),
        "e2_win_edges": int(round(stats.mean(ABM_EPI_METRICS[s]["e2_win_edges"] for s in STATES))),
        "e2_unique_hosts": int(round(stats.mean(ABM_EPI_METRICS[s]["e2_unique_hosts"] for s in STATES))),
        "e2_gi_days": f"{stats.mean(ABM_EPI_METRICS[s]['e2_gi'] for s in STATES):.2f}",
        "e2_r0_1k": f"{stats.mean(ABM_EPI_METRICS[s]['e2_r0_1k'] for s in STATES):.2f}",
        "e2_r0_10k": f"{stats.mean(ABM_EPI_METRICS[s]['e2_r0_10k'] for s in STATES):.2f}",
        "e2_re_50k": f"{stats.mean(ABM_EPI_METRICS[s]['e2_re_50k'] for s in STATES):.2f}",
        "e2_peak_7d_re": f"{stats.mean(ABM_EPI_METRICS[s]['e2_peak_re'] for s in STATES):.2f}",
        "e2_peak_tick": f"{mean_e2_tick:.1f}",
        "e2_peak_date": tick_to_date_str(mean_e2_tick),
        "e2_peak_daily": int(round(stats.mean(ABM_EPI_METRICS[s]["e2_peak_cases"] for s in STATES))),
        "e2_seeds_total": int(round(stats.mean(int(r["e2_seeds_total"]) for r in rows))),
        "e2_seeds_le_peak_exp3": int(round(stats.mean(int(r["e2_seeds_le_peak_exp3"]) for r in rows))),
        "e2_seeds_le_peak_exp3_pct": f"{stats.mean(float(r['e2_seeds_le_peak_exp3_pct']) for r in rows):.1f}",
        "e2_seeds_le_peak_exp2": int(round(stats.mean(int(r["e2_seeds_le_peak_exp2"]) for r in rows))),
        "e2_seeds_le_peak_exp2_pct": f"{stats.mean(float(r['e2_seeds_le_peak_exp2_pct']) for r in rows):.1f}",
        "e2_seeds_le_peak_exp1": int(round(stats.mean(int(r["e2_seeds_le_peak_exp1"]) for r in rows))),
        "e2_seeds_le_peak_exp1_pct": f"{stats.mean(float(r['e2_seeds_le_peak_exp1_pct']) for r in rows):.1f}",
        "reinf_all_pct": f"{stats.mean(ABM_EPI_METRICS[s]['reinf_all_pct'] for s in STATES):.1f}",
        "reinf_win_pct": f"{stats.mean(ABM_EPI_METRICS[s]['reinf_win_pct'] for s in STATES):.1f}",
        "max_inf_per_host": f"{stats.mean(ABM_EPI_METRICS[s]['max_inf'] for s in STATES):.1f}",
        "gaps_le30d_pct": f"{stats.mean(ABM_EPI_METRICS[s]['gaps_le30_pct'] for s in STATES):.2f}",
        "gaps_le60d_pct": f"{stats.mean(ABM_EPI_METRICS[s]['gaps_le60_pct'] for s in STATES):.2f}",
    })

    csv_path = out_dir / "exp3_epi_characterization.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    md_path = out_dir / "exp3_epi_characterization.md"
    with open(md_path, "w") as f:
        f.write("# Table 1: `exp3` Epidemiological Characterization Across 5 State Digital Twins\n\n")
        f.write("Generated by `scripts/utils/characterize_exp3.py` (`tickZero = 2020-11-30`, cell `0.25`, replicate `0`).\n\n")
        f.write("| State | `E1` Total | `E1` GI (d) | `E1` $R_0(1\\text{k})$ / $R_0(10\\text{k})$ | `E1` Peak Date | `E2` Total | `E2` GI (d) | `E2` $R_0(1\\text{k})$ / $R_0(10\\text{k})$ | `E2` Peak Date (Peak/day) | `E2` Seeds $\\le$ Peak (`exp3` vs `exp2`) | Window Reinf % | Gaps $\\le 30\\text{d}$ / $\\le 60\\text{d}$ |\n")
        f.write("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |\n")
        for r in rows:
            bold = "**" if r["state"] == "MEAN" else ""
            f.write(
                f"| {bold}{r['state']}{bold} | {int(r['e1_total']):,d} | {bold}{r['e1_gi_days']} d{bold} | "
                f"`{r['e1_r0_1k']}` / `{r['e1_r0_10k']}` | `{r['e1_peak_date']}` (tick {r['e1_peak_tick']}) | "
                f"{int(r['e2_total']):,d} | {bold}{r['e2_gi_days']} d{bold} | "
                f"`{r['e2_r0_1k']}` / `{r['e2_r0_10k']}` | `{r['e2_peak_date']}` ({int(r['e2_peak_daily']):,d}/d) | "
                f"**{r['e2_seeds_le_peak_exp3']} ({r['e2_seeds_le_peak_exp3_pct']}%)** vs. {r['e2_seeds_le_peak_exp2']} ({r['e2_seeds_le_peak_exp2_pct']}%) | "
                f"{r['reinf_win_pct']}% | `{r['gaps_le30d_pct']}%` / `{r['gaps_le60d_pct']}%` |\n"
            )
    return rows


# ---------------------------------------------------------------------------
# Table 2: Molecular Clock Closure (CSV + Markdown)
# ---------------------------------------------------------------------------
def build_table2_clock(bdata: dict, out_dir: Path) -> list[dict]:
    rows = []
    for st in STATES:
        sd1 = bdata["exp1"][st]["seq_div"]
        sd3 = bdata["exp3"][st]["seq_div"]
        ck1 = {(r["arm"], r["quantity"]): r for r in bdata["exp1"][st]["clock"]}
        ck3 = {(r["arm"], r["quantity"]): r for r in bdata["exp3"][st]["clock"]}

        lam1 = float(sd1["mean_substitutions"])
        lam3 = float(sd3["mean_substitutions"])
        gi1 = float(sd1["gi_mean_days"])
        gi3 = float(sd3["gi_mean_days"])

        imp1_yr = float(sd1["implied_subs_per_genome_per_year"])
        imp3_yr = float(sd3["implied_subs_per_genome_per_year"])
        imp1_site = (lam1 / gi1) * DAYS_PER_YEAR / GENOME_LEN
        imp3_site = (lam3 / gi3) * DAYS_PER_YEAR / GENOME_LEN

        t1_truth = ck1[("1S__surs", "mu_truth")]
        t3_truth = ck3[("1S__surs", "mu_truth")]
        t1_lin = ck1[("1S__surs", "mu_lineage")]
        t3_lin = ck3[("1S__surs", "mu_lineage")]
        t3_seed_w = ck3[("1S__surs", "mu_seed_trend_weighted")]
        t1_sim_surs = ck1[("1S__surs", "mu_sim")]
        t3_sim_surs = ck3[("1S__surs", "mu_sim")]
        t3_sim_lasso = ck3.get(("1S-P__lasso_greedy", "mu_sim"), t3_sim_surs)

        mu_truth1 = float(t1_truth["rate_subs_per_site_per_year"])
        mu_truth3 = float(t3_truth["rate_subs_per_site_per_year"])
        mu_lin1 = float(t1_lin["rate_subs_per_site_per_year"])
        mu_lin3 = float(t3_lin["rate_subs_per_site_per_year"])
        mu_sim1_surs = float(t1_sim_surs["rate_subs_per_site_per_year"])
        mu_sim3_surs = float(t3_sim_surs["rate_subs_per_site_per_year"])
        mu_sim3_lasso = float(t3_sim_lasso["rate_subs_per_site_per_year"])

        r2_1 = float(t1_truth["r_squared"])
        r2_3 = float(t3_truth["r_squared"])
        spath1 = mu_truth1 / imp1_site
        spath3 = mu_truth3 / imp3_site

        rows.append({
            "state": st.upper(),
            "exp1_edges": int(sd1["gi_edges"]),
            "exp3_edges": int(sd3["gi_edges"]),
            "exp1_lambda": lam1,
            "exp3_lambda": lam3,
            "exp3_frac_unchanged": float(sd3["frac_unchanged"]),
            "exp1_gi_days": gi1,
            "exp3_gi_days": gi3,
            "exp3_gi_median_iqr": f"{float(sd3['gi_median_days']):.1f} ({int(float(sd3['gi_q25_days']))}-{int(float(sd3['gi_q75_days']))})",
            "exp1_implied_per_yr": imp1_yr,
            "exp3_implied_per_yr": imp3_yr,
            "exp1_implied_per_site_yr": imp1_site,
            "exp3_implied_per_site_yr": imp3_site,
            "exp1_mu_lineage": mu_lin1,
            "exp3_mu_lineage": mu_lin3,
            "exp1_s_path": spath1,
            "exp3_s_path": spath3,
            "exp1_mu_truth": mu_truth1,
            "exp3_mu_truth": mu_truth3,
            "exp3_mu_truth_vs_sc2": mu_truth3 / REAL_SC2_CLOCK,
            "exp1_mu_sim_surs": mu_sim1_surs,
            "exp3_mu_sim_surs": mu_sim3_surs,
            "exp3_mu_sim_lasso": mu_sim3_lasso,
            "exp1_r_squared": r2_1,
            "exp3_r_squared": r2_3,
            "exp3_within_date_var_share": float(t3_lin["within_date_variance_share"]),
            "exp3_n_chains": int(t3_lin["chains"]),
            "exp3_n_chains_surviving": int(t3_seed_w["surviving_imports"]),
            "gi_unmatched_nonpos": f"{sd3['gi_parent_unmatched']}/{sd3['gi_nonpositive']}",
        })

    # Add 5-state MEAN row and Pre-Run Prediction row
    mean_row = {
        "state": "MEAN (5 States)",
        "exp1_edges": int(round(stats.mean(r["exp1_edges"] for r in rows))),
        "exp3_edges": int(round(stats.mean(r["exp3_edges"] for r in rows))),
        "exp1_lambda": stats.mean(r["exp1_lambda"] for r in rows),
        "exp3_lambda": stats.mean(r["exp3_lambda"] for r in rows),
        "exp3_frac_unchanged": stats.mean(r["exp3_frac_unchanged"] for r in rows),
        "exp1_gi_days": stats.mean(r["exp1_gi_days"] for r in rows),
        "exp3_gi_days": stats.mean(r["exp3_gi_days"] for r in rows),
        "exp3_gi_median_iqr": "4.0 (3-6)",
        "exp1_implied_per_yr": stats.mean(r["exp1_implied_per_yr"] for r in rows),
        "exp3_implied_per_yr": stats.mean(r["exp3_implied_per_yr"] for r in rows),
        "exp1_implied_per_site_yr": stats.mean(r["exp1_implied_per_site_yr"] for r in rows),
        "exp3_implied_per_site_yr": stats.mean(r["exp3_implied_per_site_yr"] for r in rows),
        "exp1_mu_lineage": stats.mean(r["exp1_mu_lineage"] for r in rows),
        "exp3_mu_lineage": stats.mean(r["exp3_mu_lineage"] for r in rows),
        "exp1_s_path": stats.mean(r["exp1_s_path"] for r in rows),
        "exp3_s_path": stats.mean(r["exp3_s_path"] for r in rows),
        "exp1_mu_truth": stats.mean(r["exp1_mu_truth"] for r in rows),
        "exp3_mu_truth": stats.mean(r["exp3_mu_truth"] for r in rows),
        "exp3_mu_truth_vs_sc2": stats.mean(r["exp3_mu_truth_vs_sc2"] for r in rows),
        "exp1_mu_sim_surs": stats.mean(r["exp1_mu_sim_surs"] for r in rows),
        "exp3_mu_sim_surs": stats.mean(r["exp3_mu_sim_surs"] for r in rows),
        "exp3_mu_sim_lasso": stats.mean(r["exp3_mu_sim_lasso"] for r in rows),
        "exp1_r_squared": stats.mean(r["exp1_r_squared"] for r in rows),
        "exp3_r_squared": stats.mean(r["exp3_r_squared"] for r in rows),
        "exp3_within_date_var_share": stats.mean(r["exp3_within_date_var_share"] for r in rows),
        "exp3_n_chains": int(round(stats.mean(r["exp3_n_chains"] for r in rows))),
        "exp3_n_chains_surviving": round(stats.mean(r["exp3_n_chains_surviving"] for r in rows), 1),
        "gi_unmatched_nonpos": "0/0",
    }
    pred_row = {
        "state": "PRE-RUN PREDICTION",
        "exp1_edges": mean_row["exp1_edges"],
        "exp3_edges": mean_row["exp3_edges"],
        "exp1_lambda": 0.3428,
        "exp3_lambda": 0.3430,
        "exp3_frac_unchanged": 0.7110,
        "exp1_gi_days": 6.87,
        "exp3_gi_days": 4.69,
        "exp3_gi_median_iqr": "4.0 (3-6)",
        "exp1_implied_per_yr": 18.21,
        "exp3_implied_per_yr": 26.60,
        "exp1_implied_per_site_yr": 6.13e-4,
        "exp3_implied_per_site_yr": 8.91e-4,
        "exp1_mu_lineage": 6.53e-4,
        "exp3_mu_lineage": 8.91e-4,
        "exp1_s_path": 0.885,
        "exp3_s_path": 0.885,
        "exp1_mu_truth": 5.42e-4,
        "exp3_mu_truth": 7.89e-4,
        "exp3_mu_truth_vs_sc2": 7.89e-4 / REAL_SC2_CLOCK,
        "exp1_mu_sim_surs": mean_row["exp1_mu_sim_surs"],
        "exp3_mu_sim_surs": 7.89e-4,
        "exp3_mu_sim_lasso": 7.89e-4,
        "exp1_r_squared": 0.314,
        "exp3_r_squared": 0.450,
        "exp3_within_date_var_share": 0.950,
        "exp3_n_chains": mean_row["exp3_n_chains"],
        "exp3_n_chains_surviving": mean_row["exp3_n_chains_surviving"],
        "gi_unmatched_nonpos": "0/0",
    }
    all_rows = rows + [mean_row, pred_row]

    csv_path = out_dir / "exp3_molecular_clock_closure.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
        w.writeheader()
        w.writerows(all_rows)

    md_path = out_dir / "exp3_molecular_clock_closure.md"
    with open(md_path, "w") as f:
        f.write("# Table 2: Molecular-Clock Decomposition (`exp1` Baseline vs. Pre-Run Prediction vs. Measured `exp3`)\n\n")
        f.write("Generated by `scripts/utils/characterize_exp3.py`. Baseline TreeTime clock (`mu_sim`) is reported on `1S__surs` (random surveillance baseline), separately from `1S-P__lasso_greedy`.\n\n")
        f.write("| State | Painted Edges (`exp3`) | $\\hat{\\lambda}$ (`exp1` $\\to$ `exp3`) | GI (`exp1` $\\to$ `exp3`) | Implied $\\lambda/\\text{GI}$ (`/yr` \\| `/site/yr`) | Within-Chain $\\mu_{\\text{lineage}}$ (`/site/yr`) | Path Factor $s_{\\text{path}}$ | Pooled Truth $\\mu_{\\text{truth}}$ (`/site/yr`, vs. `8e-4`) | TreeTime $\\mu_{\\text{sim}}$ (`1S__surs` \\| `1S-P__lasso`) | Root-to-Tip $R^2$ (`exp1` $\\to$ `exp3`) |\n")
        f.write("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |\n")
        for r in all_rows:
            is_summary = r["state"] in ("MEAN (5 States)", "PRE-RUN PREDICTION")
            b = "**" if is_summary else ""
            f.write(
                f"| {b}{r['state']}{b} | {int(r['exp3_edges']):,d} | "
                f"`{r['exp1_lambda']:.4f}` $\\to$ {b}`{r['exp3_lambda']:.4f}`{b} | "
                f"`{r['exp1_gi_days']:.2f} d` $\\to$ {b}`{r['exp3_gi_days']:.2f} d`{b} | "
                f"{b}`{r['exp3_implied_per_yr']:.2f}`{b} \\| `{r['exp3_implied_per_site_yr']:.2e}` | "
                f"`{r['exp3_mu_lineage']:.2e}` | `{r['exp3_s_path']:.3f}` | "
                f"{b}`{r['exp3_mu_truth']:.2e}` ({r['exp3_mu_truth_vs_sc2']:.2f}x){b} | "
                f"`{r['exp3_mu_sim_surs']:.2e}` \\| `{r['exp3_mu_sim_lasso']:.2e}` | "
                f"`{r['exp1_r_squared']:.3f}` $\\to$ {b}`{r['exp3_r_squared']:.3f}`{b} |\n"
            )
    return all_rows


# ---------------------------------------------------------------------------
# Tables 3 & 4: Downstream Surveillance & Phylogeography Benchmarks
# ---------------------------------------------------------------------------
SCHEME_ORDER = [
    "1S__surs",
    "1S-P__lasso_greedy",
    "1S-1_LL__lasso_greedy",
    "1S-1_LL-P__lasso_greedy",
    "4S__surs",
    "4S-P__lasso_greedy",
    "4S-4_LL__lasso_greedy",
    "4S-4_LL-P__lasso_greedy",
]


def format_recipe_name(sc: str, compact: bool = False) -> str:
    """Format canonical recipe IDs using explicit algorithm family + subgroup count
    (e.g., 'LASSO50-Greedy · 4S-P' or 'SURS · 4S') so labels scale cleanly when
    LASSO50-Stratified, LASSO1500-Greedy, or Stratified (quota) are included."""
    if "__" not in sc:
        return sc
    scenario, algo = sc.split("__", 1)
    family_map = {
        "lasso_greedy": "LASSO50-Greedy",
        "lasso_stratified": "LASSO50-Stratified",
        "stratified": "Stratified",
        "surs": "SURS",
    }
    fam = family_map.get(algo, algo)
    if compact:
        return f"{fam} ({scenario})"
    return f"{fam} · {scenario}"


def build_table3_mugration(bdata: dict, out_dir: Path) -> list[dict]:
    by_scheme: dict[str, dict[str, dict]] = {}
    by_scheme_exp1: dict[str, dict[str, dict]] = {}
    for st in STATES:
        for r in bdata["exp3"][st]["mugration"]:
            by_scheme.setdefault(r["label"], {})[st] = r
        for r in bdata["exp1"][st]["mugration"]:
            by_scheme_exp1.setdefault(r["label"], {})[st] = r

    # Rank order: top-performing LASSO50-Greedy & SURS strategies by Topological F1 / Cosine
    ranked_schemes = [
        "1S-1_LL-P__lasso_greedy",
        "4S-P__lasso_greedy",
        "4S-4_LL-P__lasso_greedy",
        "1S-P__lasso_greedy",
        "1S-1_LL__lasso_greedy",
        "1S__surs",
        "4S__surs",
        "4S-4_LL__lasso_greedy",
    ]

    rows = []
    for sc in ranked_schemes:
        if sc not in by_scheme:
            continue
        st_map = by_scheme[sc]
        st_map1 = by_scheme_exp1.get(sc, {})
        base_sc = "1S__surs" if sc.startswith("1S") else "4S__surs"
        base_map = by_scheme[base_sc]
        base_map1 = by_scheme_exp1.get(base_sc, {})

        cos_vals = [float(st_map[s]["cosine_similarity"]) for s in STATES if s in st_map]
        f1_vals = [float(st_map[s]["topological_f1"]) for s in STATES if s in st_map]
        pr_vals = [float(st_map[s]["pearson_r"]) for s in STATES if s in st_map]
        mae_vals = [float(st_map[s]["masked_mae"]) for s in STATES if s in st_map]

        cos1_vals = [float(st_map1[s]["cosine_similarity"]) for s in STATES if s in st_map1]
        f1_1_vals = [float(st_map1[s]["topological_f1"]) for s in STATES if s in st_map1]
        abs_delta_cos = [
            abs(float(st_map[s]["cosine_similarity"]) - float(st_map1[s]["cosine_similarity"]))
            for s in STATES if s in st_map and s in st_map1
        ]

        b_cos_vals = [float(base_map[s]["cosine_similarity"]) for s in STATES if s in base_map]
        b_f1_vals = [float(base_map[s]["topological_f1"]) for s in STATES if s in base_map]

        mean_cos = stats.mean(cos_vals)
        mean_f1 = stats.mean(f1_vals)
        base_cos = stats.mean(b_cos_vals)
        base_f1 = stats.mean(b_f1_vals)

        cos_gain_pct = (mean_cos / base_cos - 1.0) * 100.0 if base_cos else 0.0
        f1_gain_pct = (mean_f1 / base_f1 - 1.0) * 100.0 if base_f1 else 0.0
        cos_wins = sum(1 for s in STATES if float(st_map[s]["cosine_similarity"]) > float(base_map[s]["cosine_similarity"]))
        f1_wins = sum(1 for s in STATES if float(st_map[s]["topological_f1"]) > float(base_map[s]["topological_f1"]))
        cos_wins_10 = cos_wins + sum(
            1 for s in STATES if s in st_map1 and float(st_map1[s]["cosine_similarity"]) > float(base_map1[s]["cosine_similarity"])
        )
        f1_wins_10 = f1_wins + sum(
            1 for s in STATES if s in st_map1 and float(st_map1[s]["topological_f1"]) > float(base_map1[s]["topological_f1"])
        )

        row = {
            "scheme": sc,
            "display_name": format_recipe_name(sc, compact=False),
            "compact_name": format_recipe_name(sc, compact=True),
            "stride_baseline": base_sc,
            "mean_topological_f1": mean_f1,
            "sd_topological_f1": stats.stdev(f1_vals) if len(f1_vals) > 1 else 0.0,
            "f1_gain_vs_surs_pct": f1_gain_pct,
            "f1_wins_vs_surs": f1_wins,
            "f1_wins_10runs": f1_wins_10,
            "mean_cosine": mean_cos,
            "sd_cosine": stats.stdev(cos_vals) if len(cos_vals) > 1 else 0.0,
            "cos_gain_vs_surs_pct": cos_gain_pct,
            "cos_wins_vs_surs": cos_wins,
            "cos_wins_10runs": cos_wins_10,
            "exp1_mean_cosine": stats.mean(cos1_vals) if cos1_vals else float("nan"),
            "exp1_mean_f1": stats.mean(f1_1_vals) if f1_1_vals else float("nan"),
            "mean_abs_delta_cos_exp1_exp3": stats.mean(abs_delta_cos) if abs_delta_cos else 0.0,
            "min_cos_floor_10runs": min(cos_vals + cos1_vals) if (cos_vals + cos1_vals) else float("nan"),
            "mean_pearson_r": stats.mean(pr_vals),
            "mean_masked_mae": stats.mean(mae_vals),
        }
        for s in STATES:
            row[f"{s}_cosine"] = float(st_map[s]["cosine_similarity"])
            row[f"{s}_f1"] = float(st_map[s]["topological_f1"])
        rows.append(row)

    csv_path = out_dir / "exp3_mugration_benchmarks.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    md_path = out_dir / "exp3_mugration_benchmarks.md"
    with open(md_path, "w") as f:
        f.write("# Table 3: `exp3` Inter-County Phylogeographic Reconstruction & Cross-Experiment Resilience (`LASSO50-Greedy` vs. `SURS`)\n\n")
        f.write("Generated by `scripts/utils/characterize_exp3.py`. Ranked by 5-state mean Topological $F_1$ (flow network edge recovery). Relative gains (`%`) and state win counts (`exp3` `/5` and pooled `exp1+exp3` `/10`) are computed against the stride-matched `SURS` baseline (`SURS · 1S` or `SURS · 4S`). **Regime Drift $|\\Delta\\text{Cos}|$** measures mean per-state sensitivity when epidemic dynamics shift between `exp1` and `exp3` (lower = more regime-resilient).\n\n")
        f.write("| Strategy Recipe (`LASSO50-Greedy` vs. `SURS`) | Flow Topology $F_1$ ($\\pm$ SD) | $F_1$ vs. `SURS` (`exp3` \\| `10-Run`) | Flow Proportions Cosine ($\\pm$ SD) | Cosine vs. `SURS` (`exp3` \\| `10-Run`) | Regime Drift $|\\Delta\\text{Cos}|$ (`exp1` $\\leftrightarrow$ `exp3`) | 10-Run Cosine Floor (Min) | GA Cos/$F_1$ | MA Cos/$F_1$ | MN Cos/$F_1$ | VA Cos/$F_1$ | WA Cos/$F_1$ |\n")
        f.write("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |\n")
        for r in rows:
            is_base = r["scheme"] in ("1S__surs", "4S__surs")
            f1_vs = "Baseline" if is_base else f"**{r['f1_gain_vs_surs_pct']:+.1f}%** ({r['f1_wins_vs_surs']}/5 \\| **{r['f1_wins_10runs']}/10**)"
            cos_vs = "Baseline" if is_base else f"`{r['cos_gain_vs_surs_pct']:+.1f}%` ({r['cos_wins_vs_surs']}/5 \\| {r['cos_wins_10runs']}/10)"
            f.write(
                f"| **{r['display_name']}** (`{r['scheme']}`) | **{r['mean_topological_f1']:.3f}** ($\\pm {r['sd_topological_f1']:.3f}$) | "
                f"{f1_vs} | "
                f"**{r['mean_cosine']:.3f}** ($\\pm {r['sd_cosine']:.3f}$) | "
                f"{cos_vs} | "
                f"**{r['mean_abs_delta_cos_exp1_exp3']:.3f}** | "
                f"`{r['min_cos_floor_10runs']:.3f}` | "
                f"`{r['ga_cosine']:.3f}`/`{r['ga_f1']:.3f}` | "
                f"`{r['ma_cosine']:.3f}`/`{r['ma_f1']:.3f}` | "
                f"`{r['mn_cosine']:.3f}`/`{r['mn_f1']:.3f}` | "
                f"`{r['va_cosine']:.3f}`/`{r['va_f1']:.3f}` | "
                f"`{r['wa_cosine']:.3f}`/`{r['wa_f1']:.3f}` |\n"
            )
    return rows


def build_table4_auc(bdata: dict, out_dir: Path) -> list[dict]:
    eval_types = [
        "8_week_rolling_tree_coverage",
        "coverage_size_0",
        "coverage_size_100",
        "stride_component_coverage",
        "stride_variant_prevalence_error",
        "cumulative_infections",
        "equity_Senior_65+",
        "equity_Adult_18-49",
        "equity_Preschool_0-4",
    ]
    # Matches _HIGHER_IS_BETTER_PREFIXES in src/phylogas/benchmark/truth_runner.py:
    # Mean Reciprocal Distance (MRD = 1/|P_t| * sum 1/(d(u,S)+1)) and component
    # coverage are higher-is-better; L1 prevalence error and KL divergence are lower-is-better.
    higher_is_better_prefixes = (
        "stride_component_coverage",
        "coverage_size_",
        "8_week_rolling_tree_coverage",
        "equity_",
    )
    ranked_schemes = [
        "1S-1_LL-P__lasso_greedy",
        "4S-P__lasso_greedy",
        "4S-4_LL-P__lasso_greedy",
        "1S-P__lasso_greedy",
        "1S-1_LL__lasso_greedy",
        "1S__surs",
        "4S__surs",
        "4S-4_LL__lasso_greedy",
    ]
    # Store auc[eval_type][scheme][state] for exp3 and exp1
    auc_map: dict[str, dict[str, dict[str, float]]] = {}
    auc_map1: dict[str, dict[str, dict[str, float]]] = {}
    for st in STATES:
        for r in bdata["exp3"][st]["auc"]:
            auc_map.setdefault(r["eval_type"], {}).setdefault(r["algorithm"], {})[st] = float(r["auc"])
        for r in bdata["exp1"][st]["auc"]:
            auc_map1.setdefault(r["eval_type"], {}).setdefault(r["algorithm"], {})[st] = float(r["auc"])

    rows = []
    for sc in ranked_schemes:
        base_sc = "1S__surs" if sc.startswith("1S") else "4S__surs"
        row: dict[str, object] = {
            "scheme": sc,
            "display_name": format_recipe_name(sc, compact=False),
            "compact_name": format_recipe_name(sc, compact=True),
            "stride_baseline": base_sc,
        }
        for et in eval_types:
            vals = [auc_map.get(et, {}).get(sc, {}).get(s, float("nan")) for s in STATES]
            b_vals = [auc_map.get(et, {}).get(base_sc, {}).get(s, float("nan")) for s in STATES]
            vals = [v for v in vals if not math.isnan(v)]
            b_vals = [v for v in b_vals if not math.isnan(v)]
            m_val = stats.mean(vals) if vals else float("nan")
            b_val = stats.mean(b_vals) if b_vals else float("nan")
            row[f"{et}_mean"] = m_val
            row[f"{et}_sd"] = stats.stdev(vals) if len(vals) > 1 else 0.0
            if et.startswith(higher_is_better_prefixes):
                # Higher Mean Reciprocal Distance (MRD) or component fraction is better
                row[f"{et}_gain_vs_surs_pct"] = (m_val / b_val - 1.0) * 100.0 if b_val else 0.0
                w3 = sum(1 for s in STATES if auc_map[et][sc][s] > auc_map[et][base_sc][s])
                w1 = sum(1 for s in STATES if s in auc_map1.get(et, {}).get(sc, {}) and auc_map1[et][sc][s] > auc_map1[et][base_sc][s])
                row[f"{et}_wins_vs_surs"] = w3
                row[f"{et}_wins_10runs"] = w3 + w1
            else:
                # Lower L1/KL divergence error is better -> positive % means lower error than SURS
                row[f"{et}_gain_vs_surs_pct"] = (1.0 - m_val / b_val) * 100.0 if b_val else 0.0
                w3 = sum(1 for s in STATES if auc_map[et][sc][s] < auc_map[et][base_sc][s])
                w1 = sum(1 for s in STATES if s in auc_map1.get(et, {}).get(sc, {}) and auc_map1[et][sc][s] < auc_map1[et][base_sc][s])
                row[f"{et}_wins_vs_surs"] = w3
                row[f"{et}_wins_10runs"] = w3 + w1
        rows.append(row)

    csv_path = out_dir / "exp3_auc_truth_benchmarks.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    md_path = out_dir / "exp3_auc_truth_benchmarks.md"
    with open(md_path, "w") as f:
        f.write("# Table 4: `exp3` Epidemic Truth & Demographic Equity Coverage (`LASSO50-Greedy` vs. `SURS`)\n\n")
        f.write("Generated by `scripts/utils/characterize_exp3.py` (see [`docs/benchmarks_and_metrics.md`](../benchmarks_and_metrics.md#5-ground-truth-surveillance--phylogenetic-tree-coverage-metrics-auc_truth_rankingscsv)). For **Mean Reciprocal Distance (MRD)** tree-coverage metrics (`coverage_size_*`, `8_week_rolling_tree_coverage`, `equity_*`) and `stride_component_coverage` (Clusters Reached), **higher AUC ($\\uparrow$) = closer phylogenetic proximity / greater cluster capture** (`1/(d+1)` along the true ABM transmission graph); for `stride_variant_prevalence_error` ($L_1$) and `cumulative_infections` (KL divergence), **lower AUC ($\\downarrow$) = lower distribution error**.\n\n")
        f.write("| Strategy Recipe (`LASSO50-Greedy` vs. `SURS`) | Clusters Reached ($\\uparrow$) (`exp3` \\| `10-Run`) | Senior 65+ Equity MRD ($\\uparrow$) (`exp3` \\| `10-Run`) | Adult 18–49 Equity MRD ($\\uparrow$) (`exp3` \\| `10-Run`) | Tree Coverage $C_0$ MRD ($\\uparrow$) (`exp3` \\| `10-Run`) | 8-Wk Rolling Tree MRD ($\\uparrow$) | Variant Prevalence $L_1$ ($\\downarrow$) | Cumul. Infections KL ($\\downarrow$) |\n")
        f.write("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |\n")
        for r in rows:
            is_base = r["scheme"] in ("1S__surs", "4S__surs")
            cl_str = (
                f"**{r['stride_component_coverage_mean']:.4f}** (Base)"
                if is_base
                else f"**{r['stride_component_coverage_mean']:.4f}** (`{r['stride_component_coverage_gain_vs_surs_pct']:+.1f}%`, {r['stride_component_coverage_wins_vs_surs']}/5 \\| **{r['stride_component_coverage_wins_10runs']}/10**)"
            )
            sr_str = (
                f"**{r['equity_Senior_65+_mean']:.4f}** (Base)"
                if is_base
                else f"**{r['equity_Senior_65+_mean']:.4f}** (`{r['equity_Senior_65+_gain_vs_surs_pct']:+.1f}%`, {r['equity_Senior_65+_wins_vs_surs']}/5 \\| **{r['equity_Senior_65+_wins_10runs']}/10**)"
            )
            ad_str = (
                f"`{r['equity_Adult_18-49_mean']:.4f}` (Base)"
                if is_base
                else f"`{r['equity_Adult_18-49_mean']:.4f}` (`{r['equity_Adult_18-49_gain_vs_surs_pct']:+.1f}%`, {r['equity_Adult_18-49_wins_vs_surs']}/5 \\| **{r['equity_Adult_18-49_wins_10runs']}/10**)"
            )
            c0_str = (
                f"`{r['coverage_size_0_mean']:.4f}` (Base)"
                if is_base
                else f"`{r['coverage_size_0_mean']:.4f}` ({r['coverage_size_0_wins_vs_surs']}/5 \\| **{r['coverage_size_0_wins_10runs']}/10**)"
            )
            f.write(
                f"| **{r['display_name']}** (`{r['scheme']}`) | {cl_str} | {sr_str} | {ad_str} | {c0_str} | "
                f"`{r['8_week_rolling_tree_coverage_mean']:.4f}` | "
                f"**{r['stride_variant_prevalence_error_mean']:.4f}** | "
                f"`{r['cumulative_infections_mean']:.4f}` |\n"
            )
    return rows


# ---------------------------------------------------------------------------
# Figure 1: Epidemiological Dynamics (GI, Waning Gaps, Cohort Re, Seed Timing)
# ---------------------------------------------------------------------------
def plot_fig1_epi_dynamics(seed_dir: Path, out_dir: Path) -> Path:
    W, H = 1600, 560
    p = [svg_open(W, H)]

    # Header
    p.append(txt(32, 36, "Figure 1: exp3 Epidemiological Calibration — Generation Intervals, Cohort R0/Re, and Importation Alignment", 20, INK, weight="bold"))
    p.append(txt(32, 58, "5-state EpiHiper digital twins (GA, MA, MN, VA, WA; cell 0.25, replicate 0, tickZero = 2020-11-30). Comparing exp1, exp2, and exp3.", 13, INK2))

    # --- Panel A: Generation Interval & Waning Immunity Calibration (Left) ---
    ax_l, ax_r, ax_t, ax_b = 75, 495, 110, 450
    p.append(rect(28, 76, 495, 460, fill=BG_PANEL, stroke=GRID, rx=8))
    p.append(txt(44, 100, "A. Generation Interval (d) & Short-Gap Reinfections", 15, INK, weight="bold"))

    # Y axis: GI days from 0 to 8.0
    y_min, y_max = 0.0, 8.0
    def Y_gi(v: float) -> float:
        return ax_b - (v - y_min) / (y_max - y_min) * (ax_b - ax_t)

    for gv in [0, 2, 4, 6, 8]:
        yy = Y_gi(gv)
        p.append(line(ax_l, yy, ax_r, yy, stroke=GRID, width=1))
        p.append(txt(ax_l - 8, yy + 4, f"{gv} d", 11, INK2, anchor="end"))

    # Hart et al. (2022) reference lines (labels in open space above 5.5d on the right)
    p.append(line(ax_l, Y_gi(5.5), ax_r, Y_gi(5.5), stroke=BLUE, width=1.5, dash="5,4"))
    p.append(txt(ax_r - 6, Y_gi(5.5) - 22, "Dashed blue: Hart et al. (2022) Alpha = 5.5 d", 10, BLUE, anchor="end", weight="bold"))
    p.append(line(ax_l, Y_gi(4.7), ax_r, Y_gi(4.7), stroke=PURPLE, width=1.5, dash="5,4"))
    p.append(txt(ax_r - 6, Y_gi(5.5) - 8, "Dashed purple: Hart et al. (2022) Delta = 4.7 d", 10, PURPLE, anchor="end", weight="bold"))

    groups = [
        ("exp1 E1/E2", [6.93, 6.88, 6.76, 6.88, 6.91], ORANGE),
        ("exp2 E1", [ABM_EPI_METRICS[s]["exp2_e1_gi"] for s in STATES], AQUA),
        ("exp3 E1 (Alpha)", [ABM_EPI_METRICS[s]["e1_gi"] for s in STATES], BLUE),
        ("exp2 E2", [ABM_EPI_METRICS[s]["exp2_e2_gi"] for s in STATES], GRAY),
        ("exp3 E2 (Delta)", [ABM_EPI_METRICS[s]["e2_gi"] for s in STATES], PURPLE),
    ]
    col_w = (ax_r - ax_l) / len(groups)
    for idx, (label, vals, col) in enumerate(groups):
        cx = ax_l + (idx + 0.5) * col_w
        mean_v = stats.mean(vals)
        yy = Y_gi(mean_v)
        p.append(rect(cx - 18, yy, 36, ax_b - yy, fill=col, opacity=0.25, stroke=col, width=1.8))
        for si, v in enumerate(vals):
            jitter = (si - 2) * 5.0
            p.append(circle(cx + jitter, Y_gi(v), 4.0, fill=col))
        p.append(txt(cx, yy - 8, f"{mean_v:.2f}d", 11, col, anchor="middle", weight="bold"))
        p.append(txt(cx, ax_b + 18, label, 10, INK, anchor="middle", weight="bold"))

    # Callout box for Gamma(4, 45.75) waning immunity fix inside Panel A
    p.append(rect(44, 492, 462, 34, fill="#ffffff", stroke=GRID, rx=5))
    p.append(txt(52, 513, "Gamma(4, 45.75) Waning:  Reinf <=30d: 15.1% -> 0.03%  |  <=60d: 28.0% -> 0.56%", 10, INK, weight="bold"))

    # --- Panel B: Cohort Reproduction Numbers R0 / Re across the Wave (Center) ---
    bx_l, bx_r, bx_t, bx_b = 585, 1035, 130, 450
    p.append(rect(538, 76, 515, 460, fill=BG_PANEL, stroke=GRID, rx=8))
    p.append(txt(554, 100, "B. exp3 Cohort Reproduction Numbers (R0 & Re)", 15, INK, weight="bold"))
    p.append(txt(554, 118, "Dots = 5 states; bars = 5-state mean (E1 Alpha in blue, E2 Delta in purple)", 11, INK2))

    ry_min, ry_max = 0.0, 6.0
    def Y_re(v: float) -> float:
        return bx_b - (v - ry_min) / (ry_max - ry_min) * (bx_b - bx_t)

    for rv in [0, 1, 2, 3, 4, 5, 6]:
        yy = Y_re(rv)
        stroke_c = INK2 if rv == 1 else GRID
        dash_c = "4,3" if rv == 1 else ""
        p.append(line(bx_l, yy, bx_r, yy, stroke=stroke_c, width=1.2 if rv == 1 else 1.0, dash=dash_c))
        p.append(txt(bx_l - 8, yy + 4, f"{rv:.1f}", 11, INK2, anchor="end"))
    p.append(txt(bx_r - 4, Y_re(1.0) - 5, "Epidemic threshold (Re = 1.0)", 10, INK2, anchor="end"))

    re_metrics = [
        ("R0 (1k)", "e1_r0_1k", "e2_r0_1k"),
        ("R0 (10k)", "e1_r0_10k", "e2_r0_10k"),
        ("Re (50k)", "e1_re_50k", "e2_re_50k"),
        ("Peak 7d Re", "e1_peak_re", "e2_peak_re"),
    ]
    m_w = (bx_r - bx_l) / len(re_metrics)
    for idx, (mlabel, k1, k2) in enumerate(re_metrics):
        cx = bx_l + (idx + 0.5) * m_w
        v1 = [ABM_EPI_METRICS[s][k1] for s in STATES]
        v2 = [ABM_EPI_METRICS[s][k2] for s in STATES]
        m1, m2 = stats.mean(v1), stats.mean(v2)

        # E1 bar & dots
        x1 = cx - 20
        p.append(rect(x1 - 14, Y_re(m1), 28, bx_b - Y_re(m1), fill=BLUE, opacity=0.28, stroke=BLUE, width=1.6))
        for si, v in enumerate(v1):
            p.append(circle(x1 + (si - 2) * 3.5, Y_re(v), 3.8, fill=BLUE))
        p.append(txt(x1, Y_re(max(v1)) - 7, f"{m1:.2f}", 10, BLUE, anchor="middle", weight="bold"))

        # E2 bar & dots
        x2 = cx + 20
        p.append(rect(x2 - 14, Y_re(m2), 28, bx_b - Y_re(m2), fill=PURPLE, opacity=0.28, stroke=PURPLE, width=1.6))
        for si, v in enumerate(v2):
            p.append(circle(x2 + (si - 2) * 3.5, Y_re(v), 3.8, fill=PURPLE))
        p.append(txt(x2, Y_re(max(v2)) - 7, f"{m2:.2f}", 10, PURPLE, anchor="middle", weight="bold"))

        p.append(txt(cx, bx_b + 18, mlabel, 11, INK, anchor="middle", weight="bold"))

    p.append(rect(554, 492, 482, 34, fill="#ffffff", stroke=GRID, rx=5))
    p.append(txt(564, 513, "Malthusian growth rate r = ln(Re)/GI held invariant:  exp1 r = 0.109/d  vs.  exp3 r = 0.100/d", 10, INK, weight="bold"))

    # --- Panel C: Cumulative Empirical E2 Seed Schedule vs. Wave Peaks (Right) ---
    cx_l, cx_r, cx_t, cx_b = 1118, 1555, 130, 450
    p.append(rect(1068, 76, 504, 460, fill=BG_PANEL, stroke=GRID, rx=8))
    p.append(txt(1084, 100, "C. E2 (Delta) Importation Schedule vs. Wave Peak", 15, INK, weight="bold"))
    p.append(txt(1084, 118, "Cumulative % of empirical E2 seeds (data/seed_schedule/*-2.csv)", 11, INK2))

    t_min, t_max = 100.0, 420.0
    def X_t(t: float) -> float:
        return cx_l + (t - t_min) / (t_max - t_min) * (cx_r - cx_l)

    def Y_pct(pct: float) -> float:
        return cx_b - (pct / 100.0) * (cx_b - cx_t)

    for pv in [0, 25, 50, 75, 100]:
        yy = Y_pct(pv)
        p.append(line(cx_l, yy, cx_r, yy, stroke=GRID, width=1))
        p.append(txt(cx_l - 8, yy + 4, f"{pv}%", 11, INK2, anchor="end"))

    for tv in [120, 180, 240, 300, 360, 420]:
        xx = X_t(tv)
        p.append(line(xx, cx_t, xx, cx_b, stroke=GRID, width=1))
        p.append(txt(xx, cx_b + 15, f"t={tv}", 10, INK2, anchor="middle"))
        p.append(txt(xx, cx_b + 27, tick_to_date_str(tv)[5:], 9, INK3, anchor="middle"))

    # Vertical bands for exp2 peak (tick 206) vs exp3/exp1 peak (tick 248 / 246)
    exp2_mean_t = stats.mean(ABM_EPI_METRICS[s]["exp2_e2_peak_tick"] for s in STATES)
    exp3_mean_t = stats.mean(ABM_EPI_METRICS[s]["e2_peak_tick"] for s in STATES)
    exp1_mean_t = stats.mean(ABM_EPI_METRICS[s]["exp1_e2_peak_tick"] for s in STATES)

    p.append(line(X_t(exp2_mean_t), cx_t, X_t(exp2_mean_t), cx_b, stroke=RED, width=2.2, dash="5,4"))
    p.append(txt(X_t(exp2_mean_t) - 6, cx_t + 22, "exp2 peak (Jun 24)", 10, RED, anchor="end", weight="bold"))
    p.append(txt(X_t(exp2_mean_t) - 6, cx_t + 35, "Only 1.6% seeds in!", 10, RED, anchor="end"))

    p.append(line(X_t(exp3_mean_t), cx_t, X_t(exp3_mean_t), cx_b, stroke=PURPLE, width=2.5))
    p.append(txt(X_t(exp3_mean_t) + 6, cx_t + 22, "exp3 peak (Aug 05)", 10, PURPLE, anchor="start", weight="bold"))
    p.append(txt(X_t(exp3_mean_t) + 6, cx_t + 35, "Matches exp1 (Aug 03)", 10, ORANGE, anchor="start", weight="bold"))
    p.append(line(X_t(exp1_mean_t), cx_t, X_t(exp1_mean_t), cx_b, stroke=ORANGE, width=1.8, dash="3,3"))

    # Plot cumulative seed curves for all 5 states
    state_cols = {"ga": BLUE, "ma": ORANGE, "mn": AQUA, "va": PURPLE, "wa": YELLOW}
    for st in STATES:
        s2 = load_seed_schedule(seed_dir, st, 2)
        tot = sum(c for _, c in s2)
        cum = 0
        pts = [(X_t(100), Y_pct(0))]
        for t, c in s2:
            if t < t_min or t > t_max:
                continue
            cum += c
            pts.append((X_t(t), Y_pct(cum / tot * 100.0)))
        path_d = "M " + " L ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
        p.append(f'<path d="{path_d}" fill="none" stroke="{state_cols[st]}" stroke-width="2.0" stroke-opacity="0.85"/>')
        # Mark state's exp3 peak on its curve
        pk_t = ABM_EPI_METRICS[st]["e2_peak_tick"]
        pk_cum = sum(c for t, c in s2 if t <= pk_t) / tot * 100.0
        p.append(circle(X_t(pk_t), Y_pct(pk_cum), 4.5, fill=PURPLE))

    p.append(rect(1084, 492, 472, 34, fill="#ffffff", stroke=GRID, rx=5))
    p.append(txt(1094, 513, "Pre-peak E2 seeds restored:  36-89 seeds (exp2) -> 307-749 seeds (exp3, 4.3x-11.7x more)", 10, INK, weight="bold"))

    p.append("</svg>\n")
    svg_path = out_dir / "fig1_exp3_epi_dynamics.svg"
    svg_path.write_text("\n".join(p))
    convert_svg_to_png(svg_path, W, H)
    return svg_path


# ---------------------------------------------------------------------------
# Figure 2: Molecular Clock Closure & Two-Way Orthogonality
# ---------------------------------------------------------------------------
def plot_fig2_clock_closure(clock_rows: list[dict], out_dir: Path) -> Path:
    W, H = 1600, 560
    p = [svg_open(W, H)]
    st_rows = [r for r in clock_rows if r["state"] in ("GA", "MA", "MN", "VA", "WA")]

    p.append(txt(32, 36, "Figure 2: PhyloGAS Molecular-Clock Closure & Two-Way Orthogonality (exp1 vs. exp3)", 20, INK, weight="bold"))
    p.append(txt(32, 58, "Identity mu_truth = (lambda / GI) * s_path verified across 14,414,476 painted transmission edges in 5 states.", 13, INK2))

    # --- Panel A: Two-Way Orthogonality in (GI, lambda) Space ---
    ax_l, ax_r, ax_t, ax_b = 85, 490, 125, 460
    p.append(rect(28, 76, 485, 460, fill=BG_PANEL, stroke=GRID, rx=8))
    p.append(txt(44, 100, "A. Two-Way Orthogonality: (GI, lambda) Space", 15, INK, weight="bold"))
    p.append(txt(44, 118, "Changing lambda leaves GI fixed; changing GI leaves lambda fixed", 11, INK2))

    gi_min, gi_max = 4.0, 7.5
    lam_min, lam_max = 0.25, 0.65
    def X_gi(g: float) -> float:
        return ax_l + (g - gi_min) / (gi_max - gi_min) * (ax_r - ax_l)
    def Y_lam(l: float) -> float:
        return ax_b - (l - lam_min) / (lam_max - lam_min) * (ax_b - ax_t)

    for gv in [4.5, 5.0, 5.5, 6.0, 6.5, 7.0]:
        xx = X_gi(gv)
        p.append(line(xx, ax_t, xx, ax_b, stroke=GRID, width=1))
        p.append(txt(xx, ax_b + 16, f"{gv:.1f} d", 10, INK2, anchor="middle"))
    p.append(txt((ax_l + ax_r) / 2, ax_b + 34, "Realised Generation Interval GI (days)", 11, INK, anchor="middle", weight="bold"))

    for lv in [0.3, 0.4, 0.5, 0.6]:
        yy = Y_lam(lv)
        p.append(line(ax_l, yy, ax_r, yy, stroke=GRID, width=1))
        p.append(txt(ax_l - 8, yy + 4, f"{lv:.2f}", 10, INK2, anchor="end"))

    # Iso-clock hyperbola for lambda/GI = 9.0e-4 subs/site/yr (26.76 /yr)
    iso_pts = []
    for g_step in range(42, 74):
        g = g_step / 10.0
        l_iso = (9.0e-4 * GENOME_LEN / DAYS_PER_YEAR) * g
        if lam_min <= l_iso <= lam_max:
            iso_pts.append((X_gi(g), Y_lam(l_iso)))
    if iso_pts:
        d_str = "M " + " L ".join(f"{x:.1f},{y:.1f}" for x, y in iso_pts)
        p.append(f'<path d="{d_str}" fill="none" stroke="{AQUA}" stroke-width="1.8" stroke-dasharray="5,4"/>')
        p.append(txt(X_gi(4.6), Y_lam(0.415), "Iso-clock: lambda/GI = 9.0e-4 /site/yr", 10, AQUA, anchor="start", weight="bold"))

    # Plot the 3 clouds: (1) exp1 Poisson, (2) exp1 Lognormal, (3) exp3 Lognormal
    p1_g = stats.mean(EXP1_POISSON[s]["gi"] for s in STATES)
    p1_l = stats.mean(EXP1_POISSON[s]["lambda"] for s in STATES)
    l1_g = stats.mean(r["exp1_gi_days"] for r in st_rows)
    l1_l = stats.mean(r["exp1_lambda"] for r in st_rows)
    l3_g = stats.mean(r["exp3_gi_days"] for r in st_rows)
    l3_l = stats.mean(r["exp3_lambda"] for r in st_rows)

    # Arrows connecting centroids
    p.append(line(X_gi(p1_g), Y_lam(p1_l), X_gi(l1_g), Y_lam(l1_l), stroke=ORANGE, width=2.5))
    p.append(txt(X_gi(p1_g) - 10, Y_lam(0.49), "Step 1: Fix lambda (-41%)", 10, ORANGE, anchor="end", weight="bold"))
    p.append(txt(X_gi(p1_g) - 10, Y_lam(0.49) + 13, "GI unchanged (6.87 d)", 10, INK2, anchor="end"))

    p.append(line(X_t_gi := X_gi(l1_g), Y_lam(l1_l), X_gi(l3_g), Y_lam(l3_l), stroke=PURPLE, width=2.5))
    p.append(txt((X_t_gi + X_gi(l3_g)) / 2, Y_lam(l3_l) - 18, "Step 2: Fix GI (-31.9%)", 10, PURPLE, anchor="middle", weight="bold"))
    p.append(txt((X_t_gi + X_gi(l3_g)) / 2, Y_lam(l3_l) - 6, "lambda invariant (0.3432)", 10, INK2, anchor="middle"))

    for s in STATES:
        p.append(circle(X_gi(EXP1_POISSON[s]["gi"]), Y_lam(EXP1_POISSON[s]["lambda"]), 5.0, fill=GRAY))
    for r in st_rows:
        p.append(circle(X_gi(r["exp1_gi_days"]), Y_lam(r["exp1_lambda"]), 5.5, fill=ORANGE))
        p.append(circle(X_gi(r["exp3_gi_days"]), Y_lam(r["exp3_lambda"]), 6.0, fill=PURPLE))

    p.append(txt(X_gi(p1_g) - 12, Y_lam(p1_l) - 8, "exp1 Poisson (0.582, 6.87d)", 10, GRAY, anchor="end", weight="bold"))
    p.append(txt(X_gi(l3_g) - 6, Y_lam(l3_l) + 18, "exp3 Lognormal (0.343, 4.68d)", 10, PURPLE, anchor="start", weight="bold"))
    p.append(txt(X_gi(l1_g) + 6, Y_lam(l1_l) + 32, "exp1 Lognormal (0.343, 6.87d)", 10, ORANGE, anchor="end", weight="bold"))

    # --- Panel B: Clock Rate Hierarchy Across 5 States ---
    bx_l, bx_r, bx_t, bx_b = 585, 1085, 130, 450
    p.append(rect(528, 76, 575, 460, fill=BG_PANEL, stroke=GRID, rx=8))
    p.append(txt(544, 100, "B. Molecular Clock Rates vs. Real SARS-CoV-2 (8.0e-4 /site/yr)", 15, INK, weight="bold"))
    p.append(txt(544, 118, "Comparison across GA, MA, MN, VA, WA (1S__surs baseline for TreeTime)", 11, INK2))

    c_min, c_max = 3.5e-4, 11.0e-4
    def Y_ck(v: float) -> float:
        return bx_b - (v - c_min) / (c_max - c_min) * (bx_b - bx_t)

    for cv in [4e-4, 6e-4, 8e-4, 10e-4]:
        yy = Y_ck(cv)
        is_sc2 = abs(cv - 8e-4) < 1e-6
        p.append(line(bx_l, yy, bx_r, yy, stroke=RED if is_sc2 else GRID, width=2.0 if is_sc2 else 1.0, dash="6,4" if is_sc2 else ""))
        p.append(txt(bx_l - 8, yy + 4, f"{cv*1e4:.1f}e-4", 10, RED if is_sc2 else INK2, anchor="end", weight="bold" if is_sc2 else "normal"))
    p.append(txt(bx_l + 6, Y_ck(8e-4) - 7, "Real SARS-CoV-2 = 8.0e-4 /site/yr", 10, RED, anchor="start", weight="bold"))

    clock_series = [
        ("exp1 Truth", "exp1_mu_truth", ORANGE),
        ("exp3 Implied", "exp3_implied_per_site_yr", AQUA),
        ("exp3 Lineage", "exp3_mu_lineage", BLUE),
        ("exp3 Truth", "exp3_mu_truth", PURPLE),
        ("exp3 TreeTime", "exp3_mu_sim_surs", YELLOW),
    ]
    col_w = (bx_r - bx_l) / len(clock_series)
    for idx, (clabel, key, col) in enumerate(clock_series):
        cx = bx_l + (idx + 0.5) * col_w
        vals = [r[key] for r in st_rows]
        mv = stats.mean(vals)
        yy = Y_ck(mv)
        p.append(rect(cx - 20, yy, 40, bx_b - yy, fill=col, opacity=0.26, stroke=col, width=1.8))
        for si, r in enumerate(st_rows):
            v = r[key]
            jx = cx + (si - 2) * 5.5
            p.append(circle(jx, Y_ck(v), 4.2, fill=col))
        p.append(txt(cx, Y_ck(max(vals)) - 10, f"{mv*1e4:.2f}e-4", 10, INK, anchor="middle", weight="bold"))
        p.append(txt(cx, bx_b + 18, clabel, 10, INK, anchor="middle", weight="bold"))

    p.append(rect(544, 492, 542, 34, fill="#ffffff", stroke=GRID, rx=5))
    p.append(txt(554, 513, "Pooled Truth Clock:  5.42e-4 (0.68x, exp1) -> 7.60e-4 (0.95x, exp3; 0.99x in GA/MA/WA)", 10, INK, weight="bold"))

    # --- Panel C: Root-to-Tip Regression R^2 Jump (exp1 -> exp3) ---
    cx_l, cx_r, cx_t, cx_b = 1170, 1545, 130, 450
    p.append(rect(1118, 76, 454, 460, fill=BG_PANEL, stroke=GRID, rx=8))
    p.append(txt(1134, 100, "C. Root-to-Tip Temporal Signal (R^2)", 15, INK, weight="bold"))
    p.append(txt(1134, 118, "Faster generation turnover boosts R^2 by +61%", 11, INK2))

    r2_min, r2_max = 0.15, 0.70
    def Y_r2(v: float) -> float:
        return cx_b - (v - r2_min) / (r2_max - r2_min) * (cx_b - cx_t)

    for rv in [0.2, 0.3, 0.4, 0.5, 0.6]:
        yy = Y_r2(rv)
        p.append(line(cx_l, yy, cx_r, yy, stroke=GRID, width=1))
        p.append(txt(cx_l - 8, yy + 4, f"{rv:.1f}", 10, INK2, anchor="end"))

    st_w = (cx_r - cx_l) / len(st_rows)
    for si, r in enumerate(st_rows):
        sx = cx_l + (si + 0.5) * st_w
        y1 = Y_r2(r["exp1_r_squared"])
        y3 = Y_r2(r["exp3_r_squared"])
        p.append(line(sx, y1, sx, y3, stroke=PURPLE, width=2.5))
        p.append(circle(sx, y1, 5.0, fill=ORANGE))
        p.append(circle(sx, y3, 5.5, fill=PURPLE))
        p.append(txt(sx, y1 + 15, f"{r['exp1_r_squared']:.2f}", 9, ORANGE, anchor="middle", weight="bold"))
        p.append(txt(sx, y3 - 9, f"{r['exp3_r_squared']:.2f}", 10, PURPLE, anchor="middle", weight="bold"))
        p.append(txt(sx, cx_b + 18, r["state"], 11, INK, anchor="middle", weight="bold"))

    p.append(rect(1134, 492, 422, 34, fill="#ffffff", stroke=GRID, rx=5))
    p.append(txt(1144, 513, "5-State Mean R^2:  0.314 (exp1) -> 0.505 (exp3) | Within-date var: 95.0%", 10, INK, weight="bold"))

    p.append("</svg>\n")
    svg_path = out_dir / "fig2_exp3_clock_closure.svg"
    svg_path.write_text("\n".join(p))
    convert_svg_to_png(svg_path, W, H)
    return svg_path


# ---------------------------------------------------------------------------
# Figure 3: Surveillance & Phylogeographic Benchmarks Across 8 Strategies
# ---------------------------------------------------------------------------
def plot_fig3_surveillance(mug_rows: list[dict], auc_rows: list[dict], out_dir: Path) -> Path:
    W, H = 1600, 600
    p = [svg_open(W, H)]

    p.append(txt(32, 36, "Figure 3: Phylodynamic & Surveillance Benchmarks — LASSO50-Greedy vs. Stride-Matched SURS", 20, INK, weight="bold"))
    p.append(txt(32, 58, "LASSO50-Greedy (50 macro-subgroups) beats SURS across 10 state-regime runs (exp1+exp3) and cuts exp1<->exp3 flow regime drift |ΔCos| by up to 3.9x.", 12.5, INK2))

    # --- Panel A: Inter-County Flow Topology (F1) & Proportions (Cosine) ---
    ax_l, ax_r, ax_t, ax_b = 222, 535, 132, 485
    p.append(rect(28, 76, 765, 500, fill=BG_PANEL, stroke=GRID, rx=8))
    p.append(txt(44, 100, "A. Phylogeographic Flow Reconstruction: LASSO50-Greedy vs. SURS", 14.2, INK, weight="bold"))
    p.append(txt(44, 118, "exp3 state dots & 5-state mean bars: Topological F1 (purple) & Cosine (blue)", 10.2, INK2))
    p.append(txt(592, 118, "Mean F1 / Cos", 9.6, INK, anchor="middle", weight="bold"))
    p.append(txt(712, 118, "vs. SURS (exp3 | 10-Run)", 9.6, INK, anchor="middle", weight="bold"))

    x_min, x_max = 0.35, 0.95
    def X_m(v: float) -> float:
        return ax_l + (v - x_min) / (x_max - x_min) * (ax_r - ax_l)

    for xv in [0.40, 0.50, 0.60, 0.70, 0.80, 0.90]:
        xx = X_m(xv)
        p.append(line(xx, ax_t, xx, ax_b, stroke=GRID, width=1))
        p.append(txt(xx, ax_b + 16, f"{xv:.2f}", 10, INK2, anchor="middle"))

    row_h = (ax_b - ax_t) / len(mug_rows)
    for idx, r in enumerate(mug_rows):
        cy = ax_t + (idx + 0.5) * row_h
        is_base = r["scheme"] in ("1S__surs", "4S__surs")
        is_top = idx < 3
        if is_top:
            p.append(rect(36, cy - row_h * 0.44, 748, row_h * 0.88, fill="#eefaf6" if "LL-P" in r["scheme"] else "#f3f0ff", opacity=0.65, rx=4))
        p.append(line(ax_l, cy, ax_r, cy, stroke=GRID, width=0.8))
        lbl_col = PURPLE if is_top else (INK2 if is_base else INK)
        p.append(txt(ax_l - 8, cy + 4, str(r["display_name"]), 9.8, lbl_col, anchor="end", weight="bold" if (is_top or is_base) else "normal"))

        # State dots + 5-state mean markers
        for s in STATES:
            p.append(circle(X_m(r[f"{s}_cosine"]), cy - 5, 3.5, fill=BLUE, opacity=0.65))
            p.append(circle(X_m(r[f"{s}_f1"]), cy + 5, 3.5, fill=PURPLE, opacity=0.65))
        p.append(rect(X_m(r["mean_cosine"]) - 3, cy - 11, 6, 11, fill=BLUE, stroke="#ffffff", width=1))
        p.append(rect(X_m(r["mean_topological_f1"]) - 3, cy, 6, 11, fill=PURPLE, stroke="#ffffff", width=1))

        # Dedicated Mean F1 / Cos column (no overlap with state dots!)
        p.append(txt(592, cy + 4, f"{r['mean_topological_f1']:.3f} / {r['mean_cosine']:.3f}", 9.4, INK, anchor="middle", weight="bold"))

        # Right-side badge: F1 gain (exp3/5 | 10-run/10)
        if is_base:
            p.append(txt(712, cy + 4, f"Baseline ({r['scheme'][:2]} SURS)", 9.4, INK2, anchor="middle", weight="bold"))
        else:
            f1_g = float(r["f1_gain_vs_surs_pct"])
            f1_w = int(r["f1_wins_vs_surs"])
            f1_w10 = int(r["f1_wins_10runs"])
            cos_g = float(r["cos_gain_vs_surs_pct"])
            b_col = AQUA if (f1_g > 0 and cos_g > 0) else (PURPLE if f1_g > 0 else RED)
            win_dot = "●" if f1_w10 >= 9 else ""
            p.append(txt(712, cy + 4, f"F1 {f1_g:+.1f}% ({f1_w}/5, {f1_w10}/10{win_dot})", 9.2, b_col, anchor="middle", weight="bold"))

    p.append(rect(42, 512, 736, 52, fill="#ffffff", stroke=GRID, rx=5))
    p.append(txt(52, 531, "LASSO50-Greedy (LL-P):", 9.4, AQUA, weight="bold"))
    p.append(txt(172, 531, "1S-1_LL-P wins 10/10 runs on F1 (+5.4%) & 8/10 on Cosine (+1.1%) with highest 10-run floor (0.591).", 9.1, INK))
    p.append(txt(52, 550, "LASSO50-Greedy (P):", 9.4, PURPLE, weight="bold"))
    p.append(txt(172, 550, "4S-P & 1S-P win 9/10 runs on F1 (13/14 incl. preQC), discovering rural-urban bridges SURS misses.", 9.1, INK))

    # --- Panel B: Multi-Objective & Cross-Regime Resilience Scorecard ---
    p.append(rect(812, 76, 760, 500, fill=BG_PANEL, stroke=GRID, rx=8))
    p.append(txt(828, 100, "B. Multi-Objective & Cross-Regime Resilience Scorecard (LASSO50-Greedy vs. SURS)", 14.0, INK, weight="bold"))
    p.append(txt(828, 118, "exp3 gain (%), 10-run win record (exp1+exp3, ● >= 9/10), and exp1<->exp3 regime drift |ΔCos|", 10.2, INK2))

    # Table columns inside Panel B
    headers = [
        ("Recipe (LASSO50 / SURS)", 958, "end"),
        ("Flow F1 (10-Run)", 1052, "middle"),
        ("Clusters (10-Run)", 1192, "middle"),
        ("Senior 65+ MRD", 1334, "middle"),
        ("Regime Drift |ΔCos|", 1482, "middle"),
    ]
    for h_txt, hx, h_anc in headers:
        p.append(txt(hx, 142, h_txt, 9.8, INK, anchor=h_anc, weight="bold"))
    p.append(line(826, 150, 1558, 150, stroke=INK2, width=1.2))

    auc_by_sc = {str(r["scheme"]): r for r in auc_rows}
    base_drift = {
        "1S__surs": next(float(r["mean_abs_delta_cos_exp1_exp3"]) for r in mug_rows if r["scheme"] == "1S__surs"),
        "4S__surs": next(float(r["mean_abs_delta_cos_exp1_exp3"]) for r in mug_rows if r["scheme"] == "4S__surs"),
    }
    for idx, mr in enumerate(mug_rows):
        sc = str(mr["scheme"])
        ar = auc_by_sc[sc]
        cy = ax_t + 26 + (idx + 0.5) * ((ax_b - ax_t - 20) / len(mug_rows))
        is_base = sc in ("1S__surs", "4S__surs")
        p.append(line(826, cy + 14, 1558, cy + 14, stroke=GRID, width=0.8))
        short_lbl = sc.replace("__lasso_greedy", " (L50-G)").replace("__surs", " (SURS)")
        p.append(txt(958, cy + 4, short_lbl, 9.6, INK, anchor="end", weight="bold" if not is_base else "normal"))

        drift = float(mr["mean_abs_delta_cos_exp1_exp3"])
        if is_base:
            p.append(txt(1052, cy + 4, f"{mr['mean_topological_f1']:.3f} (Base)", 9.4, INK2, anchor="middle"))
            p.append(txt(1192, cy + 4, f"{ar['stride_component_coverage_mean']:.3f} (Base)", 9.4, INK2, anchor="middle"))
            p.append(txt(1334, cy + 4, f"{ar['equity_Senior_65+_mean']:.4f} (Base)", 9.4, INK2, anchor="middle"))
            d_col = RED if sc == "4S__surs" else ORANGE
            d_tag = "Fragile" if sc == "4S__surs" else "Base 1S"
            p.append(rect(1414, cy - 11, 136, 20, fill=d_col, opacity=0.12, rx=4))
            p.append(txt(1482, cy + 3, f"{drift:.3f} ({d_tag})", 9.4, d_col, anchor="middle", weight="bold"))
            continue

        # 1. Flow Topology F1 (exp3 % + 10-run wins)
        f1_g, f1_w10 = float(mr["f1_gain_vs_surs_pct"]), int(mr["f1_wins_10runs"])
        c1 = AQUA if f1_w10 >= 9 else (PURPLE if f1_g > 0 else RED)
        p.append(rect(989, cy - 11, 126, 20, fill=c1, opacity=0.13 if f1_g > 0 else 0.08, rx=4))
        p.append(txt(1052, cy + 3, f"{f1_g:+.1f}% ({f1_w10}/10{'●' if f1_w10>=9 else ''})", 9.5, c1, anchor="middle", weight="bold"))

        # 2. Clusters Reached (exp3 % + 10-run wins)
        cl_g = float(ar["stride_component_coverage_gain_vs_surs_pct"])
        cl_w10 = int(ar["stride_component_coverage_wins_10runs"])
        c2 = AQUA if cl_w10 >= 9 else (BLUE if cl_g > 0 else RED)
        p.append(rect(1129, cy - 11, 126, 20, fill=c2, opacity=0.13 if cl_g > 0 else 0.08, rx=4))
        p.append(txt(1192, cy + 3, f"{cl_g:+.1f}% ({cl_w10}/10{'●' if cl_w10>=9 else ''})", 9.5, c2, anchor="middle", weight="bold"))

        # 3. Senior 65+ Equity Tree Coverage MRD (exp3 % + 10-run wins, higher MRD = closer tree distance)
        sr_g = float(ar["equity_Senior_65+_gain_vs_surs_pct"])
        sr_w10 = int(ar["equity_Senior_65+_wins_10runs"])
        c3 = AQUA if sr_w10 >= 9 else (BLUE if sr_w10 >= 7 else RED)
        p.append(rect(1271, cy - 11, 126, 20, fill=c3, opacity=0.13 if sr_w10 >= 9 else 0.08, rx=4))
        p.append(txt(1334, cy + 3, f"{sr_g:+.1f}% ({sr_w10}/10{'●' if sr_w10>=9 else ''})", 9.5, c3, anchor="middle", weight="bold"))

        # 4. Regime Drift |Delta Cos| (exp1 <-> exp3)
        b_d = base_drift[str(mr["stride_baseline"])]
        ratio = b_d / drift if drift > 0 else 1.0
        c4 = AQUA if drift <= 0.030 else (BLUE if drift < b_d else ORANGE)
        p.append(rect(1414, cy - 11, 136, 20, fill=c4, opacity=0.13 if drift <= 0.035 else 0.08, rx=4))
        p.append(txt(1482, cy + 3, f"{drift:.3f} ({ratio:.1f}x stable)", 9.4, c4, anchor="middle", weight="bold"))

    p.append(rect(826, 512, 732, 52, fill="#ffffff", stroke=GRID, rx=5))
    p.append(txt(838, 531, "Cross-Regime Wins (●):", 9.8, AQUA, weight="bold"))
    p.append(txt(972, 531, "LASSO50-Greedy (P & LL-P) wins 10/10 runs on Senior 65+ & Adult Tree MRD and 9-10/10 on Flow F1 & Clusters.", 9.3, INK))
    p.append(txt(838, 550, "Regime Resilience:", 9.8, ORANGE, weight="bold"))
    p.append(txt(955, 550, "SURS · 4S drifts by |ΔCos| = 0.096 between exp1 & exp3; LASSO50-Greedy · 4S-P drifts by just 0.025 (3.9x more stable).", 9.3, INK))

    p.append("</svg>\n")
    svg_path = out_dir / "fig3_exp3_surveillance_benchmarks.svg"
    svg_path.write_text("\n".join(p))
    convert_svg_to_png(svg_path, W, H)
    return svg_path


# ---------------------------------------------------------------------------
# Figure 4: 5-State Timelines (Epicurves / Seed Schedules / Vax / Peak Shift)
# ---------------------------------------------------------------------------
def plot_fig4_state_timelines(seed_dir: Path, vax_dir: Path, incidence_cache_dir: Path | None, out_dir: Path) -> Path:
    W, H = 1600, 860
    p = [svg_open(W, H)]

    has_daily = bool(
        incidence_cache_dir
        and all((incidence_cache_dir / f"exp3_0.25_{st}_rep0_daily.csv").exists() for st in STATES)
    )
    subtitle = (
        "Daily E1 (Alpha) & E2 (Delta) exposures, empirical importation schedules, 300-day painted window, and wave peak alignment (cell 0.25)."
        if has_daily
        else "Empirical E1 & E2 importation schedules (bars), vaccination rollout, 300-day painted window, and E1/E2 wave peak alignment (exp1 vs. exp2 vs. exp3, cell 0.25)."
    )

    p.append(txt(32, 36, "Figure 4: 5-State Epidemic Timelines, Empirical Seeding Schedules, and Wave Peak Alignment (exp3 vs. exp2 & exp1)", 20, INK, weight="bold"))
    p.append(txt(32, 58, subtitle, 13, INK2))

    ax_l, ax_r = 125, 1545
    t_min, t_max = 0.0, 500.0
    def X_t(t: float) -> float:
        return ax_l + (t - t_min) / (t_max - t_min) * (ax_r - ax_l)

    panel_h = 132
    gap_h = 18
    top_start = 82

    for idx, st in enumerate(STATES):
        pt = top_start + idx * (panel_h + gap_h)
        pb = pt + panel_h
        m = ABM_EPI_METRICS[st]

        p.append(rect(28, pt, W - 56, panel_h, fill=BG_PANEL, stroke=GRID, rx=6))
        # Shaded 300-day painted window for E2 (first_v2_seed to first_v2_seed + 300)
        s2 = load_seed_schedule(seed_dir, st, 2)
        s1 = load_seed_schedule(seed_dir, st, 1)
        first_v2 = s2[0][0] if s2 else 125
        win_end = min(500, first_v2 + 300)
        p.append(rect(X_t(first_v2), pt + 6, X_t(win_end) - X_t(first_v2), panel_h - 12, fill=PURPLE, opacity=0.07))
        if idx == 0:
            p.append(txt(X_t(win_end) - 8, pt + 20, "Shaded: 300-day PhyloGAS E2 painted window", 10, PURPLE, anchor="end", weight="bold"))

        # State title & summary metrics on left
        p.append(txt(40, pt + 28, STATE_NAMES[st], 13, INK, weight="bold"))
        p.append(txt(40, pt + 48, f"E1 GI: {m['e1_gi']:.2f}d", 10, BLUE, weight="bold"))
        p.append(txt(40, pt + 64, f"E2 GI: {m['e2_gi']:.2f}d", 10, PURPLE, weight="bold"))
        p.append(txt(40, pt + 80, f"E2 Tot: {m['e2_total']/1e6:.2f}M", 10, INK2))
        p.append(txt(40, pt + 96, f"Peak: {m['e2_peak_cases']/1e3:.1f}k/d", 10, INK2))

        # Plot daily E1 & E2 incidence curves if cache CSV exists
        plot_t, plot_b = pt + 14, pb - 20
        if has_daily and incidence_cache_dir is not None:
            cache_csv = incidence_cache_dir / f"exp3_0.25_{st}_rep0_daily.csv"
            d_rows = list(csv.DictReader(open(cache_csv)))
            max_c = max(max(int(r["E1"]), int(r["E2"])) for r in d_rows) or 1
            for key, col, dash in [("E1", BLUE, "4,3"), ("E2", PURPLE, "")]:
                pts = [
                    (X_t(int(r["tick"])), plot_b - (int(r[key]) / max_c) * (plot_b - plot_t))
                    for r in d_rows
                    if int(r[key]) > 0
                ]
                if pts:
                    d_str = "M " + " L ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
                    d_attr = f' stroke-dasharray="{dash}"' if dash else ""
                    p.append(f'<path d="{d_str}" fill="none" stroke="{col}" stroke-width="2.0"{d_attr}/>')

        # Always plot E1 & E2 empirical seed schedule bars along the baseline
        max_seed = max([c for _, c in s2] + [c for _, c in s1] + [1])
        for t, c in s1:
            if t <= t_max:
                bh = (c / max_seed) * (plot_b - plot_t) * 0.65
                p.append(line(X_t(t), plot_b, X_t(t), plot_b - bh, stroke=BLUE, width=1.8, opacity=0.65))
        for t, c in s2:
            if t <= t_max:
                bh = (c / max_seed) * (plot_b - plot_t) * 0.85
                p.append(line(X_t(t), plot_b, X_t(t), plot_b - bh, stroke=PURPLE, width=1.8, opacity=0.75))

        # Vaccination cumulative rollout curve (in green)
        vax = load_vax_schedule(vax_dir, st)
        if vax:
            tot_v = sum(c for _, c in vax) or 1
            cum_v = 0
            v_pts = []
            for t, c in vax:
                if t <= t_max:
                    cum_v += c
                    v_pts.append((X_t(t), plot_b - (cum_v / tot_v) * (plot_b - plot_t) * 0.75))
            if v_pts:
                d_str = "M " + " L ".join(f"{x:.1f},{y:.1f}" for x, y in v_pts)
                p.append(f'<path d="{d_str}" fill="none" stroke="{AQUA}" stroke-width="1.6" stroke-dasharray="3,3"/>')

        # Peak markers for E1 (exp3) and E2 (exp2 vs exp1 vs exp3)
        x_e1_3 = X_t(m["e1_peak_tick"])
        p.append(line(x_e1_3, plot_t, x_e1_3, plot_b, stroke=BLUE, width=2.0))
        p.append(txt(x_e1_3 + 4, plot_t + 10, f"exp3 E1 peak: t={m['e1_peak_tick']} ({tick_to_date_str(m['e1_peak_tick'])[5:]})", 9, BLUE, weight="bold"))

        x_e2_2 = X_t(m["exp2_e2_peak_tick"])
        p.append(line(x_e2_2, plot_t, x_e2_2, plot_b, stroke=RED, width=1.8, dash="4,3"))
        p.append(txt(x_e2_2 - 4, plot_t + 24, f"exp2 E2 peak: t={m['exp2_e2_peak_tick']} (premature)", 9, RED, anchor="end", weight="bold"))

        x_e2_1 = X_t(m["exp1_e2_peak_tick"])
        p.append(line(x_e2_1, plot_t, x_e2_1, plot_b, stroke=ORANGE, width=1.8, dash="2,2"))

        x_e2_3 = X_t(m["e2_peak_tick"])
        p.append(line(x_e2_3, plot_t, x_e2_3, plot_b, stroke=PURPLE, width=2.5))
        p.append(txt(x_e2_3 + 5, plot_t + 10, f"exp3 E2 peak: t={m['e2_peak_tick']} ({tick_to_date_str(m['e2_peak_tick'])[5:]}) | exp1: t={m['exp1_e2_peak_tick']}", 9, PURPLE, weight="bold"))

        # Tick gridlines
        for tv in range(0, 501, 50):
            xx = X_t(tv)
            p.append(line(xx, plot_b, xx, plot_b + 4, stroke=INK2, width=1))
            if idx == len(STATES) - 1:
                p.append(txt(xx, pb + 14, f"Tick {tv} ({tick_to_date_str(tv)})", 10, INK2, anchor="middle"))

    p.append("</svg>\n")
    svg_path = out_dir / "fig4_exp3_state_timelines.svg"
    svg_path.write_text("\n".join(p))
    convert_svg_to_png(svg_path, W, H)
    return svg_path


def build_summary_readme(out_dir: Path) -> Path:
    """Assemble a single unified Markdown report in docs/characterization/README.md
    embedding all four SVG figures and all four Markdown tables."""
    md1 = (out_dir / "exp3_epi_characterization.md").read_text().strip()
    md2 = (out_dir / "exp3_molecular_clock_closure.md").read_text().strip()
    md3 = (out_dir / "exp3_mugration_benchmarks.md").read_text().strip()
    md4 = (out_dir / "exp3_auc_truth_benchmarks.md").read_text().strip()

    report = f"""# PhyloGAS `exp3` Multi-State Characterization & Benchmark Report

> **Backing Script:** Generated reproducibly by [`scripts/utils/characterize_exp3.py`](../../scripts/utils/characterize_exp3.py).
> **Canonical Metrics & Recipes Reference:** **[`docs/benchmarks_and_metrics.md`](../benchmarks_and_metrics.md)** (full mathematical definitions), **[`docs/clock_modes.md`](../clock_modes.md)** (clock decomposition), and **[`docs/abm_contribution.md`](../abm_contribution.md)** (two-way orthogonality identity).
> **ABM Ground-Truth Deposit (v2 / `exp3`):** [`10.5281/zenodo.23989650`](https://doi.org/10.5281/zenodo.23989650) (Concept DOI: [`10.5281/zenodo.23067669`](https://doi.org/10.5281/zenodo.23067669))
> **Downstream PhyloGAS Outputs & Benchmarks Deposit (`exp3`):** [`10.5281/zenodo.24013208`](https://doi.org/10.5281/zenodo.24013208)

---

## 1. Epidemiological & Waning Calibration (`exp1` vs. `exp2` vs. `exp3`)

![Figure 1: Epidemiological & Waning Calibration](fig1_exp3_epi_dynamics.svg)

### Per-State Importation, Vaccination, & Peak Alignment Timelines

![Figure 4: State Timelines](fig4_exp3_state_timelines.svg)

{md1}

> **Section 1 Metric Legend (see [`docs/benchmarks_and_metrics.md` §2](../benchmarks_and_metrics.md#2-abm-epidemiological--waning-characterization-metrics)):**
> - **`GI (d)`:** Realized mean generation interval $t_{{\\text{{exposed}}}}(\\text{{child}}) - t_{{\\text{{exposed}}}}(\\text{{parent}})$ across all ABM transmission edges (`5.22 d` for `E1`, `4.70 d` for `E2` in `exp3` vs. `6.87 d` in `exp1`).
> - **$R_0(1\\text{{k}})$ / $R_0(10\\text{{k}})$:** Mean secondary infections per index case across the first $1{{,}}000$ and $10{{,}}000$ wave exposures.
> - **`E2 Seeds <= Peak`:** Empirical `E2` importation events arriving on or before the simulated `E2` daily incidence peak (`13.0%` in `exp3` vs. `1.6%` in `exp2`).
> - **`Window Reinf %` & `Gaps <= 30d / <= 60d`:** Share of painted-window `E2` cases that are repeat infections on the same host, and the fraction of repeat infections occurring within $\\le 30\\text{{d}}$ (`0.03%` under `Gamma(4, 45.75)` waning in `exp3`) or $\\le 60\\text{{d}}$ (`0.56%`).

---

## 2. Molecular-Clock Closure & Two-Way Orthogonality

![Figure 2: Molecular Clock Closure](fig2_exp3_clock_closure.svg)

{md2}

> **Section 2 Metric Legend (see [`docs/benchmarks_and_metrics.md` §3](../benchmarks_and_metrics.md#3-molecular-clock-closure--sequence-divergence-metrics), [`docs/clock_modes.md`](../clock_modes.md), and [`docs/abm_contribution.md`](../abm_contribution.md)):**
> - **Closure Identity:** $\\mu_{{\\text{{truth}}}} = \\frac{{\\hat{{\\lambda}}}}{{\\text{{GI}}}} \\cdot \\frac{{365.25}}{{L}} \\cdot s_{{\\text{{path}}}}$, where $\\hat{{\\lambda}}$ is the per-transmission substitution count owned by the Genetic Painter (`0.3432`), $\\text{{GI}}$ is the generation interval owned by EpiHiper (`4.68 d`), and $s_{{\\text{{path}}}}$ (`0.844`) is the emergent path-length/coexistence factor.
> - **$\\mu_{{\\text{{lineage}}}}$ vs. $\\mu_{{\\text{{truth}}}}$ vs. $\\mu_{{\\text{{sim}}}}$:** $\\mu_{{\\text{{lineage}}}}$ is the fixed-effects within-chain root-to-tip slope; $\\mu_{{\\text{{truth}}}}$ is the pooled root-to-tip OLS slope directly on the painted FASTA (`7.60e-4` subs/site/yr, `0.95x` of the canonical SARS-CoV-2 `8.0e-4` clock); $\\mu_{{\\text{{sim}}}}$ is TreeTime's unconstrained clock estimate on the reconstructed tree.

---

## 3. Phylodynamic & Surveillance Sampling Benchmarks (`LASSO50-Greedy` vs. `SURS`)

![Figure 3: Surveillance & Mugration Benchmarks](fig3_exp3_surveillance_benchmarks.svg)

### Nomenclature Note: `LASSO50-Greedy` vs. Alternative Allocation Families

Throughout these benchmarks, **`LASSO50-Greedy`** denotes our default LASSO regression + Greedy submodular optimization allocator configured with **`50` LASSO macro-subgroups** (distinct from `LASSO50-Stratified`, fine-grained `LASSO1500-Greedy`, proportional `Stratified` quota sampling, and Simple Uniform Random Sampling `SURS`; see [`docs/benchmarks_and_metrics.md` §1](../benchmarks_and_metrics.md#1-surveillance-sampling-recipes-scenario__algorithm)). Each recipe combines an allocation algorithm with a target scenario (`<Stride>-<Target>`):
- `1S` / `4S`: 1-week vs. 4-week temporal pooling stride.
- `P`: Census population target distribution (`Pop`).
- `1_LL` / `4_LL`: Empirical clinical line-list incidence target (`Line List`).
- `1_LL-P` / `4_LL-P`: 50/50 blended line-list incidence + census population target.

### Key Insights: Cross-Regime Resilience (`exp1` vs. `exp3`) & The Complementary Strengths of `LASSO50-Greedy`

1. **Cross-Regime Resilience (`|ΔCos|` Across Epidemic Regimes) — Why `SURS` Is Regime-Fragile While `LASSO50-Greedy` Is Regime-Stable:**
   - Because Simple Uniform Random Sampling (`SURS`) passively mirrors whatever spatiotemporal distortion the epidemic + clinical ascertainment filter produces, **`SURS · 4S` suffers a massive per-state flow reconstruction swing (`mean |ΔCos| = 0.096`, max `0.109`)** when disease dynamics change from exponential waning (`exp1`, where rapid urban reinfection loops dominate 52–59% of cases) to Gamma waning (`exp3`).
   - By contrast, **`LASSO50-Greedy` actively enforces demographic and geographic macro-subgroup balance every stride**, making it **2.8x to 3.9x more resilient to epidemic regime shifts**:
     - **`LASSO50-Greedy · 4S-P`**: `mean |ΔCos| = 0.025` (**3.9x lower regime drift** than `SURS · 4S`).
     - **`LASSO50-Greedy · 1S-P`**: `mean |ΔCos| = 0.028` (5-state mean Cosine is `0.735` in `exp1` and `0.735` in `exp3` — virtually invariant).
     - **`LASSO50-Greedy · 1S-1_LL-P`**: `mean |ΔCos| = 0.029`, achieving the **highest worst-state Cosine floor (`0.591`)** across all 10 state-regime runs (`exp1` + `exp3`).
     - **`LASSO50-Greedy · 4S-4_LL-P`**: `mean |ΔCos| = 0.035` (**2.8x lower regime drift** than `SURS · 4S`, beating `SURS · 4S` on mean Cosine in both `exp1` and `exp3`).
2. **10-Run Cross-Experiment Dominance (`exp1` 5 States + `exp3` 5 States = 10 State-Regime Runs; 14 Runs Including `preQC`):**
   - **Flow Topology ($F_1$, Network Edge Recovery, $\\uparrow$):** `LASSO50-Greedy · 1S-1_LL-P` beats `SURS · 1S` in **10 / 10 runs** (**13 / 14** including `preQC`); `LASSO50-Greedy · 4S-P` and `LASSO50-Greedy · 1S-P` each win **9 / 10 runs** (**13 / 14** including `preQC`).
   - **Flow Proportions (`Cosine`, $\\uparrow$):** `LASSO50-Greedy · 1S-1_LL-P` beats `SURS · 1S` in **8 / 10 runs** (**10 / 14** including `preQC`) and ranks **#1 overall** in `exp3` on both Flow Topology ($F_1 = 0.614$) and Flow Proportions (Cosine $= 0.755$).
   - **Transmission Clusters Reached (`stride_component_coverage`, $\\uparrow$):** `LASSO50-Greedy · 1S-P` wins **10 / 10 runs** (`+4.4%` in `exp3`); `LASSO50-Greedy · 4S-P` (`+10.9%`), `4S-4_LL-P` (`+8.3%`), and `1S-1_LL-P` (`+3.2%`) each win **9 / 10 runs** (and **5 / 5 states** in `exp3`).
   - **Demographic Equity Tree Coverage (`equity_Senior_65+` & `equity_Adult_18-49` MRD, $\\uparrow$):** All four population-aware `LASSO50-Greedy` recipes (`4S-P` `+6.9%`, `1S-P` `+4.9%`, `4S-4_LL-P` `+3.7%`, and `1S-1_LL-P` `+2.9%`) win **10 / 10 runs (`5/5` in `exp1` + `5/5` in `exp3`)** on Senior 65+ and Adult 18–49 phylogenetic tree proximity, while blended `LL-P` (`1S-1_LL-P` and `4S-4_LL-P`) also wins **8 / 10 runs** on full-wave tree coverage $C_0$ MRD.
3. **Why Blended `LL-P` Outperforms Pure `P` on Cosine While Pure `P` Maximizes Network & Cluster Discovery:**
   - Nextstrain's symmetrized rate matrix $W_{{ij}} = \\frac{{1}}{{2\\mu}}(q_{{ij}}/\\pi_j + q_{{ji}}/\\pi_i)$ divides transition counts $q_{{ij}}$ by the equilibrium node share $\\pi_j$. Pure census targeting (`4S-P`) forces samples into low-incidence rural counties ($\\pi_j \\ll 1$), discovering more true county-to-county introductions (**+5.6% $F_1$, +10.9% Clusters Reached**) but amplifying $q_{{ij}}/\\pi_j$ intensity weights on sparse rural links when `exp3` has no urban reinfection trap. Blending census diversity with line-list incidence (`LL-P`) prevents rural $1/\\pi_j$ weight inflation while preserving the **5/5 state sweep** on topological flow recovery, cluster discovery, and Senior/Adult tree coverage.

{md3}

---

{md4}

> **Section 3 Metric Legend (see [`docs/benchmarks_and_metrics.md` §§4–5](../benchmarks_and_metrics.md#4-phylogeographic-flow-reconstruction-metrics-mugration_metricscsv)):**
> - **Flow Topology $F_1$ ($\\uparrow$) & Flow Cosine ($\\uparrow$):** Evaluated on the off-diagonal county-to-county transition vectors $\\mathbf{{t}}, \\mathbf{{w}} \\in \\mathbb{{R}}^{{N(N-1)}}$ ([`scoring.py`](../../src/phylogas/benchmark/scoring.py)). $F_1$ scores binary recovery of active inter-county transmission links; Cosine scores agreement of relative inter-county flow intensities.
> - **Regime Drift $|\\Delta\\text{{Cos}}|$ ($\\downarrow$) & 10-Run Cosine Floor ($\\uparrow$):** Mean per-state shift $|\\text{{Cos}}_{{s,\\text{{exp3}}}} - \\text{{Cos}}_{{s,\\text{{exp1}}}}|$ across the 5 states, and the minimum single-state Cosine across all 10 post-QC state-regime runs.
> - **Clusters Reached (`stride_component_coverage`, $\\uparrow$):** Fraction of active ground-truth imported transmission chains per stride captured by at least one sequenced genome ([`truth_runner.py`](../../src/phylogas/benchmark/truth_runner.py)).
> - **Tree Coverage & Demographic Equity MRD (`coverage_size_0`, `8_week_rolling_tree_coverage`, `equity_*`, $\\uparrow$):** Span-normalized AUC of **Mean Reciprocal Distance** $\\text{{MRD}}(P_t, S) = \\frac{{1}}{{|P_t|}}\\sum_{{u \\in P_t}} \\frac{{1}}{{d_G(u, S) + 1}}$ along the true ABM transmission graph ($1.0$ = every infection directly sequenced; higher = closer phylogenetic proximity).
> - **Variant Prevalence $L_1$ (`stride_variant_prevalence_error`, $\\downarrow$) & Cumul. Infections KL (`cumulative_infections`, $\\downarrow$):** Total variation ($L_1$) error on benchmark sub-lineage proportions and Kullback–Leibler divergence on demographic $\\times$ county infection distributions (lower = lower estimation error).
"""
    readme_path = out_dir / "README.md"
    readme_path.write_text(report)
    return readme_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate exp3 characterization tables and figures.")
    parser.add_argument("--results-dir", default="results", help="Path to active exp3 results directory")
    parser.add_argument("--exp1-archive", default="results/_archive/2026-10-16_exp1_lognormal", help="Path to archived exp1 benchmark directory")
    parser.add_argument("--seed-dir", default="data/seed_schedule", help="Path to empirical seed schedule directory")
    parser.add_argument("--vax-dir", default="data/vax_schedule", help="Path to empirical vaccination schedule directory")
    parser.add_argument("--incidence-cache", default="results/re_comparison/incidence_cache", help="Optional path to daily incidence CSV cache from scripts/calc_re.py")
    parser.add_argument("--out-dir", default="docs/characterization", help="Output directory for generated tables and figures")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    seed_dir = Path(args.seed_dir)
    vax_dir = Path(args.vax_dir)
    inc_dir = Path(args.incidence_cache) if args.incidence_cache else None

    bdata = load_benchmarks(Path(args.results_dir), Path(args.exp1_archive))

    t1 = build_table1_epi(seed_dir, out_dir)
    t2 = build_table2_clock(bdata, out_dir)
    t3 = build_table3_mugration(bdata, out_dir)
    t4 = build_table4_auc(bdata, out_dir)

    f1 = plot_fig1_epi_dynamics(seed_dir, out_dir)
    f2 = plot_fig2_clock_closure(t2, out_dir)
    f3 = plot_fig3_surveillance(t3, t4, out_dir)
    f4 = plot_fig4_state_timelines(seed_dir, vax_dir, inc_dir, out_dir)
    build_summary_readme(out_dir)

    print(f"Generated 4 CSV/Markdown tables, 4 SVG/PNG figures, and README.md in {out_dir}/:")
    for p in sorted(out_dir.iterdir()):
        print(f"  - {p.name} ({p.stat().st_size / 1024:.1f} KB)")


if __name__ == "__main__":
    main()
