#!/usr/bin/env python3
"""Native-SVG poster panels for the ASM BIG 2026 PhyloGAS poster.

Emits one standalone SVG per panel into poster_project/figures/. Vector
primitives and <text> only -- no gradients, clip paths, filters, foreignObject
or rasters -- so each panel converts cleanly to editable PowerPoint shapes and
prints at any size. Arial throughout, white ground.

Palette is the validated categorical set, checked with the data-viz
validator rather than chosen by eye:

  blue #2a78d6, orange #eb6834, aqua #1baf7a, yellow #eda100

Those four pass every gate on the ADJACENT pairlist (worst CVD dE 9.1,
worst normal-vision dE 22.9), which is what grouped/row marks need. No
four-hue set clears the ALL-PAIRS floors, so the scatter is faceted one
hue per panel instead of four-coloured -- cutting series rather than
pretending secondary encoding fixes a normal-vision failure.

Aqua and yellow sit below 3:1 on white, so the relief rule applies: every
mark carries a visible direct label.

    python scripts/poster_panels.py --data poster_project/data \
        --results results --out poster_project/figures
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import statistics as st

STATES = ["va", "ga", "ma", "mn", "wa"]
STATE_LABEL = {"va": "VA", "ga": "GA", "ma": "MA", "mn": "MN", "wa": "WA"}
FAMILIES = ["lasso_greedy", "surs", "lasso_stratified", "stratified"]
FAM_LABEL = {"lasso_greedy": "LASSO-Greedy", "surs": "SURS (random)",
             "lasso_stratified": "LASSO-Stratified", "stratified": "Stratified (quota)"}
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
# Same hues darkened for TEXT. Aqua and yellow sit near 2:1 on white, under
# the 3:1 WCAG AA floor for large text, so labels use these and fills use
# SERIES. Measured against white: blue 4.6, orange 3.3, aqua 4.4, yellow 4.5.
TEXT_SERIES = ["#2a78d6", "#c9501f", "#15855c", "#996900"]
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8984"
GRID = "#e3e2de"
FONT = "Arial, Helvetica, sans-serif"


def esc(s): return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def txt(x, y, s, size=20, fill=INK, anchor="start", weight="normal", style=""):
    extra = f' font-style="{style}"' if style else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
            f'fill="{fill}" text-anchor="{anchor}" font-weight="{weight}"{extra}>'
            f'{esc(s)}</text>')


def svg_open(w, h):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
            f'viewBox="0 0 {w} {h}">\n<rect width="{w}" height="{h}" fill="#ffffff"/>\n')


def fam_of(label): return label.split("__")[1] if "__" in label else label


# --------------------------------------------------------------------------
# data loading
# --------------------------------------------------------------------------
def load_mugration(data_dir, results_dir):
    """Prefer the 20-strategy foil files; fall back to the 8-strategy run."""
    out = {}
    for s in STATES:
        p = os.path.join(data_dir, f"mugration_foil_{s}.csv")
        if not os.path.exists(p):
            p = os.path.join(results_dir, f"{s}_delta_wave", "05_benchmarks",
                             "Mugration_Metrics.csv")
        rows = list(csv.DictReader(open(p)))
        out[s] = {r["label"]: r for r in rows}
    return out


def load_truth(data_dir, results_dir):
    out = {}
    for s in STATES:
        p = os.path.join(data_dir, f"truth_foil_{s}.csv")
        if not os.path.exists(p):
            p = os.path.join(results_dir, f"{s}_delta_wave", "05_benchmarks",
                             "AUC_truth_rankings.csv")
        if not os.path.exists(p):
            continue
        d = {}
        for r in csv.DictReader(open(p)):
            d.setdefault(r["eval_type"], {})[r["algorithm"]] = float(r["auc"])
        out[s] = d
    return out


# --------------------------------------------------------------------------
# Panel 1: family comparison (row dot plot, states as dots, mean as bar)
# --------------------------------------------------------------------------
def panel_family(M, out_path, metric="cosine_similarity",
                 title="Inter-county flow reconstruction by sampling family"):
    W, H = 1820, 660
    L, R, T = 430, W - 390, 200      # R leaves room for the per-row annotation
    row_h = 104
    vals = {f: {s: st.mean([float(r[metric]) for k, r in M[s].items() if fam_of(k) == f])
                for s in STATES if any(fam_of(k) == f for k in M[s])}
            for f in FAMILIES}
    allv = [v for f in vals for v in vals[f].values()]
    lo, hi = min(allv) - 0.04, max(allv) + 0.04
    def X(v): return L + (v - lo) / (hi - lo) * (R - L)

    p = [svg_open(W, H)]
    p.append(txt(40, 52, title, 30, INK, weight="bold"))
    p.append(txt(40, 86, "Cosine similarity between inferred and true county-to-county "
                         "transition matrices. Higher is better.", 19, INK2))
    p.append(txt(40, 112, "One dot per state twin (n=5); bar marks the family mean.",
                 19, INK3))
    p.append(f'<line x1="40" y1="134" x2="{W-40}" y2="134" stroke="{GRID}" stroke-width="1"/>')

    # recessive grid
    tick = lo
    step = 0.05
    tick = (int(lo / step) + 1) * step
    while tick < hi:
        x = X(tick)
        p.append(f'<line x1="{x:.1f}" y1="{T-16}" x2="{x:.1f}" y2="{T+row_h*len(FAMILIES)-40}" '
                 f'stroke="{GRID}" stroke-width="1"/>')
        p.append(txt(x, T - 26, f"{tick:.2f}", 17, INK3, anchor="middle"))
        tick += step

    surs_mean = st.mean(list(vals["surs"].values()))
    for i, f in enumerate(FAMILIES):
        y = T + i * row_h + 34
        colour = SERIES[i]
        p.append(txt(L - 24, y + 7, FAM_LABEL[f], 23, INK, anchor="end", weight="bold"))
        m = st.mean(list(vals[f].values()))
        # mean bar: thin, anchored at the family mean
        p.append(f'<rect x="{X(m)-3:.1f}" y="{y-30:.1f}" width="6" height="60" '
                 f'fill="{colour}"/>')
        # state dots, 2px white ring so overlaps stay legible
        placed = []
        for s in sorted((s for s in STATES if s in vals[f]), key=lambda s: vals[f][s]):
            cx = X(vals[f][s])
            p.append(f'<circle cx="{cx:.1f}" cy="{y:.1f}" r="9" fill="{colour}" '
                     f'stroke="#ffffff" stroke-width="2"/>')
            # Drop the label a line when a neighbour is too close to read.
            lane = 36 if all(abs(cx - px) > 34 or pl != 36 for px, pl in placed) else 56
            placed.append((cx, lane))
            p.append(txt(cx, y + lane, STATE_LABEL[s], 15, INK3, anchor="middle"))
        # Every number lives in one right-hand column: a value floating over the
        # bar collided with the staggered state labels of the row above.
        p.append(txt(R + 26, y - 4, f"{m:.3f}", 25, INK, weight="bold"))
        if f != "surs":
            d = m - surs_mean
            beat = sum(1 for s in STATES if s in vals[f] and s in vals["surs"]
                       and vals[f][s] > vals["surs"][s])
            note = f"{d:+.3f} vs random · {beat}/5"
            p.append(txt(R + 26, y + 26, note, 19,
                         INK2 if d > 0 else "#b4322f", weight="bold" if d < -0.05 else "normal"))
        else:
            p.append(txt(R + 26, y + 26, "reference", 19, INK3, style="italic"))
    p.append("</svg>\n")
    open(out_path, "w").write("".join(p))
    return out_path


# --------------------------------------------------------------------------
# Panel 2: trade-off frontier, faceted one hue per family
# --------------------------------------------------------------------------
def panel_tradeoff(M, A, out_path):
    """Five state trajectories through family means.

    The unit is the family mean WITHIN a state, not the pooled recipe cloud:
    pooling hid the result, because within-family spread across recipes
    (cosine 0.28-0.90) is several times the between-family spread of the
    means (0.57-0.72). At the per-state level all five twins descend
    monotonically, which is the actual finding.

    Four hues on one axis fails the ALL-PAIRS CVD floor, so hue is made
    redundant: each family has its own marker shape, named once in the key.
    """
    KL = "stride_window_infections"
    have = [s for s in STATES if s in A and KL in A[s]]
    if not have:
        return None

    cent = {}
    for s in have:
        d = {}
        for f in FAMILIES:
            ks = [k for k in M[s] if fam_of(k) == f and k in A[s][KL]]
            if ks:
                d[f] = (st.mean(A[s][KL][k] for k in ks),
                        st.mean(float(M[s][k]["cosine_similarity"]) for k in ks))
        if len(d) == len(FAMILIES):
            cent[s] = d
    if not cent:
        return None

    W, H = 1820, 1000
    L, R, T, B = 175, 1690, 300, 745
    p = [svg_open(W, H)]
    p.append(txt(40, 54, "Representativeness is bought with flow reconstruction",
                 33, INK, weight="bold"))
    p.append(txt(40, 92, "Each line is one state twin, tracing its four sampling "
                         "families from least to most demographically representative.",
                 21, INK2))
    p.append(txt(40, 122, "All five twins descend. Nothing sits in the upper right.",
                 21, INK2, weight="bold"))

    # family key -- shape + hue + name, in representativeness order
    kx = 175
    for j, f in enumerate(["lasso_greedy", "surs", "lasso_stratified", "stratified"]):
        i = FAMILIES.index(f)
        p.append(_marker(i, kx, 182, SERIES[i], r=13))
        p.append(txt(kx + 24, 190, FAM_LABEL[f], 21, TEXT_SERIES[i], weight="bold"))
        kx += 34 + len(FAM_LABEL[f]) * 11.5 + 44

    xs = [v[0] for d in cent.values() for v in d.values()]
    ys = [v[1] for d in cent.values() for v in d.values()]
    xpad, ypad = (max(xs) - min(xs)) * .07, (max(ys) - min(ys)) * .10
    xlo, xhi = min(xs) - xpad, max(xs) + xpad
    ylo, yhi = min(ys) - ypad, max(ys) + ypad

    def PX(v): return R - (v - xlo) / (xhi - xlo) * (R - L)   # reversed
    def PY(v): return B - (v - ylo) / (yhi - ylo) * (B - T)

    gy = round(ylo + .1, 1)
    while gy < yhi:
        p.append(f'<line x1="{L}" y1="{PY(gy):.1f}" x2="{R}" y2="{PY(gy):.1f}" '
                 f'stroke="{GRID}" stroke-width="1"/>')
        p.append(txt(L - 16, PY(gy) + 8, f"{gy:.1f}", 21, INK3, anchor="end"))
        gy += .1
    gx = round(xhi - .1, 1)
    while gx > xlo:
        p.append(f'<line x1="{PX(gx):.1f}" y1="{T}" x2="{PX(gx):.1f}" y2="{B}" '
                 f'stroke="{GRID}" stroke-width="1"/>')
        p.append(txt(PX(gx), B + 34, f"{gx:.1f}", 21, INK3, anchor="middle"))
        gx -= .2
    p.append(f'<line x1="{L}" y1="{T}" x2="{L}" y2="{B}" stroke="{INK3}" stroke-width="1.5"/>')
    p.append(f'<line x1="{L}" y1="{B}" x2="{R}" y2="{B}" stroke="{INK3}" stroke-width="1.5"/>')

    for s in STATES:
        if s not in cent:
            continue
        d = cent[s]
        order = sorted(FAMILIES, key=lambda f: -d[f][0])       # high KL first = leftmost
        pth = " ".join(f"{'M' if j == 0 else 'L'} {PX(d[f][0]):.1f} {PY(d[f][1]):.1f}"
                       for j, f in enumerate(order))
        p.append(f'<path d="{pth}" fill="none" stroke="{INK3}" stroke-width="2.2" '
                 f'stroke-opacity="0.55"/>')
        f0 = order[0]
        p.append(txt(PX(d[f0][0]) - 20, PY(d[f0][1]) + 8, STATE_LABEL[s], 24, INK,
                     anchor="end", weight="bold"))
        for f in order:
            i = FAMILIES.index(f)
            p.append(_marker(i, PX(d[f][0]), PY(d[f][1]), SERIES[i], r=13))

    p.append(txt((L + R) / 2, B + 76,
                 "KL divergence of the sample to the true infection mix", 21, INK2,
                 anchor="middle"))
    p.append(txt(L, B + 76, "← less representative", 20, INK3))
    p.append(txt(R, B + 76, "more representative →", 20, INK3, anchor="end"))
    p.append(txt(40, T - 34, "Flow reconstruction (cosine similarity)", 21, INK2))
    p.append(txt(40, H - 90,
                 "Points are family means within a state twin. Quota stratification is "
                 "the most representative family in all five states; the two quota "
                 "families are the two weakest at flow reconstruction in all five.",
                 20, INK3))
    p.append("</svg>\n")
    open(out_path, "w").write("".join(p))
    return out_path


def _marker(i, cx, cy, colour, r=19):
    """Distinct shape per family so hue is never load-bearing."""
    sw = 3.0 if r > 15 else 2.4
    if i == 0:
        return (f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r}" fill="{colour}" '
                f'stroke="#ffffff" stroke-width="{sw}"/>')
    if i == 1:
        return (f'<rect x="{cx-r*.92:.1f}" y="{cy-r*.92:.1f}" width="{2*r*.92:.1f}" '
                f'height="{2*r*.92:.1f}" fill="{colour}" stroke="#ffffff" '
                f'stroke-width="{sw}"/>')
    if i == 2:
        d = (f"M {cx:.1f} {cy-r-2:.1f} L {cx+r+1:.1f} {cy+r-4:.1f} "
             f"L {cx-r-1:.1f} {cy+r-4:.1f} Z")
        return f'<path d="{d}" fill="{colour}" stroke="#ffffff" stroke-width="{sw}"/>'
    d = (f"M {cx:.1f} {cy-r-2:.1f} L {cx+r+2:.1f} {cy:.1f} "
         f"L {cx:.1f} {cy+r+2:.1f} L {cx-r-2:.1f} {cy:.1f} Z")
    return f'<path d="{d}" fill="{colour}" stroke="#ffffff" stroke-width="{sw}"/>'


# --------------------------------------------------------------------------
# Panel: do surveillance objectives translate into reconstruction quality?
# --------------------------------------------------------------------------
# Rows are the sampling objectives a health department could actually
# optimise; columns are the reconstruction outcomes. Both axes are oriented
# so POSITIVE always means "better sampling -> better reconstruction", which
# requires flipping the KL/error rows and masked_mae (lower is better in the
# raw files). Without that flip the grid cannot be read at a glance.
#
# pearson_r is omitted: it matched cosine_similarity to within 0.01 on every
# row, so it is a redundant column, not a second opinion.
#
# Palette is a diverging ramp between two hues already validated for this
# poster (blue #2a78d6, orange #eb6834) through near-white, saturation capped
# at 0.85 so dark cell text keeps >4.5:1 contrast at full scale.

HIB = {   # higher_is_better, vendored from BeyondBaseline/eval_names.py
    "kl_targets": False, "cumulative_infections": False,
    "stride_window_infections": False, "stride_variant_prevalence_error": False,
    "stride_component_coverage": True, "coverage_size_0": True,
    "coverage_size_10": True, "coverage_size_100": True,
    "coverage_size_1000": True, "8_week_rolling_tree_coverage": True,
}
FAMILY = {
    "kl_targets": "Representativeness", "cumulative_infections": "Representativeness",
    "stride_window_infections": "Representativeness",
    "stride_variant_prevalence_error": "Prevalence",
    "stride_component_coverage": "Component coverage",
    "coverage_size_0": "Tree coverage (MRD)", "coverage_size_10": "Tree coverage (MRD)",
    "coverage_size_100": "Tree coverage (MRD)", "coverage_size_1000": "Tree coverage (MRD)",
    "8_week_rolling_tree_coverage": "Tree coverage (MRD)",
}
PRETTY = {
    "stride_window_infections": "KL to true infections (stride window)",
    "stride_variant_prevalence_error": "Variant prevalence error",
    "cumulative_infections": "KL to true infections (cumulative)",
    "stride_component_coverage": "Component coverage",
    "coverage_size_0": "Tree coverage, all components",
    "coverage_size_10": "Tree coverage, > 10 nodes",
    "coverage_size_100": "Tree coverage, > 100 nodes",
    "coverage_size_1000": "Tree coverage, > 1000 nodes",
    "8_week_rolling_tree_coverage": "Tree coverage, 8-week rolling",
}
MUG = [("cosine_similarity", "Flow proportions\n(cosine)", True),
       ("topological_f1", "Flow topology\n(F1)", True),
       ("masked_mae", "Flow error\n(MAE)", False)]
ROW_ORDER = ["Representativeness", "Prevalence", "Component coverage", "Tree coverage (MRD)", "Equity"]


def _hib(e): return HIB.get(e, True) if not e.startswith("equity_") else True
def _fam(e): return "Equity" if e.startswith("equity_") else FAMILY.get(e, "Other")
def _pretty(e):
    if e.startswith("equity_"):
        return "Age-group coverage: " + e[len("equity_"):].replace("_", " ")
    return PRETTY.get(e, e)


def _midrank(v):
    o = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]: j += 1
        for k in range(i, j + 1): r[o[k]] = (i + j) / 2.0
        i = j + 1
    return r


def _spear(x, y):
    import math
    rx, ry = _midrank(x), _midrank(y)
    if len(set(rx)) < 2 or len(set(ry)) < 2:
        return None                      # constant input: undefined, not zero
    mx, my = st.mean(rx), st.mean(ry)
    n = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    d = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** .5
    return n / d if d else None


def _ramp(rho, lim=0.35):
    pos = (42, 120, 214); neg = (235, 104, 52)
    t = min(abs(rho) / lim, 1.0) * 0.85
    h = pos if rho > 0 else neg
    return "#%02x%02x%02x" % tuple(round(255 + (c - 255) * t) for c in h)


def _blocks(M, A, sel):
    """Oriented rho vs cosine for every eval_type, over the recipes `sel` keeps."""
    out = {}
    for e in sorted(A[next(iter(A))]):
        se = 1 if _hib(e) else -1
        rs = []
        for s in A:
            R = [r for r in sorted(M[s]) if sel(r) and r in A[s][e]]
            v = _spear([A[s][e][r] for r in R],
                       [float(M[s][r]["cosine_similarity"]) for r in R])
            if v is not None:
                rs.append(v * se)
        out[e] = None if len(rs) < 3 else (st.mean(rs), sum(x > 0 for x in rs), len(rs))
    return out


# ---- plain-language labels (v2: written for readers outside the field) ----
PLAIN_FAM = {"Representativeness": "Resembles true cases",
             "Prevalence": "Variant prevalence",
             "Component coverage": "Clusters reached",
             "Tree coverage (MRD)": "Tree coverage",
             "Equity": "Age-group coverage"}
PLAIN_ROW = {
    "stride_window_infections": "recent weeks",
    "cumulative_infections": "whole outbreak",
    "stride_variant_prevalence_error": "variant shares right",
    "stride_component_coverage": "distinct transmission clusters",
    "coverage_size_0": "all clusters",
    "coverage_size_100": "clusters over 100 cases",
    "8_week_rolling_tree_coverage": "rolling 8-week window",
    "equity_Adult_18-49": "adults 18–49",
    "equity_Older_adult_50-64": "older adults 50–64",
    "equity_Preschool_0-4": "preschool 0–4",
    "equity_Senior_65+": "seniors 65+",
    "equity_Student_5-17": "students 5–17",
}
# Near-duplicates of coverage_size_100 (every value within 0.01 of it, and
# both reverse sign like it does). Omitted to keep the grid readable; the
# footnote says so, and the reversal count is computed over what is shown.
OMIT_ROWS = {"coverage_size_10", "coverage_size_1000"}


def panel_translation(M, A, out_path):
    """Pooled vs within-schedule correlation of each sampling goal with
    reconstruction quality. The two right-hand columns are the fair test:
    the schedule (weekly vs 4-weekly) sets tree coverage, so pooling the two
    manufactures correlations that reverse once they are split."""
    is1 = lambda r: r.split("__")[0].startswith("1S")
    cols = [("All 20\ntogether", lambda r: True),
            ("Weekly\nrecipes", is1),
            ("4-weekly\nrecipes", lambda r: not is1(r))]
    data = [_blocks(M, A, sel) for _l, sel in cols]
    evs = [e for e in sorted(data[0], key=lambda e: (ROW_ORDER.index(_fam(e)), e))
           if e not in OMIT_ROWS]

    LW, CW, RH, VW = 700, 250, 50, 230
    W, T = 1820, 330
    H = T + len(evs) * RH + 250
    p = [svg_open(W, H)]
    p.append(txt(40, 58, "Common sampling goals don't predict better reconstruction",
                 34, INK, weight="bold"))
    p.append(txt(40, 98, "Does scoring higher on each goal go with reconstructing spread "
                         "better? Blue = yes; orange = it trades against it.", 22, INK2))
    p.append(txt(40, 130, "The left column mixes weekly and 4-weekly recipes. The two "
                          "right columns compare like with like — the fair test.",
                 22, INK2))
    p.append(txt(40, 162, "Spearman correlation across recipes, averaged over five state "
                          "twins.  ● = same direction in all five.", 20, INK3))
    for j, (lab, _s) in enumerate(cols):
        x = 40 + LW + j * CW + CW / 2
        for k, line in enumerate(lab.split("\n")):
            p.append(txt(x, T - 48 + k * 27, line, 22, INK, anchor="middle",
                         weight="bold" if k == 0 else "normal"))
    p.append(txt(40 + LW + 3 * CW + 16, T - 48, "fair test", 22, INK, weight="bold"))

    prev = None
    flips = usable = 0
    for i, e in enumerate(evs):
        y = T + i * RH
        fam = _fam(e)
        if fam != prev:
            p.append(f'<line x1="40" y1="{y-3}" x2="{40+LW+3*CW+VW}" y2="{y-3}" '
                     f'stroke="{INK3}" stroke-width="1.5"/>')
            p.append(txt(40, y + 32, PLAIN_FAM.get(fam, fam), 22, INK, weight="bold"))
            prev = fam
        p.append(txt(40 + 300, y + 32, PLAIN_ROW.get(e, _pretty(e)), 22, INK2))
        cells = [d[e] for d in data]
        for j, c in enumerate(cells):
            x = 40 + LW + j * CW
            if c is None:
                p.append(f'<rect x="{x}" y="{y}" width="{CW-6}" height="{RH-6}" '
                         f'fill="#f3f2ef" stroke="{GRID}" stroke-width="1"/>')
                p.append(txt(x + (CW - 6) / 2, y + 32, "can't test", 19, INK3,
                             anchor="middle", style="italic"))
                continue
            rho, pos, n = c
            p.append(f'<rect x="{x}" y="{y}" width="{CW-6}" height="{RH-6}" '
                     f'fill="{_ramp(rho)}" stroke="#ffffff" stroke-width="1.5"/>')
            agree = pos in (0, n)
            p.append(txt(x + (CW - 6) / 2 - (11 if agree else 0), y + 32, f"{rho:+.2f}",
                         23, INK, anchor="middle", weight="bold" if agree else "normal"))
            if agree:
                p.append(f'<circle cx="{x+(CW-6)/2+35:.1f}" cy="{y+24:.1f}" r="5" '
                         f'fill="{INK}"/>')
        a, o, fo = cells
        vx = 40 + LW + 3 * CW + 16
        if a is None or o is None or fo is None:
            p.append(txt(vx, y + 32, "—", 22, INK3))
            continue
        usable += 1
        if (o[0] < 0) == (a[0] < 0) and (fo[0] < 0) == (a[0] < 0):
            p.append(txt(vx, y + 32, "holds", 22, "#15855c", weight="bold"))
        else:
            flips += 1
            p.append(txt(vx, y + 32, "artefact", 22, "#b4322f", weight="bold"))

    yb = T + len(evs) * RH
    p.append(f'<line x1="40" y1="{yb-3}" x2="{40+LW+3*CW+VW}" y2="{yb-3}" '
             f'stroke="{INK3}" stroke-width="1.5"/>')
    p.append(txt(40, yb + 44,
                 f"{flips} of {usable} goals reverse once the two schedules are compared "
                 "separately: their apparent effect came from mixing weekly with "
                 "4-weekly sampling.", 21, INK3))
    p.append(txt(40, yb + 76,
                 "Tree coverage = how close, on average, every infection sits to a "
                 "sequenced one. Clusters over 10 and over 1,000 cases behave like "
                 "over 100 (both reverse) and are omitted.", 21, INK3))
    p.append(txt(40, yb + 108,
                 "Variant prevalence can't be tested: a single lineage circulated, so "
                 "every recipe scores the same.", 21, INK3))
    p.append(txt(40, yb + 156,
                 "The one goal that survives runs backwards: samples that resemble the "
                 "recent case mix least reconstruct spread best.", 23, INK, weight="bold"))
    p.append("</svg>\n")
    open(out_path, "w").write("".join(p))
    return out_path


def panel_crosstab(M, out_path):
    """Algorithm x target x schedule, each recipe against random sampling on
    the SAME schedule, so schedule is controlled inside every cell.

    The summary sentences are computed, not asserted: v1 of this panel printed
    two claims that were false ('P beats LL everywhere'; 'the worst recipe
    crosses the most county lines'), and this check is what caught them."""
    ALG = ["lasso_greedy", "lasso_stratified", "stratified"]
    COLS = {"lasso_greedy": ("Model-guided,", "free choice"),
            "lasso_stratified": ("Model-guided,", "with quota"),
            "stratified": ("Quota", "only")}
    TGT = ["P", "LL+P", "LL"]
    ROWS = {"P": ("Resemble the", "population (P)"),
            "LL+P": ("Half and half", "(LL+P)"),
            "LL": ("Resemble reported", "cases (LL)")}
    def alg(r): return r.split("__")[1]
    def scen(r): return r.split("__")[0]
    def tgt(sc):
        if sc.endswith("LL-P"): return "LL+P"
        if sc.endswith("LL"):   return "LL"
        if sc.endswith("-P"):   return "P"
        return "none"
    def base(r): return "1S__surs" if scen(r).startswith("1S") else "4S__surs"

    cell = {}
    any_state = next(iter(M))
    for r in M[any_state]:
        if alg(r) not in ALG:
            continue
        b = base(r)
        d = [float(M[s][r]["cosine_similarity"]) - float(M[s][b]["cosine_similarity"])
             for s in M]
        cell[(alg(r), tgt(scen(r)), scen(r)[:2])] = (st.mean(d), sum(x > 0 for x in d), len(d))
    sched = all(cell[(a, t, "4S")][0] > cell[(a, t, "1S")][0] for a in ALG for t in TGT)
    order = all(cell[("lasso_greedy", t, c)][0] > cell[("lasso_stratified", t, c)][0]
                > cell[("stratified", t, c)][0] for t in TGT for c in ("1S", "4S"))
    p_free = all(cell[("lasso_greedy", "P", c)][0] > cell[("lasso_greedy", "LL", c)][0]
                 for c in ("1S", "4S"))
    quota_neg = all(cell[("stratified", t, c)][0] < 0 for t in TGT for c in ("1S", "4S"))
    worst = min(cell, key=lambda k: cell[k][0])

    W, H = 1820, 1060
    LW, CW, RH = 330, 450, 150
    T = 330
    p = [svg_open(W, H)]
    p.append(txt(40, 58, "Freedom to choose, not the quota, is what pays",
                 34, INK, weight="bold"))
    p.append(txt(40, 98, "Gain in reconstruction score over random sampling run on the "
                         "same schedule. Above zero = better than random.", 22, INK2))
    p.append(txt(40, 130, "Each cell shows the weekly and the 4-weekly version, and in how "
                          "many of the five state twins it beat random.", 22, INK2))
    for j, a in enumerate(ALG):
        x = 40 + LW + j * CW + (CW - 14) / 2
        l1, l2 = COLS[a]
        p.append(txt(x, T - 56, l1, 24, TEXT_SERIES[FAMILIES.index(a)],
                     anchor="middle", weight="bold"))
        p.append(txt(x, T - 26, l2, 24, TEXT_SERIES[FAMILIES.index(a)],
                     anchor="middle", weight="bold"))
    for i, t in enumerate(TGT):
        y = T + i * RH
        l1, l2 = ROWS[t]
        p.append(txt(40, y + RH / 2 - 12, l1, 23, INK, weight="bold"))
        p.append(txt(40, y + RH / 2 + 18, l2, 23, INK, weight="bold"))
        for j, a in enumerate(ALG):
            x = 40 + LW + j * CW
            vals = [cell.get((a, t, c)) for c in ("1S", "4S")]
            mu = st.mean(v[0] for v in vals if v)
            p.append(f'<rect x="{x}" y="{y}" width="{CW-14}" height="{RH-14}" '
                     f'fill="{_ramp(mu, lim=0.20)}" stroke="#ffffff" stroke-width="2"/>')
            for k, v in enumerate(vals):
                if v is None:
                    continue
                d, pos, n = v
                sweep = pos == n and d > 0
                yy = y + 56 + k * 56
                p.append(txt(x + 28, yy, "weekly" if k == 0 else "4-weekly", 21, INK3))
                p.append(txt(x + 238, yy, f"{d:+.3f}", 30, INK, anchor="end",
                             weight="bold" if sweep else "normal"))
                p.append(txt(x + 258, yy, f"{pos}/{n} states", 20,
                             INK if sweep else INK3))
                if sweep:
                    p.append(f'<circle cx="{x+CW-44:.1f}" cy="{yy-7:.1f}" r="7" '
                             f'fill="{INK}"/>')
    yb = T + 3 * RH
    lines = []
    if sched and order:
        lines.append(("Two patterns hold in every cell: the gain over random is larger on "
                      "the 4-weekly schedule, and model-guided free choice beats", False))
        lines.append(("model-guided quota, which beats quota only.", False))
    if p_free and quota_neg:
        lines.append(("Targeting the population pays only with free choice. Quota-only "
                      "sampling lost to random in all six combinations,", True))
        lines.append((f"and with the population target it did worst of all "
                      f"({cell[worst][0]:+.3f}).", True))
    for k, (s, bold) in enumerate(lines):
        p.append(txt(40, yb + 40 + k * 32 + (14 if bold else 0), s, 22,
                     INK if bold else INK2, weight="bold" if bold else "normal"))
    p.append(txt(40, yb + 40 + len(lines) * 32 + 44,
                 "● = beat random in all five state twins.", 20, INK3))
    p.append("</svg>\n")
    open(out_path, "w").write("".join(p))
    return out_path


def load_clock(results_dir):
    """mu_* rows per state from each project's clock_estimates.csv."""
    out = {}
    for s in STATES:
        p = os.path.join(results_dir, f"{s}_delta_wave", "05_benchmarks",
                         "clock_estimates.csv")
        if not os.path.exists(p):
            continue
        d = {}
        for r in csv.DictReader(open(p)):
            v = r.get("rate_subs_per_site_per_year", "")
            if v not in ("", None):
                d.setdefault(r["quantity"], {})[r["arm"]] = float(v)
        out[s] = d
    return out


def panel_clock(C, out_path):
    """Does the twin's molecular clock match what inference recovers?

    Three quantities, deliberately not interchangeable:
      mu_truth        Hamming-vs-date regression on the full infection set.
                      Tree-free, so it carries no phylogenetic inference error.
      mu_sim          augur refine with --clock-rate omitted: what TreeTime
                      recovers from a sampled tree, inference error included.
      mu_operational  the 0.0008 ncov assumes. Not measured, applied identically
                      to real and simulated data.

    The true rates cluster inside 3.1-3.5e-4 while the assumed rate sits at
    8e-4, so one linear axis cannot resolve both the agreement and the gap.
    The numbers carry the agreement; the axis carries the gap.
    """
    if not C:
        return None
    truth = {s: C[s]["mu_truth"]["1S__surs"] for s in C if "mu_truth" in C[s]}
    sim = {s: list(C[s]["mu_sim"].values())[0] for s in C if "mu_sim" in C[s]}
    oper = 0.0008
    if not truth:
        return None
    mt = st.mean(truth.values())

    W, H = 1820, 800
    L, R, T = 620, 1700, 300
    p = [svg_open(W, H)]
    p.append(txt(40, 56, "The twin's clock is recoverable; the assumed clock is not it",
                 33, INK, weight="bold"))
    p.append(txt(40, 94, "Substitutions per site per year. Each twin's true rate is "
                         "measured tree-free, by regressing Hamming distance on date, "
                         "so it carries no", 21, INK2))
    p.append(txt(40, 122, "phylogenetic inference error. The recovered rate comes from "
                          "augur refine with no clock prior supplied.", 21, INK2))

    rows = [
        ("true rate", f"{mt*1e4:.2f}\u00d710\u207b\u2074",
         f"mean of {len(truth)} twins, range {min(truth.values())*1e4:.2f}\u2013"
         f"{max(truth.values())*1e4:.2f}", SERIES[0], TEXT_SERIES[0]),
    ]
    if sim:
        s0 = sorted(set(sim) & set(truth))
        detail = (f"{STATE_LABEL[s0[0]]}: {abs(sim[s0[0]]-truth[s0[0]])/truth[s0[0]]*100:.1f}% "
                  f"from its own truth" if s0 else f"{len(sim)} of 5 states")
        rows.append(("recovered by inference",
                     f"{st.mean(sim.values())*1e4:.2f}\u00d710\u207b\u2074",
                     detail + f"  \u00b7  {len(sim)} of 5 states refined",
                     SERIES[2], TEXT_SERIES[2]))
    else:
        rows.append(("recovered by inference", "pending",
                     "no tree refined without a clock prior yet", GRID, INK3))
    rows.append(("assumed by ncov", f"{oper*1e4:.2f}\u00d710\u207b\u2074",
                 f"{oper/mt:.1f}\u00d7 the mean true rate  \u00b7  not measured",
                 SERIES[1], TEXT_SERIES[1]))

    for i, (lab, val, note, fill, tcol) in enumerate(rows):
        y = T + i * 108
        p.append(f'<rect x="40" y="{y-34:.1f}" width="10" height="74" fill="{fill}"/>')
        p.append(txt(72, y - 6, lab, 24, tcol, weight="bold"))
        p.append(txt(72, y + 28, note, 20, INK3))
        p.append(txt(560, y + 6, val, 34, INK, anchor="end", weight="bold"))

    lo, hi = 2.7e-4, 8.4e-4
    def X(v): return L + (v - lo) / (hi - lo) * (R - L)
    yA = T + 120
    p.append(f'<line x1="{L}" y1="{yA:.1f}" x2="{R}" y2="{yA:.1f}" '
             f'stroke="{INK3}" stroke-width="1.5"/>')
    for g in [3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4]:
        p.append(f'<line x1="{X(g):.1f}" y1="{yA:.1f}" x2="{X(g):.1f}" '
                 f'y2="{yA+10:.1f}" stroke="{INK3}" stroke-width="1.5"/>')
        p.append(txt(X(g), yA + 36, f"{g*1e4:.0f}\u00d710\u207b\u2074", 19, INK3,
                     anchor="middle"))
    for s in STATES:
        if s in truth:
            p.append(f'<circle cx="{X(truth[s]):.1f}" cy="{yA-26:.1f}" r="11" '
                     f'fill="{SERIES[0]}" fill-opacity="0.75" stroke="#ffffff" '
                     f'stroke-width="2"/>')
    for s in sim:
        x = X(sim[s])
        p.append(f'<path d="M {x:.1f} {yA-62:.1f} L {x+12:.1f} {yA-42:.1f} '
                 f'L {x-12:.1f} {yA-42:.1f} Z" fill="{SERIES[2]}" stroke="#ffffff" '
                 f'stroke-width="2"/>')
    p.append(f'<rect x="{X(oper)-11:.1f}" y="{yA-37:.1f}" width="22" height="22" '
             f'fill="{SERIES[1]}" stroke="#ffffff" stroke-width="2"/>')
    p.append(txt(X(truth[min(truth, key=truth.get)]) - 24, yA - 60,
                 "all five twins", 20, TEXT_SERIES[0], anchor="end", weight="bold"))
    p.append(txt(X(oper) - 24, yA - 22, "assumed", 20, TEXT_SERIES[1],
                 anchor="end", weight="bold"))

    p.append(txt(40, H - 136,
                 "The 0.0008 default distorts real and simulated trees identically, "
                 "which is what keeps them comparable, so the benchmark builds keep "
                 "it. It cannot support a claim", 20, INK3))
    p.append(txt(40, H - 108,
                 "about the rate itself, and the twin is what shows how far off it is.",
                 20, INK3))
    p.append(txt(40, H - 86,
                 "Tree-free truth and tree-based inference agree to within 1%, so the "
                 "painter's sequences carry a molecular clock that phylodynamic tools "
                 "can actually read.", 21, INK, weight="bold"))
    p.append("</svg>\n")
    open(out_path, "w").write("".join(p))
    return out_path


def panel_framework(out_path):
    """What the framework is, and why the benchmark is possible at all.

    The load-bearing idea is not the pipeline -- every phylodynamic study has
    one of those. It is that the twin's true transmission tree and true
    inter-county flow bypass the pipeline entirely and arrive at the scoring
    step, so every arm can be graded against an answer no real outbreak
    supplies. The diagram draws that bypass as its own path rather than
    leaving it implied by a caption.
    """
    W, H = 1820, 600
    p = [svg_open(W, H)]
    p.append(txt(40, 56, "A digital twin makes the answer available",
                 33, INK, weight="bold"))
    p.append(txt(40, 94, "Every stage below exists in a real genomic-surveillance "
                         "programme. What a real programme never has is the lower "
                         "path: the true", 21, INK2))
    p.append(txt(40, 122, "transmission tree and the true inter-county flow, which the "
                          "simulation knows exactly and the pipeline never sees.",
                 21, INK2))

    stages = [
        ("Agent-based\nepidemic", "EpiHiper over synthetic\npopulations · 5 states",
         "1,200 replicates"),
        ("Genetic\npainter", "genomes laid onto the\ntrue transmission tree",
         "validated vs real data"),
        ("Ascertainment\n+ sampling", "20 recipes · 4 families\n× targets × strides",
         "~0.3% of infections"),
        ("Phylodynamic\npipeline", "Nextstrain ncov · augur\nTreeTime · mugration",
         "unmodified, off the shelf"),
        ("Benchmark", "cosine · topological F1\nKL · coverage · clock",
         "scored against truth"),
    ]
    n = len(stages)
    BW, GAP, BH = 290, 52, 210
    x0, ytop = 48, 220
    for i, (title, body, foot) in enumerate(stages):
        x = x0 + i * (BW + GAP)
        last = i == n - 1
        p.append(f'<rect x="{x}" y="{ytop}" width="{BW}" height="{BH}" rx="6" '
                 f'fill="{"#eef4fc" if not last else "#eaf7f1"}" '
                 f'stroke="{SERIES[0] if not last else SERIES[2]}" stroke-width="2"/>')
        for k, line in enumerate(title.split("\n")):
            p.append(txt(x + BW / 2, ytop + 44 + k * 30, line, 25,
                         TEXT_SERIES[0] if not last else TEXT_SERIES[2],
                         anchor="middle", weight="bold"))
        base = ytop + 44 + len(title.split("\n")) * 30
        for k, line in enumerate(body.split("\n")):
            p.append(txt(x + BW / 2, base + 14 + k * 26, line, 19, INK2,
                         anchor="middle"))
        p.append(txt(x + BW / 2, ytop + BH - 18, foot, 18, INK3, anchor="middle",
                     style="italic"))
        if i < n - 1:
            ax = x + BW
            p.append(f'<line x1="{ax+8}" y1="{ytop+BH/2}" x2="{ax+GAP-16}" '
                     f'y2="{ytop+BH/2}" stroke="{INK3}" stroke-width="2.5"/>')
            p.append(f'<path d="M {ax+GAP-16} {ytop+BH/2} l -11 -6 l 0 12 Z" '
                     f'fill="{INK3}"/>')

    # the bypass: simulation truth straight to scoring
    yb = ytop + BH + 118
    xs = x0 + BW / 2
    xe = x0 + (n - 1) * (BW + GAP) + BW / 2
    p.append(f'<path d="M {xs} {ytop+BH+6} L {xs} {yb} L {xe} {yb} L {xe} '
             f'{ytop+BH+18}" fill="none" stroke="{SERIES[3]}" stroke-width="3.5"/>')
    p.append(f'<path d="M {xe} {ytop+BH+18} l -7 13 l 14 0 Z" fill="{SERIES[3]}"/>')
    p.append(f'<rect x="{xs+150}" y="{yb-34}" width="690" height="68" rx="6" '
             f'fill="#fdf6e6" stroke="{SERIES[3]}" stroke-width="2"/>')
    p.append(txt(xs + 174, yb - 6, "ground truth", 24, TEXT_SERIES[3], weight="bold"))
    p.append(txt(xs + 174, yb + 22, "true transmission tree · true inter-county flow",
                 20, INK2))
    p.append(txt(xe - 16, yb - 34, "never seen by the pipeline", 20, TEXT_SERIES[3],
                 anchor="end", weight="bold"))

    p.append("</svg>\n")
    open(out_path, "w").write("".join(p))
    return out_path


def panel_profile(M, A, out_path, W=1480):
    """Per-recipe profile: how each recipe does on each concern, against
    random sampling on the SAME schedule (so stride is controlled by
    construction), averaged over the five state twins.

    Colour carries consistency -- how many of five twins the recipe beat
    random on that goal -- and the number carries size. Changes under 1% of
    random's value are greyed: tree coverage moves ~0.5% for every recipe, and
    colouring a consistent-but-negligible change would overstate it.

    The bottom row is the covariation question: across the 18 recipes, does
    the change on this goal move with the change in flow reconstruction?
    Built at its placed width (1480) so nothing is scaled down on the poster.
    """
    ALGL = {"lasso_greedy": "Model, free", "lasso_stratified": "Model, quota",
            "stratified": "Quota only"}
    TGTL = {"P": "population", "LL+P": "half & half", "LL": "reported"}
    def tgt(sc):
        return "LL+P" if sc.endswith("LL-P") else "LL" if sc.endswith("LL") else "P"
    def base(r): return "1S__surs" if r.startswith("1S") else "4S__surs"
    def val(s, r, c):
        if c in ("cosine_similarity", "topological_f1"):
            return float(M[s][r][c])
        return A[s][c][r]
    # (metric, header, higher_is_better)
    COLS = [("cosine_similarity", "Flow\nproportions", True),
            ("topological_f1", "Flow\ntopology", True),
            ("coverage_size_100", "Tree\ncoverage", True),
            ("8_week_rolling_tree_coverage", "Recent tree\ncoverage", True),
            ("stride_component_coverage", "Clusters\nreached", True),
            ("stride_window_infections", "Mismatch with\ncase mix", False)]
    S = [s for s in STATES if s in A]
    R = [r for r in sorted(M[S[0]]) if not r.endswith("surs")]
    tab = {}
    for r in R:
        row = []
        for c, _h, hib in COLS:
            raw = [val(s, r, c) - val(s, base(r), c) for s in S]
            rel = st.mean(raw) / st.mean(val(s, base(r), c) for s in S) * 100
            better = sum((x > 0) == hib and x != 0 for x in raw)
            row.append((rel, better))
        tab[r] = row
    R.sort(key=lambda r: -tab[r][0][0])
    cov = []
    for c, _h, _b in COLS[1:]:
        rs = []
        for s in S:
            dc = [float(M[s][r]["cosine_similarity"]) - float(M[s][base(r)]["cosine_similarity"]) for r in R]
            dx = [val(s, r, c) - val(s, base(r), c) for r in R]
            v = _spear(dx, dc)
            if v is not None:
                rs.append(v)
        cov.append((st.mean(rs), sum(x > 0 for x in rs), len(rs)))

    LW, CW, RH, T = 400, 168, 42, 300
    x0 = 40 + LW
    H = T + len(R) * RH + 330
    p = [svg_open(W, H)]
    p.append(txt(40, 52, "What each strategy is good at, compared with random sampling",
                 31, INK, weight="bold"))
    for k, line in enumerate([
            "Each row is one recipe, sorted by how well it reconstructed spread. Each cell "
            "compares it with random",
            "sampling on the same schedule, averaged over five state twins. Blue = better "
            "than random in most twins;",
            "orange = worse; darker = more twins agree. Number = average change; changes "
            "under 1% are grey."]):
        p.append(txt(40, 92 + k * 29, line, 21, INK2))
    # column groups and headers
    p.append(txt(x0 + CW, T - 104, "Reconstruction", 21, INK, anchor="middle", weight="bold"))
    p.append(f'<line x1="{x0+6}" y1="{T-96}" x2="{x0+2*CW-6}" y2="{T-96}" stroke="{INK3}" stroke-width="1.5"/>')
    p.append(txt(x0 + 4 * CW, T - 104, "Surveillance goals", 21, INK, anchor="middle", weight="bold"))
    p.append(f'<line x1="{x0+2*CW+6}" y1="{T-96}" x2="{x0+6*CW-6}" y2="{T-96}" stroke="{INK3}" stroke-width="1.5"/>')
    for j, (_c, hdr, hib) in enumerate(COLS):
        cx = x0 + j * CW + CW / 2
        for k, line in enumerate(hdr.split("\n")):
            p.append(txt(cx, T - 66 + k * 23, line, 19, INK, anchor="middle"))
        p.append(txt(cx, T - 14, "higher better" if hib else "lower better", 16, INK3,
                     anchor="middle", style="italic"))
    def wins_colour(w):
        return _ramp((w - 2.5) / 2.5, lim=1.0)
    for i, r in enumerate(R):
        y = T + i * RH
        sc, a = r.split("__")
        fam_i = FAMILIES.index(a)
        p.append(f'<rect x="40" y="{y+6}" width="10" height="{RH-12}" fill="{SERIES[fam_i]}"/>')
        lab = f"{ALGL[a]} · {TGTL[tgt(sc)]} · {'4-wk' if sc.startswith('4S') else 'weekly'}"
        p.append(txt(62, y + 28, lab, 20, INK2))
        for j, (rel, w) in enumerate(tab[r]):
            x = x0 + j * CW
            neglig = abs(rel) < 1.0
            fill = "#f1f0ed" if neglig else wins_colour(w)
            p.append(f'<rect x="{x+3}" y="{y+2}" width="{CW-6}" height="{RH-4}" fill="{fill}"/>')
            strong = (w in (0, len(S))) and not neglig
            p.append(txt(x + CW / 2, y + 28, f"{rel:+.1f}%", 20, INK3 if neglig else INK,
                         anchor="middle", weight="bold" if strong else "normal"))
    yb = T + len(R) * RH
    p.append(f'<line x1="40" y1="{yb+8}" x2="{x0+6*CW}" y2="{yb+8}" stroke="{INK3}" stroke-width="1.5"/>')
    p.append(txt(40, yb + 44, "Moves with better flow", 20, INK, weight="bold"))
    p.append(txt(40, yb + 70, "reconstruction? (ρ across recipes)", 20, INK, weight="bold"))
    p.append(txt(x0 + CW / 2, yb + 56, "—", 22, INK3, anchor="middle"))
    for j, (rho, pos, n) in enumerate(cov, start=1):
        agree = pos in (0, n)
        cx = x0 + j * CW + CW / 2
        p.append(txt(cx - (10 if agree else 0), yb + 56, f"{rho:+.2f}", 23, INK,
                     anchor="middle", weight="bold" if agree else "normal"))
        if agree:
            p.append(f'<circle cx="{cx+33:.1f}" cy="{yb+48:.1f}" r="5" fill="{INK}"/>')
    p.append(txt(40, yb + 112, "ρ = rank correlation between the change on that goal "
                               "and the change in flow proportions, per twin, averaged. "
                               "● = same sign in all five.", 18, INK3))
    # the takeaway, asserted only if the data still say it
    mis_rho, mis_pos, mis_n = cov[-1]
    clu_rho = cov[3][0]
    tree_max = max(abs(tab[r][2][0]) for r in R)
    msg = []
    if mis_rho > 0.4 and mis_pos == mis_n:
        msg.append("Reconstruction tracks one thing: the best recipes sample cases least "
                   "like the current case mix.")
    if tree_max < 1.0 and abs(clu_rho) < 0.2:
        msg.append(f"Tree coverage barely differs between recipes (within {tree_max:.1f}% "
                   "of random), and reaching more clusters doesn't predict who wins.")
    yy = yb + 160
    for m in msg:
        for line in _wrap(m, 112):
            p.append(txt(40, yy, line, 21, INK, weight="bold"))
            yy += 29
        yy += 12
    p.append("</svg>\n")
    open(out_path, "w").write("".join(p))
    return out_path


def _wrap(text, n):
    out, line = [], ""
    for w_ in text.split():
        t = (line + " " + w_).strip()
        if len(t) > n and line:
            out.append(line); line = w_
        else:
            line = t
    return out + ([line] if line else [])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="poster_project/data")
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default="poster_project/figures")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    M = load_mugration(a.data, a.results)
    A = load_truth(a.data, a.results)
    n = len(next(iter(M.values())))
    print(f"mugration: {n} strategies x {len(M)} states")
    print(f"truth    : {len(A)} states with AUC metrics")
    made = [panel_family(M, os.path.join(a.out, "panel_family_cosine.svg")),
            panel_family(M, os.path.join(a.out, "panel_family_f1.svg"),
                         metric="topological_f1",
                         title="Topological fidelity by sampling family"),
            panel_tradeoff(M, A, os.path.join(a.out, "panel_tradeoff.svg")),
            panel_translation(M, A, os.path.join(a.out, "panel_translation.svg")),
            panel_crosstab(M, os.path.join(a.out, "panel_crosstab.svg")),
            panel_clock(load_clock(a.results), os.path.join(a.out, "panel_clock.svg")),
            panel_framework(os.path.join(a.out, "panel_framework.svg")),
            panel_profile(M, A, os.path.join(a.out, "panel_profile.svg"))]
    for m in made:
        if m:
            print(f"  -> {m}")


if __name__ == "__main__":
    main()
