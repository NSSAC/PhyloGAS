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
        p.append(txt(kx + 24, 190, FAM_LABEL[f], 21, SERIES[i], weight="bold"))
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
            panel_tradeoff(M, A, os.path.join(a.out, "panel_tradeoff.svg"))]
    for m in made:
        if m:
            print(f"  -> {m}")


if __name__ == "__main__":
    main()
