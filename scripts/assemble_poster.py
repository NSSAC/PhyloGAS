#!/usr/bin/env python3
"""Assemble the ASM BIG 2026 PhyloGAS poster: one 4800x2400 SVG (48x24 in).

v2 layout, written for readers outside phylodynamics:

  column A  the question, a three-step picture of the problem, a key to the
            strategy codes, and how a strategy is scored -- native poster
            type, not a scaled-down figure
  column B  how the digital twin works, the headline result, references
  column C  why common sampling goals mislead, a validation box for the
            synthetic genomes, conclusions

Results panels come from poster_panels.py and are embedded as nested <svg>
elements, so every mark stays a vector shape. Every number printed in the
running text is computed here from the same data files the panels use.
Previous layouts are kept under poster_project/versions/.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import statistics as st

W, H = 4800, 2400
FONT = "Arial, Helvetica, sans-serif"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8984"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
T_BLUE, T_ORANGE, T_AQUA, T_YELLOW = "#2a78d6", "#c9501f", "#15855c", "#996900"
RULE = "#d9d8d4"
STATES = ["va", "ga", "ma", "mn", "wa"]
STATE_NAME = {"va": "Virginia", "ga": "Georgia", "ma": "Massachusetts",
              "mn": "Minnesota", "wa": "Washington"}

TITLE = ("PhyloGeographic Analysis Similars: measuring the impact of "
         "public health sequencing strategies")
AUTHORS = "A. Warren, M. Wilson, G. Kaur, J. Chen, B. Lewis, M. Marathe"
AFFIL = "Biocomplexity Institute, University of Virginia, Charlottesville, VA, USA"
SITE = "nssac.github.io/vanalytica/phylogas"

BODY, LEAD = 27, 37          # running text: ~19pt at 100 px/in
WARN = []


LOGO = "poster_project/assets/uva_bi_logo.png"


def logo(p, path, x, y, w, h):
    """Place a raster logo, fitted inside the box and right-aligned.

    This is the one raster element on the poster. It is base64-embedded rather
    than linked so the SVG stays a single portable file, and it is the reason
    the poster is no longer purely vector -- which matters only for converting
    to editable PowerPoint shapes, not for printing. At the placed width the
    supplied artwork lands near 300 dpi.
    """
    import base64, struct
    raw = open(path, "rb").read()
    if raw[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"{path}: expected a PNG")
    iw, ih = struct.unpack(">II", raw[16:24])
    scale = min(w / iw, h / ih)
    dw, dh = iw * scale, ih * scale
    dx, dy = x + w - dw, y + (h - dh) / 2            # right-aligned, centred
    dpi = iw / (dw / 100.0)
    if dpi < 150:
        WARN.append(f"logo renders at {dpi:.0f} dpi (want >=150)")
    b64 = base64.b64encode(raw).decode("ascii")
    p.append(f'<image x="{dx:.1f}" y="{dy:.1f}" width="{dw:.1f}" height="{dh:.1f}" '
             f'preserveAspectRatio="xMidYMid meet" '
             f'xlink:href="data:image/png;base64,{b64}" '
             f'href="data:image/png;base64,{b64}"/>')
    return dpi


SHOT = "poster_project/assets/nextstrain_importations_screenshot.png"
CONFERENCE = ("ASM Bioinformatics, Genomics and Big Data Conference (ASM BIG)  \u00b7  "
              "Washington, DC  \u00b7  October 11\u201314, 2026")


def raster(p, path, x, y, w):
    """Embed a PNG at width w, aspect preserved. Returns the bottom edge."""
    import base64, struct
    raw = open(path, "rb").read()
    iw, ih = struct.unpack(">II", raw[16:24])
    h = w * ih / iw
    dpi = iw / (w / 100.0)
    if dpi < 150:
        WARN.append(f"{os.path.basename(path)} renders at {dpi:.0f} dpi (want >=150)")
    b64 = base64.b64encode(raw).decode("ascii")
    p.append(f'<image x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
             f'preserveAspectRatio="xMidYMid meet" xlink:href="data:image/png;base64,{b64}" '
             f'href="data:image/png;base64,{b64}"/>')
    return y + h


def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def txt(x, y, s, size=BODY, fill=INK, anchor="start", weight="normal", style=""):
    extra = f' font-style="{style}"' if style else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
            f'fill="{fill}" text-anchor="{anchor}" font-weight="{weight}"{extra}>'
            f'{esc(s)}</text>')


def wrap(text, n):
    out, line = [], ""
    for word in text.split():
        t = (line + " " + word).strip()
        if len(t) > n and line:
            out.append(line); line = word
        else:
            line = t
    return out + ([line] if line else [])


def para(p, x, y, text, n, size=BODY, lead=LEAD, fill=INK2, weight="normal"):
    for i, line in enumerate(wrap(text, n)):
        p.append(txt(x, y + i * lead, line, size, fill, weight=weight))
    return y + len(wrap(text, n)) * lead


def heading(p, x, y, s, w=360):
    p.append(txt(x, y, s, 40, INK, weight="bold"))
    p.append(f'<rect x="{x}" y="{y+14}" width="{w}" height="5" fill="{BLUE}"/>')
    return y + 66


def embed(p, path, x, y, w):
    src = open(path).read()
    m = re.match(r'<svg[^>]*width="(\d+)"[^>]*height="(\d+)"[^>]*>', src)
    pw, ph = int(m.group(1)), int(m.group(2))
    body = src[m.end():src.rindex("</svg>")]
    h = w * ph / pw
    p.append(f'<svg x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
             f'viewBox="0 0 {pw} {ph}">{body}</svg>')
    return y + h


def check(col, y, limit):
    if y > limit:
        WARN.append(f"column {col} overflows by {y-limit:.0f}px")


# --------------------------------------------------------------------------
# data for the running text
# --------------------------------------------------------------------------
def headline(data_dir):
    M = {s: {r["label"]: r for r in csv.DictReader(
         open(os.path.join(data_dir, f"mugration_foil_{s}.csv")))} for s in STATES}
    d = [float(M[s]["4S-P__lasso_greedy"]["cosine_similarity"])
         - float(M[s]["4S__surs"]["cosine_similarity"]) for s in STATES]
    m, sd = st.mean(d), st.stdev(d)
    t = m / (sd / math.sqrt(len(d)))
    return m, sum(x > 0 for x in d), t


def baselines(data_dir):
    """Absolute cosine of the two random baselines, per schedule.

    Every result on the poster is a difference from these, so the methods
    section states them: "+0.072 over random" is meaningless without knowing
    what random scored.
    """
    M = {s: {r["label"]: r for r in csv.DictReader(
         open(os.path.join(data_dir, f"mugration_foil_{s}.csv")))} for s in STATES}
    return {k: st.mean(float(M[s][v]["cosine_similarity"]) for s in STATES)
            for k, v in (("weekly", "1S__surs"), ("4-weekly", "4S__surs"))}


def sweeps(data_dir):
    """How many non-random recipes beat random in every twin, and out of how
    many. Computed, not asserted: an earlier draft of this sentence said
    "four of 20" when it is three of 18."""
    M = {s: {r["label"]: r for r in csv.DictReader(
         open(os.path.join(data_dir, f"mugration_foil_{s}.csv")))} for s in STATES}
    base = lambda r: "1S__surs" if r.startswith("1S") else "4S__surs"
    R = [r for r in sorted(M[STATES[0]]) if not r.endswith("surs")]
    n = sum(all(float(M[s][r]["cosine_similarity"])
                > float(M[s][base(r)]["cosine_similarity"]) for s in STATES) for r in R)
    return n, len(R)


def clock(results_dir):
    out = {}
    for s in STATES:
        f = os.path.join(results_dir, f"{s}_delta_wave", "05_benchmarks",
                         "clock_estimates.csv")
        if not os.path.exists(f):
            continue
        rows = list(csv.DictReader(open(f)))
        tr = [r for r in rows if r["quantity"] == "mu_truth" and r["arm"] == "1S__surs"]
        sm = [r for r in rows if r["quantity"] == "mu_sim"
              and r.get("rate_subs_per_site_per_year")]
        if tr:
            out[s] = {"truth": float(tr[0]["rate_subs_per_site_per_year"]),
                      "tips": int(float(tr[0]["tips"] or 0)),
                      "sim": float(sm[0]["rate_subs_per_site_per_year"]) if sm else None}
    return out


# --------------------------------------------------------------------------
# the three-step picture of the problem (illustrative, not data)
# --------------------------------------------------------------------------
DOTS = {"A": [(60, 80), (100, 70), (150, 85), (80, 130), (130, 125), (170, 140),
              (70, 175), (120, 175)],
        "B": [(270, 75), (320, 70), (370, 90), (285, 130), (335, 125), (380, 150),
              (300, 175), (350, 180)],
        "C": [(165, 285), (215, 280), (265, 290), (180, 330), (230, 335), (275, 340),
              (200, 380), (250, 385)]}
PTS = DOTS["A"] + DOTS["B"] + DOTS["C"]
EDGES = [(0, 1), (1, 2), (1, 4), (0, 3), (3, 6), (4, 7), (4, 5),
         (8, 9), (9, 10), (9, 12), (8, 11), (11, 14), (12, 15), (12, 13),
         (16, 17), (16, 19), (19, 22), (17, 20),      # C chain seeded from A
         (18, 21), (21, 23)]                            # C chain seeded from B
BRIDGES = [(5, 11), (7, 16), (14, 18)]          # A->B, A->C, B->C
SEQUENCED = {4, 5, 9, 11, 16, 20}               # B->C bridge left unsampled
BOXES = {"A": (30, 42, 170, 165), "B": (240, 42, 170, 165), "C": (135, 252, 170, 165)}
CENTRE = {"A": (115, 122), "B": (325, 122), "C": (220, 332)}


def counties(p, ox, oy):
    for k, (x, y, w, h) in BOXES.items():
        p.append(f'<rect x="{ox+x}" y="{oy+y}" width="{w}" height="{h}" rx="14" '
                 f'fill="#f5f4f1" stroke="{RULE}" stroke-width="2"/>')
        if k == "C":   # the two bridges into C pass above it; label inside, bottom-right
            p.append(txt(ox + x + w - 8, oy + y + h - 8, f"county {k}", 20, INK3,
                         anchor="end"))
        else:
            p.append(txt(ox + x + 6, oy + y - 9, f"county {k}", 20, INK3))


def arrow(p, x1, y1, x2, y2, colour, dashed=False, w=5):
    dx, dy = x2 - x1, y2 - y1
    L = math.hypot(dx, dy); ux, uy = dx / L, dy / L
    bx, by = x2 - ux * 18, y2 - uy * 18
    dash = ' stroke-dasharray="10 8"' if dashed else ""
    p.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{bx:.1f}" y2="{by:.1f}" '
             f'stroke="{colour}" stroke-width="{w}"{dash}/>')
    p.append(f'<path d="M {x2:.1f} {y2:.1f} L {bx-uy*10:.1f} {by+ux*10:.1f} '
             f'L {bx+uy*10:.1f} {by-ux*10:.1f} Z" fill="{colour}"/>')


def storyboard(p, x0, y0):
    FW, GAP, FH = 440, 40, 440
    titles = ["1  What really happened", "2  What gets sequenced",
              "3  What the tree can see"]
    for i, t in enumerate(titles):
        p.append(txt(x0 + i * (FW + GAP), y0, t, 28, INK, weight="bold"))
    oy = y0 + 22
    for i in range(3):
        ox = x0 + i * (FW + GAP)
        p.append(f'<rect x="{ox}" y="{oy}" width="{FW}" height="{FH}" rx="8" '
                 f'fill="#ffffff" stroke="{RULE}" stroke-width="2"/>')
        if i < 2:
            counties(p, ox, oy)
    # frame 1: every infection and transmission
    ox = x0
    for a, b in EDGES:
        (x1, y1), (x2, y2) = PTS[a], PTS[b]
        p.append(f'<line x1="{ox+x1}" y1="{oy+y1}" x2="{ox+x2}" y2="{oy+y2}" '
                 f'stroke="{INK3}" stroke-width="2.5"/>')
    for a, b in BRIDGES:
        (x1, y1), (x2, y2) = PTS[a], PTS[b]
        p.append(f'<line x1="{ox+x1}" y1="{oy+y1}" x2="{ox+x2}" y2="{oy+y2}" '
                 f'stroke="{ORANGE}" stroke-width="5"/>')
    for x, y in PTS:
        p.append(f'<circle cx="{ox+x}" cy="{oy+y}" r="10" fill="{INK2}" '
                 f'stroke="#ffffff" stroke-width="2"/>')
    # frame 2: only a few are sequenced
    ox = x0 + FW + GAP
    for a, b in EDGES + BRIDGES:
        (x1, y1), (x2, y2) = PTS[a], PTS[b]
        p.append(f'<line x1="{ox+x1}" y1="{oy+y1}" x2="{ox+x2}" y2="{oy+y2}" '
                 f'stroke="#ebeae6" stroke-width="2.5"/>')
    for i, (x, y) in enumerate(PTS):
        if i in SEQUENCED:
            p.append(f'<circle cx="{ox+x}" cy="{oy+y}" r="13" fill="{BLUE}" '
                     f'stroke="#ffffff" stroke-width="3"/>')
        else:
            p.append(f'<circle cx="{ox+x}" cy="{oy+y}" r="9" fill="#dedcd7"/>')
    # frame 3: the county-to-county links the tree recovers
    ox = x0 + 2 * (FW + GAP)
    for k, (cx, cy) in CENTRE.items():
        p.append(f'<circle cx="{ox+cx}" cy="{oy+cy}" r="58" fill="#f5f4f1" '
                 f'stroke="{RULE}" stroke-width="2"/>')
        p.append(txt(ox + cx, oy + cy + 10, k, 32, INK2, anchor="middle", weight="bold"))
    def link(a, b, found):
        (ax, ay), (bx, by) = CENTRE[a], CENTRE[b]
        L = math.hypot(bx - ax, by - ay); ux, uy = (bx - ax) / L, (by - ay) / L
        arrow(p, ox + ax + ux * 62, oy + ay + uy * 62, ox + bx - ux * 64,
              oy + by - uy * 64, BLUE if found else INK3, dashed=not found)
    link("A", "B", True); link("A", "C", True); link("B", "C", False)
    arrow(p, ox + 40, oy + 404, ox + 96, oy + 404, BLUE, w=4)
    p.append(txt(ox + 108, oy + 411, "found", 21, INK2))
    arrow(p, ox + 215, oy + 404, ox + 271, oy + 404, INK3, dashed=True, w=4)
    p.append(txt(ox + 283, oy + 411, "missed", 21, INK2))
    return oy + FH


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figures", default="poster_project/figures")
    ap.add_argument("--data", default="poster_project/data")
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default=None)
    ap.add_argument("--layout", choices=["main", "alt"], default="main",
                    help="main: the v9 poster. alt: conference line in the header, "
                         "Nextstrain screenshot in the middle column, cross-tab moved to "
                         "the top of the right column, conclusions under the screenshot.")
    a = ap.parse_args()
    if a.out is None:
        a.out = ("poster_project/poster_alt.svg" if a.layout == "alt"
                 else "poster_project/poster.svg")
    F = lambda n: os.path.join(a.figures, n)
    hd_m, hd_w, hd_t = headline(a.data)
    BASE = baselines(a.data)
    SWEEP_N, SWEEP_TOT = sweeps(a.data)
    C = clock(a.results)
    truths = [v["truth"] for v in C.values()]
    ratio = 8e-4 / st.mean(truths) if truths else None

    p = [f'<svg xmlns="http://www.w3.org/2000/svg" '
         f'xmlns:xlink="http://www.w3.org/1999/xlink" width="{W}" height="{H}" '
         f'viewBox="0 0 {W} {H}">\n<rect width="{W}" height="{H}" fill="#ffffff"/>\n']

    # ---- header ----------------------------------------------------------
    p.append(f'<rect x="0" y="0" width="{W}" height="16" fill="{BLUE}"/>')
    for i, line in enumerate(wrap(TITLE, 80)):
        p.append(txt(80, 112 + i * 70, line, 64, INK, weight="bold"))
    p.append(txt(80, 268, AUTHORS, 34, INK2))
    p.append(txt(80, 312, AFFIL, 30, INK3))
    if os.path.exists(LOGO):
        logo(p, LOGO, W - 760, 44, 680, 190)
    else:
        WARN.append(f"logo not found at {LOGO}; placeholder drawn")
        p.append(f'<rect x="{W-760}" y="44" width="680" height="190" rx="10" fill="none" '
                 f'stroke="{RULE}" stroke-width="3" stroke-dasharray="14 10"/>')
        p.append(txt(W - 420, 148, "institutional logos", 28, INK3, anchor="middle",
                     style="italic"))
    p.append(txt(W - 80, 290, SITE, 32, BLUE, anchor="end", weight="bold"))
    xD = 2560
    p.append(txt(xD, 76, "Data and code \u2014 all openly available", 28, INK, weight="bold"))
    p.append(f'<rect x="{xD}" y="88" width="560" height="4" fill="{BLUE}"/>')
    rows = [
        ("Code", "github.com/NSSAC/PhyloGAS", True),
        ("", "github.com/NSSAC/BeyondBaseline (sampling)  \u00b7  "
             "github.com/NSSAC/TwinSampler (ascertainment)", False),
        ("Simulations", "EpiHiper transmission replicates  doi:10.5281/zenodo.23067670", False),
        ("Populations", "VA 10.5281/zenodo.23088534   GA 10.5281/zenodo.23088629   "
                        "MA 10.5281/zenodo.23088723", False),
        ("", "MN 10.5281/zenodo.23088773   WA 10.5281/zenodo.23088830", False),
        ("", "synthetic populations mirrored from UVA Dataverse  doi:10.18130/V3/5LSDCY", False),
    ]
    yy = 128
    for lab, val, strong in rows:
        if lab:
            p.append(txt(xD, yy, lab, 23, INK, weight="bold"))
        p.append(txt(xD + 170, yy, val, 23, BLUE if strong else INK2,
                     weight="bold" if strong else "normal"))
        yy += 34
    p.append(f'<line x1="80" y1="356" x2="{W-80}" y2="356" stroke="{RULE}" '
             f'stroke-width="3"/>')
    TOP, BOT = 420, 2350

    # ---- column A --------------------------------------------------------
    xA, wA, nA = 80, 1400, 100
    y = heading(p, xA, TOP, "The question")
    y = para(p, xA, y, "Public-health labs can sequence only a small share of cases "
             "— in these outbreaks, fewer than 1 infection in 300. Does it matter "
             "which cases are chosen? In a real outbreak nobody can check, because the "
             "true path of spread is never known. In a simulated outbreak, it can be.",
             nA, fill=INK)
    y += 28
    y = storyboard(p, xA, y + 10)
    y += 34
    y = para(p, xA, y, "Illustration. Dots are infections, lines are transmissions, "
             "orange lines cross a county boundary. The phylogenetic tree is built only "
             "from sequenced cases, so which county-to-county links it can recover "
             "depends on which cases were picked.", nA, size=24, lead=33, fill=INK3)
    y += 30
    y = heading(p, xA, y, "What we compared")
    y = para(p, xA, y, "Each week, or every four weeks, a strategy chooses which "
             "reported cases to sequence:", nA)
    y += 10
    algos = [(BLUE, T_BLUE, "Model-guided, free choice", "LASSO-Greedy",
              "a statistical model scores each case, then cases are added one at a "
              "time to bring the sample closest to the target mix"),
             (AQUA, T_AQUA, "Model-guided, with quota", "LASSO-Stratified",
              "the same model, but filling fixed quotas per group"),
             (YELLOW, T_YELLOW, "Quota only", "Stratified",
              "fixed numbers of cases per age × race group"),
             (ORANGE, T_ORANGE, "Random", "SURS",
              "every reported case equally likely — the usual baseline")]
    for fill, tcol, name, code, desc in algos:
        p.append(f'<rect x="{xA}" y="{y-24}" width="14" height="62" fill="{fill}"/>')
        p.append(txt(xA + 30, y, name, 26, tcol, weight="bold"))
        p.append(txt(xA + 30, y + 32, code, 22, INK3))
        y2 = para(p, xA + 450, y, desc, 70, size=25, lead=33)
        y = max(y + 74, y2 + 18)
    y += 6
    y = para(p, xA, y, "… aiming for a sample that resembles:", nA)
    for code, desc in [("P \u2014 population", "everyone, tested or not; corrects for "
                                           "who gets tested"),
                       ("LL \u2014 line list", "the cases that were reported"),
                       ("LL+P", "half of each")]:
        p.append(txt(xA + 30, y, code, 26, INK, weight="bold"))
        p.append(txt(xA + 330, y, desc, 25, INK2))
        y += 38
    y += 6
    y = para(p, xA, y, "… on a schedule of:", nA)
    for code, desc in [("1S \u2014 1-week stride", "weekly, choosing from that week's cases"),
                       ("4S \u2014 4-week stride", "every four weeks, from a pool of the "
                                               "last four weeks")]:
        p.append(txt(xA + 30, y, code, 26, INK, weight="bold"))
        p.append(txt(xA + 330, y, desc, 25, INK2))
        y += 38
    y += 4
    y = para(p, xA, y, "Combinations give 20 recipes (random needs no target), each run "
             "in five state twins: VA, GA, MA, MN, WA.", nA, size=24, lead=33, fill=INK3)
    y += 30
    y = heading(p, xA, y, "How a strategy is scored")
    y = para(p, xA, y, "Nextstrain builds a tree from the sequenced cases and infers how "
             "often the virus moved between each pair of counties. We compare those "
             "shares with the true shares from the simulation using cosine similarity: "
             "1 means identical, 0 means unrelated.", nA)
    y += 2
    y = para(p, xA, y, "Every result shown is the change from random sampling on the "
             "same schedule: zero means no better than picking at random. Random "
             f"scored {BASE['weekly']:.2f} weekly and {BASE['4-weekly']:.2f} "
             f"four-weekly \u2014 only {SWEEP_N} of {SWEEP_TOT} recipes beat it in all "
             "five twins.", nA, fill=INK)
    check("A", y, BOT)

    if a.layout == "alt":
        # ---- column B (alt): framework, the twin in Nextstrain, conclusions, refs
        xB, wB = 1560, 1600
        y = embed(p, F("panel_framework.svg"), xB, TOP - 30, wB)
        y = raster(p, SHOT, xB + 60, y + 14, wB - 120)      # ~230 dpi at this width
        y = para(p, xB + 35, y + 28, "Virginia Delta-wave twin: importation locations from "
                 "a previous analysis (map, per-county z-score of introductions vs. sampling), "
                 "the Nextstrain time tree, and inferred county-to-county transmissions.",
                 128, size=21, lead=28, fill=INK3)
        y = heading(p, xB, y + 18, "Conclusions")
        have_sim = {s_: v for s_, v in C.items() if v["sim"]}
        pval = "p < 0.01" if hd_t > 4.604 else ("p < 0.05" if hd_t > 2.776 else "n.s.")
        gaps = [abs(v["sim"] - v["truth"]) / v["truth"] * 100 for v in have_sim.values()]
        concl = [
            f"Picking cases so the sample resembles the whole population, with the algorithm "
            f"free to choose from a four-week pool, reconstructed county-to-county spread "
            f"best: {hd_m:+.3f} over random, better in {hd_w} of 5 state twins (paired t = "
            f"{hd_t:.1f}, {pval}).",
            "Rigid demographic quotas lost to random in every combination tested. The same "
            "target helps or hurts depending on whether the algorithm is free to choose.",
            "Common surveillance goals did not explain the winners: tree coverage barely "
            "differed between recipes, and reaching more clusters did not predict success. "
            "The winners sampled cases least like the current case mix.",
        ]
        if gaps:
            concl.append(
                f"Nextstrain recovered the simulation's mutation rate from a single random "
                f"sample to within {max(gaps):.1f}%"
                + (f" ({len(gaps)} of 5 twins)" if len(gaps) < 5 else "")
                + ": a readable molecular clock"
                + (f", though about {ratio:.1f}\u00d7 slower than real SARS-CoV-2."
                   if ratio else "."))
        for i, c in enumerate(concl, 1):
            p.append(f'<circle cx="{xB+18}" cy="{y-9}" r="19" fill="{BLUE}"/>')
            p.append(txt(xB + 18, y, str(i), 26, "#ffffff", anchor="middle", weight="bold"))
            y = para(p, xB + 56, y, c, 118, size=25, lead=33, fill=INK)
            y += 12
        y = heading(p, xB, y + 6, "References", w=260)
        for r in ["Chen J, Hoops S, Mortveit HS, et al. EpiHiper \u2014 a high performance "
                  "computational modeling framework to support epidemic science. "
                  "PNAS Nexus. 2025;4(1):pgae557.",
                  "Hadfield J, et al. Nextstrain: real-time tracking of pathogen evolution. "
                  "Bioinformatics. 2018;34(23):4121\u20134123.",
                  "Sagulenko P, Puller V, Neher RA. TreeTime: maximum-likelihood "
                  "phylodynamic analysis. Virus Evolution. 2018;4(1):vex042."]:
            y = para(p, xB, y, r, 128, size=21, lead=27, fill=INK3)
            y += 2
        check("B", y, BOT)

        # ---- column C (alt): cross-tab on top, then the per-recipe profile
        xC, wC = 3240, 1480
        y = embed(p, F("panel_crosstab_compact.svg"), xC, TOP - 30, wC)
        y = embed(p, F("panel_profile_compact.svg"), xC, y + 20, wC)
        check("C", y, BOT)
    else:
        # ---- column B --------------------------------------------------------
        xB, wB, nB = 1560, 1600, 108
        y = embed(p, F("panel_framework.svg"), xB, TOP - 30, wB)
        y = embed(p, F("panel_crosstab.svg"), xB, y + 10, wB)
        y = heading(p, xB, y + 20, "Do the synthetic genomes behave like real ones?", w=420)
        have_sim = {s: v for s, v in C.items() if v["sim"]}
        tips = max((v["tips"] for v in C.values()), default=0)
        y = para(p, xB, y, f"Each twin's true mutation rate was fitted on {tips:,} synthetic "
                 "genomes and their dates, with no tree involved. Nextstrain, given one random "
                 "sample and no rate, estimated it from its own tree:", nB, size=25, lead=34)
        y += 6
        for s, v in have_sim.items():
            gap = abs(v["sim"] - v["truth"]) / v["truth"] * 100
            p.append(txt(xB + 20, y + 4, STATE_NAME[s], 27, INK, weight="bold"))
            p.append(txt(xB + 300, y + 4, f"true {v['truth']*1e4:.2f}\u00d710\u207b\u2074", 27, T_BLUE))
            p.append(txt(xB + 650, y + 4, f"Nextstrain {v['sim']*1e4:.2f}\u00d710\u207b\u2074", 27, T_AQUA))
            p.append(txt(xB + 1080, y + 4, f"{gap:.1f}% apart", 27, INK, weight="bold"))
            y += 42
        if not have_sim:
            WARN.append("no mu_sim rows: validation box has no comparison")
        if ratio:
            y = para(p, xB, y + 4, f"The twin's rate is about {ratio:.1f}\u00d7 slower than real "
                     "SARS-CoV-2 (\u22488\u00d710\u207b\u2074): fewer mutations give the tree less "
                     "to work with. Calibrating it upward is ongoing work.", nB, size=23,
                     lead=31, fill=INK3)
        y = heading(p, xB, y + 10, "References", w=260)
        for r in ["Chen J, Hoops S, Mortveit HS, et al. EpiHiper \u2014 a high performance "
                  "computational modeling framework to support epidemic science. "
                  "PNAS Nexus. 2025;4(1):pgae557.",
                  "Hadfield J, et al. Nextstrain: real-time tracking of pathogen evolution. "
                  "Bioinformatics. 2018;34(23):4121\u20134123.",
                  "Sagulenko P, Puller V, Neher RA. TreeTime: maximum-likelihood phylodynamic "
                  "analysis. Virus Evolution. 2018;4(1):vex042."]:
            ph = r.startswith("[ADD")
            y = para(p, xB, y, r, 128, size=21, lead=27,
                     fill="#b4322f" if ph else INK3, weight="bold" if ph else "normal")
            y += 2
        check("B", y, BOT)

        # ---- column C --------------------------------------------------------
        xC, wC = 3240, 1480
        y = embed(p, F("panel_profile.svg"), xC, TOP - 30, wC)
        y = heading(p, xC, y + 10, "Conclusions")
        have_sim = {s: v for s, v in C.items() if v["sim"]}
        pval = "p < 0.01" if hd_t > 4.604 else ("p < 0.05" if hd_t > 2.776 else "n.s.")
        gaps = [abs(v["sim"] - v["truth"]) / v["truth"] * 100 for v in have_sim.values()]
        concl = [
            f"Picking cases so the sample resembles the whole population, with the algorithm "
            f"free to choose from a four-week pool, reconstructed county-to-county spread best: "
            f"{hd_m:+.3f} over random, better in {hd_w} of 5 state twins (paired t = "
            f"{hd_t:.1f}, {pval}).",
            "Rigid demographic quotas lost to random in every combination tested. The same "
            "target helps or hurts depending on whether the algorithm is free to choose.",
            "Common surveillance goals did not explain the winners: tree coverage barely "
            "differed between recipes, and reaching more clusters did not predict success. "
            "The winners sampled cases least like the current case mix.",
        ]
        if gaps:
            concl.append(f"Nextstrain recovered the simulation's mutation rate from a single "
                         f"random sample to within {max(gaps):.1f}%"
                         + (f" ({len(gaps)} of 5 twins)" if len(gaps) < 5 else "")
                         + ", so the synthetic genomes carry a readable molecular clock.")
        for i, c in enumerate(concl, 1):
            p.append(f'<circle cx="{xC+18}" cy="{y-9}" r="19" fill="{BLUE}"/>')
            p.append(txt(xC + 18, y, str(i), 26, "#ffffff", anchor="middle", weight="bold"))
            y = para(p, xC + 56, y, c, 94, size=25, lead=33, fill=INK)
            y += 14
        check("C", y, BOT)

    p.append(f'<rect x="0" y="{H-12}" width="{W}" height="12" fill="{BLUE}"/>')
    p.append("</svg>\n")
    open(a.out, "w").write("".join(p))
    print(f"-> {a.out}")
    for w in WARN:
        print("WARNING:", w)


if __name__ == "__main__":
    main()
