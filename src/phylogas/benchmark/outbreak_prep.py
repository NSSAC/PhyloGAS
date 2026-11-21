"""Prepare focal outbreak benchmarks and orchestrate multi-method scoring.

Builds realistic, genetically pre-isolated local outbreak datasets across states
from PhyloGAS's full painted sequence sets, runs Nextstrain (``augur tree`` +
``augur refine``) to produce dated trees for ``TransPhylo``, and aggregates
phylogenetic, cluster-triage, and who-infected-whom transmission metrics across:

  - ``MicrobeTrace`` (single-linkage SNP thresholds + ``SNP + date`` baseline)
  - ``Nextstrain`` (``IQ-TREE 2`` + ``TreeTime``)
  - ``outbreaker2`` (Campbell et al. 2018)
  - ``TransPhylo`` (Didelot et al. 2017, layered on ``TreeTime`` dated tree)
  - ``JUNIPER`` (Specht et al. 2025, consensus mode ``juniper0``)
"""

from __future__ import annotations

import csv
import datetime as _dt
import json
import math
import random
import shutil
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

from . import reconstruction as recon


def _read_linelist_dates(linelist_path: str | Path | None) -> dict[str, _dt.date]:
    """Map ``alias_pid -> collection_date`` for cases ascertained in the line list."""
    out: dict[str, _dt.date] = {}
    if not linelist_path or not Path(linelist_path).exists():
        return out
    with recon._open_text(linelist_path, "rt") as fh:
        hdr = fh.readline()
        if not hdr:
            return out
        delim = "\t" if hdr.count("\t") > hdr.count(",") else ","
        cols = [c.strip() for c in hdr.rstrip("\r\n").split(delim)]
        idx = {c: i for i, c in enumerate(cols)}
        ia = idx.get("alias_pid")
        idt = idx.get("date", idx.get("collection_date"))
        if ia is None or idt is None:
            return out
        for row in csv.reader(fh, delimiter=delim):
            if len(row) <= max(ia, idt):
                continue
            al = row[ia].strip()
            ds = row[idt].strip()[:10]
            if al and len(ds) == 10 and ds[4] == "-":
                try:
                    out[al] = _dt.date.fromisoformat(ds)
                except ValueError:
                    pass
    return out


def _read_clock_truth(bench_dir: Path) -> float:
    """Read ``mu_truth`` (or ``mu_lineage``) in subs/site/year from ``clock_estimates.csv``."""
    clk_csv = bench_dir / "clock_estimates.csv"
    if clk_csv.exists():
        with open(clk_csv, "r", encoding="utf-8") as fh:
            best = None
            for r in csv.DictReader(fh):
                q = (r.get("quantity") or "").strip()
                val = r.get("rate_subs_per_site_per_year")
                if q == "mu_truth" and val:
                    try:
                        return float(val)
                    except ValueError:
                        pass
                if q == "mu_lineage" and val and best is None:
                    try:
                        best = float(val)
                    except ValueError:
                        pass
            if best is not None:
                return best
    return 6.84e-4


def _read_gi_mean(bench_dir: Path) -> float:
    seq_csv = bench_dir / "sequence_divergence.csv"
    if seq_csv.exists():
        with open(seq_csv, "r", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                v = r.get("gi_mean_days")
                if v:
                    try:
                        return float(v)
                    except ValueError:
                        pass
    return 4.673


def select_focal_outbreak_windows(
    infections_path: str | Path,
    n_clusters: int = 3,
    window_days: int = 56,
    step_days: int = 14,
    min_true_cases: int = 150,
    max_true_cases: int = 440,
    min_dominant_local_chain: int = 35,
) -> list[dict[str, Any]]:
    """Scan ``allevents_variants.csv.xz`` once and select ``n_clusters`` distinct
    county (or sub-county grid cell for populous counties like MA) 8-week windows
    with a genuine local transmission core and manageable total infection count."""
    records: list[tuple[str, str, str, int, str, str, str]] = []
    # (alias_pid, parent_alias, county_fips, exp_ord, comp_id, grid_cell, strain)
    min_ord = 10**9
    max_ord = 0

    with recon._open_text(infections_path, "rt") as fh:
        hdr = fh.readline()
        delim = "\t" if hdr.count("\t") > hdr.count(",") else ","
        cols = [c.strip() for c in hdr.rstrip("\r\n").split(delim)]
        idx = {c: i for i, c in enumerate(cols)}
        ia = idx["alias_pid"]
        ip = idx.get("alias_contact")
        icf = idx.get("county_fips", idx.get("FIPS"))
        idt = idx.get("exposure_date", idx.get("date"))
        ic = idx.get("component_id")
        ilat = idx.get("latitude", idx.get("home_latitude"))
        ilon = idx.get("longitude", idx.get("home_longitude"))
        ist = idx.get("strain")

        seen_al: set[str] = set()
        for row in csv.reader(fh, delimiter=delim):
            if len(row) <= ia:
                continue
            al = row[ia].strip()
            if not al or al in seen_al:
                continue
            seen_al.add(al)
            ds = row[idt].strip()[:10] if idt is not None and len(row) > idt else ""
            if len(ds) != 10 or ds[4] != "-":
                continue
            try:
                eo = _dt.date.fromisoformat(ds).toordinal()
            except ValueError:
                continue
            if eo < min_ord:
                min_ord = eo
            if eo > max_ord:
                max_ord = eo

            p = row[ip].strip() if ip is not None and len(row) > ip else ""
            if p in recon._NO_VAL or p == al:
                p = ""
            cf = row[icf].strip().split(".")[0] if icf is not None and len(row) > icf else "00000"
            comp = row[ic].strip() if ic is not None and len(row) > ic else ""
            st = row[ist].strip() if ist is not None and len(row) > ist else ""

            grid = cf
            if ilat is not None and ilon is not None and len(row) > max(ilat, ilon):
                try:
                    lat_g = round(float(row[ilat]), 1)
                    lon_g = round(float(row[ilon]), 1)
                    grid = f"{cf}_{lat_g:.1f}_{lon_g:.1f}"
                except ValueError:
                    pass
            records.append((al, p, cf, eo, comp, grid, st))

    def _score_grouping(use_grid: bool, lo_n: int, hi_n: int, min_dom: int) -> list[dict[str, Any]]:
        by_loc: dict[str, list[tuple[int, str, str, str, str, str]]] = defaultdict(list)
        for al, p, cf, eo, comp, grid, st in records:
            key = grid if use_grid else cf
            by_loc[key].append((eo, al, p, cf, comp, st))

        candidates: list[dict[str, Any]] = []
        for loc_key, items in by_loc.items():
            if len(items) < lo_n:
                continue
            items.sort(key=lambda x: x[0])
            t_start = min_ord + 14 if (max_ord - min_ord) >= (window_days + 14) else min_ord
            t_limit = max(t_start + window_days, max_ord + 7)
            best_for_loc: dict[str, Any] | None = None
            while t_start + window_days <= t_limit:
                t_end = t_start + window_days
                win_items = [x for x in items if t_start <= x[0] < t_end]
                n_tot = len(win_items)
                if lo_n <= n_tot <= hi_n:
                    wset = {x[1] for x in win_items}
                    ch_map: dict[str, list[str]] = defaultdict(list)
                    roots_local: list[str] = []
                    for _, al, p, _, _, _ in win_items:
                        if p in wset:
                            ch_map[p].append(al)
                        else:
                            roots_local.append(al)
                    sub_sizes: list[int] = []
                    for r0 in roots_local:
                        q = [r0]
                        for u in q:
                            q.extend(ch_map.get(u, ()))
                        sub_sizes.append(len(q))
                    sub_sizes.sort(reverse=True)
                    dom = sub_sizes[0] if sub_sizes else 0
                    if dom >= min_dom:
                        # Prefer outbreaks with a strong dominant local chain (35%..85% share)
                        # and target size near ~260 true cases (~130 at pi=0.5, ~65 at pi=0.25)
                        dom_share = dom / float(n_tot)
                        size_pen = abs(n_tot - 260) / 260.0
                        score = dom_share - 0.35 * size_pen
                        cand = {
                            "loc_key": loc_key,
                            "county_fips": win_items[0][3],
                            "use_grid": use_grid,
                            "start_ord": t_start,
                            "end_ord": t_end,
                            "window_start": _dt.date.fromordinal(t_start).isoformat(),
                            "window_end": _dt.date.fromordinal(t_end).isoformat(),
                            "n_true_window": n_tot,
                            "n_local_intros": len(roots_local),
                            "dominant_local_size": dom,
                            "top_local_sizes": sub_sizes[:5],
                            "score": score,
                            "members": [
                                {
                                    "alias_pid": al,
                                    "parent_alias": p,
                                    "county_fips": cf,
                                    "exp_ord": eo,
                                    "exposure_date": _dt.date.fromordinal(eo).isoformat(),
                                    "component_id": comp,
                                    "strain": st,
                                }
                                for (eo, al, p, cf, comp, st) in win_items
                            ],
                        }
                        if best_for_loc is None or cand["score"] > best_for_loc["score"]:
                            best_for_loc = cand
                t_start += step_days
            if best_for_loc is not None:
                candidates.append(best_for_loc)

        candidates.sort(key=lambda c: (-c["score"], c["loc_key"]))
        return candidates

    # Tier 1: Whole-county 56-day windows
    chosen: list[dict[str, Any]] = []
    used_counties: set[str] = set()
    for cand in _score_grouping(False, min_true_cases, max_true_cases, min_dominant_local_chain):
        if cand["county_fips"] not in used_counties:
            chosen.append(cand)
            used_counties.add(cand["county_fips"])
            if len(chosen) >= n_clusters:
                return chosen

    # Tier 2: Sub-county ~10km neighborhood windows (essential for 14-county MA & populous WA counties)
    for cand in _score_grouping(True, min_true_cases, max_true_cases, min_dominant_local_chain):
        if cand["loc_key"] not in {c["loc_key"] for c in chosen}:
            # Prefer distinct counties first, allow distinct sub-county cells if needed
            if cand["county_fips"] not in used_counties or len(chosen) < n_clusters:
                chosen.append(cand)
                used_counties.add(cand["county_fips"])
                if len(chosen) >= n_clusters:
                    return chosen

    # Tier 3: Relaxed bounds if still short
    for cand in _score_grouping(True, 90, 600, 20):
        if cand["loc_key"] not in {c["loc_key"] for c in chosen}:
            chosen.append(cand)
            if len(chosen) >= n_clusters:
                return chosen

    return chosen


def _isolate_focal_genetic_cluster(
    sampled_aliases: list[str],
    all_win_seqs: dict[str, str],
    dominant_comp: str,
    comp_of_win: dict[str, str],
    max_snps: int = 3,
) -> tuple[list[str], list[str]]:
    """Mimic operational outbreak pre-isolation: connect sampled sequences within
    ``<= max_snps`` SNPs, identify the focal outbreak component (preferring the
    component overlapping ``dominant_comp``), and retain any sequences within
    ``<= max_snps`` SNPs of that focal cluster (naturally keeping 0-3 SNP cousin
    re-introductions while dropping distant 10-20 SNP unrelated lineages!).

    Returns ``(isolated_sampled_aliases, compatible_true_window_aliases)`` so the
    effective local sampling rate ``pi_truth`` is measured accurately.
    """
    present = [a for a in sampled_aliases if a in all_win_seqs]
    m = len(present)
    if m <= 2:
        return present, list(all_win_seqs.keys())

    parent_uf = list(range(m))

    def _find(x: int) -> int:
        while parent_uf[x] != x:
            parent_uf[x] = parent_uf[parent_uf[x]]
            x = parent_uf[x]
        return x

    for i in range(m):
        si = all_win_seqs[present[i]]
        for j in range(i + 1, m):
            sj = all_win_seqs[present[j]]
            d = sum(1 for a, b in zip(si, sj) if a != b and a in "ACGT" and b in "ACGT")
            if d <= max_snps:
                ri, rj = _find(i), _find(j)
                if ri != rj:
                    parent_uf[max(ri, rj)] = min(ri, rj)

    groups: dict[int, list[str]] = defaultdict(list)
    for i, al in enumerate(present):
        groups[_find(i)].append(al)

    # Pick largest genetic cluster (breaking ties toward dominant_comp overlap)
    ranked = sorted(
        groups.values(),
        key=lambda g: (sum(1 for a in g if comp_of_win.get(a) == dominant_comp), len(g)),
        reverse=True,
    )
    focal_sampled = sorted(ranked[0])

    # If the focal genetic cluster has < 25 tips (e.g. a slightly wider multi-SNP burst),
    # widen threshold to 5 SNPs so we don't over-fragment a single 56-day wave
    if len(focal_sampled) < min(25, int(0.6 * m)) and max_snps < 5:
        return _isolate_focal_genetic_cluster(
            sampled_aliases, all_win_seqs, dominant_comp, comp_of_win, max_snps=5
        )

    # Determine which true infections in the window belong to this genetic cluster
    # (within <= max_snps of at least one focal_sampled sequence)
    focal_seqs = [all_win_seqs[a] for a in focal_sampled[:40]]
    compatible_true: list[str] = []
    for al, seq in all_win_seqs.items():
        if al in focal_sampled:
            compatible_true.append(al)
            continue
        for fseq in focal_seqs:
            d = sum(1 for a, b in zip(seq, fseq) if a != b and a in "ACGT" and b in "ACGT")
            if d <= max_snps:
                compatible_true.append(al)
                break

    return focal_sampled, compatible_true


def prepare_state_outbreaks(
    state: str,
    infections_path: str | Path,
    painted_fasta: str | Path,
    reference_fasta: str | Path,
    outdir: str | Path,
    linelist_path: str | Path | None = None,
    bench_dir: str | Path | None = None,
    n_clusters: int = 3,
    sampling_rates: Iterable[float] = (0.25, 0.50),
    max_cluster_snps: int = 3,
    seed: int = 2026,
) -> list[dict[str, Any]]:
    """Extract all focal outbreak datasets for ``state`` in a **single streaming
    pass** over ``painted_fasta``."""
    out_root = Path(outdir)
    out_root.mkdir(parents=True, exist_ok=True)
    bdir = Path(bench_dir) if bench_dir else Path(infections_path).parent

    mu_truth_yr = _read_clock_truth(bdir)
    gi_mean = _read_gi_mean(bdir)
    # Moment-matched Gamma GI (mean = gi_mean, sd ~ 1.48d) & sojourn (mean = 3.975d, sd = 0.65d)
    gi_sd = 1.48
    a_g = round((gi_mean / gi_sd) ** 2, 4)
    lambda_g = round(gi_mean / (gi_sd ** 2), 4)
    soj_mean, soj_sd = 3.975, 0.65
    a_s = round((soj_mean / soj_sd) ** 2, 4)
    lambda_s = round(soj_mean / (soj_sd ** 2), 4)

    print(f"[{state}] Scanning {infections_path} for {n_clusters} focal county outbreaks...")
    windows = select_focal_outbreak_windows(infections_path, n_clusters=n_clusters)
    if not windows:
        raise RuntimeError(f"No suitable focal outbreak windows found in {infections_path}")

    ll_dates = _read_linelist_dates(linelist_path)

    # Collect every alias & strain across all selected windows so we stream painted_fasta ONCE
    needed_aliases: set[str] = set()
    strain_to_alias: dict[str, str] = {}
    for w in windows:
        for m in w["members"]:
            al = m["alias_pid"]
            needed_aliases.add(al)
            if m["strain"]:
                strain_to_alias[m["strain"]] = al

    print(f"[{state}] Streaming {painted_fasta} once for {len(needed_aliases):,} candidate cases...")
    seqs_by_alias: dict[str, str] = {}
    with recon._open_text(painted_fasta, "rt") as fh:
        cur_al: str | None = None
        buf: list[str] = []
        for line in fh:
            if line.startswith(">"):
                if cur_al is not None:
                    seqs_by_alias[cur_al] = "".join(buf).upper()
                hdr_id = line[1:].strip().split()[0]
                al = strain_to_alias.get(hdr_id) or recon.alias_from_any_id(hdr_id)
                cur_al = al if (al and al in needed_aliases) else None
                buf = []
            elif cur_al is not None:
                buf.append(line.strip())
        if cur_al is not None:
            seqs_by_alias[cur_al] = "".join(buf).upper()
    print(f"[{state}] Loaded {len(seqs_by_alias):,}/{len(needed_aliases):,} painted genomes.")

    # Read Wuhan reference for augur rooting
    ref_seq = ""
    with recon._open_text(reference_fasta, "rt") as fh:
        lines = [ln.strip() for ln in fh if not ln.startswith(">")]
        ref_seq = "".join(lines).upper()

    manifests: list[dict[str, Any]] = []
    rates = list(sampling_rates)

    for c_idx, win in enumerate(windows):
        members = [m for m in win["members"] if m["alias_pid"] in seqs_by_alias]
        win_seqs = {m["alias_pid"]: seqs_by_alias[m["alias_pid"]] for m in members}
        comp_of_win = {m["alias_pid"]: m["component_id"] for m in members}
        dom_comp = Counter(comp_of_win.values()).most_common(1)[0][0]
        by_al = {m["alias_pid"]: m for m in members}

        # Pre-assign deterministic collection dates to every case in the window
        # so nested pi=0.25 vs pi=0.50 share identical case dates!
        rng_dates = random.Random(seed + 100 * (c_idx + 1))
        coll_date_of: dict[str, _dt.date] = {}
        for m in members:
            al = m["alias_pid"]
            if al in ll_dates:
                coll_date_of[al] = ll_dates[al]
            else:
                lag = max(1, int(round(rng_dates.gammavariate(a_s, 1.0 / lambda_s))))
                coll_date_of[al] = _dt.date.fromordinal(m["exp_ord"] + lag)

        # Draw nested uniform random numbers u_i in [0, 1) so pi=0.25 is a subset of pi=0.50!
        rng_draw = random.Random(seed + 1000 * (c_idx + 1))
        u_draw = {m["alias_pid"]: rng_draw.random() for m in members}

        for pi_val in rates:
            pi_tag = f"pi{int(round(pi_val * 100)):02d}"
            cluster_id = f"{state}_c{c_idx + 1}_{win['county_fips']}_{pi_tag}"
            run_dir = out_root / cluster_id
            in_dir = run_dir / "input_data"
            aug_dir = run_dir / "augur_stage"
            in_dir.mkdir(parents=True, exist_ok=True)
            aug_dir.mkdir(parents=True, exist_ok=True)

            raw_sampled = [m["alias_pid"] for m in members if u_draw[m["alias_pid"]] < pi_val]
            # Guarantee at least 25 sampled tips for JUNIPER subtree constraint
            if len(raw_sampled) < 25:
                ordered_u = sorted(members, key=lambda m: u_draw[m["alias_pid"]])
                raw_sampled = [m["alias_pid"] for m in ordered_u[: min(28, len(members))]]

            isolated, compat_true = _isolate_focal_genetic_cluster(
                raw_sampled, win_seqs, dom_comp, comp_of_win, max_snps=max_cluster_snps
            )
            isolated.sort(key=lambda al: (coll_date_of[al].isoformat(), al))

            # Assign fixed-width safe IDs S000001.. for JUNIPER & outbreaker2
            width = max(6, len(str(len(isolated))))
            lookup_rows: list[dict[str, str]] = []
            with open(in_dir / "aligned.fasta", "w", encoding="utf-8") as fa_jun, \
                 open(aug_dir / "sequences.fasta", "w", encoding="utf-8") as fa_aug:
                # Inject Wuhan reference first into augur_stage/sequences.fasta ONLY
                fa_aug.write(">Wuhan/Hu-1/2019\n")
                for k in range(0, len(ref_seq), 80):
                    fa_aug.write(ref_seq[k : k + 80] + "\n")

                for idx_s, al in enumerate(isolated, start=1):
                    sid = f"S{idx_s:0{width}d}"
                    seq = win_seqs[al]
                    fa_jun.write(f">{sid}\n")
                    fa_aug.write(f">{sid}\n")
                    for k in range(0, len(seq), 80):
                        chunk = seq[k : k + 80] + "\n"
                        fa_jun.write(chunk)
                        fa_aug.write(chunk)

                    m_info = by_al[al]
                    strain_val = m_info["strain"] or f"USA/{state.upper()}-EHip-{al}/{coll_date_of[al].year}"
                    lookup_rows.append({
                        "sample": sid,
                        "date": coll_date_of[al].isoformat(),
                        "alias_pid": al,
                        "strain": strain_val,
                        "exposure_date": m_info["exposure_date"],
                        "county_fips": m_info["county_fips"],
                        "component_id": m_info["component_id"],
                    })

            # Write ONLY `metadata.csv` inside `input_data/` (JUNIPER globs `*metadata*`!)
            with open(in_dir / "metadata.csv", "w", newline="", encoding="utf-8") as fh:
                w_csv = csv.writer(fh)
                w_csv.writerow(["sample", "date"])
                for r in lookup_rows:
                    w_csv.writerow([r["sample"], r["date"]])

            # Write `sample_lookup.csv` in `run_dir/` (NOT inside `input_data/`)
            with open(run_dir / "sample_lookup.csv", "w", newline="", encoding="utf-8") as fh:
                w_dict = csv.DictWriter(fh, fieldnames=list(lookup_rows[0].keys()))
                w_dict.writeheader()
                w_dict.writerows(lookup_rows)

            # Write `augur_stage/metadata.tsv` (with Wuhan reference)
            with open(aug_dir / "metadata.tsv", "w", newline="", encoding="utf-8") as fh:
                w_tsv = csv.writer(fh, delimiter="\t")
                w_tsv.writerow(["strain", "date", "county_fips"])
                w_tsv.writerow(["Wuhan/Hu-1/2019", "2019-12-26", "00000"])
                for r in lookup_rows:
                    w_tsv.writerow([r["sample"], r["date"], r["county_fips"]])

            # Count local introductions inside the isolated sample set using parent pointers
            iso_set = set(isolated)
            compat_set = set(compat_true)
            parent_win = {m["alias_pid"]: m["parent_alias"] for m in members}
            n_direct_true = sum(1 for al in isolated if parent_win.get(al, "") in iso_set)
            n_local_roots_compat = sum(1 for al in compat_true if parent_win.get(al, "") not in compat_set)

            pi_realized = len(isolated) / float(max(1, len(compat_true)))
            manifest = {
                "cluster_id": cluster_id,
                "state": state,
                "cluster_index": c_idx + 1,
                "loc_key": win["loc_key"],
                "county_fips": win["county_fips"],
                "window_start": win["window_start"],
                "window_end": win["window_end"],
                "pi_target": pi_val,
                "pi_truth": round(pi_realized, 6),
                "n_true_window": len(members),
                "n_true_compatible_cluster": len(compat_true),
                "n_sampled_raw": len(raw_sampled),
                "n_sampled_isolated": len(isolated),
                "n_dropped_distant_snps": len(raw_sampled) - len(isolated),
                "n_true_direct_pairs_in_sample": n_direct_true,
                "n_compatible_local_intros": n_local_roots_compat,
                "n_distinct_statewide_comps": len({r["component_id"] for r in lookup_rows}),
                "mu_truth_yr": mu_truth_yr,
                "init_mu_day": round(mu_truth_yr / 365.25, 10),
                "gi_mean_days": gi_mean,
                "a_g": a_g,
                "lambda_g": lambda_g,
                "a_s": a_s,
                "lambda_s": lambda_s,
                "run_dir": str(run_dir),
            }
            (run_dir / "outbreak_manifest.json").write_text(
                json.dumps(manifest, indent=2), encoding="utf-8"
            )
            manifests.append(manifest)
            print(
                f"  -> {cluster_id}: N={len(isolated)} tips "
                f"(from {len(compat_true)} compatible / {len(members)} window cases, "
                f"pi_eff={pi_realized:.2f}, direct_pairs={n_direct_true})"
            )

    index_path = out_root / f"outbreak_index_{state}.json"
    index_path.write_text(json.dumps(manifests, indent=2), encoding="utf-8")
    return manifests


def run_nextstrain_on_outbreak(
    run_dir: str | Path,
    ncov_dir: str | Path | None = None,
    threads: int = 4,
) -> dict[str, Any]:
    """Run ``augur tree --method iqtree`` + ``augur refine --timetree`` on a
    prepared focal outbreak directory, and export ``dated_tree_years.nwk`` for
    ``TransPhylo``."""
    from .. import nextstrain_runtime as _runtime

    rdir = Path(run_dir).resolve()
    tt_dir = rdir / "treetime"
    tt_dir.mkdir(parents=True, exist_ok=True)

    # Stage inside ncov_dir using RELATIVE paths from ncov_dir so both host
    # executions (`cwd=ncov_dir`) and container executions (`nextstrain shell`,
    # which mounts `ncov_dir` at `/nextstrain/build`) see the exact same paths!
    stage_base = None
    work_dir = rdir / "augur_stage"
    ncov_path = Path(ncov_dir).resolve() if ncov_dir and Path(ncov_dir).exists() else None
    if ncov_path is not None:
        rel_dir = Path("results") / "_outbreak_stage" / rdir.name
        stage_base = ncov_path / rel_dir
        stage_base.mkdir(parents=True, exist_ok=True)
        shutil.copy2(work_dir / "sequences.fasta", stage_base / "sequences.fasta")
        shutil.copy2(work_dir / "metadata.tsv", stage_base / "metadata.tsv")
        run_cwd = str(ncov_path)
        p_seq = str(rel_dir / "sequences.fasta")
        p_meta = str(rel_dir / "metadata.tsv")
        p_raw = str(rel_dir / "tree_raw.nwk")
        p_ref = str(rel_dir / "tree.nwk")
        p_node = str(rel_dir / "branch_lengths.json")
        host_raw = stage_base / "tree_raw.nwk"
        host_ref = stage_base / "tree.nwk"
        host_node = stage_base / "branch_lengths.json"
    else:
        run_cwd = str(work_dir)
        p_seq = "sequences.fasta"
        p_meta = "metadata.tsv"
        p_raw = "tree_raw.nwk"
        p_ref = "tree.nwk"
        p_node = "branch_lengths.json"
        host_raw = work_dir / "tree_raw.nwk"
        host_ref = work_dir / "tree.nwk"
        host_node = work_dir / "branch_lengths.json"

    t0 = time.time()
    tree_args = [
        "tree",
        "--alignment", p_seq,
        "--method", "iqtree",
        "--nthreads", str(threads),
        "--output", p_raw,
    ]
    proc1, src1 = _runtime.run("augur", tree_args, ncov_dir=ncov_path, cwd=run_cwd)
    if proc1 is None or proc1.returncode != 0 or not host_raw.exists():
        err = ((proc1.stderr or proc1.stdout) if proc1 else src1) or "augur tree failed"
        raise RuntimeError(f"augur tree failed ({src1}): {err[-600:]}")

    refine_args = [
        "refine",
        "--tree", p_raw,
        "--alignment", p_seq,
        "--metadata", p_meta,
        "--output-tree", p_ref,
        "--output-node-data", p_node,
        "--timetree",
        "--root", "Wuhan/Hu-1/2019",
        "--coalescent", "opt",
        "--date-inference", "marginal",
        "--divergence-unit", "mutations",
        "--date-confidence",
        "--no-covariance",
    ]
    proc2, src2 = _runtime.run("augur", refine_args, ncov_dir=ncov_path, cwd=run_cwd)
    if proc2 is None or proc2.returncode != 0 or not host_node.exists():
        err = ((proc2.stderr or proc2.stdout) if proc2 else src2) or "augur refine failed"
        raise RuntimeError(f"augur refine failed ({src2}): {err[-600:]}")
    elapsed = time.time() - t0

    shutil.copy2(host_raw, tt_dir / "tree_raw.nwk")
    shutil.copy2(host_ref, tt_dir / "tree.nwk")
    shutil.copy2(host_node, tt_dir / "branch_lengths.json")
    if stage_base is not None:
        shutil.rmtree(stage_base, ignore_errors=True)

    meta = recon.export_dated_newick_for_transphylo(
        tt_dir / "tree.nwk",
        tt_dir / "branch_lengths.json",
        tt_dir / "dated_tree_years.nwk",
        tt_dir / "dated_tree_meta.json",
    )
    meta["runtime_secs"] = round(elapsed, 2)
    (tt_dir / "dated_tree_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return meta


def score_outbreak_directory(
    run_dir: str | Path,
    infections_path: str | Path,
    out_csv: str | Path | None = None,
) -> list[dict[str, Any]]:
    """Score every completed reconstruction method in ``run_dir`` against the
    ABM ground truth and write ``outbreak_scores.csv``."""
    rdir = Path(run_dir)
    manifest = json.loads((rdir / "outbreak_manifest.json").read_text(encoding="utf-8"))
    lookup_csv = rdir / "sample_lookup.csv"
    id_map, coll_dates = recon.load_sample_lookup(lookup_csv)
    tip_aliases = set(coll_dates.keys())

    truth = recon.load_truth_forest(infections_path, target_aliases=tip_aliases)

    base_info = {
        "cluster_id": manifest["cluster_id"],
        "state": manifest["state"],
        "county_fips": manifest["county_fips"],
        "window_start": manifest["window_start"],
        "window_end": manifest["window_end"],
        "pi_target": manifest["pi_target"],
        "pi_truth": manifest["pi_truth"],
        "n_sampled": manifest["n_sampled_isolated"],
        "n_true_compatible": manifest["n_true_compatible_cluster"],
        "mu_truth_yr": manifest["mu_truth_yr"],
    }

    rows: list[dict[str, Any]] = []

    # 1. MicrobeTrace SNP clustering + SNP+Date donor baseline
    fasta_in = rdir / "input_data" / "aligned.fasta"
    if fasta_in.exists():
        t0 = time.time()
        snp_clusters, snp_donors = recon.score_snp_network_and_heuristic(
            fasta_in, truth, coll_dates, lookup_map=id_map
        )
        t_snp = round(time.time() - t0, 3)
        # Pick tau=1 SNP as canonical MicrobeTrace outbreak threshold (plus save full tau table)
        (rdir / "microbetrace_thresholds.json").write_text(
            json.dumps(snp_clusters, indent=2), encoding="utf-8"
        )
        tau1 = next((c for c in snp_clusters if c["snp_threshold"] == 1), snp_clusters[0] if snp_clusters else {})
        mt_row: dict[str, Any] = dict(base_info)
        mt_row.update({
            "method": "MicrobeTrace_1SNP",
            "paradigm": "distance_cluster",
            "runtime_secs": t_snp,
            "tree_chain_ari": tau1.get("chain_ari"),
            "tree_chain_f1": tau1.get("chain_f1"),
            "local_intro_ari": tau1.get("local_intro_ari"),
            "local_intro_f1": tau1.get("local_intro_f1"),
        })
        mt_row.update(recon.score_transmission(snp_donors, truth, tip_aliases))
        rows.append(mt_row)

    # 2. Nextstrain (IQ-TREE 2 + TreeTime)
    tt_dir = rdir / "treetime"
    if (tt_dir / "tree.nwk").exists() and (tt_dir / "branch_lengths.json").exists():
        ns_tree = recon.parse_newick_file(
            tt_dir / "tree.nwk",
            branch_lengths_json=tt_dir / "branch_lengths.json",
            lookup_map=id_map,
        )
        tt_meta = {}
        if (tt_dir / "dated_tree_meta.json").exists():
            tt_meta = json.loads((tt_dir / "dated_tree_meta.json").read_text(encoding="utf-8"))
        mu_est = ns_tree.meta.get("clock_rate_yr")
        ns_row = dict(base_info)
        ns_row.update({
            "method": "Nextstrain_IQTREE_TreeTime",
            "paradigm": "phylogeny_timetree",
            "runtime_secs": tt_meta.get("runtime_secs"),
            "mu_est_yr": mu_est,
            "mu_rel_err_pct": round(100.0 * abs(mu_est - manifest["mu_truth_yr"]) / manifest["mu_truth_yr"], 3)
            if mu_est is not None and manifest["mu_truth_yr"] > 0 else None,
        })
        ns_row.update(recon.score_phylogeny(ns_tree, truth, tip_collection_dates=coll_dates))
        rows.append(ns_row)

    # 3. Bayesian Transmission Tools: outbreaker2, TransPhylo, JUNIPER
    for tool_key, label in (
        ("outbreaker2", "outbreaker2"),
        ("transphylo", "TransPhylo"),
        ("juniper", "JUNIPER"),
    ):
        probs_f = rdir / f"{tool_key}_donor_probs.csv"
        nodes_f = rdir / f"{tool_key}_nodes.csv"
        params_f = rdir / f"{tool_key}_params.json"
        if not probs_f.exists() and not nodes_f.exists():
            continue

        row = dict(base_info)
        row["method"] = label
        row["paradigm"] = "transmission_tree"

        if params_f.exists():
            pmeta = json.loads(params_f.read_text(encoding="utf-8"))
            row["runtime_secs"] = pmeta.get("runtime_secs")
            mu_y = pmeta.get("mu_est_yr")
            pi_e = pmeta.get("pi_est")
            row["mu_est_yr"] = mu_y
            if mu_y is not None and manifest["mu_truth_yr"] > 0:
                row["mu_rel_err_pct"] = round(
                    100.0 * abs(mu_y - manifest["mu_truth_yr"]) / manifest["mu_truth_yr"], 3
                )
            row["pi_est"] = pi_e
            if pi_e is not None:
                row["pi_abs_err"] = round(abs(pi_e - manifest["pi_truth"]), 4)
            row["R_est"] = pmeta.get("R_est")
            row["ess_min"] = pmeta.get("ess_min")

        if nodes_f.exists():
            h_tree, point_probs = recon.parse_host_network_csv(nodes_f, lookup_csv=lookup_csv)
            row.update(recon.score_phylogeny(h_tree, truth, min_mutations_collapse=0.5, tip_collection_dates=coll_dates))
        else:
            point_probs = None

        if probs_f.exists():
            post_probs = recon.parse_donor_probs_csv(probs_f, lookup_csv=lookup_csv)
            row.update(recon.score_transmission(post_probs, truth, tip_aliases))
        elif point_probs is not None:
            row.update(recon.score_transmission(point_probs, truth, tip_aliases))

        rows.append(row)

    target_csv = Path(out_csv) if out_csv else (rdir / "outbreak_scores.csv")
    if rows:
        fieldnames: list[str] = []
        for r in rows:
            for k in r.keys():
                if k not in fieldnames:
                    fieldnames.append(k)
        target_csv.parent.mkdir(parents=True, exist_ok=True)
        with open(target_csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)
    return rows


def score_state_surveillance_trees(
    state: str,
    infections_path: str | Path,
    builds_dir: str | Path,
    out_csv: str | Path,
    juniper_runs_dir: str | Path | None = None,
) -> list[dict[str, Any]]:
    """Score existing statewide ``04_nextstrain_builds/<arm>/auspice.json`` builds
    (plus any finished ``JUNIPER`` sub-study builds and matching pruned Nextstrain
    trees on the exact same ``nova3_jun_aug`` scope) against ground truth."""
    bdir = Path(builds_dir)
    auspice_files = sorted(bdir.glob("*/auspice.json"))
    parsed_builds: list[tuple[str, recon.PhyloTree]] = []
    all_tips: set[str] = set()

    for af in auspice_files:
        arm = af.parent.name
        tr = recon.parse_auspice_json(af)
        parsed_builds.append((arm, tr))
        all_tips.update(tr.tips_map().keys())

    # Check optional Bryan JUNIPER sub-study runs (e.g. nova3_jun_aug)
    jun_entries: list[tuple[str, Path, Path]] = []
    if juniper_runs_dir and Path(juniper_runs_dir).exists():
        for sub in sorted(Path(juniper_runs_dir).glob("*")):
            if not sub.is_dir():
                continue
            nodes_f = sub / "nodes.csv"
            lkp_f = sub / "input_data" / "sample_lookup.csv"
            if nodes_f.exists() and lkp_f.exists():
                id_m, _ = recon.load_sample_lookup(lkp_f)
                all_tips.update(id_m.values())
                jun_entries.append((sub.name, nodes_f, lkp_f))

    if not all_tips:
        print(f"[{state}] No Auspice or JUNIPER trees found under {builds_dir}")
        return []

    print(f"[{state}] Loading ABM truth forest for {len(all_tips):,} surveillance tips...")
    truth = recon.load_truth_forest(infections_path, target_aliases=all_tips)

    rows: list[dict[str, Any]] = []
    for arm, tr in parsed_builds:
        print(f"[{state}] Scoring statewide Nextstrain build: {arm}...")
        sc = recon.score_phylogeny(tr, truth)
        row: dict[str, Any] = {
            "state": state,
            "scope": "statewide_full",
            "arm": arm,
            "method": "Nextstrain_IQTREE_TreeTime",
        }
        row.update(sc)
        rows.append(row)

    for jname, nodes_f, lkp_f in jun_entries:
        print(f"[{state}] Scoring regional JUNIPER build: {jname}...")
        j_tree, j_probs = recon.parse_host_network_csv(nodes_f, lookup_csv=lkp_f)
        j_tips = set(j_tree.tips_map().keys())
        j_row: dict[str, Any] = {
            "state": state,
            "scope": "nova3_jun_aug",
            "arm": jname,
            "method": "JUNIPER",
        }
        j_row.update(recon.score_phylogeny(j_tree, truth))
        j_row.update(recon.score_transmission(j_probs, truth, j_tips))
        rows.append(j_row)

        # Also score the matching Nextstrain tree pruned to the exact same tips for parity
        match_arm = "1S-P__lasso_greedy" if "lasso" in jname else "1S__surs"
        for arm, tr in parsed_builds:
            if arm == match_arm:
                pr_tr = tr.prune_to_tips(j_tips)
                ns_sub_row: dict[str, Any] = {
                    "state": state,
                    "scope": "nova3_jun_aug_pruned",
                    "arm": jname,
                    "method": "Nextstrain_IQTREE_TreeTime_pruned",
                }
                ns_sub_row.update(recon.score_phylogeny(pr_tr, truth))
                rows.append(ns_sub_row)

    if rows:
        fieldnames: list[str] = []
        for r in rows:
            for k in r.keys():
                if k not in fieldnames:
                    fieldnames.append(k)
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)
        print(f"[{state}] Wrote {len(rows)} surveillance tree score rows -> {out_csv}")
    return rows
