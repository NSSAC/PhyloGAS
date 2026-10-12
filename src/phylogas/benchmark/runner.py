"""Benchmark drivers: score sampling strategies against the simulation's truth."""

from __future__ import annotations

import glob as _glob
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from . import mugration
from .scoring import score_matrices, score_offdiagonals


# --------------------------------------------------------------------------
# shared IO
# --------------------------------------------------------------------------
def read_infections(path: str, usecols: list[str] | None = None) -> pd.DataFrame:
    """Read the ABM all-events file (the transmission ground truth).

    Uses the C parser explicitly. BeyondBaseline read this with
    ``sep=None, engine="python"``, which forces the sniffer and the pure-Python
    parser -- about 2.8x slower on a 300k-row file.

    `usecols` is a *wish list*, not a requirement: it is intersected with the
    header, so a caller may name a column that only some vintages of the table
    carry without having to probe for it first. Worth using -- the full table
    is 46 mostly-string columns and about 8.6 GB in memory at 5.3M rows, while
    the three columns a transmission edge actually needs are under 1 GB. That
    difference is what decides whether several states fit on one node.
    """
    dtypes = {
        "alias_pid": str, "alias_contact": str,
        "sim_pid": str, "pid": str, "contact_pid": str,
    }
    cols = None
    if usecols:
        header = pd.read_csv(path, nrows=0)
        have = {c.strip() for c in header.columns}
        cols = [c for c in usecols if c in have]
        missing = [c for c in usecols if c not in have]
        if not cols:
            raise SystemExit(
                f"ERROR: none of the requested columns are in {path}\n"
                f"  wanted: {', '.join(usecols)}\n"
                f"  have:   {', '.join(sorted(have))}")
        if missing:
            print(f"  (not in this table, skipped: {', '.join(missing)})")
    try:
        df = pd.read_csv(path, dtype=dtypes, usecols=cols)
    except (ValueError, pd.errors.ParserError):
        # Fall back to sniffing only if the comma-delimited read fails.
        df = pd.read_csv(path, sep=None, engine="python", dtype=dtypes, usecols=cols)
    df.columns = [c.strip() for c in df.columns]
    return df


def _pid_col(df: pd.DataFrame) -> str:
    return "alias_pid" if "alias_pid" in df.columns else "sim_pid"


def _label_from_path(p: Path) -> str:
    """Recover an algorithm label from a BeyondBaseline samples filename.

    e.g. ``run03_scenario1_LASSO-Greedy_samples.csv.xz`` -> ``scenario1/LASSO-Greedy``
    """
    name = p.name
    for suffix in (".csv.xz", ".csv.gz", ".csv"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    name = re.sub(r"_samples$", "", name)
    m = re.search(r"(scenario\d+)_(.+)$", name)
    return f"{m.group(1)}/{m.group(2)}" if m else name


# --------------------------------------------------------------------------
# phylogas benchmark mugration
# --------------------------------------------------------------------------
def benchmark_mugration(
    truth_json: str,
    samples: list[str],
    infections: str,
    out_csv: str | None = None,
    state_col: str = "county",
    duration_years: float = 1.0,
    save_matrices: str | None = None,
) -> pd.DataFrame:
    """Score one or more sampled sets against the ABM mugration truth.

    The transmission graph is read and built **once** and reused across every
    sample set, which is what keeps batch scoring cheap: graph construction is
    O(V+E) and dominates a single comparison.
    """
    print(f"Truth matrix : {truth_json}")
    master_alphabet, truth_matrix = mugration.read_traits_json(truth_json)
    print(f"  alphabet: {len(master_alphabet)} states")

    print(f"Infections   : {infections}")
    inf_df = read_infections(infections)
    print(f"  {len(inf_df):,} rows")

    pid_col = _pid_col(inf_df)
    graph = mugration.build_directed_graph(inf_df, pid_col=pid_col, contact_col="alias_contact")
    print(f"  transmission graph: {graph.number_of_nodes():,} nodes, "
          f"{graph.number_of_edges():,} edges")

    rows = []
    for spec in samples:
        for path in sorted(_glob.glob(spec)) or [spec]:
            p = Path(path)
            if not p.exists():
                print(f"  WARNING: no such samples file: {p}", file=sys.stderr)
                continue
            label = _label_from_path(p)
            sdf = pd.read_csv(p, dtype={"alias_pid": str, "sim_pid": str, "pid": str})
            s_pid = _pid_col(sdf)
            if state_col not in sdf.columns:
                print(f"  WARNING: '{state_col}' missing in {p.name}; skipping", file=sys.stderr)
                continue

            raw = sdf.set_index(s_pid)[state_col].to_dict()
            tips = {str(k).replace(".0", ""): str(v) for k, v in raw.items() if pd.notna(v)}

            res = mugration.simulate_inference_and_matrix(
                graph, tips, master_alphabet, duration_years
            )
            inferred = res["models"][state_col]["transition_matrix"] \
                if state_col in res.get("models", {}) else res["models"]["county"]["transition_matrix"]
            inferred_alpha = res["models"].get(state_col, res["models"]["county"])["alphabet"]

            scores = score_matrices(truth_matrix, master_alphabet, inferred, inferred_alpha)
            scores.update(label=label, samples_file=str(p), n_samples=len(sdf), n_tips=len(tips))
            rows.append(scores)
            print(f"  {label:34s} cos={scores['cosine_similarity']:.4f} "
                  f"F1={scores['topological_f1']:.4f} r={scores['pearson_r']:.4f}")

            if save_matrices:
                Path(save_matrices).mkdir(parents=True, exist_ok=True)
                with open(Path(save_matrices) / f"inferred_{label.replace('/', '_')}.json", "w") as fh:
                    json.dump(res, fh, indent=2)

    if not rows:
        raise SystemExit("ERROR: no sample sets were scored.")

    cols = ["label", "cosine_similarity", "topological_f1", "pearson_r", "masked_mae",
            "n_samples", "n_tips", "n_edges_truth", "n_edges_inferred", "samples_file"]
    df = pd.DataFrame(rows)[cols].sort_values("cosine_similarity", ascending=False)
    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)
        print(f"\nWrote {len(df)} rows -> {out_csv}")
    return df


# --------------------------------------------------------------------------
# realised generation interval
# --------------------------------------------------------------------------
# Candidate tick columns, best first. `exposure_tick` is the one that means
# "when this infection was infected"; `sim_tick` is the tick of the row's
# state transition, which is later by the incubation period and varies between
# rows for the same infection. Using sim_tick would measure neither the
# generation nor the serial interval, so it is only a fallback for an older
# table that lacks exposure_tick.
_GI_TICK_COLUMNS = ("exposure_tick", "sim_tick", "tick")


def _gi_tick_column(inf: pd.DataFrame) -> str | None:
    return next((c for c in _GI_TICK_COLUMNS if c in inf.columns), None)


def check_alias_convention(inf: pd.DataFrame, tick_col: str, n: int = 50_000) -> bool:
    """Confirm that `alias_pid` really is "{pid}.{exposure_tick}".

    Nothing in the generation-interval code below depends on that format -- it
    reads ticks from their own column. But `benchmark_sequence` does depend on
    it, because the painted FASTA names each genome with the alias and the only
    way to match a genome to an infection is by that string. So rather than
    leave the convention an unstated coupling between the painter, TwinSampler
    and every consumer, check it here where the tick is available both ways and
    say so when it breaks. A warning, not an error: a mismatch makes the genome
    matching suspect, not the intervals.
    """
    sub = inf if len(inf) <= n else inf.sample(n, random_state=0)
    suffix = pd.to_numeric(sub["alias_pid"].astype(str).str.rsplit(".", n=1).str[-1],
                           errors="coerce")
    own = pd.to_numeric(sub[tick_col], errors="coerce")
    bad = int((suffix != own).sum())
    if bad:
        print(f"  WARNING: alias_pid suffix disagrees with {tick_col} on "
              f"{bad:,}/{len(sub):,} checked rows. The intervals below are "
              f"unaffected, but genome matching by alias is no longer reliable.",
              file=sys.stderr)
    return bad == 0


def generation_interval(
    inf: pd.DataFrame,
    tick_col: str | None = None,
    days_per_tick: float = 1.0,
    plausible_window: int | None = 60,
) -> dict:
    """Realised infector->infectee interval, read from columns not from ids.

    This is the **generation** interval: exposure of the infector to exposure
    of the infectee. It is not the serial interval, which is onset to onset and
    would need the symptom columns. And it is the *realised* interval of this
    run, not a configured one -- EpiHiper's disease progression produces it, so
    it is a measurement of the ABM rather than a parameter of ours.

    It matters because it is the denominator of the molecular clock:
    ``subs/genome/day = mutations-per-transmission x saturation / GI``. Two of
    those three come from the painter; this one is set entirely upstream, so
    assuming it rather than measuring it silently misattributes the whole
    difference to the mutational model.

    Each row carries its own infection's tick, so the infectee's exposure needs
    no lookup at all. The infector's does, and the lookup is a join of
    ``alias_contact`` against a ``{alias_pid -> tick}`` table -- never a parse
    of the tick back out of the id string. `plausible_window` then bounds how
    far apart a joined pair may be: anything outside it is reported rather than
    dropped, since a long tail here would be a real finding and a *negative*
    interval would mean the join had matched the wrong row.
    """
    tick_col = tick_col or _gi_tick_column(inf)
    if tick_col is None or "alias_contact" not in inf.columns:
        print(f"  generation interval: skipped (no tick column; looked for "
              f"{', '.join(_GI_TICK_COLUMNS)})", file=sys.stderr)
        return {}
    if tick_col != "exposure_tick":
        print(f"  WARNING: no exposure_tick column; falling back to '{tick_col}', "
              f"which is the state-transition tick and overstates the interval.",
              file=sys.stderr)

    check_alias_convention(inf, tick_col)

    # One tick per infection. The all-events table is already one row per
    # alias_pid, but an older or filtered table can carry a row per state
    # transition, and exposure_tick is constant within an infection either way.
    own = (inf.assign(_t=pd.to_numeric(inf[tick_col], errors="coerce"))
              .dropna(subset=["_t"])
              .drop_duplicates("alias_pid")
              .set_index("alias_pid")["_t"])

    parent = inf["alias_contact"].astype(str)
    edges = ~parent.isin(["-1", "nan", "", "None"]) & inf["alias_contact"].notna()
    child_t = pd.to_numeric(inf.loc[edges, tick_col], errors="coerce")
    parent_t = parent[edges].map(own)

    gi = (child_t - parent_t).dropna()
    n_edges, n_joined = int(edges.sum()), len(gi)
    if not n_joined:
        print("  generation interval: no edge joined an infector row", file=sys.stderr)
        return {"gi_edges": n_edges, "gi_pairs_joined": 0}

    g = gi.to_numpy(dtype=float) * days_per_tick
    out = {
        "gi_edges": n_edges,
        "gi_pairs_joined": n_joined,
        "gi_parent_unmatched": n_edges - n_joined,
        "gi_mean_days": float(g.mean()),
        "gi_median_days": float(np.median(g)),
        "gi_q25_days": float(np.percentile(g, 25)),
        "gi_q75_days": float(np.percentile(g, 75)),
        "gi_min_days": float(g.min()),
        "gi_max_days": float(g.max()),
        "gi_nonpositive": int((g <= 0).sum()),
    }
    if plausible_window:
        out["gi_outside_window"] = int((g > plausible_window * days_per_tick).sum())

    print("\nRealised generation interval (infector exposure -> infectee exposure):")
    print(f"  pairs joined             {n_joined:,} of {n_edges:,} edges"
          f" ({out['gi_parent_unmatched']:,} infector row absent)")
    print(f"  mean / median            {out['gi_mean_days']:.2f} / "
          f"{out['gi_median_days']:.1f} days")
    print(f"  IQR                      {out['gi_q25_days']:.0f}-{out['gi_q75_days']:.0f}"
          f"  (range {out['gi_min_days']:.0f}-{out['gi_max_days']:.0f})")
    if out["gi_nonpositive"]:
        print(f"  WARNING: {out['gi_nonpositive']:,} pairs have an infectee exposed no "
              f"later than its infector, which means the join matched a wrong row.")
    if out.get("gi_outside_window"):
        print(f"  NOTE: {out['gi_outside_window']:,} pairs exceed the {plausible_window}-tick "
              f"plausibility window; kept, not dropped.")
    return out


# --------------------------------------------------------------------------
# phylogas benchmark sequence
# --------------------------------------------------------------------------
def benchmark_sequence(
    painted_fasta: str,
    infections: str,
    out_csv: str | None = None,
    max_pairs: int | None = 200_000,
    days_per_tick: float = 1.0,
    gi_window: int | None = 60,
) -> pd.DataFrame:
    """Measure realised parent->child divergence in the painted genomes.

    This is the check that surfaced the re-importation bug: every transmission
    edge should show a small number of substitutions. A heavy tail means
    children are being handed the wrong parent genome.

    It also reports the realised generation interval over the same edges,
    because the two together are the molecular clock: substitutions per
    transmission divided by days per transmission is subs/genome/day. Measuring
    them from one pass over one file is what keeps the numerator and the
    denominator describing the same run.
    """
    import gzip
    import lzma

    from Bio import SeqIO

    print(f"Infections : {infections}")
    # Only the edge and its tick: this rule never touches the demographic
    # columns, and loading all 46 of them costs ~8.6 GB at 5.3M rows against
    # ~0.8 GB for these three. The genome dict below is the real floor
    # (~11.6 GB for the pairs of one state), so trimming here is what lets
    # several states share a node.
    inf = read_infections(infections, usecols=["alias_pid", "alias_contact",
                                               *_GI_TICK_COLUMNS])
    if "alias_pid" not in inf.columns or "alias_contact" not in inf.columns:
        raise SystemExit("ERROR: the infections file needs alias_pid and alias_contact.")

    # Measured before the pair subsampling below, so the interval is over every
    # edge while the divergence is over a sample of them.
    gi = generation_interval(inf, days_per_tick=days_per_tick, plausible_window=gi_window)

    # Each row is one infection: alias_pid is its id and alias_contact its
    # infector's infection id -- the same "{pid}.{tick}" key the painted FASTA
    # uses -- so a transmission edge is simply that pair. (This used to rebuild
    # the edges from person ids and append the tick again, giving
    # "123.45.130"-style keys that matched no genome.)
    child = inf["alias_pid"].astype(str)
    parent = inf["alias_contact"].astype(str)
    keep = ~parent.isin(["-1", "nan", "", "None"])
    pairs = list(zip(child[keep], parent[keep]))
    print(f"  {len(pairs):,} parent-child pairs")

    if max_pairs and len(pairs) > max_pairs:
        rng = np.random.default_rng(0)
        idx = rng.choice(len(pairs), max_pairs, replace=False)
        pairs = [pairs[i] for i in idx]
        print(f"  sampling {max_pairs:,} pairs")

    # Stream the FASTA and keep only the genomes those pairs need. Holding
    # every painted genome needs ~30 kB each -- about 150 GB for a 5M-record
    # state -- which is what got these jobs killed.
    needed = {c for c, _ in pairs} | {p for _, p in pairs}
    opener = lzma.open if painted_fasta.endswith(".xz") else \
        gzip.open if painted_fasta.endswith(".gz") else open
    seqs: dict[str, str] = {}
    n_read = 0
    with opener(painted_fasta, "rt") as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            n_read += 1
            if "EHip-" not in rec.id:
                continue
            key = rec.id.split("EHip-")[1].rsplit("/", 1)[0]
            if key in needed:
                seqs[key] = str(rec.seq)
    print(f"  {n_read:,} painted sequences read, {len(seqs):,} needed for the pairs")

    dists = [
        sum(1 for a, b in zip(seqs[c], seqs[p]) if a != b)
        for c, p in pairs if c in seqs and p in seqs
    ]
    if not dists:
        sample_ids = list(needed)[:3]
        raise SystemExit(
            "ERROR: no parent-child pair had both genomes present.\n"
            f"  e.g. infection ids sought: {sample_ids}\n"
            "  Check that the painted FASTA and the infections file come from the "
            "same run.")
    print(f"  {len(dists):,} pairs scored ({len(pairs) - len(dists):,} lacked a genome, "
          f"e.g. a parent infected before the painted window)")

    d = np.asarray(dists)
    summary = {
        "pairs_scored": len(d),
        "mean_substitutions": float(d.mean()),
        "median_substitutions": float(np.median(d)),
        "max_substitutions": int(d.max()),
        "frac_unchanged": float((d == 0).mean()),
        "pairs_over_8": int((d > 8).sum()),
    }
    print("\nParent -> child divergence:")
    for k, v in summary.items():
        print(f"  {k:24s} {v}")
    if summary["pairs_over_8"]:
        print("  WARNING: pairs with >8 substitutions suggest mis-assigned parent genomes.")

    # Added after the print loop above, which would otherwise repeat the
    # interval lines generation_interval has already written.
    if gi.get("gi_mean_days"):
        # The emergent clock, assembled from the two halves rather than
        # asserted: this is what a tree-free root-to-tip regression should find.
        per_day = summary["mean_substitutions"] / gi["gi_mean_days"]
        summary["implied_subs_per_genome_per_day"] = per_day
        summary["implied_subs_per_genome_per_year"] = per_day * 365.0
        print(f"\nImplied clock = {summary['mean_substitutions']:.3f} subs/transmission "
              f"/ {gi['gi_mean_days']:.2f} d = {per_day * 365.0:.2f} subs/genome/yr")
    summary.update(gi)

    df = pd.DataFrame([summary])
    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)
        print(f"\nWrote -> {out_csv}")
    return df


# --------------------------------------------------------------------------
# phylogas benchmark compare
# --------------------------------------------------------------------------
def benchmark_compare(
    truth_json: str,
    simulated_json: str,
    augur_json: str,
    out_csv: str | None = None,
    state_col: str = "county",
) -> pd.DataFrame:
    """Compare simulated parsimony against real augur inference, same truth.

    The gap between the two rows is the cost of doing actual phylogenetic
    reconstruction from sequences: tree error, plus contamination from lineages
    that were never sampled.
    """
    master, truth = mugration.read_traits_json(truth_json)
    sim_alpha, sim_m = mugration.read_traits_json(simulated_json)
    aug_alpha, aug_m = mugration.read_traits_json(augur_json)

    rows = []
    for label, alpha, matrix in (("simulated_parsimony", sim_alpha, sim_m),
                                 ("augur_nextstrain", aug_alpha, aug_m)):
        sc = score_matrices(truth, master, matrix, alpha)
        sc["method"] = label
        rows.append(sc)
        print(f"  {label:22s} cos={sc['cosine_similarity']:.4f} "
              f"F1={sc['topological_f1']:.4f} r={sc['pearson_r']:.4f}")

    df = pd.DataFrame(rows)[
        ["method", "cosine_similarity", "topological_f1", "pearson_r", "masked_mae",
         "n_edges_truth", "n_edges_inferred"]
    ]
    delta = df.iloc[1]["cosine_similarity"] - df.iloc[0]["cosine_similarity"]
    print(f"\n  cosine delta (augur - simulated): {delta:+.4f}")
    print("  negative => sequence-based reconstruction loses signal relative to")
    print("              perfect-knowledge parsimony on the same sampled tips.")

    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)
        print(f"\nWrote -> {out_csv}")
    return df
