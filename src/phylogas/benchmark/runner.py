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
def read_infections(path: str) -> pd.DataFrame:
    """Read the ABM all-events file (the transmission ground truth).

    Uses the C parser explicitly. BeyondBaseline read this with
    ``sep=None, engine="python"``, which forces the sniffer and the pure-Python
    parser -- about 2.8x slower on a 300k-row file.
    """
    dtypes = {
        "alias_pid": str, "alias_contact": str,
        "sim_pid": str, "pid": str, "contact_pid": str,
    }
    try:
        df = pd.read_csv(path, dtype=dtypes)
    except (ValueError, pd.errors.ParserError):
        # Fall back to sniffing only if the comma-delimited read fails.
        df = pd.read_csv(path, sep=None, engine="python", dtype=dtypes)
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
# phylogas benchmark sequence
# --------------------------------------------------------------------------
def benchmark_sequence(
    painted_fasta: str,
    infections: str,
    out_csv: str | None = None,
    max_pairs: int | None = 200_000,
) -> pd.DataFrame:
    """Measure realised parent->child divergence in the painted genomes.

    This is the check that surfaced the re-importation bug: every transmission
    edge should show a small number of substitutions. A heavy tail means
    children are being handed the wrong parent genome.
    """
    import gzip
    import lzma

    from Bio import SeqIO

    print(f"Infections : {infections}")
    inf = read_infections(infections)
    pid_col = _pid_col(inf)

    tick_col = next((c for c in ("sim_tick", "exposure_tick", "tick") if c in inf.columns), None)
    if tick_col is None:
        raise SystemExit("ERROR: no tick column (sim_tick/exposure_tick/tick) in the infections file.")

    inf = inf.sort_values(tick_col, kind="stable")
    active: dict[str, str] = {}
    pairs: list[tuple[str, str]] = []
    for pid, contact, tick in zip(inf[pid_col].astype(str),
                                  inf["alias_contact"].astype(str),
                                  inf[tick_col]):
        iid = f"{pid}.{tick}"
        parent = active.get(contact)
        if contact not in ("-1", "nan", "") and parent:
            pairs.append((iid, parent))
        active[pid] = iid
    print(f"  {len(pairs):,} parent-child pairs")

    opener = lzma.open if painted_fasta.endswith(".xz") else \
        gzip.open if painted_fasta.endswith(".gz") else open
    seqs: dict[str, str] = {}
    with opener(painted_fasta, "rt") as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            if "EHip-" in rec.id:
                seqs[rec.id.split("EHip-")[1].rsplit("/", 1)[0]] = str(rec.seq)
    print(f"  {len(seqs):,} painted sequences")

    if max_pairs and len(pairs) > max_pairs:
        rng = np.random.default_rng(0)
        idx = rng.choice(len(pairs), max_pairs, replace=False)
        pairs = [pairs[i] for i in idx]
        print(f"  sampling {max_pairs:,} pairs")

    dists = [
        sum(1 for a, b in zip(seqs[c], seqs[p]) if a != b)
        for c, p in pairs if c in seqs and p in seqs
    ]
    if not dists:
        raise SystemExit("ERROR: no parent-child pair had both genomes present.")

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
