import argparse
import sys
import pandas as pd
import numpy as np
import requests
import os
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional


try:
    from pango_aliasor.aliasor import Aliasor

    def make_variant_base_map(base_variants: List[str], recombinant: bool = False) -> Dict[str, str]:
        namer = Aliasor()
        namer.enable_expansion()
        all_rules = namer.partition_focus(base_variants, recombinant=recombinant)
        lineage_base_map = {k: v for v, ks in all_rules.items() for k in ks}
        for b in base_variants:
            if b not in lineage_base_map:
                lineage_base_map[b] = b
        return lineage_base_map

except ImportError:
    print("CRITICAL ERROR: The 'pango_aliasor' library is not installed.")
    print("Please install it, e.g., using 'pip install git+ssh://git@github.com:aswarren/pango_aliasor.git'.")
    def make_variant_base_map(base_variants: List[str], recombinant: bool = False) -> Dict[str, str]:
        raise ImportError("pango_aliasor not found, function unusable.")
    PANGO_ALIASOR_AVAILABLE = False
else:
    PANGO_ALIASOR_AVAILABLE = True

try:
    from . import cluster_seeds
except ImportError:                 # run as a plain script
    import cluster_seeds

COVSPECTRUM_API_URL = 'https://lapis.cov-spectrum.org/open/v2/sample/alignedNucleotideSequences'
DEFAULT_TSV_URL = 'https://clustertracker.gi.ucsc.edu/data/hardcoded_clusters.tsv'
DEFAULT_TSV_BASENAME = 'hardcoded_clusters.tsv'

def download_file(url: str, dest_path: Path, insecure: bool = False) -> bool:
    """Download a file, optionally skipping TLS certificate verification.

    ``insecure=True`` exists because the clustertracker.gi.ucsc.edu certificate
    expired on 2025-07-02, which otherwise blocks --seed_mode entirely. It
    disables verification for this request only.
    """
    print(f"Downloading {url} to {dest_path}...")
    if insecure:
        print("  WARNING: --insecure_download set; TLS certificate verification is OFF")
        print("           for this request. The transfer is still encrypted, but the")
        print("           server's identity is NOT authenticated. Only use this with a")
        print("           host you trust, and prefer --input_file where you can.")
        try:
            import urllib3
            urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
        except Exception:
            pass
    try:
        response = requests.get(url, stream=True, verify=not insecure)
        response.raise_for_status()
        with open(dest_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print("Download complete.")
        return True
    except requests.exceptions.SSLError as e:
        print(f"Error downloading file: SSL verification failed for {url}")
        print(f"  {e}")
        print("\n  NOTE: as of 2025-07-02 the clustertracker.gi.ucsc.edu certificate is")
        print("  expired. This is an upstream problem, not a PhyloGAS one. Options:")
        print("    1. Retry with --insecure_download to bypass certificate verification.")
        print("    2. Supply the cluster TSV yourself:  --input_file <hardcoded_clusters.tsv>")
        print("       (a copy ships in data/importations/)")
        print("    3. Use bulk mode instead of --seed_mode, which queries Cov-Spectrum")
        print("       directly and is unaffected.")
        if dest_path.exists(): dest_path.unlink()
        return False
    except requests.exceptions.RequestException as e:
        print(f"Error downloading file: {e}")
        if dest_path.exists(): dest_path.unlink()
        return False

def fetch_sequences_by_strain_id(strain_ids: List[str], batch_size: int = 100) -> Optional[str]:
    """Fetches sequences for a specific list of strain IDs (for seed_mode)."""
    all_fasta_content = []
    print(f"Fetching sequences for {len(strain_ids)} strains from CovSpectrum (batch size: {batch_size})...")
    for i in range(0, len(strain_ids), batch_size):
        batch_ids = strain_ids[i:i + batch_size]
        payload = {"strain": batch_ids}
        print(f"  Fetching batch {i//batch_size + 1}/{(len(strain_ids) - 1)//batch_size + 1} ({len(batch_ids)} strains)...")
        try:
            response = requests.post(COVSPECTRUM_API_URL, json=payload, timeout=120)
            response.raise_for_status()
            fasta_data = response.text
            if not fasta_data.strip() and len(batch_ids) > 0:
                 print(f"    Warning: Batch {i//batch_size + 1} returned empty data from API.")
            all_fasta_content.append(fasta_data)
        except requests.exceptions.RequestException as e:
            print(f"  Request Error for batch {i//batch_size + 1}: {e}")
            return None
    print("Sequence fetching complete.")
    return "".join(all_fasta_content)


def fetch_sequences_by_metadata(pango: str, state: str, date_from: Optional[str], date_to: Optional[str], output_filepath: Path, include_sublineages: bool = True):
    """Fetches sequences directly from CovSpectrum based on metadata query and streams to a file.

    ``include_sublineages`` controls the trailing ``*`` in the LAPIS query.
    With it, ``B.1.617.2`` matches the whole Delta clade; without it, only
    sequences assigned to exactly that lineage. The difference is large --
    for Virginia Delta, 21,034 samples versus 836.
    """
    date_info = "all dates"
    if date_from and date_to:
        date_info = f"from {date_from} to {date_to}"

    query_lineage = f'{pango}*' if include_sublineages else pango
    scope = "including sublineages" if include_sublineages else "exact lineage only"
    print(f"\n--- Starting bulk download for {query_lineage} in {state} ({date_info}; {scope}) ---")

    params = {
        'country': 'USA',
        'division': state,
        'pangoLineage': query_lineage,
        'downloadAsFile': 'true'
    }
    # Conditionally add date parameters to the request
    if date_from:
        params['dateFrom'] = date_from
    if date_to:
        params['dateTo'] = date_to
    
    headers = {'Accept': 'text/x-fasta'}

    try:
        with requests.get(COVSPECTRUM_API_URL, params=params, headers=headers, stream=True, timeout=300) as r:
            r.raise_for_status()
            with open(output_filepath, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        print(f"Bulk download complete. Sequences saved to: {output_filepath.resolve()}")
    except requests.exceptions.RequestException as e:
        print(f"Error during bulk download for {pango}: {e}")
        # Clean up partial download
        if output_filepath.exists():
            output_filepath.unlink()

def detect_outliers_iqr(date_series: pd.Series, factor: float = 1.5) -> pd.Series:
    if date_series.empty: return pd.Series(dtype=bool, index=date_series.index)
    Q1, Q3 = date_series.quantile(0.25), date_series.quantile(0.75)
    IQR = Q3 - Q1
    if IQR == pd.Timedelta(0): return pd.Series(False, index=date_series.index, dtype=bool)
    lower_bound, upper_bound = Q1 - factor * IQR, Q3 + factor * IQR
    return (date_series < lower_bound) | (date_series > upper_bound)

def detect_outliers_zscore(date_series: pd.Series, threshold: float = 2.0) -> pd.Series:
    from scipy.stats import zscore   # only this method needs scipy
    if date_series.empty or date_series.nunique() < 2: return pd.Series(False, index=date_series.index, dtype=bool)
    numeric_dates = (date_series - date_series.min()).dt.days
    if numeric_dates.nunique() < 2: return pd.Series(False, index=date_series.index, dtype=bool)
    z_scores = zscore(numeric_dates)
    z_scores = np.nan_to_num(z_scores, nan=0.0)
    return pd.Series(np.abs(z_scores) > threshold, index=date_series.index)

def detect_outliers_chaining(date_series: pd.Series, max_gap_weeks: int = 6) -> pd.Series:
    if date_series.empty: return pd.Series(dtype=bool, index=date_series.index)
    df_proc = date_series.to_frame(name='date').copy()
    df_sorted = df_proc.sort_values(by='date').copy()
    if df_sorted['date'].nunique() <= 1: return pd.Series(False, index=date_series.index, dtype=bool)
    max_delta = pd.Timedelta(weeks=max_gap_weeks)
    df_sorted['diff_prev'] = df_sorted['date'].diff()
    df_sorted['diff_next'] = df_sorted['date'].diff(-1).abs()
    is_far_from_prev = (df_sorted['diff_prev'] > max_delta) | df_sorted['diff_prev'].isna()
    is_far_from_next = (df_sorted['diff_next'] > max_delta) | df_sorted['diff_next'].isna()
    df_sorted['is_outlier'] = is_far_from_prev & is_far_from_next
    outlier_series_sorted_index = pd.Series(df_sorted['is_outlier'].values, index=df_sorted.index)
    return outlier_series_sorted_index.reindex(date_series.index)

def _sanitize(pango: str) -> str:
    return str(pango).replace('.', '_').replace('/', '_')


def _parse_variant_starts(spec: Optional[str]) -> Dict[str, str]:
    """Parse `B.1.1.7=2020-12-01,B.1.617.2=2021-03-01` into a dict."""
    if not spec:
        return dict(cluster_seeds.VARIANT_START_DEFAULTS)
    out = {}
    for item in spec.split(','):
        item = item.strip()
        if not item:
            continue
        if '=' not in item:
            print(f"Error: --variant_start entry '{item}' is not VARIANT=YYYY-MM-DD."); exit(1)
        k, v = item.split('=', 1)
        out[k.strip()] = v.strip()
    return out


def _drop_outlier_clusters(clusters: pd.DataFrame, args) -> pd.DataFrame:
    """Filter clusters whose importation date is an outlier for their variant.

    NOTE: the schedules EpiHiper was run with were produced with NO outlier
    filtering (`--outlier_method none`), so keep it off to reproduce them.
    Filtering here removes the cluster from the schedule as well as from the
    seed set, which keeps the two consistent but changes the simulation's
    importation volume.
    """
    dates = clusters["intro_date"]
    if args.outlier_method == 'iqr':
        mask = detect_outliers_iqr(dates)
    elif args.outlier_method == 'zscore':
        mask = detect_outliers_zscore(dates)
    elif args.outlier_method == 'chaining':
        mask = detect_outliers_chaining(dates, args.chaining_max_gap_weeks)
    else:
        return clusters
    mask = mask.fillna(False)
    if not mask.any():
        print(f"  Outlier filter ({args.outlier_method}): none found")
        return clusters

    def last_dated(samples):
        dated = [d for d, _s, _a in samples if not pd.isna(d)]
        return max(dated) if dated else pd.NaT

    cand = clusters[mask].copy()
    cand["span"] = cand["samples_ordered"].map(last_dated) - cand["intro_date"]
    size = cand["sample_count"] if "sample_count" in cand.columns else pd.Series(1, index=cand.index)
    rescued = (size >= args.rescue_cluster_size) & (
        cand["span"] <= pd.Timedelta(days=args.rescue_cluster_days))
    drop_idx = cand.index[~rescued.fillna(False)]
    print(f"  Outlier filter ({args.outlier_method}): {int(mask.sum()):,} flagged, "
          f"{int(rescued.sum()):,} rescued, {len(drop_idx):,} dropped")
    return clusters.drop(index=drop_idx)


def _tick_zero(args):
    """ABM tick 0 from --tick_zero or --abm_config, or None."""
    if getattr(args, "tick_zero", None):
        return pd.to_datetime(args.tick_zero)
    if getattr(args, "abm_config", None):
        try:
            return cluster_seeds.tick_zero_from_abm_config(args.abm_config)
        except Exception as exc:
            print(f"  WARNING: could not read tickZero from {args.abm_config}: {exc}")
    return None


def _write_schedules(clusters: pd.DataFrame, output_folder_path: Path, args) -> None:
    """Write the importation schedule, and the ABM's absolute-tick form.

    The schedule is written from the clusters that actually have a sequence, so
    the day-by-day importation counts and the seed FASTA always describe the
    same set of importations.
    """
    schedule = cluster_seeds.build_schedule(clusters)
    name = f"{args.state.replace(' ', '_')}_schedule.csv"
    path = output_folder_path / name
    out = schedule.copy()
    out["date"] = pd.to_datetime(out["date"]).dt.strftime("%Y-%m-%d")
    out.to_csv(path, index=False)
    print(f"\nImportation schedule ({len(out):,} days, "
          f"{int(out['importations'].sum()):,} importations) -> {path.resolve()}")

    tick_zero = _tick_zero(args)
    if tick_zero is None:
        return
    print(f"  Absolute ticks against tick 0 = {pd.to_datetime(tick_zero).date()}:")
    for variant, frame in cluster_seeds.absolute_ticks(schedule, tick_zero).items():
        fname = (f"{args.state.replace(' ', '_')}_{_sanitize(variant)}_ticks.csv")
        frame.to_csv(output_folder_path / fname, index=False)
        print(f"    {variant:12s} ticks {frame['tick'].min()}-{frame['tick'].max()}"
              f"  -> {fname}")


def _parse_fasta_records(text: str) -> Dict[str, str]:
    """Map strain id -> FASTA record, from concatenated API responses."""
    out: Dict[str, str] = {}
    for chunk in text.split(">"):
        if not chunk.strip():
            continue
        header = chunk.split("\n", 1)[0].strip()
        out[header.split()[0] if header else ""] = ">" + chunk.rstrip("\n") + "\n"
    out.pop("", None)
    return out


def _fetch_with_fallbacks(plan, batch_size: int = 100, max_rounds: int = 5):
    """Retrieve one sequence per importation, substituting when one is missing.

    Requests go out in batches (one POST per `batch_size` strains) for
    efficiency, then the response is post-processed to see which strains came
    back. Each importation that got nothing falls back to its cluster's next
    oldest sample and is retried in the following round. A cluster that
    exhausts every sample is dropped -- and because the schedule is written
    from the surviving clusters, its importation disappears from the schedule
    too, instead of leaving the painter an off-by-one in every later pairing.

    Returns (records, dropped): a strain -> FASTA text map for the chosen
    sequences, and the rows that could not be filled.
    """
    pending = plan.copy()
    pending["attempt"] = pending["strain"]
    pending["queue"] = pending["fallbacks"].map(list)
    records: Dict[str, str] = {}
    chosen = {}
    for round_no in range(1, max_rounds + 1):
        targets = sorted(set(pending["attempt"]) - set(records))
        if not targets:
            break
        print(f"\n  Fetch round {round_no}: {len(targets):,} sequence(s)")
        text = fetch_sequences_by_strain_id(targets, batch_size=batch_size)
        if text is None:
            print("  ERROR: the API request failed; aborting fetch.")
            break
        got = _parse_fasta_records(text)
        records.update(got)
        missing_mask = ~pending["attempt"].isin(records)
        filled = pending[~missing_mask]
        for row in filled.itertuples(index=False):
            chosen[row.cluster_id] = row.attempt
        pending = pending[missing_mask].copy()
        if pending.empty:
            break
        n_missing = len(pending)
        advanced = pending["queue"].map(bool)
        if advanced.any():
            pending.loc[advanced, "attempt"] = pending.loc[advanced, "queue"].map(
                lambda q: q.pop(0))
        print(f"    {n_missing:,} importation(s) unfilled; "
              f"{int(advanced.sum()):,} have another sample to try.")
        pending = pending[advanced]
        if pending.empty:
            break
    dropped = plan[~plan["cluster_id"].isin(chosen)]
    return {cid: records[s] for cid, s in chosen.items()}, dropped


def run_seed_mode(args):
    """Schedule and seed sequences, derived together from the cluster table."""
    print("--- Running in Seed Mode ---")
    pango_lineages = [p.strip() for p in args.pango.split(',')]
    print(f"Processing for {len(pango_lineages)} Pango lineage(s): {pango_lineages}")

    if not PANGO_ALIASOR_AVAILABLE: exit(1)

    output_folder_path = Path(args.output_folder)

    input_tsv_path: Optional[Path] = None
    if args.input_file:
        input_tsv_path = Path(args.input_file)
        if not input_tsv_path.is_file(): print(f"Error: Provided input file '{input_tsv_path}' does not exist."); exit(1)
        print(f"Using provided input file: {input_tsv_path.resolve()}")
    else:
        default_tsv_in_output = output_folder_path / DEFAULT_TSV_BASENAME
        if default_tsv_in_output.is_file():
            print(f"Warning: Using existing file in output folder: '{default_tsv_in_output.resolve()}'.")
            input_tsv_path = default_tsv_in_output
        else:
            print(f"Input file not provided. Attempting to download from {DEFAULT_TSV_URL}.")
            if download_file(DEFAULT_TSV_URL, default_tsv_in_output,
                             insecure=args.insecure_download): input_tsv_path = default_tsv_in_output
            else: print(f"Error: Failed to download the default input file."); exit(1)
    if not input_tsv_path: print("Error: Could not determine input TSV file path."); exit(1)

    print(f"Loading data from {input_tsv_path}...")
    try:
        df = pd.read_csv(input_tsv_path, sep='\t', compression=('gzip' if str(input_tsv_path).endswith('.gz') else None))
    except Exception as e: print(f"Error reading TSV file '{input_tsv_path}': {e}"); exit(1)

    thresholds = _parse_variant_starts(args.variant_start)
    clusters = cluster_seeds.prepare_clusters(
        df, args.state, pango_lineages,
        include_sublineages=args.include_sublineages,
        base_map_fn=make_variant_base_map,
        thresholds=thresholds)
    if clusters.empty:
        print("No clusters found for this state and lineage set. Exiting."); exit(0)
    print(f"  {len(clusters):,} importation cluster(s) after lineage filtering")
    print(f"  Sublineages: {'INCLUDED' if args.include_sublineages else 'EXCLUDED'}")

    if args.outlier_method != 'none':
        clusters = _drop_outlier_clusters(clusters, args)
        if clusters.empty:
            print("Every cluster was filtered as an outlier. Exiting."); exit(0)

    plan = cluster_seeds.seed_plan(clusters)

    # Provenance: the strain chosen per variant, before any substitution.
    for pango in pango_lineages:
        ids = plan.loc[plan["variant"] == pango, "strain"].tolist()
        if not ids:
            continue
        name = f"{args.state.replace(' ', '_')}_{_sanitize(pango)}_seed_strains.txt"
        (output_folder_path / name).write_text("\n".join(ids))
        print(f"  {len(ids):,} seed strain(s) for {pango} -> {name}")

    if args.no_download:
        print("\n--no_download: skipping sequence retrieval.")
        _write_schedules(clusters, output_folder_path, args)
        return

    records, dropped = _fetch_with_fallbacks(plan, batch_size=args.batch_size,
                                             max_rounds=args.max_fetch_rounds)
    if len(dropped):
        banner = "!" * 78
        print(banner, file=sys.stderr)
        print(f"WARNING: {len(dropped):,} of {len(plan):,} importations have NO retrievable",
              file=sys.stderr)
        print("  sequence after exhausting every sample in their cluster. Those clusters are",
              file=sys.stderr)
        print("  DROPPED from both the seed FASTA and the schedule, so the two stay aligned.",
              file=sys.stderr)
        for row in dropped.head(10).itertuples(index=False):
            print(f"    {row.cluster_id}  {row.variant}  {pd.Timestamp(row.intro_date).date()}",
                  file=sys.stderr)
        if len(dropped) > 10:
            print(f"    ... and {len(dropped) - 10:,} more", file=sys.stderr)
        print(banner, file=sys.stderr)
        kept = set(plan["cluster_id"]) - set(dropped["cluster_id"])
        clusters = clusters[clusters["cluster_id"].isin(kept)]
        plan = plan[plan["cluster_id"].isin(kept)]

    # One FASTA per variant, written in schedule order: record i founds
    # importation i, which is what the painter assumes.
    for pango in pango_lineages:
        rows = plan[plan["variant"] == pango].sort_values(["intro_date", "cluster_id"])
        if rows.empty:
            continue
        name = f"{args.state.replace(' ', '_')}_{_sanitize(pango)}_seed_sequences.fasta"
        path = output_folder_path / name
        used = []
        with open(path, "w") as fh:
            for row in rows.itertuples(index=False):
                rec = records[row.cluster_id]
                fh.write(rec)
                used.append(rec[1:].split("\n", 1)[0].split()[0])
        print(f"  {len(rows):,} sequence(s) for {pango}, in importation order -> {path.resolve()}")

        # The manifest says which importation each record founds. The painter
        # reads it to check, day by day, that its importations line up with
        # the seeds -- a plain FASTA carries no dates, which is how a
        # scrambled pairing went unnoticed before.
        man = pd.DataFrame({
            "order": range(len(rows)),
            "strain": used,
            "cluster_id": rows["cluster_id"].values,
            "intro_date": pd.to_datetime(rows["intro_date"]).dt.strftime("%Y-%m-%d").values,
        })
        t0 = _tick_zero(args)
        if t0 is not None:
            man["tick"] = (pd.to_datetime(man["intro_date"]) - t0).dt.days
        mname = name.replace("_seed_sequences.fasta", "_seed_manifest.csv")
        man.to_csv(output_folder_path / mname, index=False)
        print(f"  seed manifest -> {mname}")

    _write_schedules(clusters, output_folder_path, args)


def run_bulk_mode(args):
    """Contains all logic for the new direct-to-CovSpectrum workflow."""
    print("--- Running in Bulk Download Mode ---")
    pango_lineages = [p.strip() for p in args.pango.split(',')]
    output_folder_path = Path(args.output_folder)

    # Validate that if one date is given, the other is too
    if (args.date_from and not args.date_to) or (not args.date_from and args.date_to):
        print("Error: If specifying a date range, both --date_from and --date_to are required.")
        exit(1)

    for pango in pango_lineages:
        pango_sanitized = pango.replace('.', '_').replace('/', '_')
        state_sanitized = args.state.replace(' ', '_')
        
        # Adjust filename based on whether dates are provided
        if args.date_from and args.date_to:
            output_filename = f"{state_sanitized}_{pango_sanitized}_{args.date_from}_{args.date_to}.fasta"
        else:
            output_filename = f"{state_sanitized}_{pango_sanitized}_all-dates.fasta"
        if not args.include_sublineages:
            output_filename = output_filename.replace(".fasta", "_exact.fasta")
            
        output_filepath = output_folder_path / output_filename
        
        fetch_sequences_by_metadata(pango, args.state, args.date_from, args.date_to,
                                    output_filepath, args.include_sublineages)

def main():
    parser = argparse.ArgumentParser(description="Prepare sequence sets from cluster data or by direct metadata query.")
    parser.add_argument("--state", required=True, type=str, help="US state to filter/query for (e.g., 'Virginia').")
    parser.add_argument("--pango", required=True, type=str, help="Comma-separated list of Pango lineages (e.g., 'B.1.1.7,B.1.617.2').")
    parser.add_argument("--output_folder", required=True, type=str, help="Path for output files.")
    parser.add_argument("--insecure_download", "--insecure-download",
                        dest="insecure_download", action='store_true',
                        help="[Seed Mode] Skip TLS certificate verification when fetching the "
                             "cluster TSV. Workaround for the expired clustertracker.gi.ucsc.edu "
                             "certificate (upstream, since 2025-07-02).")
    parser.add_argument("--include_sublineages", "--include-sublineages",
                        dest="include_sublineages", action='store_true', default=True,
                        help="Include descendant lineages, i.e. query 'B.1.617.2*' and roll "
                             "descendants up to the requested ancestor. Default: ON.")
    parser.add_argument("--no_include_sublineages", "--no-include-sublineages",
                        dest="include_sublineages", action='store_false',
                        help="Match the named lineage EXACTLY, excluding descendants. For "
                             "Virginia Delta this is 836 samples instead of 21,034.")
    parser.add_argument("--seed_mode", action='store_true', help="Enable seed-finding mode, using the cluster tracker file and outlier detection.")
    parser.add_argument("--date_from", type=str, help="[Bulk Mode] Optional start date for query (YYYY-MM-DD).")
    parser.add_argument("--date_to", type=str, help="[Bulk Mode] Optional end date for query (YYYY-MM-DD).")
    parser.add_argument("--input_file", type=str, help="[Seed Mode] Optional path to the input cluster TSV file.")
    parser.add_argument("--outlier_method", type=str, choices=['none', 'iqr', 'zscore', 'chaining'], default='none', help="[Seed Mode] Method for date outlier detection.")
    parser.add_argument("--chaining_max_gap_weeks", type=int, default=6, help="[Seed Mode] Max gap in weeks for 'chaining' outlier method.")
    parser.add_argument("--rescue_cluster_size", type=int, default=2, help="[Seed Mode] Minimum sample_count to rescue a potential outlier.")
    parser.add_argument("--rescue_cluster_days", type=int, default=365, help="[Seed Mode] Max time span within a cluster for rescue.")
    parser.add_argument("--no_download", action='store_true', help="[Seed Mode] If specified, only generate seed strain ID files and skip downloading sequences.")
    parser.add_argument("--variant_start", "--variant-start", dest="variant_start",
                        type=str, default=None,
                        help="[Seed Mode] Earliest plausible date per variant, as "
                             "VARIANT=YYYY-MM-DD[,...]. A cluster whose first sample "
                             "predates its variant's date is re-dated to its first "
                             "sample after it. Default: "
                             + ",".join(f"{k}={v}" for k, v in
                                        cluster_seeds.VARIANT_START_DEFAULTS.items()))
    parser.add_argument("--batch_size", "--batch-size", dest="batch_size", type=int,
                        default=100,
                        help="[Seed Mode] Strains per API request (default: 100).")
    parser.add_argument("--max_fetch_rounds", "--max-fetch-rounds",
                        dest="max_fetch_rounds", type=int, default=5,
                        help="[Seed Mode] How many times to retry unfilled importations "
                             "with their cluster's next oldest sample (default: 5).")
    parser.add_argument("--tick_zero", "--tick-zero", dest="tick_zero", type=str,
                        default=None,
                        help="[Seed Mode] Calendar date of the ABM's tick 0. Given, the "
                             "schedule is also written on absolute ticks, which is what "
                             f"EpiHiper seeding consumes. EpiHiper used {cluster_seeds.EPIHIPER_TICK_ZERO} "
                             "for every state, so their ticks share one calendar.")
    parser.add_argument("--abm_config", "--abm-config", dest="abm_config", type=str,
                        default=None,
                        help="[Seed Mode] Read tick 0 from an ABM config instead "
                             "(EpiHiper's config.json 'tickZero'), so the date is not "
                             "duplicated. Ignored if --tick_zero is given.")
    
    args = parser.parse_args()
    
    output_folder_path = Path(args.output_folder)
    output_folder_path.mkdir(parents=True, exist_ok=True)

    if args.seed_mode:
        run_seed_mode(args)
    else:
        run_bulk_mode(args)

    print("\nScript finished.")

if __name__ == "__main__":
    main()
