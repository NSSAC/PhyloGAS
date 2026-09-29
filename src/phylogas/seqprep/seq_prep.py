import argparse
import pandas as pd
import numpy as np
import requests
import os
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional
from scipy.stats import zscore 


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

COVSPECTRUM_API_URL = 'https://lapis.cov-spectrum.org/open/v2/sample/alignedNucleotideSequences'
DEFAULT_TSV_URL = 'https://clustertracker.gi.ucsc.edu/data/hardcoded_clusters.tsv'
DEFAULT_TSV_BASENAME = 'hardcoded_clusters.tsv'

def download_file(url: str, dest_path: Path) -> bool:
    print(f"Downloading {url} to {dest_path}...")
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        with open(dest_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print("Download complete.")
        return True
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


def fetch_sequences_by_metadata(pango: str, state: str, date_from: Optional[str], date_to: Optional[str], output_filepath: Path):
    """Fetches sequences directly from CovSpectrum based on metadata query and streams to a file."""
    date_info = "all dates"
    if date_from and date_to:
        date_info = f"from {date_from} to {date_to}"
    
    print(f"\n--- Starting bulk download for {pango} in {state} ({date_info}) ---")
    
    params = {
        'country': 'USA',
        'division': state,
        'pangoLineage': f'{pango}*',
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

def generate_schedule_file(df: pd.DataFrame, output_path: Path, state: str):
    if df.empty:
        print(f"Cannot generate schedule for {state}: No data remains.")
        return
    df_schedule = df.copy()
    min_date = df_schedule['earliest_date'].min()
    df_schedule['tick'] = (df_schedule['earliest_date'] - min_date).dt.days
    df_schedule['clusters'] = 1
    df_schedule = df_schedule.rename(columns={'earliest_date': 'date', 'pango_regularized': 'variant'})
    output_columns = ['tick', 'date', 'variant', 'clusters', 'sample_count']
    final_schedule_df = df_schedule[output_columns].sort_values(by=['tick', 'variant'])
    state_sanitized = state.replace(' ', '_')
    schedule_filename = f"{state_sanitized}_schedule.csv"
    schedule_filepath = output_path / schedule_filename
    final_schedule_df.to_csv(schedule_filepath, index=False)
    print(f"\nConsolidated importation schedule saved to: {schedule_filepath.resolve()}")


def run_seed_mode(args):
    """Contains all logic for the original seed-finding workflow."""
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
            if download_file(DEFAULT_TSV_URL, default_tsv_in_output): input_tsv_path = default_tsv_in_output
            else: print(f"Error: Failed to download the default input file."); exit(1)
    if not input_tsv_path: print("Error: Could not determine input TSV file path."); exit(1)
    
    print(f"Loading data from {input_tsv_path}...")
    try:
        df = pd.read_csv(input_tsv_path, sep='\t', compression=('gzip' if str(input_tsv_path).endswith('.gz') else None))
    except Exception as e: print(f"Error reading TSV file '{input_tsv_path}': {e}"); exit(1)
    
    if 'annotation_2' not in df.columns: print("Error: 'annotation_2' column not found."); exit(1)
    try:
        replace_map = make_variant_base_map(pango_lineages)
        df['pango_regularized'] = df['annotation_2'].map(replace_map)
        df['pango_regularized'].fillna(df['annotation_2'], inplace=True)
    except Exception as e: print(f"Error during pango_aliasor processing: {e}"); exit(1)
    
    if 'region' not in df.columns: print("Error: 'region' column not found."); exit(1)
    df_filtered_initial = df[(df['region'] == args.state) & (df['pango_regularized'].isin(pango_lineages))].copy()
    if df_filtered_initial.empty: print("No data found for state and specified Pango lineages. Exiting."); exit(0)

    processed_variant_dfs = []

    for pango in pango_lineages:
        print(f"\n{'='*20} Processing: {pango} {'='*20}")
        df_variant = df_filtered_initial[df_filtered_initial['pango_regularized'] == pango].copy()
        if df_variant.empty: continue

        for col in ['earliest_date', 'latest_date']:
             df_variant[col] = pd.to_datetime(df_variant[col], errors='coerce')
        df_variant.dropna(subset=['earliest_date', 'latest_date'], inplace=True)
        if df_variant.empty: continue
        if args.outlier_method != 'none':
            outliers_mask = pd.Series(False, index=df_variant.index, dtype=bool)
            if args.outlier_method == 'iqr': outliers_mask = detect_outliers_iqr(df_variant['earliest_date'])
            elif args.outlier_method == 'zscore': outliers_mask = detect_outliers_zscore(df_variant['earliest_date'])
            elif args.outlier_method == 'chaining': outliers_mask = detect_outliers_chaining(df_variant['earliest_date'], args.chaining_max_gap_weeks)
            
            potential_outliers = df_variant[outliers_mask]
            if not potential_outliers.empty:
                print(f"Identified {len(potential_outliers)} potential outliers. Checking for rescue candidates...")
                final_outliers_indices = []
                for idx, row in potential_outliers.iterrows():
                    internal_span = row['latest_date'] - row['earliest_date']
                    is_rescuable = (row['sample_count'] >= args.rescue_cluster_size and internal_span <= pd.Timedelta(days=args.rescue_cluster_days))
                    if not is_rescuable:
                        final_outliers_indices.append(idx)
                if final_outliers_indices:
                    df_variant.drop(final_outliers_indices, inplace=True)

        if df_variant.empty: continue
        processed_variant_dfs.append(df_variant)

        df_variant['strain_seed'] = df_variant['samples'].str.split(',', n=1, expand=True)[0].str.split('|', n=1, expand=True)[0]
        unique_strain_seeds = df_variant['strain_seed'].dropna().unique().tolist()
        if not unique_strain_seeds: continue
        
        pango_sanitized = pango.replace('.', '_').replace('/', '_')
        seed_ids_filename = f"{args.state.replace(' ', '_')}_{pango_sanitized}_seed_strains.txt"
        seed_ids_filepath = output_folder_path / seed_ids_filename
        with open(seed_ids_filepath, 'w') as f: f.write('\n'.join(unique_strain_seeds))
        print(f"Identified {len(unique_strain_seeds)} seed strains for {pango}, IDs saved to: {seed_ids_filepath.resolve()}")
        
        if not args.no_download:
            fasta_sequences = fetch_sequences_by_strain_id(unique_strain_seeds)
            if fasta_sequences and fasta_sequences.strip():
                fasta_filename = f"{args.state.replace(' ', '_')}_{pango_sanitized}_seed_sequences.fasta"
                fasta_filepath = output_folder_path / fasta_filename
                with open(fasta_filepath, 'w') as f: f.write(fasta_sequences)
                print(f"Fetched sequences for {pango} saved to: {fasta_filepath.resolve()}")
            else:
                print(f"Failed to fetch or received empty sequences for {pango}.")
        else:
            print(f"Skipping sequence download for {pango} as per --no_download flag.")

    if processed_variant_dfs:
        final_df_for_schedule = pd.concat(processed_variant_dfs, ignore_index=True)
        generate_schedule_file(final_df_for_schedule, output_folder_path, args.state)


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
            
        output_filepath = output_folder_path / output_filename
        
        fetch_sequences_by_metadata(pango, args.state, args.date_from, args.date_to, output_filepath)

def main():
    parser = argparse.ArgumentParser(description="Prepare sequence sets from cluster data or by direct metadata query.")
    parser.add_argument("--state", required=True, type=str, help="US state to filter/query for (e.g., 'Virginia').")
    parser.add_argument("--pango", required=True, type=str, help="Comma-separated list of Pango lineages (e.g., 'B.1.1.7,B.1.617.2').")
    parser.add_argument("--output_folder", required=True, type=str, help="Path for output files.")
    parser.add_argument("--seed_mode", action='store_true', help="Enable seed-finding mode, using the cluster tracker file and outlier detection.")
    parser.add_argument("--date_from", type=str, help="[Bulk Mode] Optional start date for query (YYYY-MM-DD).")
    parser.add_argument("--date_to", type=str, help="[Bulk Mode] Optional end date for query (YYYY-MM-DD).")
    parser.add_argument("--input_file", type=str, help="[Seed Mode] Optional path to the input cluster TSV file.")
    parser.add_argument("--outlier_method", type=str, choices=['none', 'iqr', 'zscore', 'chaining'], default='none', help="[Seed Mode] Method for date outlier detection.")
    parser.add_argument("--chaining_max_gap_weeks", type=int, default=6, help="[Seed Mode] Max gap in weeks for 'chaining' outlier method.")
    parser.add_argument("--rescue_cluster_size", type=int, default=2, help="[Seed Mode] Minimum sample_count to rescue a potential outlier.")
    parser.add_argument("--rescue_cluster_days", type=int, default=365, help="[Seed Mode] Max time span within a cluster for rescue.")
    parser.add_argument("--no_download", action='store_true', help="[Seed Mode] If specified, only generate seed strain ID files and skip downloading sequences.")
    
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
