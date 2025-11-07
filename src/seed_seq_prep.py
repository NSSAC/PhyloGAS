import argparse
import pandas as pd
import numpy as np
import requests
import os
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional
from scipy.stats import zscore # Added for Z-score outlier detection

# Attempt to import pango_aliasor and define make_variant_base_map
try:
    from pango_aliasor.aliasor import Aliasor

    def make_variant_base_map(base_variants: List[str], recombinant: bool = False) -> Dict[str, str]:
        """
        Regularizes variant names using pango_aliasor.
        Maps subvariants to their defined base variants.
        """
        namer = Aliasor()
        namer.enable_expansion() # As per the notebook
        all_rules = namer.partition_focus(base_variants, recombinant=recombinant)
        lineage_base_map = {k: v for v, ks in all_rules.items() for k in ks}
        
        # Ensure base variants themselves are mapped if not covered by rules
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
        response.raise_for_status()  # Raises an HTTPError for bad responses (4XX or 5XX)
        with open(dest_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print("Download complete.")
        return True
    except requests.exceptions.RequestException as e:
        print(f"Error downloading file: {e}")
        if dest_path.exists(): # Clean up partial download
            dest_path.unlink()
        return False

def fetch_sequences_from_covspectrum(strain_ids: List[str], batch_size: int = 100) -> Optional[str]:
    """
    Fetches aligned nucleotide sequences from CovSpectrum API for given strain IDs.
    Returns a single string containing all sequences in FASTA format, or None on error.
    """
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
            # Check if response is empty or indicates an issue (though API omits not found strains)
            if not fasta_data.strip() and len(batch_ids) > 0 :
                 print(f"    Warning: Batch {i//batch_size + 1} returned empty data from API, though strains were requested.")
            all_fasta_content.append(fasta_data)
        except requests.exceptions.HTTPError as e:
            print(f"  HTTP Error for batch {i//batch_size + 1}: {e}")
            print(f"  Response content: {e.response.text[:500]}...") # Show some of the error
            return None # Indicate overall failure
        except requests.exceptions.RequestException as e:
            print(f"  Request Error for batch {i//batch_size + 1}: {e}")
            return None # Indicate overall failure
        except Exception as e:
            print(f"  An unexpected error occurred during API call for batch {i//batch_size + 1}: {e}")
            return None
    print("Sequence fetching complete.")
    return "".join(all_fasta_content)

def detect_outliers_iqr(date_series: pd.Series, factor: float = 1.5) -> pd.Series:
    """
    Detects outliers in a pandas Series of datetime objects using IQR method.
    Returns a boolean Series where True indicates an outlier.
    Assumes date_series contains already parsed datetime objects and NAs are handled before.
    """
    if date_series.empty:
        return pd.Series(dtype=bool, index=date_series.index)

    Q1 = date_series.quantile(0.25)
    Q3 = date_series.quantile(0.75)
    IQR = Q3 - Q1

    # If IQR is zero (e.g., many identical dates), no outliers by this method unless factor is also 0.
    if IQR == pd.Timedelta(0):
        return pd.Series(False, index=date_series.index, dtype=bool)

    lower_bound = Q1 - factor * IQR
    upper_bound = Q3 + factor * IQR
    
    outliers_mask = (date_series < lower_bound) | (date_series > upper_bound)
    return outliers_mask

def detect_outliers_zscore(date_series: pd.Series, threshold: float = 2.0) -> pd.Series:
    """
    Detects outliers in a pandas Series of datetime objects using Z-score method.
    Returns a boolean Series where True indicates an outlier.
    Assumes date_series contains already parsed datetime objects and NAs are handled before.
    """
    if date_series.empty or date_series.nunique() < 2: # Z-score needs variance, min 2 unique values
        return pd.Series(False, index=date_series.index, dtype=bool)

    # Convert dates to days since the minimum date in the series for Z-score calculation
    numeric_dates = (date_series - date_series.min()).dt.days
    
    # Re-check for unique values after conversion to numeric
    if numeric_dates.nunique() < 2:
         return pd.Series(False, index=date_series.index, dtype=bool)

    z_scores = zscore(numeric_dates)
    z_scores = np.nan_to_num(z_scores, nan=0.0) # Handle potential NaNs if std dev is 0

    outliers_mask = np.abs(z_scores) > threshold
    return pd.Series(outliers_mask, index=date_series.index)

def detect_outliers_chaining(date_series: pd.Series, max_gap_weeks: int = 6) -> pd.Series:
    """
    Detects outliers in a pandas Series of datetime objects using a chaining method.
    An outlier is a date that is more than 'max_gap_weeks' from both its 
    preceding and succeeding date in the sorted list of dates.
    Returns a boolean Series (same index as input) where True indicates an outlier.
    Assumes date_series contains already parsed datetime objects and NAs are handled before.
    """
    if date_series.empty:
        return pd.Series(dtype=bool, index=date_series.index)

    # Work with a DataFrame to keep original indices and allow easy diff calculations
    df_proc = date_series.to_frame(name='date').copy()

    # Sort by date to calculate differences with neighbors
    df_sorted = df_proc.sort_values(by='date').copy()

    # If all dates are the same, or only one distinct date exists, no outliers by this method.
    if df_sorted['date'].nunique() <= 1:
        return pd.Series(False, index=date_series.index, dtype=bool)

    max_delta = pd.Timedelta(weeks=max_gap_weeks)

    # Difference to the previous date in sorted order
    df_sorted['diff_prev'] = df_sorted['date'].diff()
    # Difference to the next date in sorted order
    df_sorted['diff_next'] = df_sorted['date'].diff(-1).abs() # .abs() because diff(-1) is current - next

    # A date is far from previous if diff > max_delta or it's the first date (NaT diff_prev)
    is_far_from_prev = (df_sorted['diff_prev'] > max_delta) | df_sorted['diff_prev'].isna()
    # A date is far from next if diff > max_delta or it's the last date (NaT diff_next)
    is_far_from_next = (df_sorted['diff_next'] > max_delta) | df_sorted['diff_next'].isna()

    # An outlier is far from BOTH previous AND next
    df_sorted['is_outlier'] = is_far_from_prev & is_far_from_next
    
    # Create a Series from df_sorted results, indexed by its original index values (from df_proc/date_series)
    # then reindex to match the exact original input date_series's index and order.
    outlier_series_sorted_index = pd.Series(df_sorted['is_outlier'].values, index=df_sorted.index)
    return outlier_series_sorted_index.reindex(date_series.index)


def generate_schedule_file(df: pd.DataFrame, output_path: Path, state: str, pango: str):
    """Generates an importation schedule CSV from the final processed DataFrame."""
    if df.empty:
        print(f"Cannot generate schedule for {pango} in {state}: No data remains.")
        return

    df_schedule = df.copy()
    min_date = df_schedule['earliest_date'].min()
    
    # Create tick column (days since first appearance of this variant in this context)
    df_schedule['tick'] = (df_schedule['earliest_date'] - min_date).dt.days
    
    # Group by day to count clusters and sum samples
    grouped = df_schedule.groupby(['tick', 'earliest_date', 'pango_regularized'])
    
    # Count clusters (.size()) and sum samples (.sum())
    result = grouped.size().reset_index(name='clusters')
    result['sample_count'] = grouped['sample_count'].sum().values
    
    # Rename columns to match notebook output
    result = result.rename(columns={'earliest_date': 'date', 'pango_regularized': 'variant'})
    
    pango_sanitized = pango.replace('.', '_').replace('/', '_')
    state_sanitized = state.replace(' ', '_')
    schedule_filename = f"{state_sanitized}_{pango_sanitized}_schedule.csv"
    schedule_filepath = output_path / schedule_filename
    
    result.to_csv(schedule_filepath, index=False)
    print(f"Importation schedule for {pango} saved to: {schedule_filepath.resolve()}")


def main():
    parser = argparse.ArgumentParser(description="Process SARS-CoV-2 cluster data, identify seed samples, and fetch their sequences.")
    parser.add_argument("--state", required=True, type=str, help="US state to filter data for.")
    parser.add_argument("--pango", required=True, type=str, help="Comma-separated list of Pango lineages to process.")
    parser.add_argument("--input_file", type=str, help="Optional path to the input TSV file.")
    parser.add_argument("--output_folder", required=True, type=str, help="Path for output files.")
    parser.add_argument("--outlier_method", type=str, choices=['none', 'iqr', 'zscore', 'chaining'], default='none', help="Method for date outlier detection.")
    parser.add_argument("--chaining_max_gap_weeks", type=int, default=6, help="Max gap in weeks for 'chaining' outlier method.")
    parser.add_argument("--no_download", action='store_true', help="If specified, only generate seed strain ID files and skip downloading sequences.")
    parser.add_argument("--generate_schedule", action='store_true', help="If specified, generate an importation schedule CSV for each variant.")
    parser.add_argument("--rescue_cluster_size", type=int, default=2, help="Minimum sample_count above which to rescue a potential outlier cluster (default: 2).")
    parser.add_argument("--rescue_cluster_days", type=int, default=365, help="Maximum time span in days within a cluster to be considered for rescue (default: 21).")
    
    args = parser.parse_args()
    pango_lineages = [p.strip() for p in args.pango.split(',')]
    print(f"Processing for {len(pango_lineages)} Pango lineage(s): {pango_lineages}")
    if not PANGO_ALIASOR_AVAILABLE:
        print("Exiting due to missing 'pango_aliasor' library.")
        exit(1)

    # --- 1. Setup output folder ---
    output_folder_path = Path(args.output_folder)
    output_folder_path.mkdir(parents=True, exist_ok=True)
    print(f"Using output folder: {output_folder_path.resolve()}")

    # --- 2. Determine and prepare input TSV file ---
    input_tsv_path: Optional[Path] = None
    if args.input_file:
        input_tsv_path = Path(args.input_file)
        if not input_tsv_path.is_file():
            print(f"Error: Provided input file '{input_tsv_path}' does not exist.")
            exit(1)
        print(f"Using provided input file: {input_tsv_path.resolve()}")
    else:
        default_tsv_in_output = output_folder_path / DEFAULT_TSV_BASENAME
        if default_tsv_in_output.is_file():
            print(f"Warning: Input file not provided. Using existing file in output folder: '{default_tsv_in_output.resolve()}'.")
            input_tsv_path = default_tsv_in_output
        else:
            print(f"Input file not provided. Attempting to download from {DEFAULT_TSV_URL}.")
            if download_file(DEFAULT_TSV_URL, default_tsv_in_output):
                input_tsv_path = default_tsv_in_output
            else:
                print(f"Error: Failed to download the default input file. Please provide it using --input_file or ensure internet connectivity.")
                exit(1)
    
    if not input_tsv_path: # Should not happen if logic is correct
        print("Error: Could not determine input TSV file path.")
        exit(1)

    # --- 3. Load and process data ---
    print(f"Loading data from {input_tsv_path}...")
    try:
        df = pd.read_csv(input_tsv_path, sep='\t', compression=('gzip' if str(input_tsv_path).endswith('.gz') else None))
    except Exception as e: print(f"Error reading TSV file '{input_tsv_path}': {e}"); exit(1)
    print(f"Loaded {len(df)} rows from TSV.")
    
    if 'annotation_2' not in df.columns: print("Error: 'annotation_2' column not found."); exit(1)
    try:
        replace_map = make_variant_base_map(pango_lineages)
        df['pango_regularized'] = df['annotation_2'].map(replace_map)
        df['pango_regularized'].fillna(df['annotation_2'], inplace=True)
    except Exception as e: print(f"Error during pango_aliasor processing: {e}"); exit(1)
    
    if 'region' not in df.columns: print("Error: 'region' column not found."); exit(1)
    df_filtered_initial = df[(df['region'] == args.state) & (df['pango_regularized'].isin(pango_lineages))].copy()
    if df_filtered_initial.empty: print("No data found for state and specified Pango lineages. Exiting."); exit(0)

    # --- Loop through each Pango lineage for specific processing ---
    for pango in pango_lineages:
        print(f"\n{'='*20} Processing: {pango} {'='*20}")

        df_variant = df_filtered_initial[df_filtered_initial['pango_regularized'] == pango].copy()
        if df_variant.empty: print(f"No data for lineage {pango}. Skipping."); continue

        df_variant['earliest_date'] = pd.to_datetime(df_variant['earliest_date'], errors='coerce')
        df_variant.dropna(subset=['earliest_date'], inplace=True)
        if df_variant.empty: print(f"No valid 'earliest_date' entries for {pango}. Skipping."); continue

        min_date_orig, max_date_orig = df_variant['earliest_date'].min(), df_variant['earliest_date'].max()
        print(f"Original date range for {pango}: {min_date_orig.strftime('%Y-%m-%d')} to {max_date_orig.strftime('%Y-%m-%d')}")

        # --- Refined Outlier Pruning ---
        if args.outlier_method != 'none':
            outliers_mask = pd.Series(False, index=df_variant.index, dtype=bool)
            if args.outlier_method == 'iqr': outliers_mask = detect_outliers_iqr(df_variant['earliest_date'])
            elif args.outlier_method == 'zscore': outliers_mask = detect_outliers_zscore(df_variant['earliest_date'])
            elif args.outlier_method == 'chaining': outliers_mask = detect_outliers_chaining(df_variant['earliest_date'], args.chaining_max_gap_weeks)

            potential_outliers = df_variant[outliers_mask]
            if not potential_outliers.empty:
                print(f"Identified {len(potential_outliers)} potential date outliers for {pango} using {args.outlier_method} method. Checking for rescue candidates...")
                
                final_outliers_indices = []
                for idx, row in potential_outliers.iterrows():
                    is_singlet = row['sample_count'] < args.rescue_cluster_size
                    if is_singlet:
                        print(f"  - Pruning singlet cluster outlier (size {row['sample_count']}) with date {row['earliest_date'].strftime('%Y-%m-%d')}.")
                        final_outliers_indices.append(idx)
                        continue
                    
                    # Check internal time span for multi-sample clusters
                    sample_dates = [pd.to_datetime(s.split('|')[-1]) for s in row['samples'].split(',')]
                    internal_span = max(sample_dates) - min(sample_dates)
                    
                    if internal_span <= pd.Timedelta(days=args.rescue_cluster_days):
                        print(f"  - Rescuing potential outlier (size {row['sample_count']}, span {internal_span.days} days) with date {row['earliest_date'].strftime('%Y-%m-%d')}. Evidence of sustained transmission.")
                    else:
                        print(f"  - Pruning multi-sample cluster outlier (size {row['sample_count']}, span {internal_span.days} days > {args.rescue_cluster_days} days) with date {row['earliest_date'].strftime('%Y-%m-%d')}.")
                        final_outliers_indices.append(idx)
                
                if final_outliers_indices:
                    df_variant.drop(final_outliers_indices, inplace=True)
            else:
                print(f"No date outliers detected for {pango} using {args.outlier_method} method.")

        if df_variant.empty: print(f"All data for {pango} removed after outlier pruning. Skipping."); continue
        
        # --- Continue with processing the pruned df_variant ---
        min_date_new, max_date_new = df_variant['earliest_date'].min(), df_variant['earliest_date'].max()
        print(f"Date range for {pango} after pruning: {min_date_new.strftime('%Y-%m-%d')} to {max_date_new.strftime('%Y-%m-%d')}")
        
        # --- Generate Schedule (if requested) ---
        if args.generate_schedule:
            generate_schedule_file(df_variant, output_folder_path, args.state, pango)
            
        # --- Identify and Save Seed Strains ---
        df_variant['strain_seed'] = df_variant['samples'].str.split(',', n=1, expand=True)[0].str.split('|', n=1, expand=True)[0]
        unique_strain_seeds = df_variant['strain_seed'].dropna().unique().tolist()
        if not unique_strain_seeds: print(f"No seed strains identified for {pango}. Skipping."); continue
        
        pango_sanitized = pango.replace('.', '_').replace('/', '_')
        seed_ids_filename = f"{args.state.replace(' ', '_')}_{pango_sanitized}_seed_strains.txt"
        seed_ids_filepath = output_folder_path / seed_ids_filename
        with open(seed_ids_filepath, 'w') as f: f.write('\n'.join(unique_strain_seeds))
        print(f"Identified {len(unique_strain_seeds)} unique seed strains for {pango}, IDs saved to: {seed_ids_filepath.resolve()}")
        
        # --- Fetch Sequences (Conditional) ---
        if not args.no_download:
            fasta_sequences = fetch_sequences_from_covspectrum(unique_strain_seeds)
            if fasta_sequences and fasta_sequences.strip():
                fasta_filename = f"{args.state.replace(' ', '_')}_{pango_sanitized}_seed_sequences.fasta"
                fasta_filepath = output_folder_path / fasta_filename
                with open(fasta_filepath, 'w') as f: f.write(fasta_sequences)
                print(f"Fetched sequences for {pango} saved to: {fasta_filepath.resolve()}")
            else:
                print(f"Failed to fetch or received empty sequences for {pango}.")
        else:
            print(f"Skipping sequence download for {pango} as per --no_download flag.")

    print("\nScript finished for all specified Pango lineages.")

if __name__ == "__main__":
    main()