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
    print("The script cannot continue without this library.")
    # Define a dummy function so the script can be parsed, but exit if called.
    def make_variant_base_map(base_variants: List[str], recombinant: bool = False) -> Dict[str, str]:
        raise ImportError("pango_aliasor not found, function unusable.")
    # It's better to exit early in main if this is the case.
    PANGO_ALIASOR_AVAILABLE = False
else:
    PANGO_ALIASOR_AVAILABLE = True


# --- Constants ---
COVSPECTRUM_API_URL = 'https://lapis.cov-spectrum.org/open/v2/sample/alignedNucleotideSequences'
DEFAULT_TSV_URL = 'https://clustertracker.gi.ucsc.edu/data/hardcoded_clusters.tsv'
DEFAULT_TSV_BASENAME = 'hardcoded_clusters.tsv'


def download_file(url: str, dest_path: Path) -> bool:
    """Downloads a file from a URL to a destination path."""
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
            response = requests.post(COVSPECTRUM_API_URL, json=payload, timeout=120) # Increased timeout
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

def main():
    parser = argparse.ArgumentParser(description="Process SARS-CoV-2 cluster data, identify seed samples, and fetch their sequences.")
    parser.add_argument("--state", required=True, type=str, help="US state to filter data for (e.g., 'Washington', 'California'). Corresponds to 'region' column.")
    parser.add_argument("--pango", required=True, type=str, help="Pango lineage to filter data for (e.g., 'B.1.1.7', 'BA.1'). This will be matched against regularized lineage names.")
    parser.add_argument("--input_file", type=str, help="Optional path to the input TSV file. If not provided, the script will attempt to use/download 'hardcoded_clusters.tsv'.")
    parser.add_argument("--output_folder", required=True, type=str, help="Path to the folder where output files (seed IDs, FASTA sequences, downloaded TSV) will be saved.")
    parser.add_argument("--outlier_method", type=str, choices=['none', 'iqr', 'zscore', 'chaining'], default='none', help="Method for 'earliest_date' outlier detection: 'none' (default), 'iqr', 'zscore', or 'chaining'.")
    parser.add_argument("--chaining_max_gap_weeks", type=int, default=6, help="Max gap in weeks for 'chaining' outlier method (default: 6).")
    
    args = parser.parse_args()

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
        # The notebook file had .gz, the URL is plain .tsv
        compression = 'gzip' if str(input_tsv_path).endswith('.gz') else None
        df = pd.read_csv(input_tsv_path, sep='\t', compression=compression)
    except Exception as e:
        print(f"Error reading TSV file '{input_tsv_path}': {e}")
        exit(1)
    
    print(f"Loaded {len(df)} rows from TSV.")

    # --- 4. Pango Lineage Regularization ---
    print("Regularizing Pango lineages using 'annotation_2' column...")
    if 'annotation_2' not in df.columns:
        print(f"Error: 'annotation_2' column not found in the input TSV. Available columns: {df.columns.tolist()}")
        exit(1)

    extended_target_list = [args.pango] # Regularize based on the target pango
    try:
        replace_map = make_variant_base_map(extended_target_list)
    except Exception as e:
        print(f"Error during pango_aliasor processing: {e}")
        print("This might be due to an issue with the pango_aliasor library or its data.")
        exit(1)

    df['pango_regularized'] = df['annotation_2'].map(replace_map)
    df['pango_regularized'].fillna(df['annotation_2'], inplace=True) # Keep original if not in map
    print(f"Pango lineage regularization complete. New column 'pango_regularized' created.")

    # --- 5. Filter data by state and pango lineage ---
    print(f"Filtering data for state: '{args.state}' and Pango lineage (regularized): '{args.pango}'...")
    if 'region' not in df.columns:
        print(f"Error: 'region' column not found in the input TSV. Needed for state filtering. Available columns: {df.columns.tolist()}")
        exit(1)
        
    original_row_count = len(df)
    df_filtered = df[(df['region'] == args.state) & (df['pango_regularized'] == args.pango)].copy() # Use .copy()

    print(f"Filtered from {original_row_count} to {len(df_filtered)} rows based on state and Pango lineage.")
    if df_filtered.empty:
        print(f"No data found for state '{args.state}' and Pango lineage '{args.pango}'. Exiting.")
        exit(0)

    # --- 6. Date Handling and Outlier Pruning ---
    print("Processing 'earliest_date' column...")
    if 'earliest_date' not in df_filtered.columns:
        print(f"Error: 'earliest_date' column not found. Available columns: {df_filtered.columns.tolist()}")
        exit(1)

    df_filtered['earliest_date'] = pd.to_datetime(df_filtered['earliest_date'], errors='coerce')
    
    rows_before_nat_drop = len(df_filtered)
    df_filtered.dropna(subset=['earliest_date'], inplace=True)
    if len(df_filtered) < rows_before_nat_drop:
        print(f"Dropped {rows_before_nat_drop - len(df_filtered)} rows due to unparseable 'earliest_date' values.")

    if df_filtered.empty:
        print("No valid 'earliest_date' entries after initial parsing and NaT drop. Exiting.")
        exit(0)

    min_date_orig = df_filtered['earliest_date'].min()
    max_date_orig = df_filtered['earliest_date'].max()
    print(f"Original date range for filtered data: {min_date_orig.strftime('%Y-%m-%d')} to {max_date_orig.strftime('%Y-%m-%d')}")

    outliers_mask = pd.Series(False, index=df_filtered.index, dtype=bool) # Initialize to no outliers

    if args.outlier_method == 'iqr':
        print("Using IQR method for outlier detection (factor=1.5).")
        outliers_mask = detect_outliers_iqr(df_filtered['earliest_date'], factor=1.5)
    elif args.outlier_method == 'zscore':
        print("Using Z-score method for outlier detection (threshold=2.0).")
        outliers_mask = detect_outliers_zscore(df_filtered['earliest_date'], threshold=2.0)
    elif args.outlier_method == 'chaining':
        print(f"Using chaining method for outlier detection (max_gap_weeks={args.chaining_max_gap_weeks}).")
        outliers_mask = detect_outliers_chaining(df_filtered['earliest_date'], max_gap_weeks=args.chaining_max_gap_weeks)
    elif args.outlier_method == 'none':
        print("No outlier removal method selected.")

    if args.outlier_method != 'none' and outliers_mask.any():
        outliers_df = df_filtered[outliers_mask]
        print(f"Identified {len(outliers_df)} date outliers using {args.outlier_method} method.")
        
        # Report bounds for IQR method (based on distribution before these outliers were removed)
        if args.outlier_method == 'iqr':
            Q1_orig_dist = df_filtered['earliest_date'].quantile(0.25)
            Q3_orig_dist = df_filtered['earliest_date'].quantile(0.75)
            if pd.notna(Q1_orig_dist) and pd.notna(Q3_orig_dist): # Check if quantiles are valid
                IQR_orig_dist = Q3_orig_dist - Q1_orig_dist
                if IQR_orig_dist > pd.Timedelta(0): 
                    lower_bound_applied = Q1_orig_dist - 1.5 * IQR_orig_dist
                    upper_bound_applied = Q3_orig_dist + 1.5 * IQR_orig_dist
                    print(f"  (IQR method decision bounds: {lower_bound_applied.strftime('%Y-%m-%d')} to {upper_bound_applied.strftime('%Y-%m-%d')})")

        outlier_dates_str = ", ".join(sorted(outliers_df['earliest_date'].dt.strftime('%Y-%m-%d').unique()))
        print(f"Outlier dates being pruned: {outlier_dates_str}")
        df_filtered = df_filtered[~outliers_mask].copy() # Use .copy()
    elif args.outlier_method != 'none': # A method was chosen, but no outliers found
        print(f"No date outliers detected using {args.outlier_method} method.")

    if df_filtered.empty: # Check if all data was pruned
        print("All data removed after outlier pruning. Exiting.")
        exit(0)
    
    min_date_new = df_filtered['earliest_date'].min()
    max_date_new = df_filtered['earliest_date'].max()
    print(f"Date range after potential pruning: {min_date_new.strftime('%Y-%m-%d')} to {max_date_new.strftime('%Y-%m-%d')}")

    # --- 7. Identify Seed Samples and Save IDs ---
    print("Identifying seed strains from 'samples' column...")
    if 'samples' not in df_filtered.columns:
        print(f"Error: 'samples' column not found. Available columns: {df_filtered.columns.tolist()}")
        exit(1)

    # Extract strain_seed (first part of the first sample entry)
    try:
        df_filtered['strain_seed'] = df_filtered['samples'].str.split(',', n=1, expand=True)[0].str.split('|', n=1, expand=True)[0]
    except Exception as e:
        print(f"Error processing 'samples' column to extract strain_seed: {e}")
        print("Ensure the 'samples' column format is as expected (e.g., 'strain|genbank|date,...').")
        exit(1)
        
    unique_strain_seeds = df_filtered['strain_seed'].dropna().unique().tolist()

    if not unique_strain_seeds:
        print("No seed strains identified after processing. Exiting.")
        exit(0)
    
    print(f"Identified {len(unique_strain_seeds)} unique seed strains.")

    seed_ids_filename = f"{args.state.replace(' ', '_')}_{args.pango.replace('.', '_')}_seed_strains.txt"
    seed_ids_filepath = output_folder_path / seed_ids_filename
    
    with open(seed_ids_filepath, 'w') as f:
        for strain_id in unique_strain_seeds:
            f.write(f"{strain_id}\n")
    print(f"Seed strain IDs saved to: {seed_ids_filepath.resolve()}")

    # --- 8. Fetch Sequences from CovSpectrum ---
    fasta_sequences = fetch_sequences_from_covspectrum(unique_strain_seeds)

    if fasta_sequences is None:
        print("Failed to fetch sequences from CovSpectrum. Output FASTA file will not be created.")
    elif not fasta_sequences.strip():
        print("Warning: Fetched sequences from CovSpectrum are empty. Output FASTA file will be empty or not created.")
    else:
        fasta_filename = f"{args.state.replace(' ', '_')}_{args.pango.replace('.', '_')}_seed_sequences.fasta"
        fasta_filepath = output_folder_path / fasta_filename
        with open(fasta_filepath, 'w') as f:
            f.write(fasta_sequences)
        print(f"Fetched sequences saved to: {fasta_filepath.resolve()}")

    print("Script finished.")

if __name__ == "__main__":
    main()
