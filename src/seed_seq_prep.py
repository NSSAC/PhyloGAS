import argparse
import pandas as pd
import numpy as np
import requests
import os
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional

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
# Comprehensive list from the notebook for robust aliasing
PANGO_TARGET_LIST_FOR_ALIASING = [
    'B.1.1.7', 'B.1.617.2', 'BA.1', 'BA.2', 'BA.4', 'BA.5', 'XBB', 'XBB.1.5', 
    'XBB.1.16', 'XBB.1.9', 'BA.2.86', 'KP.3', 'LP.8', 'XEC'
]

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
            # Depending on policy, could stop all or just skip this batch.
            # For now, let's try to continue with other batches if possible, but signal overall failure.
            return None # Indicate overall failure
        except requests.exceptions.RequestException as e:
            print(f"  Request Error for batch {i//batch_size + 1}: {e}")
            return None # Indicate overall failure
        except Exception as e:
            print(f"  An unexpected error occurred during API call for batch {i//batch_size + 1}: {e}")
            return None

    print("Sequence fetching complete.")
    return "".join(all_fasta_content)

def main():
    parser = argparse.ArgumentParser(description="Process SARS-CoV-2 cluster data, identify seed samples, and fetch their sequences.")
    parser.add_argument("--state", required=True, type=str, help="US state to filter data for (e.g., 'Washington', 'California'). Corresponds to 'region' column.")
    parser.add_argument("--pango", required=True, type=str, help="Pango lineage to filter data for (e.g., 'B.1.1.7', 'BA.1'). This will be matched against regularized lineage names.")
    parser.add_argument("--input_file", type=str, help="Optional path to the input TSV file. If not provided, the script will attempt to use/download 'hardcoded_clusters.tsv'.")
    parser.add_argument("--output_folder", required=True, type=str, help="Path to the folder where output files (seed IDs, FASTA sequences, downloaded TSV) will be saved.")
    
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

    # Ensure the user-provided pango lineage is part of the aliasing process
    # to handle cases where it might be a direct base variant itself.
    #extended_target_list = list(set(PANGO_TARGET_LIST_FOR_ALIASING + [args.pango]))
    extended_target_list = [args.pango]
    
    try:
        replace_map = make_variant_base_map(extended_target_list)
    except Exception as e:
        print(f"Error during pango_aliasor processing: {e}")
        print("This might be due to an issue with the pango_aliasor library or its data.")
        exit(1)

    df['pango_regularized'] = df['annotation_2'].map(replace_map)
    # For lineages not in replace_map (e.g. not covered by target_list or its sublineages),
    # pango_regularized will be NaN. We can fill them with original annotation_2 if needed,
    # or rely on filtering to remove them. Let's fill with original.
    df['pango_regularized'].fillna(df['annotation_2'], inplace=True)
    print(f"Pango lineage regularization complete. New column 'pango_regularized' created.")

    # --- 5. Filter data by state and pango lineage ---
    print(f"Filtering data for state: '{args.state}' and Pango lineage (regularized): '{args.pango}'...")
    if 'region' not in df.columns:
        print(f"Error: 'region' column not found in the input TSV. Needed for state filtering. Available columns: {df.columns.tolist()}")
        exit(1)
        
    original_row_count = len(df)
    df_filtered = df[(df['region'] == args.state) & (df['pango_regularized'] == args.pango)].copy() # Use .copy() to avoid SettingWithCopyWarning

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
    
    # Drop rows where date conversion failed
    rows_before_nat_drop = len(df_filtered)
    df_filtered.dropna(subset=['earliest_date'], inplace=True)
    if len(df_filtered) < rows_before_nat_drop:
        print(f"Dropped {rows_before_nat_drop - len(df_filtered)} rows due to unparseable 'earliest_date' values.")

    if df_filtered.empty:
        print("No valid 'earliest_date' entries after initial parsing and NaT drop. Exiting.")
        exit(0)

    # Outlier detection using IQR
    min_date_orig = df_filtered['earliest_date'].min()
    max_date_orig = df_filtered['earliest_date'].max()
    print(f"Original date range for filtered data: {min_date_orig.strftime('%Y-%m-%d')} to {max_date_orig.strftime('%Y-%m-%d')}")

    Q1 = df_filtered['earliest_date'].quantile(0.25)
    Q3 = df_filtered['earliest_date'].quantile(0.75)
    IQR = Q3 - Q1
    
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    outliers_mask = (df_filtered['earliest_date'] < lower_bound) | (df_filtered['earliest_date'] > upper_bound)
    outliers_df = df_filtered[outliers_mask]

    if not outliers_df.empty:
        print(f"Identified {len(outliers_df)} date outliers outside range [{lower_bound.strftime('%Y-%m-%d')}, {upper_bound.strftime('%Y-%m-%d')}].")
        outlier_dates_str = ", ".join(sorted(outliers_df['earliest_date'].dt.strftime('%Y-%m-%d').unique()))
        print(f"Outlier dates being pruned: {outlier_dates_str}")
        df_filtered = df_filtered[~outliers_mask].copy() # Use .copy()
    else:
        print("No date outliers detected based on IQR method.")

    if df_filtered.empty:
        print("All data removed after outlier pruning. Exiting.")
        exit(0)
    
    min_date_new = df_filtered['earliest_date'].min()
    max_date_new = df_filtered['earliest_date'].max()
    print(f"Date range after pruning: {min_date_new.strftime('%Y-%m-%d')} to {max_date_new.strftime('%Y-%m-%d')}")

    # --- 7. Identify Seed Samples and Save IDs ---
    print("Identifying seed strains from 'samples' column...")
    if 'samples' not in df_filtered.columns:
        print(f"Error: 'samples' column not found. Needed to extract seed strains. Available columns: {df_filtered.columns.tolist()}")
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
