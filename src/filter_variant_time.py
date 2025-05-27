import argparse
import pandas as pd
import numpy as np
from scipy.stats import zscore
from typing import List, Dict, Optional
import us
import os
import subprocess
from pathlib import Path
import requests

# --- Functions imported/adapted from seed_seq_prep.py ---
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
    print("Please install it, e.g., using 'pip install pango-aliasor'.")
    print("The script cannot continue without this library.")
    def make_variant_base_map(base_variants: List[str], recombinant: bool = False) -> Dict[str, str]:
        raise ImportError("pango_aliasor not found, function unusable.")
    PANGO_ALIASOR_AVAILABLE = False
else:
    PANGO_ALIASOR_AVAILABLE = True

def detect_outliers_iqr(date_series: pd.Series, factor: float = 1.5) -> pd.Series:
    if date_series.empty:
        return pd.Series(dtype=bool, index=date_series.index)
    Q1 = date_series.quantile(0.25)
    Q3 = date_series.quantile(0.75)
    IQR = Q3 - Q1
    if IQR == pd.Timedelta(0):
        return pd.Series(False, index=date_series.index, dtype=bool)
    lower_bound = Q1 - factor * IQR
    upper_bound = Q3 + factor * IQR
    outliers_mask = (date_series < lower_bound) | (date_series > upper_bound)
    return outliers_mask

def detect_outliers_zscore(date_series: pd.Series, threshold: float = 2.0) -> pd.Series:
    if date_series.empty or date_series.nunique() < 2:
        return pd.Series(False, index=date_series.index, dtype=bool)
    numeric_dates = (date_series - date_series.min()).dt.days
    if numeric_dates.nunique() < 2:
         return pd.Series(False, index=date_series.index, dtype=bool)
    z_scores = zscore(numeric_dates)
    z_scores = np.nan_to_num(z_scores, nan=0.0)
    outliers_mask = np.abs(z_scores) > threshold
    return pd.Series(outliers_mask, index=date_series.index)

def detect_outliers_chaining(date_series: pd.Series, max_gap_weeks: int = 6) -> pd.Series:
    if date_series.empty:
        return pd.Series(dtype=bool, index=date_series.index)
    df_proc = date_series.to_frame(name='date').copy()
    df_sorted = df_proc.sort_values(by='date').copy()
    if df_sorted['date'].nunique() <= 1:
        return pd.Series(False, index=date_series.index, dtype=bool)
    max_delta = pd.Timedelta(weeks=max_gap_weeks)
    df_sorted['diff_prev'] = df_sorted['date'].diff()
    df_sorted['diff_next'] = df_sorted['date'].diff(-1).abs()
    is_far_from_prev = (df_sorted['diff_prev'] > max_delta) | df_sorted['diff_prev'].isna()
    is_far_from_next = (df_sorted['diff_next'] > max_delta) | df_sorted['diff_next'].isna()
    df_sorted['is_outlier'] = is_far_from_prev & is_far_from_next
    outlier_series_sorted_index = pd.Series(df_sorted['is_outlier'].values, index=df_sorted.index)
    return outlier_series_sorted_index.reindex(date_series.index)
# --- End of imported functions ---

DEFAULT_METADATA_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.metadata.tsv.gz"
DEFAULT_METADATA_FILENAME = "public-latest.metadata.tsv.gz"
DEFAULT_MAT_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.ShUShER.pb.gz"
DEFAULT_MAT_FILENAME = "public-latest.all.masked.ShUShER.pb.gz"

# Define the target schema based on the input_tsv example provided by the user
# 'region' is added because id_to_state creates it, and simulated_tsv might have it.
TARGET_METADATA_COLUMNS = [
    'strain', 'genbank_accession', 'date', 'country', 'host',
    'completeness', 'length', 'Nextstrain_clade', 'pangolin_lineage',
    'Nextstrain_clade_usher', 'pango_lineage_usher', 'region'
]


def download_file_if_needed(url: str, dest_folder: Path, filename: str) -> Optional[Path]:
    """Downloads a file if it doesn't exist in the destination folder."""
    dest_path = dest_folder / filename
    if dest_path.exists():
        print(f"Using existing file: {dest_path}")
        return dest_path
    
    print(f"Downloading {url} to {dest_path}...")
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        with open(dest_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print("Download complete.")
        return dest_path
    except requests.exceptions.RequestException as e:
        print(f"Error downloading {filename}: {e}")
        if dest_path.exists():
            dest_path.unlink(missing_ok=True)
        return None

def id_to_state(strain: str, country: str) -> str:
    """
    Parse the strain column to extract the state for samples from the USA.
    """
    if country == "USA":
        try:
            state_abbreviation = strain.split("/")[1].split("-")[0].upper()
            state = us.states.lookup(state_abbreviation)
            if state is not None:
                return state.name
        except IndexError:
            pass
    return "" # Return empty string or np.nan if not USA or not parsable


def process_tsv(input_file: Path,
                pango_lineage: str,
                outlier_method: str = 'chaining',
                chaining_max_gap_weeks: int = 6,
                iqr_factor: float = 1.5,
                zscore_threshold: float = 2.0) -> Optional[pd.DataFrame]:
    """
    Process the TSV file to filter rows based on pango lineage and date range (after outlier removal).
    """
    print(f"Processing TSV: {input_file} for Pango lineage: {pango_lineage}")
    try:
        # Assuming 'date' is the correct column name in public-latest.metadata.tsv.gz
        df = pd.read_csv(input_file, sep='\t', compression='gzip', low_memory=False)
        df['date'] = pd.to_datetime(df['date'], errors='coerce')
        df.dropna(subset=['date'], inplace=True)
    except Exception as e:
        print(f"Error reading or processing input TSV {input_file}: {e}")
        return None
    
    if 'pango_lineage_usher' not in df.columns:
        print(f"Error: 'pango_lineage_usher' column not found in {input_file}. This column is required for Pango regularization.")
        return None

    # Pango Lineage Regularization
    if not PANGO_ALIASOR_AVAILABLE: return None # Guard against missing library
    replace_map = make_variant_base_map([pango_lineage])
    df['regularized'] = df['pango_lineage_usher'].map(replace_map)
    
    lineage_df = df[df['regularized'] == pango_lineage].copy()
    if lineage_df.empty:
        print(f"No rows found for Pango lineage '{pango_lineage}' after regularization.")
        return None
    print(f"Found {len(lineage_df)} rows for Pango lineage '{pango_lineage}'.")

    # Outlier Detection
    print(f"Applying '{outlier_method}' outlier detection for dates...")
    valid_dates_for_outlier = lineage_df['date'].dropna()
    outliers_mask = pd.Series(False, index=valid_dates_for_outlier.index) # Default to no outliers

    if valid_dates_for_outlier.empty:
        print("No valid dates to perform outlier detection on.")
    elif outlier_method == 'iqr':
        outliers_mask = detect_outliers_iqr(valid_dates_for_outlier, factor=iqr_factor)
    elif outlier_method == 'zscore':
        outliers_mask = detect_outliers_zscore(valid_dates_for_outlier, threshold=zscore_threshold)
    elif outlier_method == 'chaining':
        outliers_mask = detect_outliers_chaining(valid_dates_for_outlier, max_gap_weeks=chaining_max_gap_weeks)
    
    if outliers_mask.any():
        num_outliers = outliers_mask.sum()
        print(f"Identified {num_outliers} date outliers using {outlier_method} method. Pruning them.")
        # Apply mask to lineage_df by aligning indices
        lineage_df = lineage_df.loc[valid_dates_for_outlier[~outliers_mask].index]
    elif outlier_method != 'none':
        print(f"No date outliers detected using {outlier_method} method.")

    if lineage_df.empty:
        print("No data remaining after outlier removal.")
        return None

    earliest_date = lineage_df['date'].min()
    latest_date = lineage_df['date'].max()
    
    if pd.isna(earliest_date) or pd.isna(latest_date):
        print("No valid earliest/latest dates found after outlier processing.")
        return None
    
    print(f"Effective date range for '{pango_lineage}': {earliest_date.strftime('%Y-%m-%d')} to {latest_date.strftime('%Y-%m-%d')}")
    
    # Filter rows within two weeks of the established earliest and latest dates
    two_weeks = pd.Timedelta(weeks=2)
    time_window_filtered_df = lineage_df[
        (lineage_df['date'] >= earliest_date - two_weeks) &
        (lineage_df['date'] <= latest_date + two_weeks)
    ].copy()
    
    print(f"Filtered to {len(time_window_filtered_df)} rows within +/- 2 weeks of effective date range.")
    return time_window_filtered_df


def create_reductions(df: pd.DataFrame,
                      output_path: Path,
                      base_name_prefix: str,
                      num_reductions: int,
                      target_state: str,
                      input_mat_path: Path,
                      ablate: bool):
    """
    Create reductions of the table, save metadata TSVs, and run matUtils extract.
    """
    # Ensure 'region' column exists for target_state filtering
    # df['region'] is created outside if not present, using id_to_state
    if 'region' not in df.columns:
        print("Warning: 'region' column not present in DataFrame. Creating it for USA samples.")
        df['region'] = df.apply(lambda x: id_to_state(x['strain'], x['country']), axis=1)


    target_mask = (df['region'] == target_state) & (df['regularized'] == df['regularized'].iloc[0] if not df.empty else True) # Ensure we only remove target pango
    
    try:
        abbr_state = us.states.lookup(target_state).abbr
    except:
        print(f"Warning: Could not find abbreviation for state '{target_state}'.Using full name for filenames.")
        abbr_state = target_state.replace(" ", "_")

    print(f"Original number of rows for target Pango lineage in {target_state}: {target_mask.sum()}")

    reduction_iterations = []
    # Iteration 0: "complete" reduction (all target state samples removed)
    reduction_iterations.append({"factor_val": "complete", "i_val": 0})

    if ablate:
        for i in range(1, num_reductions + 1): # 2-fold, 4-fold etc.
            reduction_iterations.append({"factor_val": 2**i, "i_val": i})
            
    for item in reduction_iterations:
        reduction_factor = item["factor_val"]
        i = item["i_val"]

        current_df_for_reduction = df.copy()
        
        if reduction_factor == "complete":
            # Remove all target state samples for the target Pango lineage
            reduced_df = current_df_for_reduction[~target_mask]
            label = "complete"
            num_removed = target_mask.sum()
        else: # Numerical reduction factor (2**i)
            target_state_samples = current_df_for_reduction[target_mask]
            other_samples = current_df_for_reduction[~target_mask]
            
            reduced_target_state_samples = target_state_samples.iloc[::int(reduction_factor), :]
            reduced_df = pd.concat([reduced_target_state_samples, other_samples])
            label = str(reduction_factor)
            num_removed = len(target_state_samples) - len(reduced_target_state_samples)

        print(f"Reduction '{label}': {num_removed} rows of target state/lineage removed. Resulting in {len(reduced_df)} total rows.")

        # Ensure all TARGET_METADATA_COLUMNS exist, fill with NaN if not
        for col in TARGET_METADATA_COLUMNS:
            if col not in reduced_df.columns:
                reduced_df[col] = np.nan
        final_output_df = reduced_df[TARGET_METADATA_COLUMNS].sort_values(by='date')


        # --- Save Metadata TSV ---
        output_filename_base = f"{base_name_prefix}_{abbr_state}_reduction_{label}"
        metadata_output_file = output_path / f"{output_filename_base}.tsv"
        final_output_df.to_csv(metadata_output_file, sep='\t', index=False)
        print(f"Metadata for reduction '{label}' written to {metadata_output_file}")

        # --- Prepare for and Run matUtils extract ---
        ids_to_keep = final_output_df['strain'].dropna().unique().tolist()
        if not ids_to_keep:
            print(f"No strain IDs to keep for reduction '{label}'. Skipping matUtils extract.")
            continue

        ids_file = output_path / f"{output_filename_base}_ids.txt"
        with open(ids_file, 'w') as f_ids:
            f_ids.write('\n'.join(ids_to_keep))
        
        output_mat_name = f"{output_filename_base}.pb"
        
        matutils_cmd = [
            "matUtils", "extract",
            "--input-mat", str(input_mat_path),
            "--output-directory", str(output_path.resolve()),
            "--write-mat", output_mat_name,
            "--samples", str(ids_file.resolve())
        ]
        
        print(f"Running matUtils extract for reduction '{label}':")
        print(f"  Command: {' '.join(matutils_cmd)}")
        
        try:
            result = subprocess.run(matutils_cmd, capture_output=True, text=True, check=True)
            print(f"  matUtils STDOUT:\n{result.stdout}")
            if result.stderr:
                print(f"  matUtils STDERR:\n{result.stderr}")
            print(f"  matUtils extract completed. Output tree: {output_path / output_mat_name}")
        except subprocess.CalledProcessError as e:
            print(f"  Error running matUtils extract for reduction '{label}':")
            print(f"  Return code: {e.returncode}")
            print(f"  STDOUT:\n{e.stdout}")
            print(f"  STDERR:\n{e.stderr}")
        except FileNotFoundError:
            print("Error: matUtils command not found. Please ensure it is installed and in your PATH.")
            # Optionally, decide if the script should terminate or continue without matUtils
            break # Stop further matUtils attempts if it's not found

def align_simulated_df(sim_df: pd.DataFrame) -> pd.DataFrame:
    """Aligns simulated DataFrame to the TARGET_METADATA_COLUMNS schema."""
    aligned_df = pd.DataFrame()
    for col in TARGET_METADATA_COLUMNS:
        if col in sim_df.columns:
            aligned_df[col] = sim_df[col]
        else:
            # Handle specific known mappings if any, or default to NaN
            if col == 'pango_lineage_usher': # Example: if sim has 'pangolin_lineage', maybe use that
                aligned_df[col] = sim_df['pangolin_lineage'] if 'pangolin_lineage' in sim_df.columns else np.nan
            else:
                aligned_df[col] = np.nan
    
    # Ensure date is datetime
    if 'date' in aligned_df.columns:
        aligned_df['date'] = pd.to_datetime(aligned_df['date'], errors='coerce')
    return aligned_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter a metadata TSV based on pango lineage and date ranges, "
                                                 "perform reductions, and extract corresponding subtrees using matUtils.")
    parser.add_argument("output_path", type=str, help="Path to the output folder for all generated files.")
    parser.add_argument("pango_lineage", help="Pango lineage to filter by (e.g., BA.1, B.1.617.2).")
    parser.add_argument("target_state", help="Target US state for reduction (e.g., Virginia, California).")
    
    parser.add_argument("--input_file", type=str, default=None,
                        help="Optional path to the input metadata TSV file (e.g., public-latest.metadata.tsv.gz). "
                             "If not provided, attempts to download to output_path.")
    parser.add_argument("--input_mat", type=str, default=None,
                        help="Optional path to the input MAT protobuf file (e.g., public-latest.all.masked.ShUShER.pb.gz). "
                             "If not provided, attempts to download to output_path.")
    parser.add_argument("--simulated_tsv", type=str, default=None,
                        help="Optional path to a TSV file with simulated data to concatenate.")
    
    parser.add_argument("--reductions", type=int, default=3,
                        help="Number of 2-fold reduction steps if --ablate is used (e.g., 3 means 2-fold, 4-fold, 8-fold). Default: 3.")
    parser.add_argument("--ablate", action='store_true',
                        help="Perform N-fold ablations of target state samples in addition to complete removal.")

    parser.add_argument("--outlier_method", type=str, choices=['none', 'iqr', 'zscore', 'chaining'], default='chaining',
                        help="Method for 'earliest_date' outlier detection (default: chaining).")
    parser.add_argument("--chaining_max_gap_weeks", type=int, default=6,
                        help="Max gap in weeks for 'chaining' outlier method (default: 6).")
    parser.add_argument("--iqr_factor", type=float, default=1.5, help="IQR factor for 'iqr' outlier method.")
    parser.add_argument("--zscore_threshold", type=float, default=2.0, help="Z-score threshold for 'zscore' outlier method.")

    args = parser.parse_args()

    if not PANGO_ALIASOR_AVAILABLE:
        print("Exiting due to missing 'pango_aliasor' library.")
        exit(1)

    output_folder = Path(args.output_path)
    output_folder.mkdir(parents=True, exist_ok=True)

    # --- Handle Input Metadata File ---
    metadata_file_path: Optional[Path]
    if args.input_file:
        metadata_file_path = Path(args.input_file)
        if not metadata_file_path.exists():
            print(f"Error: Provided input_file '{metadata_file_path}' does not exist.")
            exit(1)
    else:
        metadata_file_path = download_file_if_needed(DEFAULT_METADATA_URL, output_folder, DEFAULT_METADATA_FILENAME)
    
    if not metadata_file_path or not metadata_file_path.exists():
        print("Error: Could not obtain input metadata file. Exiting.")
        exit(1)

    # --- Handle Input MAT File ---
    input_mat_path: Optional[Path]
    if args.input_mat:
        input_mat_path = Path(args.input_mat)
        if not input_mat_path.exists():
            print(f"Error: Provided input_mat '{input_mat_path}' does not exist.")
            exit(1)
    else:
        input_mat_path = download_file_if_needed(DEFAULT_MAT_URL, output_folder, DEFAULT_MAT_FILENAME)

    if not input_mat_path or not input_mat_path.exists():
        print("Error: Could not obtain input MAT file. Exiting.")
        exit(1)

    # --- Process main metadata ---
    filtered_df = process_tsv(metadata_file_path,
                              args.pango_lineage,
                              args.outlier_method,
                              args.chaining_max_gap_weeks,
                              args.iqr_factor,
                              args.zscore_threshold)

    if filtered_df is None or filtered_df.empty:
        print("No data after initial processing and filtering. Exiting.")
        exit(0)

    # Add 'region' column using id_to_state for main df
    print("Mapping strain IDs to states for USA entries...")
    filtered_df['region'] = filtered_df.apply(lambda x: id_to_state(x['strain'], x['country']), axis=1)
    
    # --- Handle Simulated TSV ---
    if args.simulated_tsv:
        sim_tsv_path = Path(args.simulated_tsv)
        if sim_tsv_path.exists():
            print(f"Loading simulated data from: {sim_tsv_path}")
            try:
                sim_df = pd.read_csv(sim_tsv_path, sep='\t', low_memory=False)
                aligned_sim_df = align_sim_ulated_df(sim_df)
                
                # Before concat, ensure consistent dtypes for common problematic columns like 'date'
                if 'date' in filtered_df.columns:
                     filtered_df['date'] = pd.to_datetime(filtered_df['date'], errors='coerce')
                if 'date' in aligned_sim_df.columns:
                    aligned_sim_df['date'] = pd.to_datetime(aligned_sim_df['date'], errors='coerce')
                
                # Also ensure 'regularized' column for simulated data if needed for downstream logic
                # For now, assume simulated data doesn't need Pango regularization for filtering purposes here
                # or it's already implicitly handled by not being part of target_mask for removal.
                # If it *should* be regularized for some reason, that logic needs to be added.
                if 'regularized' not in aligned_sim_df.columns and 'regularized' in filtered_df.columns:
                    aligned_sim_df['regularized'] = np.nan # Or derive if possible


                filtered_df = pd.concat([filtered_df, aligned_sim_df], ignore_index=True)
                filtered_df.sort_values(by='date', inplace=True) # Sort after concat
                print(f"Concatenated simulated data. Total rows now: {len(filtered_df)}")
            except Exception as e:
                print(f"Error processing simulated TSV {sim_tsv_path}: {e}")
        else:
            print(f"Warning: Simulated TSV file not found: {sim_tsv_path}")

    # --- Create Reductions ---
    if args.target_state:
        # Sanitize pango_lineage and target_state for use in filenames
        sanitized_pango = args.pango_lineage.replace("/", "_").replace(".", "")
        sanitized_state = args.target_state.replace(" ", "_")
        base_name_prefix = f"{sanitized_state}_{sanitized_pango}"
        
        create_reductions(filtered_df,
                          output_folder,
                          base_name_prefix,
                          args.reductions,
                          args.target_state,
                          input_mat_path,
                          args.ablate)
    else:
        print("No target_state provided. Skipping reductions and matUtils steps.")
        # Optionally, save the globally filtered_df and extract one MAT for it if desired.
        # For now, follows the old script's conditional execution.

    print("Script finished.")