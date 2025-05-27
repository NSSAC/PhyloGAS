import argparse
import pandas as pd
import numpy as np
# from scipy.stats import zscore # zscore is used by imported outlier functions
import us
import os
import subprocess
from pathlib import Path
import requests

# --- Attempt to import functions from seed_seq_prep.py ---
# Define placeholders for functions and flags
make_variant_base_map = None
detect_outliers_iqr = None
detect_outliers_zscore = None
detect_outliers_chaining = None
SSP_PANGO_ALIASOR_AVAILABLE = False  # Flag for pango_aliasor availability via seed_seq_prep
SEED_SEQ_PREP_IMPORTS_SUCCESSFUL = False

try:
    from seed_seq_prep import (
        make_variant_base_map as ssp_make_variant_base_map_imported,
        detect_outliers_iqr as ssp_detect_outliers_iqr_imported,
        detect_outliers_zscore as ssp_detect_outliers_zscore_imported,
        detect_outliers_chaining as ssp_detect_outliers_chaining_imported,
        PANGO_ALIASOR_AVAILABLE as ssp_pango_available_flag
    )
    # Assign imported functions to module-level names
    make_variant_base_map = ssp_make_variant_base_map_imported
    detect_outliers_iqr = ssp_detect_outliers_iqr_imported
    detect_outliers_zscore = ssp_detect_outliers_zscore_imported
    detect_outliers_chaining = ssp_detect_outliers_chaining_imported
    SSP_PANGO_ALIASOR_AVAILABLE = ssp_pango_available_flag
    SEED_SEQ_PREP_IMPORTS_SUCCESSFUL = True
    print("Successfully imported helper functions from seed_seq_prep.py")
except ImportError as e:
    print(f"ERROR: Could not import required components from seed_seq_prep.py: {e}")
    print("Ensure seed_seq_prep.py is in the same directory or accessible via PYTHONPATH,")
    print("and that its dependencies (like pango_aliasor) are installed.")
    # Functions remain None, SSP_PANGO_ALIASOR_AVAILABLE remains False

# --- Constants ---
DEFAULT_METADATA_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.metadata.tsv.gz"
DEFAULT_METADATA_FILENAME = "public-latest.metadata.tsv.gz"
DEFAULT_MAT_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.ShUShER.pb.gz"
DEFAULT_MAT_FILENAME = "public-latest.all.masked.ShUShER.pb.gz"

TARGET_METADATA_COLUMNS = [
    'strain', 'genbank_accession', 'date', 'country', 'host',
    'completeness', 'length', 'Nextstrain_clade', 'pangolin_lineage',
    'Nextstrain_clade_usher', 'pango_lineage_usher', 'region'
]

def download_file_if_needed(url: str, dest_folder: Path, filename: str) -> Optional[Path]:
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
    if country == "USA":
        try:
            parts = strain.split("/")
            if len(parts) > 1:
                state_abbreviation = parts[1].split("-")[0].upper()
                state = us.states.lookup(state_abbreviation)
                if state is not None:
                    return state.name
        except IndexError:
            pass # Strain format might be unexpected
        except Exception: # Catch other potential errors from us.states.lookup
            pass
    return ""


def process_tsv(input_file: Path,
                pango_lineage: str,
                outlier_method: str = 'chaining',
                chaining_max_gap_weeks: int = 6,
                iqr_factor: float = 1.5,
                zscore_threshold: float = 2.0) -> Optional[pd.DataFrame]:
    print(f"Processing TSV: {input_file} for Pango lineage: {pango_lineage}")
    try:
        df = pd.read_csv(input_file, sep='\t', compression='gzip', low_memory=False)
        df['date'] = pd.to_datetime(df['date'], errors='coerce')
        df.dropna(subset=['date'], inplace=True)
    except Exception as e:
        print(f"Error reading or processing input TSV {input_file}: {e}")
        return None
    
    if 'pango_lineage_usher' not in df.columns:
        print(f"Error: 'pango_lineage_usher' column not found in {input_file}. This column is required for Pango regularization.")
        return None

    # Pango Lineage Regularization (uses imported function)
    try:
        replace_map = make_variant_base_map([pango_lineage]) # type: ignore
    except Exception as e: # Catch potential errors from make_variant_base_map if pango_aliasor fails
        print(f"Error during Pango lineage regularization: {e}")
        return None

    df['regularized'] = df['pango_lineage_usher'].map(replace_map)
    
    lineage_df = df[df['regularized'] == pango_lineage].copy()
    if lineage_df.empty:
        print(f"No rows found for Pango lineage '{pango_lineage}' after regularization.")
        return None
    print(f"Found {len(lineage_df)} rows for Pango lineage '{pango_lineage}'.")

    print(f"Applying '{outlier_method}' outlier detection for dates...")
    valid_dates_for_outlier = lineage_df['date'].dropna()
    outliers_mask = pd.Series(False, index=valid_dates_for_outlier.index) 

    if valid_dates_for_outlier.empty:
        print("No valid dates to perform outlier detection on.")
    elif outlier_method == 'iqr':
        outliers_mask = detect_outliers_iqr(valid_dates_for_outlier, factor=iqr_factor) # type: ignore
    elif outlier_method == 'zscore':
        outliers_mask = detect_outliers_zscore(valid_dates_for_outlier, threshold=zscore_threshold) # type: ignore
    elif outlier_method == 'chaining':
        outliers_mask = detect_outliers_chaining(valid_dates_for_outlier, max_gap_weeks=chaining_max_gap_weeks) # type: ignore
    
    if outliers_mask.any():
        num_outliers = outliers_mask.sum()
        print(f"Identified {num_outliers} date outliers using {outlier_method} method. Pruning them.")
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
                      num_reduction_steps: int, # Renamed from num_reductions for clarity
                      target_state: str,
                      input_mat_path: Path,
                      ablate: bool) -> Optional[Path]:
    """
    Create reductions of the table, save metadata TSVs, and run matUtils extract.
    Returns the path to the "complete" reduction .pb file if created.
    """
    complete_pb_file_path = None
    if 'region' not in df.columns:
        print("Warning: 'region' column not present. Creating it for USA samples.")
        df['region'] = df.apply(lambda x: id_to_state(x['strain'], x['country']), axis=1)

    target_pango = df['regularized'].iloc[0] if not df.empty and 'regularized' in df.columns and not df['regularized'].empty else None
    if target_pango is None:
        print("Warning: Could not determine target Pango lineage from 'regularized' column for precise masking.")
        target_mask = (df['region'] == target_state) # Less precise mask
    else:
        target_mask = (df['region'] == target_state) & (df['regularized'] == target_pango)
    
    try:
        abbr_state = us.states.lookup(target_state).abbr
    except Exception:
        print(f"Warning: Could not find abbreviation for state '{target_state}'. Using full name for filenames.")
        abbr_state = target_state.replace(" ", "_")

    print(f"Original number of rows for target Pango lineage in {target_state}: {target_mask.sum()}")

    reduction_iterations = []
    reduction_iterations.append({"factor_val": "complete", "label_suffix": "complete"})

    if ablate:
        for i in range(1, num_reduction_steps + 1):
            reduction_iterations.append({"factor_val": 2**i, "label_suffix": str(2**i)})
            
    for item in reduction_iterations:
        reduction_factor_val = item["factor_val"]
        label_suffix = item["label_suffix"]
        current_df_for_reduction = df.copy()
        
        if reduction_factor_val == "complete":
            reduced_df = current_df_for_reduction[~target_mask]
            num_removed = target_mask.sum()
        else: # Numerical reduction factor
            target_state_samples = current_df_for_reduction[target_mask]
            other_samples = current_df_for_reduction[~target_mask]
            
            reduced_target_state_samples = target_state_samples.iloc[::int(reduction_factor_val)]
            reduced_df = pd.concat([reduced_target_state_samples, other_samples])
            num_removed = len(target_state_samples) - len(reduced_target_state_samples)

        print(f"Reduction '{label_suffix}': {num_removed} rows of target state/lineage removed. Resulting in {len(reduced_df)} total rows.")

        for col in TARGET_METADATA_COLUMNS:
            if col not in reduced_df.columns:
                reduced_df[col] = np.nan
        final_output_df = reduced_df[TARGET_METADATA_COLUMNS].sort_values(by='date')

        output_filename_base = f"{base_name_prefix}_{abbr_state}_reduction_{label_suffix}"
        metadata_output_file = output_path / f"{output_filename_base}.tsv"
        final_output_df.to_csv(metadata_output_file, sep='\t', index=False)
        print(f"Metadata for reduction '{label_suffix}' written to {metadata_output_file}")

        ids_to_keep = final_output_df['strain'].dropna().unique().tolist()
        if not ids_to_keep:
            print(f"No strain IDs to keep for reduction '{label_suffix}'. Skipping matUtils extract.")
            continue

        ids_file = output_path / f"{output_filename_base}_ids.txt"
        with open(ids_file, 'w') as f_ids:
            f_ids.write('\n'.join(ids_to_keep))
        
        output_mat_name = f"{output_filename_base}.pb"
        current_pb_path = output_path / output_mat_name
        
        matutils_cmd = [
            "matUtils", "extract",
            "--input-mat", str(input_mat_path.resolve()),
            "--output-directory", str(output_path.resolve()),
            "--write-mat", output_mat_name, # This is just the basename
            "--samples", str(ids_file.resolve())
        ]
        
        print(f"Running matUtils extract for reduction '{label_suffix}':")
        print(f"  Command: {' '.join(matutils_cmd)}")
        
        try:
            result = subprocess.run(matutils_cmd, capture_output=True, text=True, check=True)
            print(f"  matUtils STDOUT:\n{result.stdout}")
            if result.stderr: print(f"  matUtils STDERR:\n{result.stderr}")
            print(f"  matUtils extract completed. Output tree: {current_pb_path}")
            if reduction_factor_val == "complete":
                complete_pb_file_path = current_pb_path
        except subprocess.CalledProcessError as e:
            print(f"  Error running matUtils extract for reduction '{label_suffix}':")
            print(f"  Return code: {e.returncode}\n  STDOUT:\n{e.stdout}\n  STDERR:\n{e.stderr}")
        except FileNotFoundError:
            print("Error: matUtils command not found. Please ensure it is installed and in your PATH.")
            return None # Stop if matUtils is not found
    return complete_pb_file_path


def align_simulated_df(sim_df: pd.DataFrame) -> pd.DataFrame:
    aligned_df = pd.DataFrame()
    # Map simulated_tsv columns to target_tsv columns
    # This is a basic mapping, can be made more sophisticated
    sim_to_target_map = {
        'strain': 'strain',
        'date': 'date',
        'country': 'country',
        'region': 'region', # if 'division' is state name in sim_df
        'pangolin_lineage': 'pangolin_lineage', # Will also be used for pango_lineage_usher if not present
        # Add other direct mappings if they exist
    }
    
    # First, create columns based on direct mapping
    for sim_col, target_col in sim_to_target_map.items():
        if sim_col in sim_df.columns:
            aligned_df[target_col] = sim_df[sim_col]

    # Ensure all TARGET_METADATA_COLUMNS exist in aligned_df
    for col in TARGET_METADATA_COLUMNS:
        if col not in aligned_df.columns:
            # Special handling for pango_lineage_usher if not directly mapped
            if col == 'pango_lineage_usher' and 'pangolin_lineage' in aligned_df.columns:
                aligned_df[col] = aligned_df['pangolin_lineage']
            else:
                aligned_df[col] = np.nan # Default to NaN
    
    if 'date' in aligned_df.columns:
        aligned_df['date'] = pd.to_datetime(aligned_df['date'], errors='coerce')
    return aligned_df[TARGET_METADATA_COLUMNS] # Ensure correct order and selection


def run_simulation_placement(sim_fasta_path: Path, 
                             complete_pb_path: Path, 
                             output_folder: Path,
                             base_name_prefix_for_sim: str,
                             abbr_state_for_sim: str):
    print("\n--- Running simulation placement steps ---")
    
    # Command 1: faToVcf
    sim_fasta_stem = sim_fasta_path.stem
    output_vcf_filename = f"{sim_fasta_stem}.vcf"
    output_vcf_path = output_folder / output_vcf_filename
    
    fatovcf_cmd = ["faToVcf", str(sim_fasta_path.resolve()), str(output_vcf_path.resolve())]
    print(f"Running faToVcf:\n  Command: {' '.join(fatovcf_cmd)}")
    try:
        result_fatovcf = subprocess.run(fatovcf_cmd, capture_output=True, text=True, check=True)
        print(f"  faToVcf STDOUT:\n{result_fatovcf.stdout}")
        if result_fatovcf.stderr: print(f"  faToVcf STDERR:\n{result_fatovcf.stderr}")
        print(f"  faToVcf completed. Output VCF: {output_vcf_path}")
    except subprocess.CalledProcessError as e:
        print(f"  Error running faToVcf:\n  Return code: {e.returncode}\n  STDOUT:\n{e.stdout}\n  STDERR:\n{e.stderr}")
        return
    except FileNotFoundError:
        print("Error: faToVcf command not found. Please ensure it is installed and in your PATH.")
        return

    # Command 2: usher
    # Construct a descriptive output name for the usher output
    usher_output_pb_filename = f"{base_name_prefix_for_sim}_{abbr_state_for_sim}_simulated_{sim_fasta_stem}.pb"
    usher_output_pb_path = output_folder / usher_output_pb_filename
    
    usher_cmd = [
        "usher",
        "-i", str(complete_pb_path.resolve()),
        "-v", str(output_vcf_path.resolve()),
        "-o", str(usher_output_pb_path.resolve()),
        "-T", "30" # Number of threads
    ]
    print(f"Running usher for simulation placement:\n  Command: {' '.join(usher_cmd)}")
    try:
        result_usher = subprocess.run(usher_cmd, capture_output=True, text=True, check=True)
        print(f"  usher STDOUT:\n{result_usher.stdout}")
        if result_usher.stderr: print(f"  usher STDERR:\n{result_usher.stderr}")
        print(f"  usher completed. Output tree with simulated sequences: {usher_output_pb_path}")
    except subprocess.CalledProcessError as e:
        print(f"  Error running usher:\n  Return code: {e.returncode}\n  STDOUT:\n{e.stdout}\n  STDERR:\n{e.stderr}")
    except FileNotFoundError:
        print("Error: usher command not found. Please ensure it is installed and in your PATH.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter metadata, perform reductions, and optionally place simulated sequences.")
    parser.add_argument("output_path", type=str, help="Path to the output folder for all generated files.")
    parser.add_argument("pango_lineage", help="Pango lineage to filter by.")
    parser.add_argument("target_state", help="Target US state for reduction.")
    
    parser.add_argument("--input_file", type=str, default=None, help="Optional input metadata TSV (gzipped). Downloads if not given.")
    parser.add_argument("--input_mat", type=str, default=None, help="Optional input MAT protobuf (gzipped). Downloads if not given.")
    parser.add_argument("--simulated_tsv", type=str, default=None, help="Optional TSV with simulated data to concatenate.")
    parser.add_argument("--sim_fasta", type=str, default=None, help="Optional FASTA file of simulated sequences for placement.")
    
    parser.add_argument("--reductions", type=int, default=3, dest="num_reduction_steps", help="Number of N-fold reduction steps if --ablate is used (default: 3).")
    parser.add_argument("--ablate", action='store_true', help="Perform N-fold ablations in addition to complete removal.")

    parser.add_argument("--outlier_method", type=str, choices=['none', 'iqr', 'zscore', 'chaining'], default='chaining', help="Date outlier detection method (default: chaining).")
    parser.add_argument("--chaining_max_gap_weeks", type=int, default=6, help="Max gap for 'chaining' outlier method (default: 6 weeks).")
    parser.add_argument("--iqr_factor", type=float, default=1.5, help="IQR factor for 'iqr' outlier method.")
    parser.add_argument("--zscore_threshold", type=float, default=2.0, help="Z-score threshold for 'zscore' outlier method.")

    args = parser.parse_args()

    if not SEED_SEQ_PREP_IMPORTS_SUCCESSFUL:
        print("Exiting due to failure in importing critical components from seed_seq_prep.py.")
        exit(1)
    if not SSP_PANGO_ALIASOR_AVAILABLE: # Check the flag imported from seed_seq_prep
        print("Pango_aliasor library is not available (as per seed_seq_prep.py). Regularization will fail. Exiting.")
        exit(1)

    output_folder = Path(args.output_path)
    output_folder.mkdir(parents=True, exist_ok=True)

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

    filtered_df = process_tsv(metadata_file_path,
                              args.pango_lineage,
                              args.outlier_method,
                              args.chaining_max_gap_weeks,
                              args.iqr_factor,
                              args.zscore_threshold)

    if filtered_df is None or filtered_df.empty:
        print("No data after initial processing and filtering. Exiting.")
        exit(0)

    print("Mapping strain IDs to states for USA entries...")
    filtered_df['region'] = filtered_df.apply(lambda x: id_to_state(x['strain'], x['country'] if 'country' in x and pd.notna(x['country']) else ""), axis=1)
    
    if args.simulated_tsv:
        sim_tsv_path = Path(args.simulated_tsv)
        if sim_tsv_path.exists():
            print(f"Loading simulated data from: {sim_tsv_path}")
            try:
                sim_df = pd.read_csv(sim_tsv_path, sep='\t', low_memory=False)
                # Add 'region' to sim_df if 'division' exists and country is USA (example)
                if 'division' in sim_df.columns and 'country' in sim_df.columns:
                     sim_df['region'] = sim_df.apply(lambda x: x['division'] if x['country'] == 'USA' else '', axis=1)

                aligned_sim_df = align_simulated_df(sim_df)
                
                if 'date' in filtered_df.columns: filtered_df['date'] = pd.to_datetime(filtered_df['date'], errors='coerce')
                if 'date' in aligned_sim_df.columns: aligned_sim_df['date'] = pd.to_datetime(aligned_sim_df['date'], errors='coerce')
                
                if 'regularized' not in aligned_sim_df.columns and 'regularized' in filtered_df.columns:
                    aligned_sim_df['regularized'] = np.nan

                filtered_df = pd.concat([filtered_df, aligned_sim_df], ignore_index=True)
                filtered_df.sort_values(by='date', inplace=True)
                print(f"Concatenated simulated data. Total rows now: {len(filtered_df)}")
            except Exception as e:
                print(f"Error processing simulated TSV {sim_tsv_path}: {e}")
        else:
            print(f"Warning: Simulated TSV file not found: {sim_tsv_path}")

    complete_pb_tree_path = None
    base_name_prefix_for_sim = "" # To be used by sim placement if needed
    abbr_state_for_sim = ""

    if args.target_state:
        sanitized_pango = args.pango_lineage.replace("/", "_").replace(".", "")
        sanitized_state_for_filename = args.target_state.replace(" ", "_")
        base_name_prefix_for_sim = f"{sanitized_state_for_filename}_{sanitized_pango}"
        try:
            abbr_state_for_sim = us.states.lookup(args.target_state).abbr
        except:
            abbr_state_for_sim = sanitized_state_for_filename
        
        complete_pb_tree_path = create_reductions(filtered_df,
                                                  output_folder,
                                                  base_name_prefix_for_sim,
                                                  args.num_reduction_steps,
                                                  args.target_state,
                                                  input_mat_path,
                                                  args.ablate)
    else:
        print("No target_state provided. Skipping reductions and matUtils steps.")

    if args.sim_fasta and complete_pb_tree_path:
        sim_fasta_file = Path(args.sim_fasta)
        if sim_fasta_file.exists():
            run_simulation_placement(sim_fasta_file, 
                                     complete_pb_tree_path, 
                                     output_folder,
                                     base_name_prefix_for_sim,
                                     abbr_state_for_sim)
        else:
            print(f"Error: --sim_fasta file not found: {args.sim_fasta}")
    elif args.sim_fasta and not complete_pb_tree_path:
        print("Warning: --sim_fasta provided, but the 'complete' reduction tree was not generated. Skipping simulation placement.")

    print("Script finished.")