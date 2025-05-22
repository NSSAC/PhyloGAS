

"""

awarren

date:  06 oct 2023 (actually wrote much earlier).

Purpose:
Take a collection of genomic sequences, and do one of two analyses, or both.
Each analysis is one execution of the code, even if you do both analyses---in this
case, you run the code twice.


Analysis 1:
Determine the Shannon entropy (in the form of a threshold) for each column (i.e., location, slot)
in the sequences.

Analysis 2:
Given one sequence, determine a second (following) sequence by possibly modifying
one or more "slots" from the first sequence.
This modification is based on the Shannon entropies calculated in analysis 1.
So you can chain a list of sequences and their changes this way.

"""


from matplotlib import pyplot as plt
import matplotlib as mpl
import networkx as nx
from locale import currency
from multiprocessing import current_process
from Bio import AlignIO
import pandas as pd
import subprocess
import os
import numpy as np
import glob
import sys
from Bio import SeqIO
import random
import time
import lzma
import json
import gzip # Added for aligned_to_df

import argparse
from itertools import islice, cycle
import math # Added for rate limiting

# ====================================
# Constants.
__version__ = '0.0.12' # Incremented version
# Analysis types
ENTROPY_ANALYSIS="entropy_analysis"
GEN_SEQUENCE_ANALYSIS="generate_sequence_analysis"
BOTH="both"

# Compression types
XZ="xz"
PARQUET="parquet" # not currently supported

# Includes ambiguous nucleotides - moved to global scope for broader use
LETTERS = np.array(
    ["A", "C", "G", "T", "N", "R", "K", "S", "Y", "M", "W", "B", "H", "D", "V"]
)
# ====================================
def getClas():
    """
    Read in command line arguments.
    :return:
    """
    parser = argparse.ArgumentParser()

    # For both analyses.
    parser.add_argument('--analysis_type', type=str, dest='analysis_type', required=True,
                        choices=[ENTROPY_ANALYSIS, GEN_SEQUENCE_ANALYSIS, BOTH],
                        help='Type of analysis to run.')
    parser.add_argument("--threshold_file", type=str,dest="threshold_file",required=True, help="file containing column threshold values.")
    parser.add_argument("--base_threshold_df", type=str,dest="base_threshold_df",required=True, help="base name of files containing threshold dfs (expects .npy extension for prob_matrix).") # Clarified help
    parser.add_argument("--align_fasta", type=str, default=None, nargs='?', dest="align_fasta", required=False, help="path to alignment file in FASTA format")
    parser.add_argument("--seed_fasta", type=str, default=None, nargs='?', dest="seed_fasta", required=False, help="path to seed file in FASTA format; defaults to align_fasta if not set")
    parser.add_argument("--random_number_seed", type=int, dest="random_number_seed", required=True, help="if < 0, then random assignment")

    # For genomic sequences analysis.
    parser.add_argument("--start_date", default="2021-05-31", dest="start_date", required=False, type=str, help="simulation alignment to date")
    parser.add_argument("--input_graph_csv", type=str,dest="input_graph_csv",required=False, help="directed graph file; nodes are infections.")
    parser.add_argument("--output_prefix", default="syn_gen", type=str, dest="output_prefix", required=False, help="prefix for output file name (for fasta and metadata files)")
    paint_group = parser.add_mutually_exclusive_group(required=False)
    paint_group.add_argument("--input_graph_painted_state", type=str, dest="input_graph_painted_state", default="var1E", help="Infection state that gets painted")
    paint_group.add_argument("--input_graph_painted_prefix", type=str, dest="input_graph_painted_prefix", default=None, help="Prefix for infection states that get painted")
    parser.add_argument("--proportional", default=True, action="store_true", dest="proportional", required=False, help="use proportional letter choices")
    parser.add_argument("--neutral", default=False, action="store_false", dest="proportional", required=False, help="use neutral letter choices")
    parser.add_argument("--poor", default=False, action="store_true", dest="poor", required=False, help="use poor mutational model")
    parser.add_argument("--limit", default=None, type=int, dest="limit", required=False, help="maximum number of items to process")
    tick_group = parser.add_mutually_exclusive_group(required=False)
    tick_group.add_argument("--end_tick", type=int, dest="end_tick",
                            help="End tick for processing the input_graph_csv. Mutually exclusive with --num_ticks.")
    tick_group.add_argument("--num_ticks", type=int, dest="num_ticks",
                            help="Number of ticks to process from the start of the simulation or from --start_tick if provided. Mutually exclusive with --end_tick.")
    parser.add_argument("--start_tick", type=int, dest="start_tick", required=False, default=0,
                        help="Start tick for processing the input_graph_csv. Defaults to 0.")

    parser.add_argument("--reference", default=None, type=str, dest="reference", required=False, help="add reference sequence to the output")
    parser.add_argument("--compression", default=None, type=str, dest="compression_type", required=False, help="add compression method -- None, xz, or parquet",
                       choices=["None",XZ])
    parser.add_argument("--persontrait_file", default=None, type=str, dest="persontrait_file", required=False, help="the full path to the persontrait data file with additional data")
    parser.add_argument("--add_metadata", default=None, type=str, dest="add_metadata", required=False, help="the columns (comma-delimited) from the persontrait_file to include in the metadata output")
    parser.add_argument("--location", default='{"country":"USA","division":"Virginia","divisionAbbr":"VA","region":"North America"}', type=str, dest="location", required=False, help="the location data for the infection record")
    parser.add_argument("--reference_location", default='{"country":"China","division":"Wuhan","divisionAbbr":"Hu","region":"Asia","date":"2019-12-26"}', type=str, dest="reference_location", required=False, help="the location data for the reference infection record")
    
    # START ADDED ARGUMENTS FOR RATE LIMITING
    parser.add_argument("--rate_limit", action="store_true", default=False, dest="rate_limit",
                        help="Enable rate-limiting of mutations based on within-host dynamics and iSNV paper.")
    parser.add_argument("--initial_viral_load", type=float, default=10.0, dest="initial_viral_load",
                        help="Initial viral load for rate-limiting model (relevant if --rate_limit is used).")
    # END ADDED ARGUMENTS

    parser.add_argument('--version', action='version', version=f'genetic_painter {__version__}')
    args = parser.parse_args()

    if (args.align_fasta == None):
        if (args.analysis_type != GEN_SEQUENCE_ANALYSIS):
            print("  Error.")
            print("  args.align_fasta has value None, which is not allowed.")
            parser.print_help()
            print("  Terminate.")    
            sys.exit(0)
        elif (args.seed_fasta == None):
            print("  Error.")
            print("  Either args.align_fasta or args.seed_fasta must be set.")
            parser.print_help()
            print("  Terminate.")    
            sys.exit(0)
            
    if args.rate_limit and args.initial_viral_load <= 0:
        parser.error("--initial_viral_load must be positive if --rate_limit is used.")

    return args


# ====================================
def write_output_entropy(args, thresh, prob_matrix, entropy_values):

    # Filename and base filename.
    threshold_file = args.threshold_file
    base_threshold_df_path = args.base_threshold_df # This is the base path, .npy will be added by np.save
    entropy_file = base_threshold_df_path + "_entropy.csv"

    # Write the entropy values to file.
    try:
        fh_out = open(entropy_file,"w")
    except:
        print("   Error")
        print("   Trying to open the output file, where entropy values are to be written.")
        print("   This failed.")
        print("   File name: ", entropy_file)
        print("   Terminate.")
        exit(1)
    # write out all the entropy values joined by newline
    fh_out.write("\n".join([str(e) for e in entropy_values]))
    fh_out.close()

    # Write the thresholds to file.
    try:
        fh_out = open(threshold_file,"w")
    except:
        print("   Error")
        print("   Trying to open the output file, where thresholds are to be written.")
        print("   This failed.")
        print("   File name: ", threshold_file)
        print("   Terminate.")
        exit(1)

    for ithresh in thresh:
        fh_out.write(str(ithresh) + "\n")

    fh_out.close()

    try:
        # np.save will add .npy extension if not present in base_threshold_df_path
        np.save(base_threshold_df_path, prob_matrix, allow_pickle=False)
        print(f"  Probability matrix saved to {base_threshold_df_path}.npy")
    except:
        print("   Error")
        print("   Trying to write to a probablity matrix to npy file.")
        print("   This failed.")
        print("   File name base: ", base_threshold_df_path)
        print("   Terminate.")
        exit(1)

    return


# ====================================
def load_thresholds_and_dfs(args):

    # Filenames to write things to.
    threshold_file = args.threshold_file
    # base_threshold_df is the base path, .npy is assumed by np.load if not present
    prob_matrix_file = args.base_threshold_df 
    if not prob_matrix_file.endswith('.npy'):
        prob_matrix_file += '.npy'


    # Output lists.
    thresh=list()
    # thresh_detail=list() # This seems to be unused if prob_matrix is loaded directly

    # Read thresholds from file.
    try:
        fh_in = open(threshold_file,"r")
    except:
        print("   Error")
        print("   Trying to open the output file, where thresholds are to be read.")
        print("   This failed.")
        print("   File name: ", threshold_file)
        print("   Terminate.")
        exit(1)

    for aline in fh_in:
        sline = aline.strip()
        if (len(sline)==0 or sline[0]=="#"):
            continue
        ithresh=(float)(sline)
        thresh.append(ithresh)
    fh_in.close()

    try:
        prob_matrix = np.load(prob_matrix_file, allow_pickle=False)
        print(f"  Probability matrix loaded from {prob_matrix_file}")
    except FileNotFoundError:
        print(f"   Error: Probability matrix file not found: {prob_matrix_file}")
        print("   Ensure you have run the entropy_analysis first or provided the correct path.")
        print("   Terminate.")
        exit(1)
    except Exception as e:
        print(f"   Error loading probability matrix from {prob_matrix_file}: {e}")
        print("   Terminate.")
        exit(1)
        
    return thresh, prob_matrix


# ====================================
def main():

    args = getClas()


    # Seed random numbers.
    # If number is < 0, then using random seeding.
    if args.random_number_seed >= 0:
        random.seed(args.random_number_seed)
        np.random.seed(args.random_number_seed)

    analysis_type = args.analysis_type

    if analysis_type == BOTH or analysis_type==ENTROPY_ANALYSIS:
        # Compute the shannon entropies for the colummns of a
        # group of sequences.
        print("  \n\n --- doing entropy calculations --- \n\n")
        compute_entropy(args)

    if analysis_type == BOTH or analysis_type==GEN_SEQUENCE_ANALYSIS:
        # Determine perturbations in a series of sequences.
        print("  \n\n --- generating sequences --- \n\n") # Added print statement
        generate_sequences(args)


    return

def aligned_to_df(align_file):
     # read in alignment to pandas dataframe
    print('reading alignment file into pandas dataframe.....')
    if align_file.endswith('.gz'):
        open_func = gzip.open
    elif align_file.endswith('.xz'):
        open_func = lzma.open
    else:
        open_func = open

    with open_func(align_file, 'rt') as file:
        align = AlignIO.read(file, 'fasta')
    
    # name = [] # Not used
    # description = [] # Not used
    # for record in align:
    #     name.append(record.name)
    #     description.append(record.description)
    align_list_of_lists = [list(str(record.seq)) for record in align]
    align2 = pd.DataFrame(align_list_of_lists) # More direct conversion
    print(f"  Alignment dimensions: {align2.shape}")
    return align2

def df_to_entropy(align2):
    # create threshold list, each column threshold included
    print('calculating entropy and getting the threshold...')
    thresh = []
    thresh_detail_dfs = [] # Renamed from thresh_detail to avoid confusion with list of Series
    entropy_values = []
    for i in range(len(align2.columns)):
        df1 = align2.iloc[:, i].value_counts(normalize=True)

        # If we have N (unknown) then use the frequencies of others only
        # (distribute N equally between the other letters and gaps).
        # special case: If we have only N or N and gap only, then distribute N equally among ACTG.
        
        # Create a copy to avoid SettingWithCopyWarning if df1 is a slice
        df1_copy = df1.copy()

        if "-" in df1_copy.index:
            if len(df1_copy) == 1: # Only "-"
                value = 0.25
                # Replace df1_copy with a new Series
                df1_copy = pd.Series([value, value, value, value], index=["A", "C", "G", "T"])
            else:
                # Drop "-" and renormalize
                df1_copy = df1_copy.drop(index="-")
                if not df1_copy.empty: # Check if anything is left after dropping "-"
                    df1_copy = df1_copy / df1_copy.sum()
                else: # If only "-" was present with other non-ACGTN chars that also got dropped
                      # or if it was all N and -, this case might need refinement based on desired behavior
                    value = 0.25
                    df1_copy = pd.Series([value, value, value, value], index=["A", "C", "G", "T"])


        df1_copy.name = i # Set name for pd.concat
        thresh_detail_dfs.append(df1_copy)

        e_act, thresh_val = column_entropy_thresh(df1_copy) # Use df1_copy
        thresh.append(thresh_val)
        entropy_values.append(e_act)
        
    # Create a matrix of probabilities for each letter at each position,
    # fill missing with 0 and reindex to maintain order
    prob_matrix_df = pd.concat(thresh_detail_dfs, axis=1).reindex(LETTERS).fillna(0)
    prob_matrix = prob_matrix_df.values # Convert to numpy array
    return thresh, prob_matrix, entropy_values


# ====================================
def compute_entropy(args):
    """
    Compute the entropy for each column (i.e., each position)
    of a collection of genomic sequences.
    :param args:  the CLAs.
    :return:
    """

    # start_date = args.start_date  # start of Delta strain - UNUSED in this function

    # read in alignment to pandas dataframe
    align2 = aligned_to_df(args.align_fasta)
    thresh, prob_matrix, entropy_values = df_to_entropy(align2) # prob_matrix is now numpy array

    write_output_entropy(args, thresh, prob_matrix, entropy_values)

    return

# ====================================
def create_aug_metadata_dict(metadata_cols, pid, pid_df=None, standardized_replacements=None):
    temp_dict = {}
    if len(metadata_cols) > 0:
        for col in metadata_cols:
            # Standardize pid_df access, ensuring pid_df is a Series for a single pid
            # or handle cases where pid_df might be None or col not present
            if standardized_replacements is not None:
                real_col = standardized_replacements[col]
            else:
                real_col = col
            col_value = "NA" # Default
            if pid_df is not None:
                if real_col == "pid": # pid itself is not usually in persontrait_df by that name
                    col_value = pid
                elif real_col in pid_df.index: # Check if column exists for this pid_df (Series)
                    val_from_df = pid_df[real_col]
                    if col in ["sex", "gender"]:
                        if pd.isna(val_from_df): col_value = "NA"
                        elif val_from_df == 1: col_value = "male"
                        elif val_from_df == 2: col_value = "female"
                        else: col_value = str(val_from_df) # Or "unknown"
                    else:
                        col_value = str(val_from_df) if not pd.isna(val_from_df) else "NA"
            temp_dict[col] = col_value # Use [] for assignment
    return temp_dict


# ====================================
# START HELPER FUNCTIONS FOR RATE LIMITING
def calculate_replication_cycles(initial_virions, target_early_population, burst_size):
    """Calculates the number of replication cycles to reach target_early_population."""
    if initial_virions <= 0 or target_early_population <= 0 or burst_size <= 1:
        return 1 # Avoid math errors, assume at least 1 cycle if inputs are problematic
    if initial_virions >= target_early_population:
        return 1 # Already at or above target, assume 1 cycle for potential mutations
    
    # Formula: target = initial * (burst_size ^ cycles)
    # cycles = log_burst_size(target / initial)
    try:
        cycles = math.log(target_early_population / initial_virions, burst_size)
        return math.ceil(cycles)
    except ValueError: # e.g. log of zero or negative
        return 1

def sample_power_law(min_frequency, max_frequency, alpha=2.0):
    """Samples a frequency from a power-law distribution P(x) ~ x^-alpha."""
    # Using inverse transform sampling: F(x) = (x^(1-alpha) - min^(1-alpha)) / (max^(1-alpha) - min^(1-alpha))
    # Solve for x: x = [(F(x) * (max^(1-alpha) - min^(1-alpha))) + min^(1-alpha)] ^ (1/(1-alpha))
    # F(x) is u (random number from 0 to 1)
    u = random.random()
    # Handle alpha = 1 case separately if needed, but typical iSNV alpha is ~2
    if alpha == 1.0: # Avoid division by zero if 1-alpha is zero
        # P(x) ~ 1/x. CDF is (ln(x) - ln(min)) / (ln(max) - ln(min))
        # x = exp(u * (ln(max) - ln(min)) + ln(min))
        # x = exp(u*ln(max/min) + ln(min)) = min * (max/min)^u
        return min_frequency * ((max_frequency / min_frequency) ** u)

    # Normal case for alpha != 1
    # Numerator for the exponent term
    term_min = min_frequency**(1.0 - alpha)
    term_max = max_frequency**(1.0 - alpha)
    
    sampled_value = (u * (term_max - term_min) + term_min)**(1.0 / (1.0 - alpha))
    return max(min_frequency, min(max_frequency, sampled_value)) # Ensure bounds

# END HELPER FUNCTIONS FOR RATE LIMITING

# ====================================
def generate_sequences(args):

    output_file_prefix = args.output_prefix

    if args.compression_type is None or args.compression_type == "None":
        fasta_to_write = output_file_prefix + ".sequences.fasta"
        metadata_file_to_write = output_file_prefix + ".metadata.tsv"
    elif args.compression_type == XZ:
        fasta_to_write = output_file_prefix + ".sequences.fasta.xz"
        metadata_file_to_write = output_file_prefix + ".metadata.tsv.xz"
    else: 
        print("   Warning: Unsupported compression type. Continuing with no compression.")
        fasta_to_write = output_file_prefix + ".sequences.fasta"
        metadata_file_to_write = output_file_prefix + ".metadata.tsv"

    augment_metadata = False
    if args.persontrait_file and args.add_metadata:
        augment_metadata = True
        aug_metadata_columns = args.add_metadata.split(",")
        try:
            persontrait_df = pd.read_csv(args.persontrait_file).set_index("pid")
        except FileNotFoundError:
            print(f"  Error: persontrait_file {args.persontrait_file} not found. Cannot augment metadata.")
            augment_metadata = False # Turn off augmentation
        except KeyError: # 'pid' not in columns
            print(f"  Error: 'pid' column not found in {args.persontrait_file}. Cannot augment metadata.")
            augment_metadata = False
    elif args.persontrait_file or args.add_metadata:
        print("   Info: persontrait_file and add_metadata must BOTH be provided to augment metadata. Not augmenting.")


    use_poor_mut_model = args.poor
    # use_proportional = args.proportional # This is used to select letters_to_use
    seq_limit = args.limit
    input_graph_csv = args.input_graph_csv
    start_date = args.start_date

    thresh, prob_matrix = load_thresholds_and_dfs(args)
    
    print('reading in the network data....')
    df = pd.read_csv(input_graph_csv)

    if args.compression_type == XZ:
        seq_file = lzma.open(fasta_to_write, 'wb')
        metadata_file = lzma.open(metadata_file_to_write, 'wb')
    else:
        seq_file = open(fasta_to_write, 'w')
        metadata_file = open(metadata_file_to_write, 'w')

    line_keys=["virus","region","country","division","divisionExposure","date","strain"]
    custom_sim_keys = ["sim_pid", "sim_tick"]
    line_keys += custom_sim_keys
    meta_line = "\t".join(line_keys)

    if augment_metadata:
        # Standardize known column renames
        standardized_aug_cols = []
        standardized_replacements = {}
        for col in aug_metadata_columns:
            if col == "gender": 
                standardized_aug_cols.append("sex")
                standardized_replacements["sex"]="gender"
            elif col == "home_latitude": 
                standardized_aug_cols.append("latitude")
                standardized_replacements["latitude"] = "home_latitude"
            elif col == "home_longitude": 
                standardized_aug_cols.append("longitude")
                standardized_replacements["longitude"] = "home_longitude"
            else: 
                standardized_aug_cols.append(col)
                standardized_replacements[col] = col
        
        aug_metadata_str = "\t".join(standardized_aug_cols)
        meta_line += "\t" + aug_metadata_str

    meta_line += "\n"

    if args.compression_type == XZ:
        metadata_file.write(meta_line.encode())
    else:
        metadata_file.write(meta_line)

    ref_location_dict = json.loads(args.reference_location)
    if args.reference is not None:
        align_ref = AlignIO.read(args.reference, "fasta")
        infection = InfectionRecord()
        country_ref=ref_location_dict['country'] # Use distinct names for clarity
        division_ref=ref_location_dict['division']
        divisionAbbr_ref=ref_location_dict['divisionAbbr']
        region_ref=ref_location_dict['region']
        date_ref=ref_location_dict['date']
        #f"{division_ref}-{divisionAbbr_ref}-1/{date_ref.split('-')[0]}""
        infection.fromEpihiper("ncov", region_ref, country_ref, division_ref, division_ref, date_ref, align_ref[0].id)

        aug_metadata_dict_ref = {} # Initialize for reference
        if augment_metadata:
             # For reference, PID is typically not applicable unless you have specific metadata for it
            aug_metadata_dict_ref = create_aug_metadata_dict(standardized_aug_cols, pid="reference_strain", standardized_replacements=standardized_replacements) # Pass standardized
        
        add_to_fasta(str(align_ref[0].seq), infection, seq_file, args.compression_type)
        write_metadata(metadata_file, infection, line_keys, args.compression_type, 
                       aug_metadata_columns=standardized_aug_cols if augment_metadata else None, # Pass standardized
                       aug_metadata_dict=aug_metadata_dict_ref if augment_metadata else None)


    location_dict = json.loads(args.location)
    country = location_dict["country"]
    division = location_dict["division"]
    divisionAbbr = location_dict["divisionAbbr"]
    region = location_dict["region"]

    loop_counter=0
    infection_counter = 0 # Moved initialization here

    paint_this = (lambda state: state == args.input_graph_painted_state)
    if args.input_graph_painted_prefix:
        paint_this = (lambda state: state.startswith(args.input_graph_painted_prefix))

    transitions_to_paint = df[df["exit_state"].map(paint_this)]
    seed_transitions_mask = transitions_to_paint["contact_pid"] == -1 # Use mask for efficiency

    seed_df = transitions_to_paint[seed_transitions_mask] # Renamed from 'seed' to 'seed_df'
    
    if args.seed_fasta == None:
        align_seed_records = list(AlignIO.read(args.align_fasta, 'fasta')) # Read once
    else:
        align_seed_records = list(AlignIO.read(args.seed_fasta, 'fasta')) # Read once
    
    # Assign seed sequences
    seed_pids = seed_df["pid"].tolist()
    N = len(seed_pids)
    M = len(align_seed_records)
    
    current_sequences = {}
    if M == 0 and N > 0:
        print("Error: No sequences in seed FASTA file, but seed transitions exist. Cannot proceed.")
        sys.exit(1)

    temp_seed_seqs = {}
    seed_seq_dict = {}
    if N <= M:
        for i, (pid_val, tick) in enumerate(zip(seed_df["pid"], seed_df["tick"])):
            temp_seed_seqs[pid_val] = np.array(list(align_seed_records[i].seq))
            seed_seq_dict[f"{pid_val}.{tick}"] = align_seed_records[i]
    else:  # set temp_seed_seqs to first N sequences
        for i, (pid_val, tick) in enumerate(zip(seed_df["pid"], seed_df["tick"])):
            if i < M:
                temp_seed_seqs[pid_val] = np.array(list(align_seed_records[i].seq))
                seed_seq_dict[f"{pid_val}.{tick}"] = align_seed_records[i]
            else:
                break  
    #else: # N > M, cycle through align_seed_records
    #    align_str_list = [np.array(list(r.seq)) for r in align_seed_records]
    #    cycled_seqs = list(islice(cycle(align_str_list), N))
    #    for i, pid_val in enumerate(seed_pids):
    #        temp_seed_seqs[pid_val] = cycled_seqs[i]
            
    current_sequences.update(temp_seed_seqs)


    # uncomment to keep seed sequences out
    #transitions_to_paint_df = transitions_to_paint[~seed_transitions_mask][ # Renamed
    #    ["pid", "contact_pid", "tick"] # Removed exit_state as it's already filtered
    #]
    transitions_to_paint_df = transitions_to_paint[["pid", "contact_pid", "tick"]]

 # --- START: TICK-BASED FILTERING ---
    print(f"  Initial number of transitions to paint: {len(transitions_to_paint_df)}")

    # Apply start_tick
    # The default for args.start_tick is 0, so this filter will always be applied.
    # If you want to truly skip start_tick unless specified, default should be None and check for it.
    # Assuming default=0 means we always filter from at least tick 0.
    if args.start_tick > 0: # Only filter if start_tick is greater than the absolute minimum
        transitions_to_paint_df = transitions_to_paint_df[transitions_to_paint_df["tick"] >= args.start_tick]
        print(f"  After applying --start_tick {args.start_tick}: {len(transitions_to_paint_df)} transitions")

    # Apply end_tick or num_ticks (mutually exclusive due to argparse group)
    if args.end_tick is not None:
        # Filter by end_tick (inclusive of end_tick)
        transitions_to_paint_df = transitions_to_paint_df[transitions_to_paint_df["tick"] <= args.end_tick]
        print(f"  After applying --end_tick {args.end_tick}: {len(transitions_to_paint_df)} transitions")
    elif args.num_ticks is not None:
        # Determine the effective end tick based on num_ticks and start_tick
        # The starting point for num_ticks is args.start_tick if specified, 
        # otherwise it's the minimum tick present in the (potentially already start_tick filtered) data.
        
        if not transitions_to_paint_df.empty:
            # If args.start_tick was used to filter, current_min_tick_for_num_ticks is effectively args.start_tick
            # If args.start_tick was 0 (or less than actual min tick), then use the actual min tick in the current df
            current_min_tick_for_num_ticks = transitions_to_paint_df["tick"].min()
            # If start_tick was specified and is > current_min_tick, num_ticks should start from start_tick
            effective_start_for_num_ticks = max(current_min_tick_for_num_ticks, args.start_tick)
            
            calculated_end_tick = effective_start_for_num_ticks + args.num_ticks - 1 # -1 because num_ticks includes the start_tick itself

            transitions_to_paint_df = transitions_to_paint_df[transitions_to_paint_df["tick"] <= calculated_end_tick]
            print(f"  After applying --num_ticks {args.num_ticks} (from effective start {effective_start_for_num_ticks}, calculated end: {calculated_end_tick}): {len(transitions_to_paint_df)} transitions")
        else:
            print(f"  DataFrame empty before applying --num_ticks, no further filtering.")
            
    # --- END: TICK-BASED FILTERING ---

    
    if args.limit and args.limit > 0 : # Check if args.limit is set and positive
        if len(transitions_to_paint_df) > args.limit:
            transitions_to_paint_df = transitions_to_paint_df.iloc[:args.limit]
            print(f"  After applying --limit {args.limit}: {len(transitions_to_paint_df)} transitions")
        else:
            print(f"  Number of transitions ({len(transitions_to_paint_df)}) is already within --limit {args.limit}.")

    if not transitions_to_paint_df.empty:
        transitions_to_paint_df["date"] = pd.to_datetime(start_date) + transitions_to_paint_df["tick"].map(pd.offsets.Day)
    elif not current_sequences: # No seeds initialized AND no transitions from other sources
        print("  No seed sequences initialized and no transitions to process. Exiting.")
        seq_file.close()
        metadata_file.close()
        print("Done generating sequences (no work performed).")
        return
    else: # Seeds might be initialized, but no subsequent transitions in the filtered range
        print("  No transitions to process after filtering (seeds may have been initialized).")
        # The script might still write out the initial seed sequences if any were processed
        # Or it might just end if the loop below doesn't run.
        # Decide if you want to write just seeds if no transmissions. For now, it will proceed.
    
    thresh_np = np.array(thresh) # Convert list to numpy array for determine_change

    # prob_matrix is loaded by load_thresholds_and_dfs and should always be available here.
    # If prob_matrix loading failed, the script would have exited earlier.
    if prob_matrix.size == 0: # Should not happen if load_thresholds_and_dfs succeeded
        print("Error: prob_matrix is empty. Cannot determine sequence length. Exiting.")
        # Close files if open
        if 'seq_file' in locals() and not seq_file.closed: seq_file.close()
        if 'metadata_file' in locals() and not metadata_file.closed: metadata_file.close()
        sys.exit(1)
        
    example_sequence_length = prob_matrix.shape[1]
    print(f"  Using sequence length from prob_matrix: {example_sequence_length}")

    # Optional: Sanity check if seeds exist and match this length
    if current_sequences:
        first_seed_seq_len = len(next(iter(current_sequences.values())))
        if first_seed_seq_len != example_sequence_length:
            print(f"  Warning: Seed sequence length ({first_seed_seq_len}) does not match "
                  f"prob_matrix length ({example_sequence_length}). Proceeding with prob_matrix length.")
            # Potentially raise an error here if this is considered a critical mismatch:
            # sys.exit("Critical error: Seed sequence length and prob_matrix length mismatch.")


    if not args.proportional:
        letters_to_use = np.array(["A", "C", "G", "T"])
        n_letters = len(letters_to_use)
        # Create a prob_matrix where each of A, C, G, T has 0.25 probability, others 0
        neutral_prob_matrix = np.zeros((len(LETTERS), example_sequence_length))
        acgt_indices = [np.where(LETTERS == L)[0][0] for L in letters_to_use]
        neutral_prob_matrix[acgt_indices, :] = 1.0 / n_letters
        cumulative_probs_matrix = np.cumsum(neutral_prob_matrix, axis=0)

    else:
        letters_to_use = LETTERS # This is already a np.array
        cumulative_probs_matrix = np.cumsum(prob_matrix, axis=0) # prob_matrix is from loaded data

    assert example_sequence_length == cumulative_probs_matrix.shape[1], "Sequence length must match columns in probability matrix"
    assert len(LETTERS) == cumulative_probs_matrix.shape[0], "Number of global LETTERS must match rows in cumulative probability matrix"

    # --- Rate Limiting constants (can be made CLI args later) ---
    mutation_rate_per_cycle = 3.40e-6
    peak_viral_load = 1e6  # Example peak viral load
    # Fraction of peak viral load to define "early" phase, for calculating replication cycles
    rt_early_population_threshold = 0.01 
    # Probability a mutation occurring in "early" phase (defined by cycles) becomes major
    rt_early_mutation_probability = 0.8 
    # Min/max burst size for sampling
    min_burst_size = 10
    max_burst_size = 1000


    for _, pid, contact_pid, date_obj, tick in transitions_to_paint_df[ # Use date_obj to avoid name clash
        ["pid", "contact_pid", "date", "tick"] 
    ].itertuples():
        loop_counter += 1
        infection_id=f"{pid}.{tick}"
        cur_strain_id=f"{country}/{divisionAbbr}-EHip-{infection_id}/{date_obj.year}"

        if loop_counter % 1000 == 0:
            print(f"    Processed {loop_counter} graph edges; Decorated {infection_counter} infections.")

        #check infection_id if seed sequence
        seed_fasta = seed_seq_dict.get(infection_id, None)
            
        if seed_fasta is None:
            if contact_pid not in current_sequences:
                print(f"Warning: contact_pid {contact_pid} not found in current_sequences. Skipping mutation for pid {pid}.")
                # Optionally, assign a default/random sequence or skip adding this pid to fasta
                continue 
                
            seq_to_change_arr = current_sequences[contact_pid] # This is a NumPy array of chars

        if seed_fasta is not None: #seed sequence doesn't change
            new_seq_arr = np.array(list(seed_fasta.seq))
            cur_strain_id = seed_fasta.id
        elif use_poor_mut_model: # This overrides other mutation models
            new_seq_str = poor_mut_model(seq_to_change_arr) # poor_mut_model expects array, returns string
            new_seq_arr = np.array(list(new_seq_str))
        elif args.rate_limit:
            # --- Rate Limiting Logic ---
            burst_size = random.randint(min_burst_size, max_burst_size)
            effective_initial_load = args.initial_viral_load
            
            replication_cycles_for_early_phase = calculate_replication_cycles(
                effective_initial_load,
                rt_early_population_threshold * peak_viral_load,
                burst_size
            )
            replication_cycles_for_early_phase = max(1, replication_cycles_for_early_phase)

            num_potential_mutations = np.random.poisson(
                mutation_rate_per_cycle * len(seq_to_change_arr) * replication_cycles_for_early_phase
            )
            num_potential_mutations = min(num_potential_mutations, len(seq_to_change_arr))

            final_change_mask = np.zeros(len(seq_to_change_arr), dtype=bool)

            if num_potential_mutations > 0:
                # --- Weighted Site Selection ---
                # Weights: Higher for more entropy (lower threshold)
                # Ensure weights are non-negative
                site_weights = np.maximum(0.0, 1.0 - (thresh_np / 100.0)) 
                
                # Do not allow mutating existing gaps
                gap_indices = np.where(seq_to_change_arr == '-')[0]
                site_weights[gap_indices] = 0.0

                # Normalize weights (optional, but good practice for probabilities)
                # If sum_weights is 0, it means no sites are mutable, handle this.
                sum_weights = np.sum(site_weights)
                if sum_weights > 0:
                    # Using random.choices for weighted sampling.
                    # Note: random.choices samples WITH replacement by default.
                    # If num_potential_mutations is small relative to genome length,
                    # duplicates are rare. If concerned, could sample more and take unique,
                    # or implement a more complex weighted sampling without replacement.
                    potential_mutation_indices = random.choices(
                        population=range(len(seq_to_change_arr)), 
                        weights=site_weights, 
                        k=num_potential_mutations
                    )
                    
                    for site_idx in potential_mutation_indices:
                        # Fixation probability for this "successful" iSNV to become major
                        if random.random() < rt_early_mutation_probability:
                            # The check for gap 'if seq_to_change_arr[site_idx] != '-':' is now implicitly
                            # handled by site_weights[gap_indices] = 0.0, as such sites
                            # should not be chosen by random.choices if their weight is 0.
                            # However, a direct check before assignment is still a safeguard.
                            if seq_to_change_arr[site_idx] != '-': # Safeguard
                                final_change_mask[site_idx] = True
                # --- End Weighted Site Selection ---
            
            new_seq_arr = weighted_change(
                seq_to_change_arr, final_change_mask, cumulative_probs_matrix, letters=letters_to_use
            )

        else: # Original model (not poor, not rate-limited)
            change_mask = determine_change(thresh_np, seq_to_change_arr)
            new_seq_arr = weighted_change(
                seq_to_change_arr, change_mask, cumulative_probs_matrix, letters=letters_to_use
            )
        
        current_sequences[pid] = new_seq_arr # Store the array
        new_seq_str = "".join(new_seq_arr.tolist()) # Convert to string for FASTA

        infection = InfectionRecord()
        infection.fromEpihiper(
            "ncov",
            region, country, division, division, # Assuming divisionExposure is same as division
            date_obj.strftime("%Y-%m-%d"),
            cur_strain_id
        )
        infection.populate_sim_details(pid, tick)

        aug_metadata_dict_current = {} # Initialize for current infection
        if augment_metadata:
            try:
                pid_df_series = persontrait_df.loc[pid] # This should be a Series
                
                aug_metadata_dict_current = create_aug_metadata_dict(
                    standardized_aug_cols, pid, pid_df_series, standardized_replacements # Pass standardized
                )
            except KeyError: # pid not in persontrait_df
                 aug_metadata_dict_current = create_aug_metadata_dict(standardized_aug_cols, pid, standardized_replacements=standardized_replacements) # Will fill with NA


        add_to_fasta(new_seq_str, infection, seq_file, args.compression_type)
        write_metadata(
            metadata_file, infection, line_keys, args.compression_type,
            aug_metadata_columns=standardized_aug_cols if augment_metadata else None, # Pass standardized
            aug_metadata_dict=aug_metadata_dict_current if augment_metadata else None
        )
        infection_counter += 1

    seq_file.close()
    metadata_file.close()

    print("Done generating sequences.")
    return


# ====================================
def column_entropy_thresh(freq_df): # freq_df is a pandas Series
    e_act = 0
    # For Shannon entropy, typically log base 2 is used for bits, or ln for nats.
    # The formula for max entropy E_max = -log(1/N) = log(N) where N is alphabet size.
    # If using all LETTERS, N = len(LETTERS). If ACGTN, N=5.
    # The original code implies N=5 (A,C,G,T, and implicitly N or something else making up the 5th category for p_xm)
    # Let's stick to the paper's likely intention or common practice. If it's DNA/RNA, N=4 (or 5 with N).
    # The provided freq_df here is *after* filtering out '-', so it contains actual characters.
    
    alphabet_size_for_max_entropy = 4 # Assuming ACGT for max entropy reference point
    # If freq_df is empty or sums to zero, handle to avoid division by zero or NaN
    if freq_df.empty or freq_df.sum() == 0:
        return 0, 0 # Or some other default for no information/all gaps

    for p_xi in freq_df: # Iterate over values (frequencies)
        if p_xi > 0: # log(0) is undefined
             e_act -= p_xi * np.log(p_xi) # Using natural log (nats)

    # Max entropy for an alphabet of size N is log(N)
    # The original code used p_xm = 1/5.0 ... e_max += p_xm * np.log(p_xm) which is -log(5)
    # This implies comparison to a 5-symbol alphabet.
    # If we only consider ACGT for e_max, then it's -log(4).
    # Let's keep the original e_max logic for consistency unless specified otherwise.
    # This e_max calculation is a bit unusual if freq_df can have more/less than 5 symbols.
    # A more standard H_max = log(len(freq_df.index)) if all symbols in freq_df are equally likely
    # Or H_max = log(alphabet_size_for_max_entropy)
    
    # Replicating original e_max:
    # e_max_val = -np.log(5.0) # This seems to be the intended reference max entropy
    # A more standard approach: if the alphabet is ACGT, max entropy is log(4). If ACGTN, log(5).
    # If freq_df contains only ACGT, then using log(5) as max might be strange.
    # Let's use log of the number of unique characters in the column, or a fixed alphabet like ACGTN.
    
    # For consistency with the original (1 - e_act/e_max) * 100:
    # e_max needs to be negative if e_act is negative (as calculated from p*log(p))
    # So, if p_xm = 1/N, e_max_contrib = (1/N) * log(1/N). Sum N times: N * (1/N) * log(1/N) = log(1/N) = -log(N)
    
    # Consider alphabet size for max entropy. If it's ACGT, then 4. If ACGTN, then 5.
    # The original code used '5' implicitly in p_xm = 1/float(5).
    ref_alphabet_size = 5 
    e_max_val = -np.log(ref_alphabet_size) # e.g., for ACGTN, all equally likely

    if e_max_val == 0: # Avoid division by zero if, somehow, e_max_val is 0
        thresh = 0
    else:
        # Normalized entropy: H_norm = e_act / log(num_symbols_in_col)
        # The formula used: (1 - (e_act / e_max_val)) * 100
        # If e_act is close to e_max_val (high diversity), ratio is ~1, thresh is ~0.
        # If e_act is close to 0 (low diversity, one symbol dominates), ratio is ~0, thresh is ~100.
        # This seems correct: high threshold means high consistency (low entropy).
        thresh = (1 - (e_act / e_max_val)) * 100
    
    if np.isnan(thresh):
        thresh = 0 # If only one symbol in column, e_act can be 0. If e_max_val is also 0 (e.g. 1 symbol alphabet), NaN.

    return e_act, max(0, min(100, thresh)) # Clamp threshold 0-100


# ====================================
def determine_change(thresh_array, seq_to_change_array): # Now takes np arrays
    comparison_values = np.random.randint(0, 100, len(thresh_array))
    # random_selection is True if val > threshold (i.e. more random than consistent column allows mutation)
    random_selection = comparison_values > thresh_array 
    # If the position is a gap, then don't change it
    return random_selection & (seq_to_change_array != "-")


# determine_change_old is not used.

# LETTERS is now global

def weighted_change(sequence_array, change_mask, cumulative_prob_matrix, letters=LETTERS):
    # sequence_array is a numpy array of characters
    # change_mask is a boolean numpy array
    # cumulative_prob_matrix is pre-calculated
    
    # Ensure sequence_array is a copy if it's going to be modified directly
    # and the original is needed elsewhere. Here, it's modified and returned.
    output_sequence_array = sequence_array.copy()

    # Only generate random values and find letters for positions that need to change
    num_to_change = np.sum(change_mask)
    if num_to_change == 0:
        return output_sequence_array

    # Get indices of positions to change
    change_indices = np.where(change_mask)[0]

    # Generate random numbers only for these positions
    random_values_for_change = np.random.rand(num_to_change)

    # Select relevant columns from cumulative_prob_matrix
    relevant_cum_probs = cumulative_prob_matrix[:, change_indices]

    # Determine the index of the letter for each position to change
    # (random_values_for_change[np.newaxis, :] < relevant_cum_probs) broadcasts random_values
    # .argmax(axis=0) finds the first True, which corresponds to the letter index
    letter_indices_for_change = (random_values_for_change[np.newaxis, :] < relevant_cum_probs).argmax(axis=0)
    
    new_letters_for_change = letters[letter_indices_for_change]

    # Apply changes
    output_sequence_array[change_indices] = new_letters_for_change
    
    return output_sequence_array


# thresh_detail_to_prob_matrix is not used if prob_matrix loaded from .npy

# commit_change is not used.
# intermediate_mut_model is not used.


def poor_mut_model(sequence_array): # Takes numpy array
    new_seq_list = [] # Build as list then join
    for nucleotide_char in sequence_array: # Iterate over chars in array
        # Original logic: change_val = np.random.randint(1, 2) -- this is always 1. So always try to change.
        # Assuming intent was 50% chance to change:
        if random.random() < 0.5: # 50% chance to enter this block
            # Original new_nucleotide = np.random.randint(1, 5) maps to ACGT, 5 was '-'
            # Let's use ACGT directly for simplicity.
            # Not clear if '-' was intended. If so, np.random.choice(['A','C','G','T','-'])
            new_nucleotide_char = random.choice(['A', 'C', 'G', 'T'])
            new_seq_list.append(new_nucleotide_char)
        else:
            new_seq_list.append(nucleotide_char) # Keep original

    return "".join(new_seq_list) # Returns string


class InfectionRecord:
    def __init__(self):
        # Initialize all known keys to None or a sensible default
        self.inf_dict = {
            "virus": None, "age": None, "country": None, "countryExposure": None,
            "date": None, "dateSubmitted": None, "died": None, "division": None,
            "divisionExposure": None, "fullyVaccinated": None, "strain": None,
            "gisaidClade": None, "gisaidEpiIsl": None, "hospitalized": None,
            "host": "Homo sapiens", # Default host
            "location": None, "month": None, 
            "nextcladePangoLineage": None, "nextstrainClade": None,
            "originatingLab": None, "pangoLineage": None, "region": None,
            "regionExposure": None, "samplingStrategy": None, "sex": None,
            "sraAccession": None, "strainold": None, "submittingLab": None, "year": None,
            "sim_pid": None, "sim_tick": None 
        }
    def populate_sim_details(self, sim_pid, sim_tick):
        self.inf_dict["sim_pid"] = sim_pid
        self.inf_dict["sim_tick"] = sim_tick

    def get(self,key,default=None):
        # Ensure that if default is None, we actually return None string if not present,
        # as TSV expects string representations or empty strings.
        val = self.inf_dict.get(key)
        if val is None:
            return str(default) if default is not None else "" # Return empty string for TSV if None
        return str(val) # Ensure string output for TSV

    def fromEpihiper(self, virus, region, country, division, divisionExposure, date, strain):
        self.inf_dict["virus"] = virus
        self.inf_dict["region"] = region # Note: GISAID often uses region for continent
        self.inf_dict["country"] = country
        self.inf_dict["division"] = division
        self.inf_dict["divisionExposure"] = divisionExposure
        self.inf_dict["date"] = date
        self.inf_dict["strain"] = strain
        if date:
            try:
                year, month, _ = date.split('-')
                self.inf_dict["year"] = year
                self.inf_dict["month"] = month
            except ValueError: # Date not in YYYY-MM-DD
                pass


def write_metadata(metadata_file, infection_record, line_keys, compression_type, aug_metadata_columns=None, aug_metadata_dict=None):
    """
    Writes metadata to a file. aug_metadata_columns should be the standardized list.
    """
    # Build list of values corresponding to line_keys first
    values = [infection_record.get(key) for key in line_keys]

    if aug_metadata_columns and aug_metadata_dict:
        for col_key in aug_metadata_columns: # Iterate over standardized keys
            # Get value from aug_metadata_dict; handle if key might be missing (shouldn't if create_aug_metadata_dict is robust)
            values.append(str(aug_metadata_dict.get(col_key, ""))) # "" if missing
    
    meta_line = "\t".join(values) + "\n"

    if compression_type == XZ:
        metadata_file.write(meta_line.encode('utf-8')) # Specify encoding
    else:
        metadata_file.write(meta_line)


def add_to_fasta(seq_str, infection_record, seq_file, compression_type): 
    # seq should be a string here
    seq_line = ">" + str(infection_record.get("strain")) + "\n" + seq_str + "\n"
    if compression_type == XZ:
        seq_file.write(seq_line.encode('utf-8')) # Specify encoding
    else:
        seq_file.write(seq_line)


# find_seq is not used.
# dfs_edges_with_ticks is not used.


if __name__ == '__main__':
    begin_time = time.time()
    main()
    end_time = time.time()
    time_s=end_time-begin_time
    time_hr=(float)(time_s)/3600.0
    print(f"   Execution time (s): {time_s:.2f}, (hr): {time_hr:.2f}")
    print("   --- good termination ---")

