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

import argparse
from itertools import islice, cycle

# ====================================
# Constants.
__version__ = '0.0.11'
# Analysis types
ENTROPY_ANALYSIS="entropy_analysis"
GEN_SEQUENCE_ANALYSIS="generate_sequence_analysis"
BOTH="both"

# Compression types
XZ="xz"
PARQUET="parquet" # not currently supported

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
    parser.add_argument("--base_threshold_df", type=str,dest="base_threshold_df",required=True, help="base name of files containing threshold dfs.")
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
    parser.add_argument("--limit", default=16521, type=int, dest="limit", required=False, help="maximum number of items to process")
    parser.add_argument("--reference", default=None, type=str, dest="reference", required=False, help="add reference sequence to the output")
    parser.add_argument("--compression", default=None, type=str, dest="compression_type", required=False, help="add compression method -- None, xz, or parquet",
                       choices=["None",XZ])
    parser.add_argument("--persontrait_file", default=None, type=str, dest="persontrait_file", required=False, help="the full path to the persontrait data file with additional data")
    parser.add_argument("--add_metadata", default=None, type=str, dest="add_metadata", required=False, help="the columns (comma-delimited) from the persontrait_file to include in the metadata output")
    parser.add_argument("--location", default='{"country":"USA","division":"Virginia","divisionAbbr":"VA","region":"North America"}', type=str, dest="location", required=False, help="the location data for the infection record")
    parser.add_argument("--reference_location", default='{"country":"China","division":"Wuhan","divisionAbbr":"Hu","region":"Asia","date":"2019-12-26"}', type=str, dest="reference_location", required=False, help="the location data for the reference infection record")
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

    # print(args.output_prefix)
    # print(args.proportional)
    # print(args.poor)

    return args


# ====================================
def write_output_entropy(args, thresh, prob_matrix, entropy_values):

    # Filename and base filename.
    threshold_file = args.threshold_file
    base_threshold_df = args.base_threshold_df
    entropy_file = base_threshold_df + "_entropy.csv"

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
        np.save(base_threshold_df, prob_matrix, allow_pickle=False)
    except:
        print("   Error")
        print("   Trying to write to a probablity matrix to npy file.")
        print("   This failed.")
        print("   File name base: ", base_threshold_df)
        print("   Terminate.")
        exit(1)

    return


# ====================================
def load_thresholds_and_dfs(args):

    # Filenames to write things to.
    threshold_file = args.threshold_file
    base_threshold_df = args.base_threshold_df

    # Output lists.
    thresh=list()
    thresh_detail=list()

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

    prob_matrix = np.load(args.base_threshold_df, allow_pickle=False)
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
    name = []
    description = []
    for record in align:
        name.append(record.name)
        description.append(record.description)
    align2 = pd.DataFrame(align)
    return align2

def df_to_entropy(align2):
    # create threshold list, each column threshold included
    print('calculating entropy and getting the threshold...')
    thresh = []
    thresh_detail=[]
    entropy_values = []
    for i in range(len(align2.columns)):
        df1 = align2.iloc[:, i].value_counts(normalize=True)

        # If we have N (unknown) then use the frequencies of others only
        # (distribute N equally between the other letters and gaps).
        # special case: If we have only N or N and gap only, then distribute N equally among ACTG.
        if "-" in df1.index:
            if len(df1) == 1:
                # This should be inconsequential, as we keep gaps as is
                value = 0.25
                df1 = pd.Series(
                    [value, value, value, value], index=["A", "C", "G", "T"]
                )
            else:
                # Don't use gaps for entropy calculation
                df1 = df1.drop(index="-")
                # Normalize the frequencies so they sum to 1.
                df1 = df1 / df1.sum()
        df1.name = i
        thresh_detail.append(df1)

        e_act, thresh_val = column_entropy_thresh(df1)
        thresh.append(thresh_val)
        entropy_values.append(e_act)
    # Create a matrix of probabilities for each letter at each position,
    # fill missing with 0 and reindex to maintain order
    prob_matrix = pd.concat(thresh_detail, axis=1).reindex(LETTERS).fillna(0).values
    return thresh, prob_matrix, entropy_values


# ====================================
def compute_entropy(args):
    """
    Compute the entropy for each column (i.e., each position)
    of a collection of genomic sequences.
    :param args:  the CLAs.
    :return:
    """

    start_date = args.start_date  # start of Delta strain

    # read in alignment to pandas dataframe
    align2 = aligned_to_df(args.align_fasta)
    thresh, thresh_detail, entropy_values = df_to_entropy(align2)


    ##########################################################

    write_output_entropy(args, thresh, thresh_detail, entropy_values)

    return

# ====================================
def create_aug_metadata_dict(metadata_cols, pid, pid_df=None):
    # Creates a dictionary of the desired metadata items that can be posted to
    # add_to_fasta
    temp_dict = {}

    if len(metadata_cols) > 0:
        for col in metadata_cols:
            if pid_df is None:
                temp_dict.update({col: "NA"})
            else:
                if col in ["sex", "gender"]:
                    if pid_df[col] == 1:
                        col_value = "male"
                    else:
                        col_value = "female"
                elif col == "pid":
                    col_value = pid
                else:
                    col_value = pid_df[col]
                temp_dict.update({col: col_value})
    return temp_dict


# ====================================
def generate_sequences(args):

    # Set these values to run the good or poor mutational model
    output_file_prefix = args.output_prefix

    # Check for special compression types
    if args.compression_type is None or args.compression_type == "None":
        fasta_to_write = output_file_prefix + ".sequences.fasta"
        metadata_file_to_write = output_file_prefix + ".metadata.tsv"
    elif args.compression_type == XZ:
        fasta_to_write = output_file_prefix + ".sequences.fasta.xz"
        metadata_file_to_write = output_file_prefix + ".metadata.tsv.xz"
    else: 
        # compression type is invalid
        print("   Warning")
        print("   Unsupported compression type specified.")
        print("   Supported compression types are None and xz")
        print("   Continuing with no compression.")
        fasta_to_write = output_file_prefix + ".sequences.fasta"
        metadata_file_to_write = output_file_prefix + ".metadata.tsv"

    # Check to see if persontrait_file is defined -- if so, augmenting metadata
    if args.persontrait_file is None and args.add_metadata is None:
        augment_metadata = False
    elif args.persontrait_file is None or args.add_metadata is None:
        # We need both the persontrait file and the columns to augment metadata
        # If either one is missing, then we don't add metadata
        print("   Info")
        print("   persontrait_file and add_metadata must be used together")
        print("   Since one is missing, not augmenting metadata")
        augment_metadata = False
    else:
        # augmenting metadata
        augment_metadata = True
        aug_metadata_columns = args.add_metadata.split(",")
        persontrait_df = pd.read_csv(args.persontrait_file).set_index("pid")

    use_poor_mut_model = args.poor
    use_proportional = args.proportional
    seq_limit = args.limit
    input_graph_csv = args.input_graph_csv
    start_date = args.start_date

    # Load the dataframe for each threshold.
    thresh, prob_matrix = load_thresholds_and_dfs(args)
    # prob_matrix = thresh_detail_to_prob_matrix(thresh_detail)
    # Read in network data
    print('reading in the network data....')
    df = pd.read_csv(input_graph_csv)

    ##########################################################

    # creating .fasta and .tsv files to append
    if args.compression_type == XZ:
        seq_file = lzma.open(fasta_to_write, 'wb')
        metadata_file = lzma.open(metadata_file_to_write, 'wb')
    else:
        seq_file = open(fasta_to_write, 'w')
        metadata_file = open(metadata_file_to_write, 'w')

    # Prep metadata TSV file with required column names:
    # https://docs.nextstrain.org/projects/ncov/en/latest/guides/data-prep/local-data.html#required-metadata

    # need two data structures
    # one list to just append to to generate MSA of all sequences as we go
    # another dict that maps node to it's most recent sequence

    if args.seed_fasta == None:
        align = AlignIO.read(args.align_fasta, 'fasta')
    else:
        align = AlignIO.read(args.seed_fasta, 'fasta')

    # This assumes the data read into df is in chronological order (ascending
    # according to tick)
    current_sequences = {}
    i = 0
    max_seed_value_index = len(align) - 1
    print(max_seed_value_index)
    strain_id = 0  # initialize strain_id for fasta, for now just increment a value, in the future could use node pid but would have to append to it in the case of multiple infections for a given pid

    sequences_mutated = 0

    line_keys=["virus","region","country","division","divisionExposure","date","strain"]
    meta_line = "\t".join(line_keys)

    if augment_metadata:
        # Add persontrait "add_metadata" columns
        aug_metadata_str = "\t".join(aug_metadata_columns).replace("gender","sex").replace("home_latitude","latitude").replace("home_longitude","longitude")
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
        country=ref_location_dict['country']
        division=ref_location_dict['division']
        divisionAbbr=ref_location_dict['divisionAbbr']
        region=ref_location_dict['region']
        date=ref_location_dict['date']
        infection.fromEpihiper("ncov", region, country, division, division, date, f"{division}-{divisionAbbr}-1/{date.split('-')[0]}")

        if augment_metadata:
            aug_metadata_dict = create_aug_metadata_dict(aug_metadata_columns,pid=-1)
            add_to_fasta(
                str(align_ref[0].seq), infection, seq_file, args.compression_type
            )
            write_metadata(metadata_file, infection, line_keys, args.compression_type, aug_metadata_columns, aug_metadata_dict)
        else:
            add_to_fasta(
                str(align_ref[0].seq), infection, seq_file, args.compression_type
            )
            write_metadata(metadata_file, infection, line_keys, args.compression_type)

    location_dict = json.loads(args.location)
    country = location_dict["country"]
    division = location_dict["division"]
    divisionAbbr = location_dict["divisionAbbr"]
    region = location_dict["region"]

    loop_counter=0

    def check_prefix(state):
        return state.startswith(args.input_graph_painted_prefix)
    def check_state(state):
        return state == args.input_graph_painted_state
    paint_this = check_state
    if args.input_graph_painted_prefix:
        paint_this = check_prefix

    transitions_to_paint = df[df["exit_state"].map(paint_this)]
    seed_transitions = transitions_to_paint["contact_pid"] == -1

    # Assign seed sequences
    seed = transitions_to_paint[seed_transitions]
    seed = seed["pid"].to_frame()
    # For each pid add sequence in align if keep adding form the begining
    # if the pid list is larger than align
    N = len(seed)
    M = len(align)
    if N <= M:
        seed["seq"] = [np.array(list(r.seq)) for r in align[:N]]
    else:
        align_str = [np.array(list(r.seq)) for r in align]
        seed["seq"] = list(islice(cycle(align_str), N))
    current_sequences = seed[["pid", "seq"]].set_index("pid")["seq"].to_dict()

    # Work on on remaining transitions
    transitions_to_paint = transitions_to_paint[~seed_transitions][
        ["pid", "contact_pid", "tick"]
    ]
    # Get only the seq_limit transitions
    if seq_limit:
        transitions_to_paint = transitions_to_paint.iloc[:seq_limit]

    # Convert tick to date relative to start date
    transitions_to_paint["date"] = transitions_to_paint["tick"].map(pd.offsets.Day)
    transitions_to_paint["date"] = (
        pd.to_datetime(start_date) + transitions_to_paint["date"]
    )
    thresh = np.array(thresh)
    sequence = next(current_sequences.values().__iter__())
    if not args.proportional:
        letters_to_use = np.array(["A", "C", "G", "T"])
        n_letters = len(letters_to_use)
        prob_matrix = np.repeat(
            np.array([1 / n_letters] * n_letters)[:, np.newaxis],
            repeats=len(sequence),
            axis=1,
        )
    else:
        letters_to_use = LETTERS
    cumulative_probs_matrix = np.cumsum(prob_matrix, axis=0)
    # Validate
    assert (
        len(sequence) == cumulative_probs_matrix.shape[1]
    ), "Sequence length must match the number of columns in the probability matrix"
    assert (
        len(letters_to_use) == cumulative_probs_matrix.shape[0]
    ), "Number of letters must match the number of rows in the probability matrix"

    for _, pid, contact_pid, date in transitions_to_paint[
        ["pid", "contact_pid", "date"]
    ].itertuples():
        loop_counter += 1
        if loop_counter % 1000 == 0:
            # print("    number of graph edges processed:  ",loop_counter)
            print("    number of infections decorated:  ", strain_id)
            # Get the parent's sequence, mutate it, and append result to fasta
            # print('Mutating sequence, adding to fasta.....')
        seq_to_change = current_sequences[contact_pid]
        change = determine_change(thresh, seq_to_change)
        # print(pid)
        if use_poor_mut_model:
            new_seq = poor_mut_model(seq_to_change)
        else:
            # new_seq = weight_change(seq_to_change, change, thresh_detail, use_proportional)
            # new_seq = weighted_change_np(seq_to_change, change, prob_matrix)
            new_seq = weighted_change(
                seq_to_change, change, cumulative_probs_matrix, letters=letters_to_use
            )
        current_sequences[pid] = new_seq

        new_seq = "".join(new_seq.tolist())

        infection = InfectionRecord()
        # get country, division, divisionAbbr, region from json parameter args.location
        infection.fromEpihiper(
            "ncov",
            region,
            country,
            division,
            division,
            date.strftime("%Y-%m-%d"),
            f"{country}/{divisionAbbr}-EHip-{strain_id}/{date.year}",
        )

        # Get augmented values for pid
        if augment_metadata:
            pid_df = persontrait_df.loc[pid]
            aug_metadata_dict = create_aug_metadata_dict(
                aug_metadata_columns, pid, pid_df
            )

            add_to_fasta(new_seq, infection, seq_file, args.compression_type)
            write_metadata(
                metadata_file,
                infection,
                line_keys,
                args.compression_type,
                aug_metadata_columns,
                aug_metadata_dict,
            )
        else:
            add_to_fasta(new_seq, infection, seq_file, args.compression_type)
            write_metadata(metadata_file, infection, line_keys, args.compression_type)
        strain_id += 1
    seq_file.close()
    metadata_file.close()

    print("Done")

    # calculates threshold for nucleotide change based on shannon
    # column entropy

    return


# ====================================
def column_entropy_thresh(freq_df):
    e_act = 0
    e_max = 0
    for i in range(len(freq_df.index)):
        p_xi = freq_df.iloc[i]
        e_act += p_xi * np.log(p_xi)

        p_xm = 1 / float(5)
        e_max += p_xm * np.log(p_xm)

    thresh = (1 - (e_act / e_max)) * 100
    # print(thresh)

    if np.isnan(thresh):
        thresh = 0

    return e_act, thresh


# ====================================
# Determine nucleotide change
def determine_change(thresh, seq_to_change):
    comparison_values = np.random.randint(0, 100, len(thresh))
    random_selection = comparison_values > thresh
    # If the position is a gap, then don't change it
    return random_selection & (seq_to_change != "-")


def determine_change_old(thresh):
    change = []
    for threshold in thresh:
        val = np.random.randint(0, 100)
        # the threshold is a measure of the consistency of the column
        # if the consistency is high, there should be less chance to change it
        # if the consistency is low it should be easier to change
        if val > threshold:
            change.append(True)
        else:
            change.append(False)
    return change


# Includes ambiguous nucleotides
LETTERS = np.array(
    ["A", "C", "G", "T", "N", "R", "K", "S", "Y", "M", "W", "B", "H", "D", "V"]
)


def weighted_change(sequence, change, cumulative_prob_matrix, letters=LETTERS):
    # Generate random numbers for each column
    random_values = np.random.rand(cumulative_prob_matrix.shape[1])

    # Compute cumulative probabilities for each column

    # Determine the index of the letter for each position
    # If our random number is higher than the cumulative probability, we take the next index
    letter_indices = (random_values[np.newaxis, :] < cumulative_prob_matrix).argmax(
        axis=0
    )
    new_sequence = letters[letter_indices]
    # Apply change based on change mask
    sequence = sequence.copy()
    sequence[change] = new_sequence[change]
    return sequence


def thresh_detail_to_prob_matrix(thresh_detail):
    # TODO: proportional
    letters = LETTERS
    prob_matrix = np.zeros((len(letters), len(thresh_detail)))
    filtered_thresh_detail = [
        df[df["letter"].isin(letters)].copy() for df in thresh_detail
    ]

    for col_index, df in enumerate(filtered_thresh_detail):
        df["change_value"] = df["change_value"].astype(float)
        total_count = df["change_value"].sum()
        for letter_index, letter in enumerate(letters):
            if letter in df["letter"].values:
                prob_matrix[letter_index, col_index] = (
                    df.loc[df["letter"] == letter, "change_value"].values[0]
                    / total_count
                )
            else:
                prob_matrix[letter_index, col_index] = 0

    return prob_matrix


def commit_change(index, change):
    new_seq = []
    for (nucleotide, change_val) in zip(index, change):
        if change_val:
            new_nucleotide = np.random.randint(1, 5)
            if new_nucleotide == 1:
                new_nucleotide = 'A'
                new_seq.append(new_nucleotide)
            if new_nucleotide == 2:
                new_nucleotide = 'T'
                new_seq.append(new_nucleotide)
            if new_nucleotide == 3:
                new_nucleotide = 'G'
                new_seq.append(new_nucleotide)
            if new_nucleotide == 4:
                new_nucleotide = 'C'
                new_seq.append(new_nucleotide)
            if new_nucleotide == 5:
                new_nucleotide = '-'
                new_seq.append(new_nucleotide)
        else:
            new_nucleotide = nucleotide
            new_seq.append(new_nucleotide)

    new_seq = ''.join(new_seq)
    return new_seq


def intermediate_mut_model(index, change):
    new_seq = []
    for (nucleotide, change_val) in zip(index, change):
        if change_val:
            if nucleotide == 'A':
                new_nucleotide = 'T'
                new_seq.append(new_nucleotide)
            if nucleotide == 'T':
                new_nucleotide = 'G'
                new_seq.append(new_nucleotide)
            if nucleotide == 'G':
                new_nucleotide = 'C'
                new_seq.append(new_nucleotide)
            if nucleotide == 'C':
                new_nucleotide = 'A'
                new_seq.append(new_nucleotide)
            else:
                new_seq.append(nucleotide)
        else:
            new_nucleotide = nucleotide
            new_seq.append(new_nucleotide)

    new_seq = ''.join(new_seq)
    return new_seq


def poor_mut_model(sequence):
    new_seq = []
    for nucleotide in sequence:
        change_val = np.random.randint(1, 2)
        if change_val == 1:
            new_nucleotide = np.random.randint(1, 5)
            if new_nucleotide == 1:
                new_nucleotide = 'A'
                new_seq.append(new_nucleotide)
            if new_nucleotide == 2:
                new_nucleotide = 'T'
                new_seq.append(new_nucleotide)
            if new_nucleotide == 3:
                new_nucleotide = 'G'
                new_seq.append(new_nucleotide)
            if new_nucleotide == 4:
                new_nucleotide = 'C'
                new_seq.append(new_nucleotide)
            if new_nucleotide == 5:
                new_nucleotide = '-'
                new_seq.append(new_nucleotide)
        else:
            new_nucleotide = nucleotide
            new_seq.append(new_nucleotide)

    new_seq = ''.join(new_seq)
    return new_seq


class InfectionRecord:
    def __init__(self):
        self.inf_dict = {
            "virus": None,
            "age": None,
            "country": None,
            "countryExposure": None,
            "date": None,
            "dateSubmitted": None,
            "died": None,
            "division": None,
            "divisionExposure": None,
            "fullyVaccinated": None,
            "strain": None,
            "gisaidCloade": None,
            "gisaidEpiIsl": None,
            "hospitalized": None,
            "host": None,
            "location": None,
            "month": None,
            "nextcladePangoLineage": None,
            "nextstrainClade": None,
            "originatingLab": None,
            "pangoLineage": None,
            "region": None,
            "regionExposure": None,
            "samplingStrategy": None,
            "sex": None,
            "sraAccession": None,
            "strainold": None,
            "submittingLab": None,
            "year":None
        }
    def get(self,key,default=None):
        return self.inf_dict.get(key,default)

    def fromEpihiper(self, virus, region, country, division, divisionExposure, date, strain):
        self.inf_dict["virus"] = virus
        self.inf_dict["region"] = country
        self.inf_dict["country"] = country
        self.inf_dict["division"] = division
        self.inf_dict["divisionExposure"] = divisionExposure
        self.inf_dict["date"] = date
        self.inf_dict["strain"] = strain

# Model the GISAID / Nextstrain metadata file for now
# https://docs.nextstrain.org/projects/ncov/en/latest/guides/data-prep/local-data.html
# virus,age,country,countryExposure,date,dateSubmitted,died,division,divisionExposure,fullyVaccinated,strain,gisaidClade,gisaidEpiIsl,hospitalized,host,location,month,nextcladePangoLineage,nextstrainClade,originatingLab,pangoLineage,region,regionExposure,samplingStrategy,sex,sraAccession,strainold,submittingLab,year
# ncov,,USA,USA,2021-09-20,2021-10-11,,Virginia,Virginia,,OK455686,,EPI_ISL_5088839,,Homo sapiens,,9,,21J,,AY.122,North America,North America,,,,USA/VA-CDC-LC0291093/2021,,2021

def write_metadata(metadata_file, infection, line_keys, compression_type, aug_metadata_columns=None, aug_metadata_dict=None):
    """
    Writes metadata to a file with optional compression and additional metadata columns.

    Parameters:
    metadata_file (file object): The file object to write the metadata to.
    infection (dict): A dictionary containing infection data.
    line_keys (list): A list of keys to extract from the infection dictionary.
    compression_type (str): The type of compression to use (e.g., 'XZ' for XZ compression).
    aug_metadata_columns (list, optional): A list of additional metadata columns to include. Defaults to None.
    aug_metadata_dict (dict, optional): A dictionary containing additional metadata values. Defaults to None.

    Returns:
    None
    """

    if aug_metadata_columns is None:
        meta_line = "\t".join([infection.get(key) for key in line_keys]) + "\n"
    else:
        meta_line = "\t".join([infection.get(key) for key in line_keys])
        for col in aug_metadata_columns:
            meta_line += "\t" + str(aug_metadata_dict[col])
        meta_line += "\n"

    if compression_type == XZ:
        metadata_file.write(meta_line.encode())
    else:
        metadata_file.write(meta_line)


def add_to_fasta(seq, infection, seq_file, compression_type): 
    seq_line = ">" + str(infection.get("strain")) + "\n" + seq + "\n"
    if compression_type == XZ:
        seq_file.write(seq_line.encode())
    else:
        seq_file.write(seq_line)


def find_seq(node):
    # dgaulton: going to parse this each time, need to modify this for
    # re-infections, could be two entries
    records = list(SeqIO.parse(fasta_to_write, "fasta"))
    rec = 0
    for i in range(0, len(records)):
        print(records[i].name)
        print(node)
        if records[i].name == str(node):
            rec = records[i]
    return rec


def dfs_edges_with_ticks(G, source=None, depth_limit=None):
    if source is None:
        # edges for all components
        nodes = G
    else:
        # edges for components with source
        nodes = [source]
    visited = set()
    if depth_limit is None:
        depth_limit = len(G)

    for start in nodes:
        tick = -1

        if start in visited:
            continue
        # modify visitied to include tick - TODO make sure this is treated as a
        # string and not as an int
        visited.add(start + "_" + tick)
        stack = [(start, depth_limit, iter(G[start]))]
        while stack:
            parent, depth_now, children = stack[-1]
            try:
                # need to confirm if iter for children in MultiDiGraph will
                # treat the node all at once or if will go edge by edge between
                # the neighbor - might need to flatten the results of
                # iter(G[start]) to treat a new edge with an old neighbor as
                # new - that would occur above at stack =
                child = next(children)
                if child not in visited:
                    yield parent, child
                    visited.add(child + "_" + tick)
                    if depth_now > 1:
                        stack.append((child, depth_now - 1, iter(G[child])))
            except StopIteration:
                stack.pop()


if __name__ == '__main__':
    begin_time = time.time()
    main()
    end_time = time.time()
    time_s=end_time-begin_time
    time_hr=(float)(time_s)/3600.0
    print("   Execution time (s), (hr): ",str(time_s),str(time_hr))
    print("   --- good termination ---")
