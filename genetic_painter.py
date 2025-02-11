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

# ====================================
# Constants.
__version__ = '0.0.11'
# Analysis types
ENTROPY_ANALYSIS="entropy_analysis"
GEN_SEQUENCE_ANALYSIS="generate_sequence_analysis"
BOTH="both"

# Whether to write all threshold DFs to individual files (the previous method)
# or write all threshold DFs to one file.
INDIVIDUAL_FILES="individual_files"
ALL_IN_ONE_FILE="all_in_one_file"

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
    # parser.add_argument("align_fasta", type=str, default=None, nargs='?', help="path to alignment file in FASTA format")
    parser.add_argument("--align_fasta", type=str, default=None, nargs='?', dest="align_fasta", required=False, help="path to alignment file in FASTA format")
    parser.add_argument("--seed_fasta", type=str, default=None, nargs='?', dest="seed_fasta", required=False, help="path to seed file in FASTA format; defaults to align_fasta if not set")
    parser.add_argument("--threshold_df_num_files", type=str, dest="threshold_df_num_files", required=False,
                        choices=[INDIVIDUAL_FILES, ALL_IN_ONE_FILE],
                        help="whether all threshold DFs get written to one file or individual files.",default="ALL_IN_ONE_FILE")
    parser.add_argument("--random_number_seed", type=int, dest="random_number_seed", required=True, help="if < 0, then random assignment")

    # For genomic sequences analysis.
    parser.add_argument("--start_date", default="2021-05-31", dest="start_date", required=False, type=str, help="simulation alignment to date")
    parser.add_argument("--input_graph_csv", type=str,dest="input_graph_csv",required=False, help="directed graph file; nodes are infections.")
    parser.add_argument("--input_graph_painted_state", type=str,dest="input_graph_painted_state",default="var1E",required=False, help="Exit state that gets painted")
    parser.add_argument("--output_prefix", default="syn_gen", type=str, dest="output_prefix", required=False, help="prefix for output file name (for fasta and metadata files)")
    parser.add_argument("--proportional", default=True, action="store_true", dest="proportional", required=False, help="use proportional letter choices")
    parser.add_argument("--neutral", default=False, action="store_false", dest="proportional", required=False, help="use neutral letter choices")
    parser.add_argument("--poor", default=False, action="store_true", dest="poor", required=False, help="use poor mutational model")
    parser.add_argument("--limit", default=16521, type=int, dest="limit", required=False, help="maximum number of items to process")
    parser.add_argument("--reference", default=None, type=str, dest="reference", required=False, help="add reference sequence to the output")
    parser.add_argument("--compression", default=None, type=str, dest="compression_type", required=False, help="add compression method -- None, xz, or parquet",
                       choices=["None",XZ])
    parser.add_argument("--persontrait_file", default=None, type=str, dest="persontrait_file", required=False, help="the full path to the persontrait data file with additional data")
    parser.add_argument("--add_metadata", default=None, type=str, dest="add_metadata", required=False, help="the columns (comma-delimited) from the persontrait_file to include in the metadata output")
    #country="USA" division="Virginia" divisionAbbr="VA" region="North America"
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
def write_output_entropy(args, thresh, thresh_detail):

    # Filename and base filename.
    threshold_file = args.threshold_file
    base_threshold_df = args.base_threshold_df

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

    # Write out a DF to file; one DF for each threshold above.
    # ... or ...
    # put all DFs in one file.
    if (args.threshold_df_num_files==INDIVIDUAL_FILES):
        size_detail = len(thresh_detail)
        for itime in range(0, size_detail):
            df_the = thresh_detail[itime]
            filename = base_threshold_df+"_"+str(itime)+".csv"

            try:
                df_the.to_csv(filename)
            except:
                print("   Error")
                print("   Trying to write to a CSV file using a DF, where threshold DFs are to be written.")
                print("   This failed.")
                print("   CSV file name: ", filename)
                print("   Terminate.")
                exit(1)
    else:
        size_detail = len(thresh_detail)
        filename = base_threshold_df + ".csv"
        for itime in range(0, size_detail):
            df_the = thresh_detail[itime]
            if itime==0:
                fh_out = open(filename, "w")
            else:
                fh_out = open(filename, "a")
            fh_out.write("+-------------\n")
            fh_out.close()

            try:
                df_the.to_csv(filename,mode="a")
            except:
                print("   Error")
                print("   Trying to write to a CSV file using a DF, where threshold DFs are to be written.")
                print("   This failed.")
                print("   CSV file name: ", filename)
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

    # Read in the len(thresh) number of dataframes; put into list.
    size01 = len(thresh)
    if (args.threshold_df_num_files==INDIVIDUAL_FILES):
        for itime in range(0,size01):
            filename = base_threshold_df+"_"+str(itime)+".csv"
            df_one = pd.read_csv(filename)
            thresh_detail.append(df_one)
    else:
        # All data in one file.
        filename = base_threshold_df + ".csv"
        # create an Empty DataFrame object
        df_one = pd.DataFrame(data=None, columns=['letter','change_value'])
        fh_in = open(filename,"r")
        # Read the first line just to get rid of it.
        dash_string = fh_in.readline()
        for aline in fh_in:
            sline = aline.strip()
            if sline[0] == "+":
                # Found next entry, so stop entering into this DF.
                # Add this DF to list.
                thresh_detail.append(df_one)
                # Create an Empty DataFrame object
                df_one = pd.DataFrame(data=None, columns=['letter', 'change_value'])
            else:
                tokens = sline.split(",")
                # df_one.append([tokens[0], tokens[1]])
                df_one.loc[len(df_one)] = [tokens[0], tokens[1]]
        # The last DF needs to be added to list.
        thresh_detail.append(df_one)

        # for line in finalText.splitlines():
        #     print(line)
        #     m = re.findall(r'\w+', line)
        #     print(m)
        #     matches = re.findall(r'\w+', line)
        #     df.loc[len(df)] = [matches[1], matches[6]]
        #     df.loc[len(df)] = [matches[9], matches[14]]

    return thresh, thresh_detail


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

    if analysis_type == BOTH and args.threshold_df_num_files == ALL_IN_ONE_FILE:
        # Have to put CSV extension on the file with threshold DFs, in this case.
        # This is so the filename is well-formed when opening to read contents.
        print("  \n\n --- doing next sequence calculations --- \n\n")
        args.base_threshold_df = args.base_threshold_df + ".csv"

    if analysis_type == BOTH or analysis_type==GEN_SEQUENCE_ANALYSIS:
        # Determine perturbations in a series of sequences.
        generate_sequences(args)


    return

# ====================================
def compute_entropy(args):
    """
    Compute the entropy for each column (i.e., each position)
    of a collection of genomic sequences.
    :param args:  the CLAs.
    :return:
    """

    start_date = args.start_date  # start of Delta strain

    # cjk:  move this to the sequence-generating method.
    # # Set these values to run the good or poor mutational model
    # # output_file_prefix = "test_new_metadata.good_mut_model"
    # output_file_prefix = args.output_prefix
    # fasta_to_write = output_file_prefix + ".sequences.fasta"
    # metadata_file_to_write = output_file_prefix + ".metadata.tsv"
    # use_poor_mut_model = args.poor
    # use_proportional = args.proportional
    # seq_limit = args.limit

    ##########################################################
    # read in alignment to pandas dataframe
    print('reading alignment file into pandas dataframe.....')
    align = AlignIO.read(args.align_fasta, 'fasta')
    name = []
    description = []
    for record in align:
        name.append(record.name)
        description.append(record.description)
    align2 = pd.DataFrame(align)

    # create threshold list, each column threshold included
    print('calculating entropy and setting the threshold...')
    thresh = []
    thresh_detail=[]
    df = pd.DataFrame()
    for i in range(len(align2.columns)):
        df1 = align2.iloc[:, i].value_counts(normalize=True)
        #df1 = df1.divide(len(align2.index))
        thresh_detail.append(df1)
        thresh.append(column_entropy_thresh(df1))

    ##########################################################

    write_output_entropy(args, thresh, thresh_detail)

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

    # Load thresholds into list.
    # Load the dataframe for each threshold.
    thresh, thresh_detail = load_thresholds_and_dfs(args)


    # Read in network data
    print('reading in the network data....')
    # df = pd.read_csv('contact_network_va_delta.csv')
    df = pd.read_csv(input_graph_csv)


    # Create connection dataframe (pid, and contact_pid columns)
    connections1 = df.loc[:, "pid"]  # pid
    connections2 = df.loc[:, "contact_pid"]  # contact_pid
    # cjk comment out; not used.
    # connections = pd.concat([connections1, connections2],
    #                         axis=1, ignore_index=False)

    # Create id dataframe (with tick and exit_state columns)
    id1 = df.loc[:, "tick"]  # tick
    id2 = df.loc[:, "exit_state"]  # exit_state

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

    # kuhlman.  Need to read in fasta file again to get length.
    if args.seed_fasta == None:
        align = AlignIO.read(args.align_fasta, 'fasta')
    else:
        align = AlignIO.read(args.seed_fasta, 'fasta')
    
    align2 = pd.DataFrame(align)


    # This assumes the data read into df is in chronological order (ascending
    # according to tick)
    current_sequences = {}
    i = 0
    max_seed_value_index = len(align2) - 1


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
        align = AlignIO.read(args.reference, 'fasta')
        infection = InfectionRecord()
        country=ref_location_dict['country']
        division=ref_location_dict['division']
        divisionAbbr=ref_location_dict['divisionAbbr']
        region=ref_location_dict['region']
        date=ref_location_dict['date']
        infection.fromEpihiper("ncov", region, country, division, division, date, f"{division}-{divisionAbbr}-1/{date.split('-')[0]}")
        
        if augment_metadata:
            aug_metadata_dict = create_aug_metadata_dict(aug_metadata_columns,pid=-1)
            add_to_fasta(str(align[0].seq), infection, seq_file, metadata_file, line_keys, args.compression_type, aug_metadata_columns, aug_metadata_dict)
        else:
            add_to_fasta(str(align[0].seq), infection, seq_file, metadata_file, line_keys, args.compression_type)
    
    location_dict = json.loads(args.location)
    country = location_dict["country"]
    division = location_dict["division"]
    divisionAbbr = location_dict["divisionAbbr"]
    region = location_dict["region"]

    loop_counter=0
    painted_exit_state = args.input_graph_painted_state
    for pid, contact_pid, tick, exit_state in zip(
            connections1, connections2, id1, id2):
        # kuhlman:  to give indication of progress.
        loop_counter += 1
        if loop_counter%1000 == 0:
            #print("    number of graph edges processed:  ",loop_counter)
            print("    number of infections decorated:  ",strain_id)
        if exit_state == painted_exit_state and strain_id < seq_limit:
#            print(tick)
#            print(contact_pid)
#            print(pid)
#            print(exit_state)
            if contact_pid == -1:  # seed case
                # grab a new real sequence
                #print('Adding seed seq to .fasta ........')
                index = align2.iloc[i].values.tolist()
                index = ''.join(index)
                if i == max_seed_value_index:  # wrap to the beginning of seed sequences if we've run out
                    i = 0
                else:
                    i += 1

                # update the node's current sequence
                current_sequences[pid] = index

            else:
                # Get the parent's sequence, mutate it, and append result to fasta
                #print('Mutating sequence, adding to fasta.....')
                seq_to_change = current_sequences[contact_pid]
                change = determine_change(thresh)
                #print(pid)
                if (use_poor_mut_model):
                    new_seq = poor_mut_model(seq_to_change)
                else:
                    new_seq = weight_change(seq_to_change, change, thresh_detail, use_proportional)

                date = pd.to_datetime(start_date) + pd.DateOffset(days=tick)
                infection = InfectionRecord()
                #get country, division, divisionAbbr, region from json parameter args.location
                infection.fromEpihiper("ncov", region, country, division, division, date.strftime("%Y-%m-%d"), f"{country}/{divisionAbbr}-EHip-{strain_id}/{date.year}")

                # Get augmented values for pid
                if augment_metadata:
                    pid_df = persontrait_df.loc[pid]
                    aug_metadata_dict = create_aug_metadata_dict(aug_metadata_columns,pid,pid_df)

                if augment_metadata:
                    add_to_fasta(new_seq, infection, seq_file, metadata_file, line_keys, args.compression_type, aug_metadata_columns, aug_metadata_dict)
                else:
                    add_to_fasta(new_seq, infection, seq_file, metadata_file, line_keys, args.compression_type)

                strain_id += 1

                current_sequences[pid] = new_seq

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

        #p_xm = 1 / len(freq_df.index)
        p_xm = 1 / float(5)
        e_max += p_xm * np.log(p_xm)

    thresh = (1 - (e_act / e_max)) * 100
    #print(thresh)

    if np.isnan(thresh):
        thresh = 0

    return thresh


# ====================================
# Determine nulceotide change
def determine_change(thresh):
    change = []
    for threshold in thresh:
        val = np.random.randint(0, 100)
        #the threshold is a measure of the consistency of the column
        #if the consistency is high, there should be less chance to change it
        #if the consistency is low it should be easier to change
        if val > threshold:
            change.append(True)
        else:
            change.append(False)
    return change


def weight_change(index, change, letter_odds, proportional=False):
    new_seq = []

    for (nucleotide, change_val, odds_val) in zip(index, change, letter_odds):
        if change_val and nucleotide != '-':
            letter_list=[]
            weight_list=[]
            if proportional:
                #python 3.7 order guaranteed but just in case
                for key, item in odds_val.iterrows():
                    if item["change_value"] == "proportion" or item["letter"] not in set(['A','C','G','T']):
                        continue
                    letter_list.append(item["letter"])
                    weight_list.append(float(item["change_value"]))
                new_nucleotide = random.choices(letter_list, weights=weight_list,k=1)[0]
            else:
                #set the letters to equal weight except for the gap symbol.
                letter_list=["A","C","T","G"]
                gap_odds=odds_val.get("-",0)
                weight_list = [(1-gap_odds)/float(len(letter_list))] * len(letter_list) #make four copies
                weight_list.append(gap_odds)
                letter_list.append("-")
                new_nucleotide = random.choices(letter_list, weights=weight_list,k=1)[0]
        else:
            new_nucleotide = nucleotide
        new_seq.append(new_nucleotide)
    new_seq = ''.join(new_seq)
    return new_seq


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

#Model the GISAID / Nextstrain metadata file for now
#https://docs.nextstrain.org/projects/ncov/en/latest/guides/data-prep/local-data.html
#virus,age,country,countryExposure,date,dateSubmitted,died,division,divisionExposure,fullyVaccinated,strain,gisaidClade,gisaidEpiIsl,hospitalized,host,location,month,nextcladePangoLineage,nextstrainClade,originatingLab,pangoLineage,region,regionExposure,samplingStrategy,sex,sraAccession,strainold,submittingLab,year
#ncov,,USA,USA,2021-09-20,2021-10-11,,Virginia,Virginia,,OK455686,,EPI_ISL_5088839,,Homo sapiens,,9,,21J,,AY.122,North America,North America,,,,USA/VA-CDC-LC0291093/2021,,2021
def add_to_fasta(seq, infection, seq_file, metadata_file, line_keys, compression_type, aug_metadata_columns=None, aug_metadata_dict=None):
    # seq_file = open(fasta_to_write, "a")
    # metadata_file = open(metadata_file_to_write, "a")
    if aug_metadata_columns is None:
        meta_line="\t".join([infection.get(key) for key in line_keys])+"\n"
    else:
        meta_line="\t".join([infection.get(key) for key in line_keys])
        for col in aug_metadata_columns:
            meta_line += "\t" + str(aug_metadata_dict[col])
        meta_line += "\n"

    if compression_type == XZ:
        seq_line=">" + str(infection.get("strain")) + "\n" + seq + "\n"
        seq_file.write(seq_line.encode())
        # meta_line="\t".join([infection.get(key) for key in line_keys])+"\n"
        metadata_file.write(meta_line.encode())
    else:
        seq_file.write(">" + str(infection.get("strain")) + "\n" + seq + "\n")
        metadata_file.write(meta_line)
    # seq_file.close()
    # metadata_file.close()


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
                # tick = child.
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

