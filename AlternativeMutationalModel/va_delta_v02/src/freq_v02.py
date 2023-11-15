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


import argparse

# ====================================
# Constants.

# Analysis types
ENTROPY_ANALYSIS="entropy_analysis"
GEN_SEQUENCE_ANALYSIS="generate_sequence_analysis"
BOTH="both"

# Whether to write all threshold DFs to individual files (the previous method)
# or write all threshold DFs to one file.
INDIVIDUAL_FILES="individual_files"
ALL_IN_ONE_FILE="all_in_one_file"

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
    parser.add_argument("--align_fasta", type=str, default=None, nargs='?', dest="align_fasta", required=True, help="path to alignment file in FASTA format")
    parser.add_argument("--threshold_df_num_files", type=str, dest="threshold_df_num_files", required=True,
                        choices=[INDIVIDUAL_FILES, ALL_IN_ONE_FILE],
                        help="whether all threshold DFs get written to one file or individual files.")

    # For genomic sequences analysis.
    parser.add_argument("--start_date", default="2021-05-31", dest="start_date", required=False, type=str, help="simulation alignment to date")
    parser.add_argument("--input_graph_csv", type=str,dest="input_graph_csv",required=False, help="directed graph file; nodes are genomic sequences.")
    parser.add_argument("--output_prefix", default="syn_gen", type=str, dest="output_prefix", required=False, help="prefix for output file name (for fasta and metadata files)")
    parser.add_argument("--proportional", default=False, action="store_true", dest="proportional", required=False, help="use proportional letter choices")
    parser.add_argument("--poor", default=False, action="store_true", dest="poor", required=False, help="use poor mutational model")
    parser.add_argument("--limit", default=16521, type=int, dest="limit", required=False, help="maximum number of items to process")
    parser.add_argument("--reference", default=None, type=str, dest="reference", required=False, help="add reference sequence to the output")



    args = parser.parse_args()
    if (args.align_fasta == None):
        print("  Error.")
        print("  args.align_fasta has value None, which is not allowed.")
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
            fh_out.write("--------------\n")
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
    for itime in range(0,size01):
        filename = base_threshold_df+"_"+str(itime)+".csv"
        df_one = pd.read_csv(filename)
        thresh_detail.append(df_one)


    return thresh, thresh_detail


# ====================================
def main():


    args = getClas()

    analysis_type = args.analysis_type

    if analysis_type == BOTH or analysis_type==ENTROPY_ANALYSIS:
        # Compute the shannon entropies for the colummns of a
        # group of sequences.
        compute_entropy(args)

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
def generate_sequences(args):


    # Set these values to run the good or poor mutational model
    output_file_prefix = args.output_prefix
    fasta_to_write = output_file_prefix + ".sequences.fasta"
    metadata_file_to_write = output_file_prefix + ".metadata.tsv"
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
    connections1 = df.iloc[:, 1]  # pid
    connections2 = df.iloc[:, 3]  # contact_pid
    # cjk comment out; not used.
    # connections = pd.concat([connections1, connections2],
    #                         axis=1, ignore_index=False)

    # Create id dataframe (with tick and exit_state columns)
    id1 = df.iloc[:, 0]  # tick
    id2 = df.iloc[:, 2]  # exit_state
    # cjk:   ids = pd.concat([id1, id2], axis=1)  # dgaulton: looks like this is never used


    # Assigning ticks to correct nodes for coloring based on time
    timelist = dict(zip(connections1, id1))
    # dgaulton: where did this number come from? Is this from 20 tick cutoff?
    # - looks like picking out that item and moving it
    # cjk:  see comment immediately above by dgualton:  should this be hardcoded, or an input?
    pos = list(timelist.keys()).index(476724)
    items = list(timelist.items())
    items.insert(pos, (-1, 0))
    timelist = dict(items)
    colors_list = list(timelist.values())


    # draw directed graph network
    if False:
        print('creating network digraph .......')
        G = nx.from_pandas_edgelist(
            df,
            'contact_pid',
            'pid',
            create_using=nx.MultiDiGraph(),
            edge_key='tick')  # tick will be edge[2]

    ##########################################################

    # creating .fasta and .tsv files to append
    seq_file = open(fasta_to_write, 'w')
    metadata_file = open(metadata_file_to_write, 'w')

    # Prep metadata TSV file with required column names:
    # https://docs.nextstrain.org/projects/ncov/en/latest/guides/data-prep/local-data.html#required-metadata

    # need two data structures
    # one list to just append to to generate MSA of all sequences as we go
    # another dict that maps node to it's most recent sequence

    # kuhlman.  Need to read in fasta file again to get length.
    align = AlignIO.read(args.align_fasta, 'fasta')
    align2 = pd.DataFrame(align)


    # This assumes the data read into df is in chronological order (ascending
    # according to tick)
    current_sequences = {}
    i = 0
    max_seed_value_index = len(align2) - 1


    strain_id = 0  # initialize strain_id for fasta, for now just increment a value, in the future could use node pid but would have to append to it in the case of multiple infections for a given pid

    sequences_mutated = 0
    
    line_keys=["virus","region","country","division","divisionExposure","date","strain"]
    metadata_file.write("\t".join(line_keys)+"\n")

    if args.reference != None:
        align = AlignIO.read(args.reference, 'fasta')
        infection = InfectionRecord()
        country="China"
        division="Wuhan"
        divisionAbbr="Hu"
        region="Asia"
        date="2019-12-26"
        infection.fromEpihiper("ncov", region, country, division, division, date, "Wuhan-Hu-1/2019")
        add_to_fasta(str(align[0].seq), infection, seq_file, metadata_file, line_keys)

    loop_counter=0
    for pid, contact_pid, tick, exit_state in zip(
            connections1, connections2, id1, id2):
        # kuhlman:  to give indication of progress.
        loop_counter += 1
        if loop_counter%1000 == 0:
            print("    number of graph edges processed:  ",loop_counter)
        if exit_state == "var1E" and strain_id < seq_limit:
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
                country="USA"
                division="Virginia"
                divisionAbbr="VA"
                region="North America"
                infection.fromEpihiper("ncov", region, country, division, division, date.strftime("%Y-%m-%d"), f"{country}/{divisionAbbr}-EHip-{strain_id}/{date.year}")
                add_to_fasta(new_seq, infection, seq_file, metadata_file, line_keys)

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
        if change_val:
            letter_list=[]
            weight_list=[]
            if proportional:
                #python 3.7 order guaranteed but just in case
                for item in odds_val.items(): letter_list.append(item[0]), weight_list.append(item[1])
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
def add_to_fasta(seq, infection, seq_file, metadata_file, line_keys):
    # seq_file = open(fasta_to_write, "a")
    # metadata_file = open(metadata_file_to_write, "a")
    seq_file.write(">" + str(infection.get("strain")) + "\n" + seq + "\n")
    metadata_file.write("\t".join([infection.get(key) for key in line_keys])+"\n")
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

