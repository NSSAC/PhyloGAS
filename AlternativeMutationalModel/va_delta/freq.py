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


import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_prefix", default="syn_gen", type=str, help="prefix for output file name")
    parser.add_argument("--proportional", default=False, action="store_true", help="use proportional letter choices")
    parser.add_argument("--poor", default=False, action="store_true", help="use poor mutational model")
    parser.add_argument("--limit", default=16521, type=int, help="maximum number of items to process")
    parser.add_argument("--start_date", default="2021-05-31", type=str, help="simulation alignment to date")
    parser.add_argument("--reference", default=None, type=str, help="add reference sequence to the output")
    parser.add_argument("align_fasta", type=str, default=None, nargs='?', help="path to alignment file in FASTA format")
    

    args = parser.parse_args()
    if (args.align_fasta == None):
        parser.print_help()
        sys.exit(0)

    print(args.output_prefix)
    print(args.proportional)
    print(args.poor)
    
    start_date = args.start_date  # start of Delta strain

    # Set these values to run the good or poor mutational model
    # output_file_prefix = "test_new_metadata.good_mut_model"
    output_file_prefix = args.output_prefix
    fasta_to_write = output_file_prefix + ".sequences.fasta"
    metadata_file_to_write = output_file_prefix + ".metadata.tsv"
    use_poor_mut_model = args.poor
    use_proportional = args.proportional
    seq_limit = args.limit

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

    # Read in network data
    print('reading in the network data....')
    df = pd.read_csv('contact_network_va_delta.csv')


    # Create connection dataframe (pid, and contact_pid columns)
    connections1 = df.iloc[:, 1]  # pid
    connections2 = df.iloc[:, 3]  # contact_pid
    connections = pd.concat([connections1, connections2],
                            axis=1, ignore_index=False)

    # Create id dataframe (with tick and exit_state columns)
    id1 = df.iloc[:, 0]  # tick
    id2 = df.iloc[:, 2]  # exit_state
    ids = pd.concat([id1, id2], axis=1)  # dgaulton: looks like this is never used


    # Assigning ticks to correct nodes for coloring based on time
    timelist = dict(zip(connections1, id1))
    # dgaulton: where did this number come from? Is this from 20 tick cutoff?
    # - looks like picking out that item and moving it
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
        infection.fromEpiHiper("ncov", region, country, division, division, date, "Wuhan-Hu-1/2019")
        add_to_fasta(align.seq, infection, seq_file, metadata_file, line_keys)

    for pid, contact_pid, tick, exit_state in zip(
            connections1, connections2, id1, id2):
        if exit_state == "var1E" and strain_id < seq_limit:
#            print(tick)
#            print(contact_pid)
#            print(pid)
            print(exit_state)
            if contact_pid == -1:  # seed case
                # grab a new real sequence
                print('Adding seed seq to .fasta ........')
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
                print('Mutating sequence, adding to fasta.....')
                seq_to_change = current_sequences[contact_pid]
                change = determine_change(thresh)
                print(pid)
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
                infection.fromEpiHiper("ncov", region, country, division, division, date, f"{country}/{divisionAbbr}-EHip-{strain_id}/{date.year}")
                add_to_fasta(new_seq, infection, seq_file, metadata_file, line_keys)

                strain_id += 1

                current_sequences[pid] = new_seq

    seq_file.close()
    metadata_file.close()

    print("Done")

    # calculates threshold for nucleotide change based on shannon
    # column entropy


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
    def fromEpihiper(self, virus, region, country, division, divisionExposure, date, strain):
        self.my_dict["virus"] = virus
        self.my_dict["region"] = country
        self.my_dict["country"] = country
        self.my_dict["division"] = division
        self.my_dict["divisionExposure"] = divisionExposure
        self.my_dict["date"] = date
        self.my_dict["strain"] = strain

#Model the GISAID / Nextstrain metadata file for now
#https://docs.nextstrain.org/projects/ncov/en/latest/guides/data-prep/local-data.html
#virus,age,country,countryExposure,date,dateSubmitted,died,division,divisionExposure,fullyVaccinated,strain,gisaidClade,gisaidEpiIsl,hospitalized,host,location,month,nextcladePangoLineage,nextstrainClade,originatingLab,pangoLineage,region,regionExposure,samplingStrategy,sex,sraAccession,strainold,submittingLab,year
#ncov,,USA,USA,2021-09-20,2021-10-11,,Virginia,Virginia,,OK455686,,EPI_ISL_5088839,,Homo sapiens,,9,,21J,,AY.122,North America,North America,,,,USA/VA-CDC-LC0291093/2021,,2021
def add_to_fasta(seq, infection, seq_file, metadata_file, line_keys):
    # seq_file = open(fasta_to_write, "a")
    # metadata_file = open(metadata_file_to_write, "a")
    seq_file.write(">" + str(strain_id) + "\n" + seq + "\n")
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
    main()
