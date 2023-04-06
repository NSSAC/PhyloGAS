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

# Set these values to run the good or poor mutational model
# output_file_prefix = "test_new_metadata.good_mut_model"
output_file_prefix = sys.argv[1] if sys.argv[1] is not None else "simulation_output"
fasta_to_write = output_file_prefix + ".sequences.fasta"
metadata_file_to_write = output_file_prefix + ".metadata.tsv"
use_poor_mut_model = True if len(
    sys.argv) > 2 and sys.argv[2] == "-p" else False
seq_limit = 16521

# calculates threshold for nucleotide change based on shannon
# column entropy


def column_entropy_thresh(freq_df):
    e_act = 0
    e_max = 0
    for i in range(len(freq_df.index)):
        p_xi = freq_df.iloc[i]
        e_act += p_xi * np.log(p_xi)

        p_xm = 1 / len(freq_df.index)
        e_max += p_xm * np.log(p_xm)

    thresh = (1 - (e_act / e_max)) * 100
    print(thresh)

    if np.isnan(thresh):
        thresh = 0

    return thresh

# Determine nulceotide change


def determine_change(thresh):
    change = []
    for threshold in thresh:
        # import pdb; pdb.set_trace()
        val = np.random.randint(0, 100)
        if val < threshold:
            change.append('Yes')
        else:
            change.append('No')
    return change


def commit_change(index, change):
    new_seq = []
    for (nucleotide, change_val) in zip(index, change):
        if change_val == 'Yes':
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
        if change_val == 'Yes':
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


def add_to_fasta(seq, strain_id, date):
    # seq_file = open(fasta_to_write, "a")
    # metadata_file = open(metadata_file_to_write, "a")
    seq_file.write(">" + str(strain_id) + "\n" + seq + "\n")
    metadata_file.write(
        str(strain_id) +
        "\t" +
        str(date) +
        "\tncov\tNorth America\n")
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


##########################################################
# read in alignment to pandas dataframe
print('reading alignment file into pandas dataframe.....')
align = AlignIO.read('delta.fa', 'fasta')
name = []
description = []
for record in align:
    name.append(record.name)
    description.append(record.description)
align2 = pd.DataFrame(align)

# create threshold list, each column threshold included
print('calculating entropy and setting the threshold...')
thresh = []
df = pd.DataFrame()
for i in range(len(align2.columns)):
    import pdb
    pdb.set_trace()
    df1 = align2.iloc[:, i].value_counts()
    df1 = df1.divide(len(align2.index))
    thresh.append(column_entropy_thresh(df1) * 100)

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
metadata_file.write("strain\tdate\tvirus\tregion\n")

# need two data structures
# one list to just append to to generate MSA of all sequences as we go
# another dict that maps node to it's most recent sequence

# This assumes the data read into df is in chronological order (ascending
# according to tick)
current_sequences = {}
i = 0
max_seed_value_index = len(align2) - 1

start_date = "2021-05-31"  # start of Delta strain

strain_id = 0  # initialize strain_id for fasta, for now just increment a value, in the future could use node pid but would have to append to it in the case of multiple infections for a given pid

sequences_mutated = 0

for pid, contact_pid, tick, exit_state in zip(
        connections1, connections2, id1, id2):
    if exit_state == "var1E" and strain_id < seq_limit:
        print(tick)
        print(contact_pid)
        print(pid)
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
                new_seq = commit_change(seq_to_change, change)

            date = pd.to_datetime(start_date) + pd.DateOffset(days=tick)

            add_to_fasta(new_seq, strain_id, date)

            strain_id += 1

            current_sequences[pid] = new_seq

seq_file.close()
metadata_file.close()

print("Done")
