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
import re
import json
from Bio import bgzf
import gzip # Added for aligned_to_df

import argparse
from itertools import islice, cycle
import math # Added for rate limiting

# Works both as an installed package (`phylogas.painter.mutational_models`) and
# when this file is run directly from its own directory, which the existing
# run.03.vadelta.* scripts still do.
try:
    from .mutational_models import registry as model_registry
    from .. import ids as _ids
except ImportError:
    # Run directly from this directory: mutational_models/ is a subdirectory,
    # but ids.py lives one level up in the package, so put it on the path.
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from mutational_models import registry as model_registry
    import ids as _ids

# Force line-buffered stdout/stderr so progress is visible in SLURM logs as it
# happens. Without this, stdout is block-buffered (4-8 KB) when redirected to a
# file and a long-running job appears to hang with an empty log.
try:
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)
except AttributeError:  # pragma: no cover - Python < 3.7
    pass

# ====================================
# Constants.
__version__ = '0.0.12' # Incremented version
# Analysis types
ENTROPY_ANALYSIS="entropy_analysis"
GEN_SEQUENCE_ANALYSIS="generate_sequence"
BOTH="both"

# Compression types
XZ="xz"
BGZF="bgzf"
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
    parser.add_argument("--compression", default=None, type=str, dest="compression_type", required=False, 
                        help="add compression method -- None, xz, bgzf, or parquet",
                        choices=["None", XZ, BGZF]) # Added BGZF here
    parser.add_argument("--compression_level", default=1, type=int, dest="compression_level", required=False,
                        choices=range(0, 10), metavar="[0-9]",
                        help="Compression level. For xz this is the LZMA preset, for bgzf the zlib level. "
                             "Default 1. NOTE: the library default (6) is ~40x slower than 1 on this "
                             "workload for <5%% size benefit and is the usual cause of multi-hour runs.")
    parser.add_argument("--compression_threads", default=0, type=int, dest="compression_threads", required=False,
                        help="If >0 and --compression xz, pipe output through the external multi-threaded "
                             "`xz -T<n>` binary instead of Python's single-threaded lzma module. "
                             "0 (default) uses in-process lzma.")
    parser.add_argument("--persontrait_file", default=None, type=str, dest="persontrait_file", required=False, help="the full path to the persontrait data file with additional data")
    parser.add_argument("--add_metadata", default=None, type=str, dest="add_metadata", required=False, help="the columns (comma-delimited) from the persontrait_file to include in the metadata output")
    parser.add_argument("--location", default='{"country":"USA","division":"Virginia","divisionAbbr":"VA","region":"North America"}', type=str, dest="location", required=False, help="the location data for the infection record")
    parser.add_argument("--reference_location", default='{"country":"China","division":"Wuhan","divisionAbbr":"Hu","region":"Asia","date":"2019-12-26"}', type=str, dest="reference_location", required=False, help="the location data for the reference infection record")
    
    # START ADDED ARGUMENTS FOR RATE LIMITING
    parser.add_argument("--linelist_filter", "--linelist-filter", dest="linelist_filter",
                        action="append", default=None, metavar="[LABEL=]FILE",
                        help="Also write a FASTA/metadata pair containing only the "
                             "infections named in FILE. Repeatable, so one run can emit "
                             "the full tree plus any number of sampled subsets (e.g. one "
                             "per sampling strategy). FILE may be a TwinSampler linelist "
                             "or a BeyondBaseline samples file; identifiers are read from "
                             "alias_pid, infection_id, strain, sim_pid or pid. Without "
                             "this flag a single unrestricted output is written, as before.")
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
def _open_text_maybe_compressed(path):
    """Open a text file that may be plain, .gz or .xz."""
    p = str(path)
    if p.endswith(".xz"):
        return lzma.open(p, "rt")
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    return open(p, "r")


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
        
    return np.array(thresh), prob_matrix


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
        # Get raw counts, including Ns, gaps, and ambiguity codes
        raw_counts = align2.iloc[:, i].value_counts()
        total_seqs = raw_counts.sum()
        
        # 1. Filter down to ONLY canonical bases
        canonical_bases = ['A', 'C', 'G', 'T']
        acgt_counts = raw_counts[raw_counts.index.isin(canonical_bases)]
        
        # 2. Calculate the "Valid Fraction" multiplier
        # What percentage of the column is a valid, canonical base?
        valid_fraction = acgt_counts.sum() / total_seqs if total_seqs > 0 else 0
        
        # 3. Normalize the ACGT counts to calculate pure biological entropy
        if acgt_counts.sum() > 0:
            acgt_probs = acgt_counts / acgt_counts.sum()
            e_act, raw_thresh = column_entropy_thresh(acgt_probs)
        else:
            # If there are NO canonical bases (e.g., all Ns or gaps)
            e_act, raw_thresh = 0, 100.0 
            # Provide a uniform fallback distribution for the prob matrix
            acgt_probs = pd.Series([0.25, 0.25, 0.25, 0.25], index=["A", "C", "G", "T"])
            
        # 4. Scale the weight by the valid fraction
        # Since weight = 1.0 - (thresh/100), we calculate the base weight, apply the 
        # valid fraction penalty, and then convert it back into a threshold format.
        base_weight = max(0.0, 1.0 - (raw_thresh / 100.0))
        penalized_weight = base_weight * valid_fraction
        final_thresh = (1.0 - penalized_weight) * 100.0
        
        # Store the clean, ACGT-only probability distribution
        acgt_probs.name = i
        thresh_detail_dfs.append(acgt_probs)
        thresh.append(final_thresh)
        entropy_values.append(e_act)
        
    # Create a matrix of probabilities for each letter at each position,
    # fill missing with 0 and reindex to maintain order.
    # Because we only passed ACGT probs, all Ns, gaps, etc. will become exactly 0.0!
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
class _ExternalCompressor:
    """Writes through an external multi-threaded compressor (e.g. `xz -T8`).

    Python's `lzma` module is single-threaded and releases the GIL only in
    coarse chunks, so it serialises the whole generation loop. Handing the
    bytes to a separate `xz` process lets compression overlap with sequence
    generation and use multiple cores.
    """

    def __init__(self, path, argv):
        self._fh = open(path, "wb")
        self._proc = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=self._fh)
        self.closed = False

    def write(self, data):
        self._proc.stdin.write(data)

    def close(self):
        if self.closed:
            return
        self.closed = True
        self._proc.stdin.close()
        rc = self._proc.wait()
        self._fh.close()
        if rc != 0:
            raise RuntimeError(f"external compressor exited with status {rc}")


def _open_output_writers(args, fasta_to_write, metadata_file_to_write):
    """Open the FASTA and metadata handles honouring compression type/level.

    Compression is by far the dominant cost of the generation loop. Measured on
    29,903 bp SARS-CoV-2 records:

        xz preset 6 (lzma default)   ~243 records/s
        xz preset 1                ~9,000 records/s
        xz preset 0               ~10,500 records/s
        bgzf/gzip level 6            ~186 records/s
        bgzf/gzip level 1          ~3,100 records/s
        uncompressed              ~56,000 records/s

    At 5.35M records the library default (preset 6) alone costs ~6 hours of
    wall time, which is why the job appears to hang. Level 1 costs ~10 minutes
    for well under 5% extra output size on this highly redundant data.
    """
    level = args.compression_level
    threads = getattr(args, "compression_threads", 0)

    if args.compression_type == XZ:
        if threads > 0:
            xz_argv_base = ["xz", f"-T{threads}", f"-{level}", "-c"]
            print(f"  Compression: external `{' '.join(xz_argv_base)}` (multi-threaded)")
            seq_file = _ExternalCompressor(fasta_to_write, xz_argv_base)
            metadata_file = _ExternalCompressor(metadata_file_to_write, xz_argv_base)
        else:
            print(f"  Compression: in-process lzma, preset={level}")
            seq_file = lzma.open(fasta_to_write, 'wb', preset=level)
            metadata_file = lzma.open(metadata_file_to_write, 'wb', preset=level)
    elif args.compression_type == BGZF:
        print(f"  Compression: bgzf/gzip, level={level}")
        # Use Biopython's BGZF writer for the FASTA
        seq_file = bgzf.BgzfWriter(fasta_to_write, 'wb', compresslevel=level)
        # Standard gzip is fine for the metadata TSV
        metadata_file = gzip.open(metadata_file_to_write, 'wb', compresslevel=level)
    else:
        seq_file = open(fasta_to_write, 'w')
        metadata_file = open(metadata_file_to_write, 'w')

    return seq_file, metadata_file


# ====================================
class _OutputSet:
    """One FASTA + metadata pair, optionally restricted to a set of infections.

    The painter always walks the whole transmission tree -- a child's genome is
    derived from its parent's, so nothing can be skipped during computation.
    These objects decide only what reaches disk.

    With no --linelist-filter, a single unrestricted set is created and
    behaviour is exactly as before. Each --linelist-filter file adds another
    set that writes only the infections named in it, so one pass can emit the
    full tree plus any number of sampled subsets.
    """

    def __init__(self, label, prefix, args, keys=None, source=None):
        self.label = label
        self.prefix = prefix
        self.keys = keys            # None => write everything
        self.source = source
        self.written = 0
        fasta, meta = _output_paths(prefix, args.compression_type)
        self.fasta_path, self.meta_path = fasta, meta
        self.seq_file, self.metadata_file = _open_output_writers(args, fasta, meta)

    def wants(self, infection_id):
        return self.keys is None or infection_id in self.keys

    def close(self):
        for fh in (self.seq_file, self.metadata_file):
            try:
                fh.close()
            except Exception:
                pass


def _output_paths(prefix, compression_type):
    if compression_type is None or compression_type == "None":
        return prefix + ".sequences.fasta", prefix + ".metadata.tsv"
    if compression_type == XZ:
        return prefix + ".sequences.fasta.xz", prefix + ".metadata.tsv.xz"
    if compression_type == BGZF:
        return prefix + ".sequences.fasta.gz", prefix + ".metadata.tsv.gz"
    return prefix + ".sequences.fasta", prefix + ".metadata.tsv"


# Columns that can identify an infection in a filter file, in priority order.
# alias_pid / infection_id are "{pid}.{tick}" and match exactly. strain is the
# painter's own ID and is parsed back. A bare pid matches every infection of
# that person, which is coarser but sometimes what you have.
_FILTER_KEY_COLUMNS = ("alias_pid", "infection_id", "strain", "sim_pid", "pid")


def _read_filter_keys(path):
    """Read a linelist/samples file and return (infection_ids, pids, column).

    Accepts anything pandas can read, compressed or not. Returns exact
    infection ids where available, plus a fallback set of bare pids.
    """
    df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    col = next((c for c in _FILTER_KEY_COLUMNS if c in df.columns), None)
    if col is None:
        sys.exit(
            f"  Error: {path} has none of the recognised identifier columns.\n"
            f"         Looked for: {', '.join(_FILTER_KEY_COLUMNS)}\n"
            f"         Found: {', '.join(list(df.columns)[:10])}"
        )

    vals = df[col].dropna().astype(str).str.strip()
    infection_ids, pids = set(), set()
    if col == "strain":
        # USA/VA-EHip-{pid}.{tick}/{year}
        for v in vals:
            if "EHip-" in v:
                infection_ids.add(v.split("EHip-")[1].rsplit("/", 1)[0])
    elif col in ("alias_pid", "infection_id"):
        for v in vals:
            (infection_ids if "." in v else pids).add(v)
    else:
        pids.update(v.replace(".0", "") for v in vals)
    return infection_ids, pids, col


def _build_output_sets(args):
    """Create the output sets implied by --output_prefix and --linelist-filter."""
    sets = [_OutputSet("all", args.output_prefix, args)]

    for spec in (args.linelist_filter or []):
        # "label=path" or just "path" (label derived from the filename)
        if "=" in spec and not os.path.exists(spec):
            label, path = spec.split("=", 1)
        else:
            path = spec
            label = re.sub(r"(\.csv|\.tsv)?(\.xz|\.gz)?$", "", os.path.basename(path))
            label = re.sub(r"_samples$", "", label)
        if not os.path.exists(path):
            sys.exit(f"  Error: --linelist-filter file not found: {path}")

        ids, pids, col = _read_filter_keys(path)
        if not ids and not pids:
            print(f"  Warning: {path} yielded no identifiers; skipping", file=sys.stderr)
            continue
        keys = ids if ids else None
        s = _OutputSet(label, f"{args.output_prefix}.{label}", args, keys=keys, source=path)
        s.pid_keys = pids if not ids else None
        sets.append(s)
        n = len(ids) if ids else len(pids)
        print(f"  Filter '{label}': {n:,} identifiers from column '{col}' ({path})")
    return sets


# ====================================
def generate_sequences(args):

    output_file_prefix = args.output_prefix

    # Create the output directory if the prefix names one that does not exist
    # yet. Without this the run dies only after loading the network and the
    # seed FASTA, which can be several minutes of wasted work.
    _out_parent = os.path.dirname(os.path.abspath(output_file_prefix))
    if _out_parent:
        os.makedirs(_out_parent, exist_ok=True)

    augment_metadata = False
    # These are passed to create_infection_record unconditionally, so they must
    # exist even when metadata augmentation is off or fails to load.
    persontrait_df = None
    aug_metadata_columns = []
    standardized_aug_cols = []
    standardized_replacements = {}
    if args.persontrait_file and args.add_metadata:
        augment_metadata = True
        aug_metadata_columns = args.add_metadata.split(",")

        # This file must be the DERIVED demographics table, not the raw
        # EpiHiper persontrait file. The v2.4.0 persontrait opens with a JSON
        # schema line and lacks county/home_latitude/home_longitude/latino
        # entirely. Build the right file with:
        #     phylogas build-demographics --persontrait ... --person ... \
        #         --household ... --residence ... --fips ... --out ...
        #
        # Failures here are FATAL rather than a silent downgrade: previously a
        # bad file disabled augmentation and produced a complete-looking run
        # whose demographic columns were all empty.
        try:
            with _open_text_maybe_compressed(args.persontrait_file) as _fh:
                _first = _fh.readline().lstrip()
            _skip = 1 if _first.startswith("{") else 0
            if _skip:
                print(f"  Note: {args.persontrait_file} begins with an EpiHiper JSON schema "
                      f"line; skipping it.")
                print( "        (This looks like a RAW persontrait file. The painter expects the "
                       "derived")
                print( "         demographics table - see `phylogas build-demographics`.)")
            persontrait_df = pd.read_csv(args.persontrait_file, skiprows=_skip).set_index("pid")
        except FileNotFoundError:
            sys.exit(f"  Error: persontrait_file {args.persontrait_file} not found.")
        except KeyError:
            sys.exit(
                f"  Error: no 'pid' column in {args.persontrait_file}.\n"
                f"         Build the demographics table with `phylogas build-demographics`."
            )

        # Verify the requested columns actually exist before painting millions
        # of records. The painter renames a few on the way in.
        _alias = {"gender": "sex", "home_latitude": "latitude", "home_longitude": "longitude"}
        _have = set(persontrait_df.columns)
        _missing = [c for c in aug_metadata_columns
                    if c not in _have and _alias.get(c, c) not in _have]
        if _missing:
            sys.exit(
                f"  Error: --add_metadata requested columns that are not in\n"
                f"         {args.persontrait_file}:\n"
                f"           missing : {', '.join(_missing)}\n"
                f"           present : {', '.join(sorted(_have)[:12])}\n"
                f"         If this is a raw EpiHiper persontrait file, build the derived\n"
                f"         demographics table first:\n"
                f"           phylogas build-demographics --persontrait <pt> --person <p> \\\n"
                f"               --household <hh> --residence <rl> --fips <fips> --out <out>"
            )
    elif args.persontrait_file or args.add_metadata:
        print("   Info: persontrait_file and add_metadata must BOTH be provided to augment metadata. Not augmenting.")

    use_poor_mut_model = args.poor
    # use_proportional = args.proportional # This is used to select letters_to_use
    seq_limit = args.limit
    input_graph_csv = args.input_graph_csv
    start_date = args.start_date

    thresh, prob_matrix = load_thresholds_and_dfs(args)

    print('reading in the network data....')
    begin_time = time.time()
    # Only these columns are ever used. Reading 5 of 5 columns over ~91M rows
    # wastes several GB; `location_id` in particular is never referenced.
    # exit_state is read as a dictionary/category: it has ~135 distinct values
    # over 91M rows, so categorical encoding is ~50x smaller than object dtype.
    usecols = ["tick", "pid", "exit_state", "contact_pid"]
    try:
        df = pd.read_csv(
            input_graph_csv,
            engine="pyarrow",
            usecols=usecols,
            dtype={"tick": "int32", "pid": "int64",
                   "contact_pid": "int64", "exit_state": "category"},
        )
    except (ImportError, ValueError, TypeError) as exc:
        print(f"pyarrow read failed ({type(exc).__name__}), falling back to the C engine")
        df = pd.read_csv(
            input_graph_csv,
            usecols=usecols,
            dtype={"tick": "int32", "pid": "int64",
                   "contact_pid": "int64", "exit_state": "category"},
        )
    end_time = time.time()
    time_s = end_time - begin_time
    print(f"Done. Time: {time_s:.2f} s  rows={len(df):,}  "
          f"mem={df.memory_usage(deep=True).sum() / 1e9:.2f} GB")

    output_sets = _build_output_sets(args)
    if len(output_sets) > 1:
        print(f"  Writing {len(output_sets)} output sets "
              f"({', '.join(o.label for o in output_sets)})")
    # The unrestricted set keeps the historical variable names so the header
    # and reference-record code below is unchanged.
    seq_file = output_sets[0].seq_file
    metadata_file = output_sets[0].metadata_file

    line_keys=["virus","region","country","division","divisionExposure","date","strain","real_strain"]
    # alias_pid == infection_id == "{pid}.{tick}". Emitting it directly means
    # downstream joins against TwinSampler/BeyondBaseline files are a key match
    # rather than parsing it back out of the strain ID.
    custom_sim_keys = ["sim_pid", "sim_tick", "alias_pid"]
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

    # Header goes to every output set.
    for _o in output_sets:
        if args.compression_type in [XZ, BGZF]:
            _o.metadata_file.write(meta_line.encode('utf-8'))
        else:
            _o.metadata_file.write(meta_line)

    ref_location_dict = json.loads(args.reference_location)
    if args.reference is not None:
        align_ref = AlignIO.read(args.reference, "fasta")
        infection = InfectionRecord()
        country_ref=ref_location_dict['country'] # Use distinct names for clarity
        division_ref=ref_location_dict['division']
        divisionAbbr_ref=ref_location_dict['divisionAbbr']
        region_ref=ref_location_dict['region']
        date_ref=ref_location_dict['date']
        # f"{division_ref}-{divisionAbbr_ref}-1/{date_ref.split('-')[0]}""
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

    # Select the painted rows. `exit_state` is categorical, so resolve the
    # predicate against the ~135 distinct categories and then use a vectorised
    # isin() rather than calling a Python lambda once per row. On 91M rows this
    # is ~26x faster (28s -> ~1s).
    states = df["exit_state"]
    if hasattr(states, "cat"):
        categories = states.cat.categories
        if args.input_graph_painted_prefix:
            wanted = [c for c in categories if c.startswith(args.input_graph_painted_prefix)]
        else:
            wanted = [c for c in categories if c == args.input_graph_painted_state]
        print(f"  Painting exit_states: {sorted(wanted)}")
        transitions_to_paint = df[states.isin(wanted)]
    else:
        if args.input_graph_painted_prefix:
            mask = states.str.startswith(args.input_graph_painted_prefix)
        else:
            mask = states == args.input_graph_painted_state
        transitions_to_paint = df[mask]
    seed_transitions_mask = transitions_to_paint["contact_pid"] == -1 # Use mask for efficiency

    seed_df = transitions_to_paint[
        seed_transitions_mask
    ].copy()  # Renamed from 'seed' to 'seed_df'

    # Seed sequences do not have to be aligned to each other, and the file is
    # commonly distributed gzip/xz compressed, so use SeqIO with a
    # compression-aware handle rather than AlignIO on a bare path.
    seed_path = args.align_fasta if args.seed_fasta is None else args.seed_fasta
    if seed_path.endswith('.gz'):
        _seed_open = gzip.open
    elif seed_path.endswith('.xz'):
        _seed_open = lzma.open
    else:
        _seed_open = open
    with _seed_open(seed_path, 'rt') as _seed_handle:
        align_seed_records = list(SeqIO.parse(_seed_handle, 'fasta'))  # Read once
    print(f"  Loaded {len(align_seed_records):,} seed sequences from {seed_path}")

    # Assign seed sequences
    seed_pids = seed_df["pid"].tolist()
    N = len(seed_pids)
    M = len(align_seed_records)

    current_sequences = {}
    active_infections = {}
    if M == 0 and N > 0:
        print("Error: No sequences in seed FASTA file, but seed transitions exist. Cannot proceed.")
        sys.exit(1)

    temp_seed_seqs = {}
    seed_seq_dict = {}
    # Sequences are stored as 1-byte ASCII ('S1') rather than numpy's default
    # 4-byte UCS-4 ('<U1'). This cuts resident memory per genome from 117 KB to
    # 29 KB and makes the FASTA conversion ~190x faster (tobytes() vs
    # ''.join(tolist())).
    #
    # NOTE: only the genome is pre-allocated here. active_infections is
    # deliberately NOT populated: an importation must not become the pid's
    # "current" infection until the main loop reaches its tick. Setting it here
    # made re-importations resolve to the wrong genome (see the main loop).
    for i, (pid_val, tick) in enumerate(zip(seed_df["pid"], seed_df["tick"])):
        if i >= M:
            break
        infection_id = _ids.alias_pid(pid_val, tick)
        temp_seed_seqs[infection_id] = np.frombuffer(
            str(align_seed_records[i].seq).encode('ascii'), dtype='S1'
        )
        seed_seq_dict[infection_id] = align_seed_records[i]

    # Loudly report any importation that could not be given a genome. Each one
    # is silently dropped from the output AND orphans every downstream
    # transmission chain it would have founded, so this is a data-loss event
    # rather than a cosmetic warning.
    dropped_importations = max(0, N - M)
    if N > M:
        dropped = seed_df.iloc[M:]
        n_dropped = dropped_importations
        first_tick = int(dropped["tick"].min())
        last_tick = int(dropped["tick"].max())
        kept_last_tick = int(seed_df.iloc[M - 1]["tick"]) if M > 0 else None
        banner = "!" * 78
        print(banner, file=sys.stderr)
        print(f"WARNING: {n_dropped:,} of {N:,} importations have NO seed sequence "
              f"and will be DROPPED.", file=sys.stderr)
        print(f"  Seed FASTA supplied : {M:,} sequences ({seed_path})", file=sys.stderr)
        print(f"  Importations needed : {N:,}", file=sys.stderr)
        print(f"  Shortfall           : {n_dropped:,} ({n_dropped / N * 100:.1f}% of importations)",
              file=sys.stderr)
        print(f"  Dropped ticks       : {first_tick}-{last_tick}"
              + (f" (seeds run out after tick {kept_last_tick})" if kept_last_tick is not None else ""),
              file=sys.stderr)
        print("  Consequence: these importations are absent from the output, and every",
              file=sys.stderr)
        print("  transmission descending from them is skipped (see the 'Skipped' count at",
              file=sys.stderr)
        print("  the end of the run). Importations are consumed in tick order, so the loss",
              file=sys.stderr)
        print("  is concentrated at the END of the simulated period.", file=sys.stderr)
        print("  Fix: supply at least as many seed sequences as importations via --seed_fasta.",
              file=sys.stderr)
        print(banner, file=sys.stderr)

    current_sequences.update(temp_seed_seqs)

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

    # The tick -> date mapping only has as many distinct values as there are
    # ticks (300 here), so build a small lookup table instead of constructing
    # millions of Timestamp objects. The previous per-row .apply() cost ~62s on
    # 8M rows and is replaced by a dict lookup in the main loop.
    base_date_for_conversion = pd.to_datetime(args.start_date)
    tick_to_date = {}
    tick_to_datestr = {}
    tick_to_year = {}
    if not transitions_to_paint_df.empty or not seed_df.empty:
        all_ticks = set()
        if not transitions_to_paint_df.empty:
            all_ticks.update(transitions_to_paint_df["tick"].unique().tolist())
        if not seed_df.empty:
            all_ticks.update(seed_df["tick"].unique().tolist())
        for t in all_ticks:
            d = base_date_for_conversion + pd.Timedelta(days=(int(t) - args.start_tick))
            tick_to_date[int(t)] = d
            tick_to_datestr[int(t)] = d.strftime("%Y-%m-%d")
            tick_to_year[int(t)] = d.year

    if not transitions_to_paint_df.empty:
        pass
    elif not current_sequences: # No seeds initialized AND no transitions from other sources
        print("  No seed sequences initialized and no transitions to process. Exiting.")
        for _o in output_sets:
            _o.close()
        print("Done generating sequences (no work performed).")
        return
    else: # Seeds might be initialized, but no subsequent transitions in the filtered range
        print("  No transitions to process after filtering (seeds may have been initialized).")
        # The script might still write out the initial seed sequences if any were processed
        # Or it might just end if the loop below doesn't run.
        # Decide if you want to write just seeds if no transmissions. For now, it will proceed.

    # prob_matrix is loaded by load_thresholds_and_dfs and should always be available here.
    # If prob_matrix loading failed, the script would have exited earlier.
    if prob_matrix.size == 0: # Should not happen if load_thresholds_and_dfs succeeded
        print("Error: prob_matrix is empty. Cannot determine sequence length. Exiting.")
        # Close files if open
        for _o in locals().get("output_sets", []):
            _o.close()
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

    if use_poor_mut_model:
        mutational_model = model_registry["poor"]()
    elif args.rate_limit:
        mutational_model = model_registry["rate_limited"](args.initial_viral_load, thresh, prob_matrix, LETTERS)
    else:
        mutational_model = model_registry["simple"](thresh, cumulative_probs_matrix, LETTERS)
    print(f"Using {mutational_model} for mutations")

    missing_contact_count = 0
    reimportation_count = 0
    total_transmissions = len(transitions_to_paint_df)
    loop_begin_time = time.time()
    # Report roughly every 1% of the work, clamped to a sane range, so the log
    # shows progress promptly on small runs without flooding on large ones.
    progress_interval = min(100000, max(1000, total_transmissions // 100))
    print(f"  Painting {total_transmissions:,} infection events "
          f"(progress every {progress_interval:,})")

    # IMPORTANT: seeds and transmissions are walked in a SINGLE chronological
    # pass.
    #
    # Previously the seed rows (contact_pid == -1) were all assigned up front
    # and emitted in their own loop, then the transmission loop ran separately.
    # That silently broke re-importations: a pid can acquire an imported genome
    # at tick T2 *after* already having been infected by contact at T1 < T2.
    # Pre-seeding set active_infections[pid] = "pid.T2", but the transmission
    # loop then overwrote it with "pid.T1" when it processed the earlier event,
    # so anyone that pid infected after T2 inherited the *T1* genome instead of
    # the newly imported one - producing children that look wildly divergent
    # from their recorded parent.
    #
    # Walking every painted row in tick order and applying seeds at the moment
    # they occur keeps active_infections consistent with simulation time.
    # EpiHiper emits rows in non-decreasing tick order; sort defensively with a
    # stable kind so ties keep their original within-tick order.
    ordered_events = transitions_to_paint_df.sort_values(
        "tick", kind="stable"
    )[["pid", "contact_pid", "tick"]]

    for _, pid, contact_pid, tick in ordered_events.itertuples():
        infection_id = _ids.alias_pid(pid, tick)
        is_seed = contact_pid == -1

        if is_seed:
            # An importation: the genome comes from the seed FASTA, not from a
            # parent. Only emit it if a seed sequence was actually allocated.
            new_sequence = current_sequences.get(infection_id)
            if new_sequence is None:
                continue
            seed_fasta = seed_seq_dict.get(infection_id)
            if active_infections.get(pid) is not None:
                reimportation_count += 1
        else:
            seed_fasta = None
            contact_infection_id = active_infections.get(contact_pid)
            if contact_infection_id is None:
                # Fallback/Safety: the infector has no tracked genome. With a
                # prefix filter that captures a whole variant this should be 0;
                # a non-zero count means the painted subgraph is not closed.
                # Counted rather than printed per-event: an unthrottled print
                # here is itself a major slowdown and log-size problem.
                missing_contact_count += 1
                continue
            parent_sequence = current_sequences.get(contact_infection_id)
            if parent_sequence is None:
                missing_contact_count += 1
                continue

            new_sequence = mutational_model.mutate(parent_sequence)

        # A pid can be infected more than once (1.06M of 4.11M pids in the
        # example data, up to 7 times). Its previous infection's genome can
        # never be referenced again, so drop it to bound memory.
        previous_infection_id = active_infections.get(pid)
        if previous_infection_id is not None and previous_infection_id != infection_id:
            current_sequences.pop(previous_infection_id, None)

        current_sequences[infection_id] = new_sequence
        active_infections[pid] = infection_id

        create_infection_record(
            new_sequence,
            pid,
            tick,
            tick_to_datestr[int(tick)],
            tick_to_year[int(tick)],
            infection_id,  # infection_id
            seed_fasta,    # non-None only for importations
            country,
            region,
            division,
            divisionAbbr,
            augment_metadata,
            persontrait_df,
            standardized_aug_cols,
            standardized_replacements,
            output_sets,
            line_keys,
            args.compression_type,
        )
        loop_counter += 1
        if loop_counter % progress_interval == 0:
            elapsed = time.time() - loop_begin_time
            rate = loop_counter / elapsed if elapsed > 0 else 0.0
            remaining = total_transmissions - loop_counter
            eta_h = (remaining / rate / 3600.0) if rate > 0 else float('nan')
            msg = (f"    {loop_counter:,}/{total_transmissions:,} edges "
                   f"({100.0 * loop_counter / total_transmissions:5.1f}%)  "
                   f"{rate:,.0f} rec/s  elapsed {elapsed / 60:6.1f} min  ETA {eta_h:5.2f} h  "
                   f"live_seqs={len(current_sequences):,}")
            if hasattr(mutational_model, "allowed_to_mutate"):
                msg += f"  mutated={mutational_model.allowed_to_mutate:,}"
            if missing_contact_count:
                msg += f"  skipped={missing_contact_count:,}"
            print(msg)
        infection_counter += 1

    for _o in output_sets:
        _o.close()

    total_elapsed = time.time() - loop_begin_time
    print(f"  Wrote {loop_counter:,} records in {total_elapsed / 60:.1f} min "
          f"({loop_counter / total_elapsed if total_elapsed > 0 else 0:,.0f} rec/s)")
    if len(output_sets) > 1:
        print("  Output sets:")
        for _o in output_sets:
            print(f"    {_o.label:20s} {_o.written:>9,} records  {_o.fasta_path}")
    if reimportation_count:
        print(f"  Handled {reimportation_count:,} re-importations (a pid receiving an "
              f"imported genome after an earlier infection).")
    if dropped_importations:
        print(f"  Dropped {dropped_importations:,} importations that had no seed sequence "
              f"(seed FASTA supplied {M:,} for {N:,} importations).", file=sys.stderr)
    if missing_contact_count:
        print(f"  Skipped {missing_contact_count:,} transmissions whose infector had no "
              f"tracked genome.", file=sys.stderr)
        if dropped_importations:
            print(f"    This is expected here: the {dropped_importations:,} dropped "
                  f"importations above orphaned their descendant chains.", file=sys.stderr)
        else:
            print(f"    Expect 0 when the painted states form a closed subgraph; a non-zero "
                  f"value means some infectors were filtered out by the painted-state "
                  f"selection.", file=sys.stderr)
    print("Done generating sequences.")
    return


def process_transmission(
    infection_id, tick, contact_infection_id, seed_seq_dict, current_sequences, mutational_model
):
    """Kept for API compatibility; the main loop now inlines this so it can
    also evict the infector's superseded genome."""
    seq_to_change_arr = current_sequences[contact_infection_id]
    new_seq_arr = mutational_model.mutate(seq_to_change_arr)

    current_sequences[infection_id] = new_seq_arr  # Store the array
    return new_seq_arr


def create_infection_record(
    sequence,
    pid,
    tick,
    date_str,
    date_year,
    infection_id,
    seed_fasta,
    country,
    region,
    division,
    divisionAbbr,
    augment_metadata,
    persontrait_df,
    standardized_aug_cols,
    standardized_replacements,
    output_sets,
    line_keys,
    compression_type,
):
    # 'S1' arrays convert via tobytes() in ~2.7 us; the old
    # "".join(sequence.tolist()) path on '<U1' took ~486 us per record.
    if sequence.dtype.kind == 'S':
        new_seq_str = sequence.tobytes().decode('ascii')
    else:
        new_seq_str = "".join(sequence.tolist())
    cur_strain_id = _ids.strain_id_from_alias(country, divisionAbbr, infection_id, date_year)
    infection = InfectionRecord()
    infection.fromEpihiper(
        "ncov",
        region,
        country,
        division,
        division,  # Assuming divisionExposure is same as division
        date_str,
        cur_strain_id,
    )
    fasta_metadata=None
    if seed_fasta is not None:  # seed sequence doesn't change
        fasta_metadata = seed_fasta.id
    infection.populate_sim_details(pid, tick, fasta_metadata)
    infection.inf_dict["alias_pid"] = infection_id

    aug_metadata_dict_current = {}  # Initialize for current infection
    if augment_metadata:
        try:
            pid_df_series = persontrait_df.loc[pid]  # This should be a Series

            aug_metadata_dict_current = create_aug_metadata_dict(
                standardized_aug_cols,
                pid,
                pid_df_series,
                standardized_replacements,  # Pass standardized
            )
        except KeyError:  # pid not in persontrait_df
            aug_metadata_dict_current = create_aug_metadata_dict(
                standardized_aug_cols,
                pid,
                standardized_replacements=standardized_replacements,
            )  # Will fill with NA

    for _o in output_sets:
        if not _o.wants(infection_id):
            # A filter keyed on bare pids matches any infection of that person.
            pid_keys = getattr(_o, "pid_keys", None)
            if not (pid_keys and str(pid) in pid_keys):
                continue
        add_to_fasta(new_seq_str, infection, _o.seq_file, compression_type, fasta_metadata)
        write_metadata(
            _o.metadata_file,
            infection,
            line_keys,
            compression_type,
            aug_metadata_columns=(
                standardized_aug_cols if augment_metadata else None
            ),
            aug_metadata_dict=aug_metadata_dict_current if augment_metadata else None,
        )
        _o.written += 1


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
        return 0, 100 # Default to max conservation (100) if no data

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
    # Since the input freq_df is now strictly filtered to ACGT, the reference size is 4.
    ref_alphabet_size = 4 
    e_max_val = -np.log(1/float(ref_alphabet_size)) # np.log(4)

    if e_max_val == 0: # Avoid division by zero
        thresh = 100
    else:
        # Normalized entropy: H_norm = e_act / log(num_symbols_in_col)
        # The formula used: (1 - (e_act / e_max_val)) * 100
        # If e_act is close to e_max_val (high diversity), ratio is ~1, thresh is ~0.
        # If e_act is close to 0 (low diversity, one symbol dominates), ratio is ~0, thresh is ~100.
        # This seems correct: high threshold means high consistency (low entropy).
        thresh = (1 - (e_act / e_max_val)) * 100

    if np.isnan(thresh):
        thresh = 100

    return e_act, max(0, min(100, thresh)) # Clamp threshold 0-100


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
            "sim_pid": None, "sim_tick": None, "real_strain": None,
            "alias_pid": None
        }

    def populate_sim_details(self, sim_pid, sim_tick, real_strain=None):
        self.inf_dict["sim_pid"] = sim_pid
        self.inf_dict["sim_tick"] = sim_tick
        self.inf_dict["real_strain"] = real_strain

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

    if compression_type in [XZ, BGZF]:
        metadata_file.write(meta_line.encode('utf-8')) # Specify encoding
    else:
        metadata_file.write(meta_line)


def add_to_fasta(seq_str, infection_record, seq_file, compression_type, metadata=None): 
    # seq should be a string here
    if metadata is not None:
        seq_line = ">" + str(infection_record.get("strain")) + f" {metadata}" + "\n" + seq_str + "\n"
    else:
        seq_line = ">" + str(infection_record.get("strain")) + "\n" + seq_str + "\n"
    if compression_type in [XZ, BGZF]:
        seq_file.write(seq_line.encode('utf-8')) # Specify encoding
    else:
        seq_file.write(seq_line)


if __name__ == '__main__':
    begin_time = time.time()
    main()
    end_time = time.time()
    time_s=end_time-begin_time
    time_hr=(float)(time_s)/3600.0
    print(f"   Execution time (s): {time_s:.2f}, (hr): {time_hr:.2f}")
    print("   --- good termination ---")
