#!/usr/bin/env python3

import argparse
import sys
import gzip
import lzma
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(
        description="Subset a FASTA file (uncompressed, .gz, or .xz) using the 'strain' column from a metadata CSV."
    )
    parser.add_argument(
        "-m", "--metadata", 
        required=True, 
        help="Path to the input metadata CSV file."
    )
    parser.add_argument(
        "-f", "--fasta", 
        required=True, 
        help="Path to the input FASTA file (can end in .fasta, .fasta.gz, or .fasta.xz)."
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="Path to save the output subset FASTA file (compression determined by .gz or .xz extension)."
    )
    return parser.parse_args()

def smart_open(filepath, mode):
    """
    Dynamically opens files using the correct library based on their extension.
    'mode' should be 'rt' (read text) or 'wt' (write text).
    """
    if filepath.endswith('.xz'):
        return lzma.open(filepath, mode)
    elif filepath.endswith('.gz'):
        return gzip.open(filepath, mode)
    else:
        # Standard uncompressed file (wants 'r' or 'w' without the 't' for standard open, 
        # though 'rt'/'wt' is valid in Python 3, we'll strip the 't' just to be perfectly standard)
        return open(filepath, mode.replace('t', ''))

def main():
    args = parse_args()

    # 1. Read the metadata and work out which sequences to keep.
    #
    # Accepts either a painter metadata file (has `strain`) or a
    # TwinSampler/BeyondBaseline linelist or samples file (has `alias_pid`).
    # The painter writes strain IDs as USA/VA-EHip-{pid}.{tick}/{year}, and now
    # also emits alias_pid directly, so both can be matched without the caller
    # having to reformat anything.
    print(f"Reading metadata from: {args.metadata}")
    try:
        df = pd.read_csv(args.metadata, low_memory=False, dtype=str)
        df.columns = [c.strip() for c in df.columns]
    except Exception as e:
        print(f"Error reading the metadata CSV: {e}", file=sys.stderr)
        sys.exit(1)

    ID_COLUMNS = ("strain", "alias_pid", "infection_id", "sim_pid", "pid")
    col = next((c for c in ID_COLUMNS if c in df.columns), None)
    if col is None:
        print(f"Error: {args.metadata} has none of the recognised identifier columns.\n"
              f"       Looked for: {', '.join(ID_COLUMNS)}\n"
              f"       Found: {', '.join(list(df.columns)[:10])}", file=sys.stderr)
        sys.exit(1)

    values = set(df[col].dropna().astype(str).str.strip().unique())
    target_strains = values if col == "strain" else set()
    target_aliases = set() if col == "strain" else values
    print(f"Found {len(values)} unique identifiers in column '{col}'.")

    # 2. Stream through the input FASTA and write to the output FASTA
    print(f"Streaming input FASTA file: {args.fasta}")
    print(f"Writing extracted sequences to: {args.output}")
    
    found_strains = set()
    
    try:
        # Use our smart_open helper for both input ('rt') and output ('wt')
        with smart_open(args.fasta, 'rt') as f_in, smart_open(args.output, 'wt') as f_out:
            keep_sequence = False
            
            for line in f_in:
                if line.startswith('>'):
                    # Extract the strain ID (up to the first space or newline)
                    header_id = line[1:].strip().split()[0]

                    # Match on the full strain ID, or on the alias_pid embedded
                    # in it (USA/VA-EHip-{pid}.{tick}/{year}).
                    hit = header_id in target_strains
                    if not hit and target_aliases and "EHip-" in header_id:
                        alias = header_id.split("EHip-")[1].rsplit("/", 1)[0]
                        hit = alias in target_aliases or alias.split(".")[0] in target_aliases
                    if hit:
                        keep_sequence = True
                        found_strains.add(header_id)
                        f_out.write(line)
                    else:
                        keep_sequence = False
                elif keep_sequence:
                    # If the switch is ON, write the sequence lines
                    f_out.write(line)
                    
    except Exception as e:
        print(f"Error processing the FASTA files: {e}", file=sys.stderr)
        sys.exit(1)

    # 3. Report the results and anything that could not be matched.
    print("\n--- Summary ---")
    print(f"Requested identifiers : {len(values)}")
    print(f"Successfully extracted: {len(found_strains)} sequences.")

    n_missing = len(values) - len(found_strains)
    if n_missing > 0:
        print(f"Missing sequences     : {n_missing}")
        if target_strains:
            missing = sorted(target_strains - found_strains)[:20]
            print("\nWARNING: listed in the metadata but absent from the FASTA:")
            for m in missing:
                print(f"  - {m}")
            if n_missing > len(missing):
                print(f"  ... and {n_missing - len(missing)} more")
        else:
            print("\nWARNING: some identifiers did not match any sequence header.")
    else:
        print("\nSuccess: every requested identifier was found.")

if __name__ == "__main__":
    main()
