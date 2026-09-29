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

    # 1. Read the metadata using pandas and extract the target strains
    print(f"Reading metadata from: {args.metadata}")
    try:
        # pandas automatically handles uncompressed, .gz, or .xz CSV files
        df = pd.read_csv(args.metadata, low_memory=False)
        if 'strain' not in df.columns:
            print("Error: The metadata file must contain a 'strain' column.", file=sys.stderr)
            sys.exit(1)
            
        # Get unique strains to avoid duplicating effort
        target_strains = set(df['strain'].dropna().astype(str).unique())
        print(f"Found {len(target_strains)} unique strains in the metadata.")
        
    except Exception as e:
        print(f"Error reading the metadata CSV: {e}", file=sys.stderr)
        sys.exit(1)

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
                    
                    if header_id in target_strains:
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

    # 3. Report the results and missing sequences
    missing_strains = target_strains - found_strains
    
    print("\n--- Summary ---")
    print(f"Successfully extracted: {len(found_strains)} sequences.")
    print(f"Missing sequences: {len(missing_strains)}")

    if missing_strains:
        print("\nWARNING: The following strains were listed in the metadata but missing from the FASTA file:")
        for missing in sorted(list(missing_strains)):
            print(f"  - {missing}")
    else:
        print("\nSuccess: All metadata strains were successfully found in the FASTA file!")

if __name__ == "__main__":
    main()
