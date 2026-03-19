#!/usr/bin/env python3

import argparse
import sys
import gzip
import pandas as pd
import pyfastx

def parse_args():
    parser = argparse.ArgumentParser(
        description="Subset a BGZF compressed FASTA file using the 'strain' column from a metadata CSV."
    )
    parser.add_argument(
        "-m", "--metadata", 
        required=True, 
        help="Path to the input metadata CSV file."
    )
    parser.add_argument(
        "-f", "--fasta", 
        required=True, 
        help="Path to the input compressed FASTA file."
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="Path to save the output subset FASTA file (use .gz extension to compress)."
    )
    return parser.parse_args()

def main():
    args = parse_args()

    # 1. Read the metadata using pandas and extract the target strains
    print(f"Reading metadata from: {args.metadata}")
    try:
        df = pd.read_csv(args.metadata)
        if 'strain' not in df.columns:
            print("Error: The metadata file must contain a 'strain' column.", file=sys.stderr)
            sys.exit(1)
            
        # Get unique strains to avoid duplicating effort
        target_strains = set(df['strain'].dropna().astype(str).unique())
        print(f"Found {len(target_strains)} unique strains in the metadata.")
        
    except Exception as e:
        print(f"Error reading the metadata CSV: {e}", file=sys.stderr)
        sys.exit(1)

    # 2. Open the input FASTA file and prepare the output file
    print(f"Indexing/Opening input FASTA file: {args.fasta}")
    found_count = 0
    missing_strains =[]

    try:
        # pyfastx natively handles bgzf/gzip compression
        fasta = pyfastx.Fasta(args.fasta)
        
        print(f"Extracting sequences to: {args.output}")
        
        # Smart open: use gzip if the output filename ends with .gz
        open_func = gzip.open if args.output.endswith('.gz') else open
        mode = 'wt' if args.output.endswith('.gz') else 'w'
        
        with open_func(args.output, mode) as out_fasta:
            for strain in target_strains:
                # Check if the strain exists in the FASTA index
                if strain in fasta:
                    # Retrieve the sequence record
                    record = fasta[strain]
                    # Write to the output file in standard FASTA format
                    out_fasta.write(f">{record.name}\n{record.seq}\n")
                    found_count += 1
                else:
                    missing_strains.append(strain)
                    
    except Exception as e:
        print(f"Error processing the FASTA files: {e}", file=sys.stderr)
        sys.exit(1)

    # 3. Report the results and missing sequences
    print("\n--- Summary ---")
    print(f"Successfully extracted: {found_count} sequences.")
    print(f"Missing sequences: {len(missing_strains)}")

    if missing_strains:
        print("\nWARNING: The following strains were listed in the metadata but missing from the FASTA file:")
        for missing in sorted(missing_strains):
            print(f"  - {missing}")
    else:
        print("\nSuccess: All metadata strains were successfully found in the FASTA file!")

if __name__ == "__main__":
    main()
