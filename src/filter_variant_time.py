import argparse
import pandas as pd
import numpy as np
from scipy.stats import zscore
from typing import List, Dict
#https://github.com/aswarren/pango_aliasor
from pango_aliasor.aliasor import Aliasor
import us
import os

def make_variant_base_map(base_variants: List[str]) -> Dict[str, str]:
    """
    Create a mapping of variants to their base variants.
    """
    namer = Aliasor()
    namer.enable_expansion()
    all_rules = namer.partition_focus(base_variants)
    lineage_base_map = {k: v for v, ks in all_rules.items() for k in ks}
    
    for b in base_variants:
        if b not in lineage_base_map:
            lineage_base_map[b] = b
            
    return lineage_base_map

def filter_outliers(dates: pd.Series) -> pd.Timestamp:
    """
    Filter outliers from a series of dates using z-score.
    """
    dates = pd.to_datetime(dates, errors='coerce').dropna()
    if dates.empty:
        return None
    z_scores = zscore((dates - dates.min()).dt.days)
    non_outliers = dates[np.abs(z_scores) <= 2]
    return non_outliers if not non_outliers.empty else None

def process_tsv(input_file: str, pango_lineage: str):
    """
    Process the TSV file to filter rows based on the earliest and latest dates
    for a given pango lineage.
    """
    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t', parse_dates=['date'], infer_datetime_format=True)
    
    # Create the 'regularized' column
    replace_map = make_variant_base_map([pango_lineage])
    df['regularized'] = df['pango_lineage_usher'].map(replace_map)
    
    # Filter rows matching the pango lineage
    lineage_df = df[df['regularized'] == pango_lineage]
    
    # Find the earliest and latest non-outlier dates
    accepted_dates = filter_outliers(lineage_df['date'])
    if accepted_dates is None:
        print("No valid dates found for the specified lineage.")
        return
    earliest_date = accepted_dates.min()
    latest_date = accepted_dates.max()
    
    if earliest_date is None or latest_date is None:
        print("No valid dates found for the specified lineage.")
        return
    else:
        print(f"Earliest date: {earliest_date}, Latest date: {latest_date}")
    
    # Filter rows within two weeks of the earliest and latest dates
    earliest_date = pd.to_datetime(earliest_date)
    latest_date = pd.to_datetime(latest_date)
    two_weeks = pd.Timedelta(weeks=2)
    filtered_df = lineage_df[
        (pd.to_datetime(lineage_df['date'], errors='coerce') >= earliest_date - two_weeks) &
        (pd.to_datetime(lineage_df['date'], errors='coerce') <= latest_date + two_weeks)
    ]
    
    return filtered_df

def create_reductions(df: pd.DataFrame, output_base: str, num_reductions: int, target_state: str, mode: str):
    """
    Create multiple 2-fold reductions of the table and save them to files,
    applying reductions only to rows matching the target state.
    """
    df['region'] = df[df["country"] == "USA"].apply(lambda x: id_to_state(x['strain'], x['country']), axis=1)
    
    # Create a mask for the target state
    target_mask = df['region'] == target_state
    #get state abbreviation
    abbr_state = us.states.lookup(target_state).abbr
    print(f"Original number of {target_state} rows: {target_mask.sum()}")
    
    for i in range(0, num_reductions + 1):
        reduction_factor = 2 ** i
        
        if reduction_factor == 0:#code for no in-state samples
            combined_df = df[~target_mask]
            reduction_factor ="complete"
        else:
            # Apply reduction only to rows matching the target state
            reduced_target_df = df[target_mask].iloc[::reduction_factor, :]
            # Print the number of rows after reduction
            print(f"Reduction {reduction_factor}: {reduced_target_df.shape[0]} rows")
            
            # Combine reduced target rows with the rest of the rows
            combined_df = pd.concat([reduced_target_df, df[~target_mask]])
        
        if mode == "metadata":
        # Save the combined DataFrame to a file
            reduced_output_file = output_base + f'_{abbr_state}_reduction_{reduction_factor}.tsv'
            combined_df.to_csv(reduced_output_file, sep='\t', index=False)
            print(f"Reduction {reduction_factor} written to {reduced_output_file}")
        elif mode == "script":
            # Create a bash script to extract the reduced target rows
            #set reduced output file to be the protobuf file *.pb
            reduced_output_file = output_base + f'_{abbr_state}_reduction_{reduction_factor}.pb'
            reduced_output_folder = os.path.dirname(reduced_output_file)
            reduced_output_name = os.path.basename(reduced_output_file)
            reduced_ids_file = output_base + f'_{abbr_state}_reduction_{reduction_factor}_ids.txt'
            script_file = output_base + f'_{abbr_state}_reduction_{reduction_factor}.sh'
            with open(script_file, 'w') as f, open(reduced_ids_file, 'w') as id_file:
                # Write the reduced target IDs to the file
                id_file.write('\n'.join(combined_df['strain'].tolist()))
                f.write("#!/bin/bash\n")
                f.write("if [ $# -lt 1 ]; then\n")
                f.write("  echo 'Usage: $0 <input-mat>'\n")
                f.write("  exit 1\n")
                f.write("fi\n")
                f.write("\n")
                f.write(f"echo 'Running matUtils extract with the following command:'\n")
                f.write(f"echo 'matUtils extract --input-mat $1 --write-mat {reduced_output_name} --output-directory {reduced_output_folder} --samples {reduced_ids_file}'\n")
                f.write("\n")
                f.write(f"matUtils extract --input-mat $1 --write-mat {reduced_output_name} --output-directory {reduced_output_folder} --samples {reduced_ids_file}\n")
            print(f"Script for reduction {reduction_factor} written to {script_file}")


def id_to_state(strain: str, country: str) -> str:
    """
    Parse the strain column to extract the state for samples from the USA using the `us` package.

    Args:
        strain (str): The strain identifier (e.g., "USA/CA-12345/2021").
        country (str): The country of the sample (e.g., "USA").

    Returns:
        str: The full state name if found, otherwise an empty string.
    """
    try:
        # Extract the state abbreviation from the strain
        state_abbreviation = strain.split("/")[1].split("-")[0].upper()
        # Use the `us` package to look up the state
        state = us.states.lookup(state_abbreviation)
        if state is not None:
            return state.name#.replace(" ", "_")  # Return the full state name with underscores
    except IndexError:
        # Handle cases where the strain format is unexpected
        pass
    return ""


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter a TSV table based on pango lineage and date ranges.")
    parser.add_argument("input_file", help="Path to the input TSV file.")
    parser.add_argument("output_path", help="Path to the output TSV file.")
    parser.add_argument("pango_lineage", help="Pango lineage to filter by.")
    parser.add_argument("--reductions", type=int, default=5, help="Number of 2-fold reductions to create.")
    #reduction target state
    parser.add_argument("--target_state", type=str, default=None, help="Target state for reduction. e.g. Virginia")
    #two modes one to filter and produce metadata, two produce the bash script with matutils extract and id list
    parser.add_argument("--mode", type=str, default="script", choices=["metadata", "script"], help="Mode of operation: 'filter' or 'script'.")
    args = parser.parse_args()
    #combine output path with input file name minus extension
    output_base = args.output_path + os.path.splitext(os.path.basename(args.input_file))[0]+"_" + args.pango_lineage
    args = parser.parse_args()
    #filter by time / pango lineage
    filtered_df=process_tsv(args.input_file, args.pango_lineage)
    if args.target_state != None and args.reductions > 0:
        #create a region column
        create_reductions(filtered_df, output_base, args.reductions, args.target_state, args.mode)


