import argparse
import pandas as pd
import numpy as np
from scipy.stats import zscore
from typing import List, Dict
#https://github.com/aswarren/pango_aliasor
from pango_aliasor.aliasor import Aliasor

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

def process_tsv(input_file: str, output_file: str, pango_lineage: str):
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

def create_reductions(df: pd.DataFrame, output_file: str, num_reductions: int):
    """
    Create multiple 2-fold reductions of the table and save them to files.
    """
    for i in range(1, num_reductions + 1):
        reduction_factor = 2 ** i
        reduced_df = df.iloc[::reduction_factor, :]  # Take every nth row based on the reduction factor
        reduced_output_file = output_file.replace('.tsv', f'_reduction_{reduction_factor}.tsv')
        reduced_df.to_csv(reduced_output_file, sep='\t', index=False)
        print(f"Reduction {reduction_factor} written to {reduced_output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter a TSV table based on pango lineage and date ranges.")
    parser.add_argument("input_file", help="Path to the input TSV file.")
    parser.add_argument("output_file", help="Path to the output TSV file.")
    parser.add_argument("pango_lineage", help="Pango lineage to filter by.")
    parser.add_argument("--reductions", type=int, default=0, help="Number of 2-fold reductions to create.")

    args = parser.parse_args()
    
    filtered_df=process_tsv(args.input_file, args.output_file, args.pango_lineage)
    # Write the filtered table to the output file
    filtered_df.to_csv(output_file, sep='\t', index=False)
    print(f"Filtered table written to {output_file}")
    # Create reductions if specified
    if args.reductions > 0:
        create_reductions(filtered_df, args.output_file, args.reductions)