import pandas as pd
import numpy as np
import argparse
from Bio import AlignIO



#line_keys=["virus","region","country","division","divisionExposure","date","strain"]

def sample_metadata_uniform(metadata_df, num_samples, time_range_days=7):
    """
    Samples metadata uniformly at random within a specified time range.
    
    Parameters:
    metadata_df (pd.DataFrame): The metadata DataFrame.
    num_samples (int): The number of samples to draw.
    time_range_days (int): The time range in days for sampling. Default is 7 days.
    
    Returns:
    pd.DataFrame: The sampled metadata DataFrame.
    """
    # Convert date column to datetime
    metadata_df['date'] = pd.to_datetime(metadata_df['date'])
    
    # Get the minimum and maximum dates
    min_date = metadata_df['date'].min()
    max_date = metadata_df['date'].max()
    
    # Generate random start dates within the range
    start_dates = pd.to_datetime(np.random.choice(pd.date_range(min_date, max_date - pd.Timedelta(days=time_range_days)), num_samples))
    
    sampled_metadata = pd.DataFrame()
    
    for start_date in start_dates:
        end_date = start_date + pd.Timedelta(days=time_range_days)
        sample = metadata_df[(metadata_df['date'] >= start_date) & (metadata_df['date'] <= end_date)]
        sampled_metadata = pd.concat([sampled_metadata, sample.sample(n=1)])
    
    return sampled_metadata


def sample_metadata(metadata_df, strategy, num_samples, time_range_days=7):
    """
    Samples metadata using the specified strategy.
    
    Parameters:
    metadata_df (pd.DataFrame): The metadata DataFrame.
    strategy (str): The sampling strategy to use.
    num_samples (int): The number of samples to draw.
    time_range_days (int): The time range in days for sampling. Default is 7 days.
    
    Returns:
    pd.DataFrame: The sampled metadata DataFrame.
    """
    if strategy == 'uniform':
        return sample_metadata_uniform(metadata_df, num_samples, time_range_days)
    else:
        raise ValueError(f"Unknown sampling strategy: {strategy}")
    

def filter_fasta_by_metadata(metadata_df, fasta_path, output_path):
    """
    Filters a FASTA file to include only sequences with IDs present in the metadata DataFrame.
    
    Parameters:
    metadata_df (pd.DataFrame): The metadata DataFrame containing the IDs to filter by.
    fasta_path (str): The path to the input FASTA file.
    output_path (str): The path to the output filtered FASTA file.
    """
    ids_to_keep = set(metadata_df['strain'])
    with open(output_path, 'w') as output_handle:
        for record in AlignIO.read(fasta_path, "fasta"):
            if record.id in ids_to_keep:
                AlignIO.write(record, output_handle, "fasta")

# Example usage:
# metadata_df = read_metadata('path_to_metadata.csv')
# sampled_metadata = sample_metadata(metadata_df, 'uniform', 10, 7)
def read_metadata(file_path):
    """
    Reads metadata from a TSV file into a pandas DataFrame.
    
    Parameters:
    file_path (str): The path to the metadata TSV file.
    
    Returns:
    pd.DataFrame: The metadata DataFrame.
    """
    return pd.read_csv(file_path, sep='\t')

def main():
    parser = argparse.ArgumentParser(description='Sample metadata from a TSV file.')
    parser.add_argument('metadata_path', type=str, help='The path to the metadata TSV file.')
    parser.add_argument('--num_samples', type=int, default=10, help='The number of samples to draw.')
    parser.add_argument('--time_range_days', type=int, default=7, help='The time range in days for sampling.')
    parser.add_argument('--strategy', type=str, default='uniform', help='The sampling strategy to use.')
    
    args = parser.parse_args()
    
    metadata_df = read_metadata(args.file_path)
    sampled_metadata = sample_metadata(metadata_df, args.strategy, args.num_samples, args.time_range_days)
    
    parser.add_argument('fasta_path', type=str, help='The path to the input FASTA file.')
    parser.add_argument('output_path', type=str, help='The path to the output directory.')

    args = parser.parse_args()
    
    metadata_df = read_metadata(args.metadata_path)
    sampled_metadata = sample_metadata(metadata_df, args.strategy, args.num_samples, args.time_range_days)
    
    # Generate output file names
    metadata_base_name = args.metadata_path.split('/')[-1].replace('.tsv', '')
    fasta_base_name = args.fasta_path.split('/')[-1].replace('.fasta', '')
    fasta_out_name = f"{fasta_base_name}_samples_{args.num_samples}_days_{args.time_range_days}"
    meta_out_name = f"{metadata_base_name}_samples_{args.num_samples}_days_{args.time_range_days}"

    
    output_metadata_path = f"{args.output_path}/{meta_out_name}.tsv"
    output_fasta_path = f"{args.output_path}/{fasta_out_name}.fasta"
    
    # Write sampled metadata to file
    sampled_metadata.to_csv(output_metadata_path, sep='\t', index=False)
    
    # Filter FASTA file by sampled metadata
    filter_fasta_by_metadata(sampled_metadata, args.fasta_path, output_fasta_path)

if __name__ == '__main__':
    main()