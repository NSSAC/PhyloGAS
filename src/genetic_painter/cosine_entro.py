import argparse
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
from genetic_painter import aligned_to_df, df_to_entropy
import gzip
import lzma

def read_entropy_profile(file_path):
    if file_path.endswith('.gz'):
        open_func = gzip.open
    elif file_path.endswith('.xz'):
        open_func = lzma.open
    else:
        open_func = open

    with open_func(file_path, 'rt') as file:
        entropy_values = [float(line.strip()) for line in file]
    return np.array(entropy_values)

def read_entropy_msa(file_path):

    df = aligned_to_df(file_path)
    thresh, thresh_detail, entropy_values = df_to_entropy(df)
    return np.array(entropy_values)

def main():
    parser = argparse.ArgumentParser(description="Calculate cosine similarity between two entropy vectors.")
    group1 = parser.add_mutually_exclusive_group(required=True)
    group1.add_argument('--entropy_profile_1', type=str, help="Path to the first entropy profile file.")
    group1.add_argument('--entropy_msa_1', type=str, help="Path to the first MSA file.")
    
    group2 = parser.add_mutually_exclusive_group(required=True)
    group2.add_argument('--entropy_profile_2', type=str, help="Path to the second entropy profile file.")
    group2.add_argument('--entropy_msa_2', type=str, help="Path to the second MSA file.")
    
    if len(vars(parser.parse_args())) == 0:
        parser.print_help()
        parser.exit()
    
    args = parser.parse_args()
    
    if args.entropy_profile_1:
        entropy_vector_1 = read_entropy_profile(args.entropy_profile_1)
    else:
        entropy_vector_1 = read_entropy_msa(args.entropy_msa_1)
    
    if args.entropy_profile_2:
        entropy_vector_2 = read_entropy_profile(args.entropy_profile_2)
    else:
        entropy_vector_2 = read_entropy_msa(args.entropy_msa_2)
    
    similarity = cosine_similarity(entropy_vector_1.reshape(1, -1), entropy_vector_2.reshape(1, -1))
    print(f"Cosine similarity: {similarity[0][0]}")

if __name__ == "__main__":
    main()
