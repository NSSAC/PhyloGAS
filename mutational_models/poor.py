from .base import _MutationalModel
import numpy as np
import random

class PoorMutationalModel(_MutationalModel):
    name = 'poor'
    
    def mutate(self, sequence):
        new_seq_list = [] # Build as list then join
        for nucleotide_char in sequence: # Iterate over chars in array
            # Original logic: change_val = np.random.randint(1, 2) -- this is always 1. So always try to change.
            # Assuming intent was 50% chance to change:
            if random.random() < 0.5: # 50% chance to enter this block
                # Original new_nucleotide = np.random.randint(1, 5) maps to ACGT, 5 was '-'
                # Let's use ACGT directly for simplicity.
                # Not clear if '-' was intended. If so, np.random.choice(['A','C','G','T','-'])
                new_nucleotide_char = random.choice(['A', 'C', 'G', 'T'])
                new_seq_list.append(new_nucleotide_char)
            else:
                new_seq_list.append(nucleotide_char) # Keep original

        return np.array("".join(new_seq_list))