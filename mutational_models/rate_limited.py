import sys
from .base import _MutationalModel
import random
import numpy as np
import math

class RateLimitedMutationalModel(_MutationalModel):
    name = 'rate_limited'

    # --- Rate Limiting constants ---
    mutation_rate_per_cycle = 3.40e-6
    peak_viral_load = 1e6  # Example peak viral load
    # Fraction of peak viral load to define "early" phase, for calculating replication cycles
    rt_early_population_threshold = 0.01 
    # Probability a mutation occurring in "early" phase (defined by cycles) becomes major
    rt_early_mutation_probability = 0.8 
    # Min/max burst size for sampling
    min_burst_size = 10
    max_burst_size = 1000
    allowed_to_mutate = 0

    def __init__(self, initial_viral_load, thresholds, prob_matrix, letters_to_use):
        self.initial_viral_load = initial_viral_load
        self.thresholds = thresholds
        self.prob_matrix = prob_matrix
        self.letters_to_use = letters_to_use

    def mutate(self, sequence):
        # --- Rate Limiting Logic ---
        burst_size = random.randint(self.min_burst_size, self.max_burst_size)
        effective_initial_load = self.initial_viral_load

        replication_cycles_for_early_phase = self.calculate_replication_cycles(
            effective_initial_load,
            self.rt_early_population_threshold * self.peak_viral_load,
            burst_size
        )
        replication_cycles_for_early_phase = max(1, replication_cycles_for_early_phase)
        N = len(sequence)
        num_potential_mutations = np.random.poisson(
            self.mutation_rate_per_cycle * N * replication_cycles_for_early_phase
        )
        num_potential_mutations = min(num_potential_mutations, N)

        final_change_mask = np.zeros(N, dtype=bool)
    
        if num_potential_mutations > 0:
            self.allowed_to_mutate+=1
            # --- Weighted Site Selection ---
            # Weights: Higher for more entropy (lower threshold)
            # Ensure weights are non-negative
            site_weights = np.maximum(0.0, 1.0 - (self.thresholds / 100.0)) 

            # Do not allow mutating existing gaps
            gap_indices = np.where(sequence == '-')[0]
            site_weights[gap_indices] = 0.0

            # Normalize weights (optional, but good practice for probabilities)
            # If sum_weights is 0, it means no sites are mutable, handle this.
            sum_weights = np.sum(site_weights)
            if sum_weights > 0:
                # Using random.choices for weighted sampling.
                # Note: random.choices samples WITH replacement by default.
                # If num_potential_mutations is small relative to genome length,
                # duplicates are rare. If concerned, could sample more and take unique,
                # or implement a more complex weighted sampling without replacement.
                potential_mutation_indices = random.choices(
                    population=range(N), 
                    weights=site_weights, 
                    k=num_potential_mutations
                )

                for site_idx in potential_mutation_indices:
                    # Fixation probability for this "successful" iSNV to become major
                    if random.random() < self.rt_early_mutation_probability:
                        # The check for gap 'if seq_to_change_arr[site_idx] != '-':' is now implicitly
                        # handled by site_weights[gap_indices] = 0.0, as such sites
                        # should not be chosen by random.choices if their weight is 0.
                        # However, a direct check before assignment is still a safeguard.
                        if sequence[site_idx] != '-': # Safeguard
                            final_change_mask[site_idx] = True
            # --- End Weighted Site Selection ---
            else:
                print("\nFATAL ERROR: All site weights are zero. No mutations possible.", file=sys.stderr)
                print("This is likely because the threshold file contains all 100s due to a conserved MSA.", file=sys.stderr)
                sys.exit(1)
                
            new_seq_arr = self.weighted_change(sequence, final_change_mask)
        else:
            new_seq_arr = sequence.copy() # No mutations, return original sequence as array for consistency
        return new_seq_arr
    
    def weighted_change(self, sequence_array, change_mask):
        """
        Forces divergence to canonical ACGT bases using the raw probability matrix.
        """
        output_sequence_array = sequence_array.copy()
        num_to_change = np.sum(change_mask)
        if num_to_change == 0:
            return output_sequence_array

        canonical_letters = np.array(['A', 'C', 'G', 'T'])
        canonical_indices = np.array([np.where(self.letters_to_use == char)[0][0] for char in canonical_letters])
        canonical_char_to_idx = {char: i for i, char in enumerate(canonical_letters)}

        change_indices = np.where(change_mask)[0]
        original_letters_at_change_sites = sequence_array[change_indices]
        new_letters = np.empty_like(original_letters_at_change_sites)

        for i, site_idx in enumerate(change_indices):
            original_letter = original_letters_at_change_sites[i]
            full_site_probs = self.prob_matrix[:, site_idx]
            acgt_probs = full_site_probs[canonical_indices].copy()

            original_letter_canonical_idx = canonical_char_to_idx.get(original_letter)
            if original_letter_canonical_idx is not None:
                acgt_probs[original_letter_canonical_idx] = 0.0

            sum_acgt_probs = np.sum(acgt_probs)
            if sum_acgt_probs > 0:
                normalized_acgt_probs = acgt_probs / sum_acgt_probs
                new_letters[i] = random.choices(canonical_letters, weights=normalized_acgt_probs, k=1)[0]
            else:
                new_letters[i] = original_letter

        output_sequence_array[change_indices] = new_letters
        return output_sequence_array
    
    def calculate_replication_cycles(self, initial_virions, target_early_population, burst_size):
        """Calculates the number of replication cycles to reach target_early_population."""
        if initial_virions <= 0 or target_early_population <= 0 or burst_size <= 1:
            return 1 # Avoid math errors, assume at least 1 cycle if inputs are problematic
        if initial_virions >= target_early_population:
            return 1 # Already at or above target, assume 1 cycle for potential mutations
        
        # Formula: target = initial * (burst_size ^ cycles)
        # cycles = log_burst_size(target / initial)
        try:
            cycles = math.log(target_early_population / initial_virions, burst_size)
            return math.ceil(cycles)
        except ValueError: # e.g. log of zero or negative
            return 1

    def sample_power_law(self, min_frequency, max_frequency, alpha=2.0):
        """Samples a frequency from a power-law distribution P(x) ~ x^-alpha."""
        # Using inverse transform sampling: F(x) = (x^(1-alpha) - min^(1-alpha)) / (max^(1-alpha) - min^(1-alpha))
        # Solve for x: x = [(F(x) * (max^(1-alpha) - min^(1-alpha))) + min^(1-alpha)] ^ (1/(1-alpha))
        # F(x) is u (random number from 0 to 1)
        u = random.random()
        # Handle alpha = 1 case separately if needed, but typical iSNV alpha is ~2
        if alpha == 1.0: # Avoid division by zero if 1-alpha is zero
            # P(x) ~ 1/x. CDF is (ln(x) - ln(min)) / (ln(max) - ln(min))
            # x = exp(u * (ln(max) - ln(min)) + ln(min))
            # x = exp(u*ln(max/min) + ln(min)) = min * (max/min)^u
            return min_frequency * ((max_frequency / min_frequency) ** u)

        # Normal case for alpha != 1
        # Numerator for the exponent term
        term_min = min_frequency**(1.0 - alpha)
        term_max = max_frequency**(1.0 - alpha)
        sampled_value = (u * (term_max - term_min) + term_min)**(1.0 / (1.0 - alpha))
        return max(min_frequency, min(max_frequency, sampled_value)) # Ensure bounds