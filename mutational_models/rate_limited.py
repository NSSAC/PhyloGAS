from .simple import SimpleMutationalModel
import random
import numpy as np
import math

class RateLimitedMutationalModel(SimpleMutationalModel):

    name = 'rate_limited'

    # --- Rate Limiting constants (can be made CLI args later) ---
    mutation_rate_per_cycle = 3.40e-6
    peak_viral_load = 1e6  # Example peak viral load
    # Fraction of peak viral load to define "early" phase, for calculating replication cycles
    rt_early_population_threshold = 0.01 
    # Probability a mutation occurring in "early" phase (defined by cycles) becomes major
    rt_early_mutation_probability = 0.8 
    # Min/max burst size for sampling
    min_burst_size = 10
    max_burst_size = 1000

    def __init__(self, initial_viral_load, thresholds, cumulative_probs_matrix, letters_to_use):
        self.initial_viral_load = initial_viral_load
        super().__init__(thresholds, cumulative_probs_matrix, letters_to_use)

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

        new_seq_arr = self.weighted_change(
            sequence, final_change_mask
        )
        return new_seq_arr
    
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