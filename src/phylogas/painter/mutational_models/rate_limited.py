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

    def __init__(self, initial_viral_load, thresholds, prob_matrix, letters_to_use):
        self.initial_viral_load = initial_viral_load
        self.thresholds = np.asarray(thresholds, dtype=np.float64)
        self.prob_matrix = prob_matrix
        self.letters_to_use = np.asarray(letters_to_use)
        self.allowed_to_mutate = 0

        N = self.thresholds.shape[0]

        # ---- Precomputed site-selection table -------------------------------
        # The thresholds are fixed for the whole run, so the per-site weights,
        # their normalised CDF and the list of mutable sites can all be built
        # once instead of being rebuilt on every transmission.
        #
        # The previous implementation called
        #     random.choices(population=range(N), weights=site_weights, k=n)
        # which rebuilds a 29,903-entry cumulative-weight list on *every* call
        # (~5.2 ms). Sampling from a precomputed CDF with np.searchsorted is
        # ~3.7 us, a ~1400x speedup on the site-selection step alone.
        site_weights = np.maximum(0.0, 1.0 - (self.thresholds / 100.0))
        total_weight = site_weights.sum()
        if total_weight <= 0:
            print("\nFATAL ERROR: All site weights are zero. No mutations possible.", file=sys.stderr)
            print("This is likely because the threshold file contains all 100s due to a conserved MSA.", file=sys.stderr)
            sys.exit(1)

        self.site_weights = site_weights
        self._site_cdf = np.cumsum(site_weights)
        self._site_cdf /= self._site_cdf[-1]

        # ---- Precomputed substitution table ---------------------------------
        # For each site, the ACGT probabilities are fixed. Precompute them so
        # weighted_change does a couple of array lookups rather than slicing
        # and renormalising the probability matrix per mutated site.
        canonical = np.array(['A', 'C', 'G', 'T'])
        self._canonical = canonical
        self._canonical_bytes = canonical.astype('S1')
        canonical_indices = np.array(
            [np.where(self.letters_to_use == char)[0][0] for char in canonical]
        )
        # acgt_probs[site, base] -> probability of `base` at `site`
        self._acgt_probs = np.ascontiguousarray(
            np.asarray(prob_matrix)[canonical_indices, :].T
        )
        # Map a base character to its canonical column index, via a 256-entry
        # byte lookup table so we can vectorise over the mutated sites.
        self._byte_to_canonical = np.full(256, -1, dtype=np.int8)
        for i, ch in enumerate(canonical):
            self._byte_to_canonical[ord(ch)] = i
            self._byte_to_canonical[ord(ch.lower())] = i

        # Derive this model's Generator from the already-seeded legacy global
        # RNG so that --random_number_seed still makes runs reproducible.
        self._rng = np.random.default_rng(np.random.randint(0, 2**32 - 1))
        # Poisson rate is constant apart from the burst-size-dependent cycle
        # count, so cache the per-cycle term.
        self._lambda_per_cycle = self.mutation_rate_per_cycle * N

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
        num_potential_mutations = self._rng.poisson(
            self._lambda_per_cycle * replication_cycles_for_early_phase
        )
        num_potential_mutations = min(num_potential_mutations, N)

        # Fast path: with the default parameters ~82% of transmissions draw
        # zero mutations. Return the parent array by reference so the child
        # shares its memory (see weighted_change for the same optimisation).
        if num_potential_mutations == 0:
            return sequence

        self.allowed_to_mutate += 1

        # --- Weighted Site Selection (precomputed CDF) ---
        # Apply the fixation probability up front: instead of drawing k sites
        # and rejecting each with probability (1 - p), draw a Binomial number
        # of survivors and sample only those. Statistically identical, but
        # avoids the per-candidate Python loop.
        num_fixed = self._rng.binomial(
            num_potential_mutations, self.rt_early_mutation_probability
        )
        if num_fixed == 0:
            return sequence

        picks = np.searchsorted(self._site_cdf, self._rng.random(num_fixed))
        np.clip(picks, 0, N - 1, out=picks)
        change_indices = np.unique(picks)

        return self.weighted_change(sequence, change_indices)

    def weighted_change(self, sequence_array, change_indices):
        """
        Forces divergence to canonical ACGT bases using the raw probability matrix.

        `change_indices` is an array of site indices to mutate. A boolean mask
        of length N is also accepted for backwards compatibility.
        """
        change_indices = np.asarray(change_indices)
        if change_indices.dtype == bool:
            change_indices = np.where(change_indices)[0]

        if change_indices.size == 0:
            # MEMORY FIX: Return the original array reference! Do not copy.
            # This allows millions of identical infections to share the same RAM.
            return sequence_array

        # Gaps are never mutated. Sites are drawn from the entropy-derived CDF
        # which is independent of the sequence, so filter them here.
        originals = sequence_array[change_indices]
        if sequence_array.dtype.kind == 'S':
            keep = originals != b'-'
            original_codes = self._byte_to_canonical[
                np.frombuffer(originals.tobytes(), dtype=np.uint8)
            ]
        else:
            keep = originals != '-'
            original_codes = self._byte_to_canonical[
                originals.astype('S1').view(np.uint8)
            ]

        if not keep.all():
            change_indices = change_indices[keep]
            original_codes = original_codes[keep]
            if change_indices.size == 0:
                return sequence_array

        # Probabilities of the four canonical bases at each mutated site.
        probs = self._acgt_probs[change_indices].copy()

        # Zero out the current base so a mutation always changes the sequence,
        # then renormalise.
        rows = np.arange(probs.shape[0])
        valid = original_codes >= 0
        probs[rows[valid], original_codes[valid]] = 0.0

        totals = probs.sum(axis=1)
        usable = totals > 0
        if not usable.any():
            return sequence_array

        # Vectorised categorical sampling: one uniform draw per site compared
        # against the row-wise CDF.
        cdf = np.cumsum(probs[usable], axis=1)
        cdf /= cdf[:, -1:]
        draws = self._rng.random((cdf.shape[0], 1))
        choice = (draws > cdf).sum(axis=1)

        # Only allocate a copy once we know a real change is happening.
        output_sequence_array = sequence_array.copy()
        target_indices = change_indices[usable]
        if sequence_array.dtype.kind == 'S':
            output_sequence_array[target_indices] = self._canonical_bytes[choice]
        else:
            output_sequence_array[target_indices] = self._canonical[choice]
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
