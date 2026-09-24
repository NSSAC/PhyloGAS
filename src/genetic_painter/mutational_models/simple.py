from .base import _MutationalModel
import numpy as np

class SimpleMutationalModel(_MutationalModel):
    name = 'simple'

    def __init__(self, thresholds, cumulative_probs_matrix, letters_to_use):
        self.thresholds = thresholds
        self.cumulative_probs_matrix = cumulative_probs_matrix
        self.letters_to_use = np.asarray(letters_to_use)
        # Sequences are stored as 1-byte ASCII ('S1'); keep a matching copy of
        # the alphabet so assignments do not silently coerce or raise.
        self._letters_bytes = self.letters_to_use.astype('S1')

    def mutate(self, sequence):
        change_mask = self.determine_change(self.thresholds, sequence)
        new_seq = self.weighted_change(sequence, change_mask)
        return new_seq

    def determine_change(self, thresh_array, seq_to_change_array):
        comparison_values = np.random.randint(0, 100, len(thresh_array))
        # random_selection is True if val > threshold (i.e. more random than consistent column allows mutation)
        random_selection = comparison_values > thresh_array 
        # If the position is a gap, then don't change it. The gap literal has to
        # match the array dtype, otherwise numpy compares bytes against str and
        # returns all-False, which would let gaps mutate.
        gap = b'-' if seq_to_change_array.dtype.kind == 'S' else '-'
        return random_selection & (seq_to_change_array != gap)

    def weighted_change(self, sequence_array, change_mask):
        # sequence_array is a numpy array of characters
        # change_mask is a boolean numpy array
        # cumulative_prob_matrix is pre-calculated

        # Only generate random values and find letters for positions that need to change
        num_to_change = np.sum(change_mask)
        if num_to_change == 0:
            # Return the parent array by reference rather than copying. Most
            # transmissions produce no change, so this lets identical
            # infections share memory instead of each holding a full genome.
            return sequence_array

        # Only allocate once we know something actually changes.
        output_sequence_array = sequence_array.copy()

        change_indices = np.where(change_mask)[0]

        # Generate random numbers only for these positions
        random_values_for_change = np.random.rand(num_to_change)

        # Select relevant columns from cumulative_prob_matrix
        relevant_cum_probs = self.cumulative_probs_matrix[:, change_indices]

        # Determine the index of the letter for each position to change
        # (random_values_for_change[np.newaxis, :] < relevant_cum_probs) broadcasts random_values
        # .argmax(axis=0) finds the first True, which corresponds to the letter index
        letter_indices_for_change = (random_values_for_change[np.newaxis, :] < relevant_cum_probs).argmax(axis=0)

        if sequence_array.dtype.kind == 'S':
            new_letters_for_change = self._letters_bytes[letter_indices_for_change]
        else:
            new_letters_for_change = self.letters_to_use[letter_indices_for_change]

        # Apply changes
        output_sequence_array[change_indices] = new_letters_for_change

        return output_sequence_array