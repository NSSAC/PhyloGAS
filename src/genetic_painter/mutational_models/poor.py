from .base import _MutationalModel
import numpy as np


class PoorMutationalModel(_MutationalModel):
    name = 'poor'

    # Deliberately unrealistic null model: every site has an independent 50%
    # chance of being replaced by a uniformly random canonical base.
    _ACGT = np.array(['A', 'C', 'G', 'T'])
    _ACGT_BYTES = _ACGT.astype('S1')

    def mutate(self, sequence):
        n = len(sequence)
        change_mask = np.random.random(n) < 0.5
        num_to_change = int(change_mask.sum())
        if num_to_change == 0:
            return sequence

        # Previously this built a Python list per site and returned
        # np.array("".join(...)), a 0-d string array, which broke the caller's
        # per-base indexing and FASTA conversion. Stay in numpy and preserve
        # the input dtype instead.
        output = sequence.copy()
        picks = np.random.randint(0, 4, num_to_change)
        alphabet = self._ACGT_BYTES if sequence.dtype.kind == 'S' else self._ACGT
        output[change_mask] = alphabet[picks]
        return output
