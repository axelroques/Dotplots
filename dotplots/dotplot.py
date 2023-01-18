
from .dot import Dot

from collections import Counter
import numpy as np


class Dotplot:

    def __init__(self, sequences, freq_threshold) -> None:
        """
        Inputs:
            - sequences = list of scanpaths
            - freq_threshold = frequency threshold on the 
            inverse of the occurrences
        """

        # Store input parameter
        self.freq_threshold = freq_threshold

        # Concatenate sequences
        self.sequence, self.boundaries = self._concatenate_sequences(
            sequences
        )
        print('self.sequence =', self.sequence)
        print('self.boundaries =', self.boundaries)

        # Build and filter dotplot representation
        self.M = self._build(self.sequence, freq_threshold)
        # print('self.M =', self.M)

    def process(self):
        """
        Compute linear regressions in each submatrices (corresponding
        to pairs of scanpaths) to find AoI patterns.
        """

        return

    @staticmethod
    def _concatenate_sequences(sequences):
        """
        Concatenate input scanpaths.
        Also returns a list with the different boundaries of the
        future submatrices (corresponding to all pairwise comparison 
        of scanpaths).
        """

        # Initialize lists
        sequence = []
        boundaries = []

        i_boundary_start = 0
        n = len(sequences)
        for i in range(n):

            # Add sequence to the list
            sequence += sequences[i]

            # Compute boundaries
            j_boundary_start = len(sequences[i])
            for j in range(i+1, n):

                i_boundary_end = i_boundary_start + len(sequences[i])
                j_boundary_end = j_boundary_start + len(sequences[j])
                boundaries.append(
                    (
                        slice(i_boundary_start, i_boundary_end),
                        slice(j_boundary_start, j_boundary_end)
                    )
                )

                j_boundary_start += len(sequences[j])

            i_boundary_start += len(sequences[i])

        return sequence, boundaries

    @staticmethod
    def _build(sequence, freq_threshold):
        """
        Build dotplot. 
        Also includes the inverse-frequency filter step.
        """

        # Count occurrences of each AoI
        occurrences = Counter(sequence)

        # Generator object
        generator = (
            (i, j)
            for i in range(len(sequence))
            for j in range(i, len(sequence))
        )

        # Matrix generation
        n = len(sequence)
        M = np.empty((n, n), dtype=object)
        for i, j in generator:

            # Only concerned with AoI matches
            aoi_1 = sequence[i]
            aoi_2 = sequence[j]
            if aoi_1 == aoi_2:

                # Create Dot object
                dot = Dot(i, j, aoi_1)

                # Add inverse frequency information
                dot.compute_frequency(occurrences[aoi_1], freq_threshold)

                # Store Dot object
                M[i, j] = dot

        return M

    def __repr__(self) -> str:
        return f'{self.M}'
