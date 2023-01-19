
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

    def process(self):
        """
        Compute linear regressions in each submatrices (corresponding
        to pairs of scanpaths) to find AoI patterns.
        """

        self.patterns = []

        # Iterate over each pair of scanpaths
        for boundary in self.boundaries:

            # Retrieve appropriate submatrix
            M = self.M[boundary]
            print('Boundary =', boundary)
            print('submatrix =', M)

            # List significant dots
            dots = M[M != None]
            sgnfct_dots = [
                dot for dot in dots if dot.is_significant == True
            ]
            print('all dots =', dots)
            print('significant dots =', sgnfct_dots)

            # Linear regression
            self.patterns += self._linear_regression(sgnfct_dots)

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

    @staticmethod
    def _linear_regression(dots, patterns=[]):
        """
        Linear regression analysis pipeline to find common gaze
        patterns.
        """

        # Get x, y arrays
        x = np.array([dot.j for dot in dots])
        y = np.array([dot.i for dot in dots])
        print('x =', x)
        print('y =', y)

        # Linear regression
        f, r_squared = Dotplot._linear_fit(x, y)
        print('r_squared =', r_squared)

        # If linear fit is poor, return the current patterns
        if r_squared < 0.5:
            return patterns

        # Otherwise, compute the distance between all significant
        # dots and the linear fit
        # *** Here it the simple vertical distance that is computed,
        # perhaps this should be the orthogonal distance to the fit?
        distance = np.abs(y - f(x))
        print('distance =', distance)

        # Identify dots within \mu + \sigma of the fit
        closest_dots_indices = np.where(
            distance < (distance.mean() + distance.std())
        )[0]
        print('closest_dots_indices =', closest_dots_indices)

        closest_dots = [dots[index] for index in closest_dots_indices]
        print('closest_dots =', closest_dots)

        # Do another linear fit with the closest dots
        x = np.array([dot.j for dot in closest_dots])
        y = np.array([dot.i for dot in closest_dots])
        f, _ = Dotplot._linear_fit(x, y)

        # Compute distance and select only the closest dots once again
        distance = np.abs(y - f(x))
        closest_dots_indices = np.where(
            distance < (distance.mean() + distance.std())
        )[0]
        closest_dots = [closest_dots[index] for index in closest_dots_indices]
        print('new closest_dots_indices =', closest_dots_indices)
        print('new closest_dots =', closest_dots)

        # Add closest_dots to the list of patterns
        patterns += closest_dots

        # Repeat process after excluding the closest dots
        new_dots = [
            dot for dot in dots if dot not in closest_dots
        ]
        print('dots =', dots)
        print('new dots =', new_dots)

        return Dotplot._linear_regression(new_dots, patterns=patterns)

    @ staticmethod
    def _linear_fit(x, y):
        """
        Helper function for the linear fit.
        """

        # Fit
        coeffs = np.polyfit(x, y, deg=1)

        # Create function with fit coeffs
        f = np.poly1d(coeffs)

        # r-squared computation
        correlation = np.corrcoef(x, y)[0, 1]
        r_squared = correlation ** 2

        return f, r_squared

    def __repr__(self) -> str:
        return f'{self.M}'
