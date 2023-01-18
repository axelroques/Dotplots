
class Dot:

    def __init__(self, i, j, aoi) -> None:

        self.i = i
        self.j = j
        self.aoi = aoi

        # By default, the Dot is initialized as 'insignificant'
        self.is_significant = False

    def compute_frequency(self, count, freq_threshold):
        """
        Add the inverse frequency value.
        """

        self.inv_freq = 1 / count

        # If inverse frequency is greater than the threshold value
        # the dot is considered significant
        if self.inv_freq >= freq_threshold:
            self.is_significant = True

        return

    def __repr__(self) -> str:
        return '⚫'
