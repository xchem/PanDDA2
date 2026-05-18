import re

from ..interfaces import *

class FilterIgnoreDatasets:
    def __init__(self, ignore_string: str):
        self.ignore_string = ignore_string
        if ignore_string is not None:
            self.datasets_to_exclude = [x for x in self.ignore_string.split(",")]
        else:
            self.datasets_to_exclude = None

    def exclude(self, dtag):
        if self.datasets_to_exclude is None:
            return True
        else:
            if dtag in self.datasets_to_exclude:
                return False
            else:
                return True


    def __call__(self, datasets: Dict[str, DatasetInterface]):

        good_rfree_dtags = filter(
            lambda dtag: self.exclude(dtag),
            datasets,
        )

        new_datasets = {dtag: datasets[dtag] for dtag in good_rfree_dtags}

        return new_datasets

    def description(self):
        return f"Filtered because in ignore set!"