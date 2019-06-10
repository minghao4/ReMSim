# paired_end_read.py
from typing import Tuple

"""
TODO: [Docstring]:: write docstrings for read classes
"""


class SeqRead:
    # TODO: [class] write sequence read class and associated methods.
    def __init__(self, start, end, seq) -> None:
        ...


class PairedEndRead:
    # TODO: [class]:: write paired end read class and associated methods.
    def __init__(self, read1_start: int, read1_end: int, read2_start: int, read2_end: int) -> None:
        self.read1: Tuple[int, int] = (read1_start, read1_end)
        self.read2: Tuple[int, int] = ()
