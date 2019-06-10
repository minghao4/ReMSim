# simulator.py

import csv
from pathlib import Path
from typing import Dict, List

from Bio import SeqIO

"""
TODO: [Doctring]:: add docstring for simulator class.
"""


class Simulator:
    def __init__(self, chrom_seq: str, chrom_met_calls: Dict[int, float]) -> None:
        """
        Constructor.
        """
        self.sequence: str = chrom_seq
        self.met_calls: Dict[int, float] = chrom_met_calls
        self.reads: List[str] = []  # list of read objects; string type as placeholder

    def sim_reads() -> None:
        # TODO: [Sim]:: write sim_reads method
        # while x < # of reads
        # random seed datetime
        # sim start/end read_1 (150, or user specified)
        # draw the gap from normal distribution (mean 300, sd 25?)
        # discard if whole insert passes end of chromosome and redo
        # read 2 start =
        # sim end read_2 (150, or user specified)
        # sim_met()
        # add to self.reads
        ...

    def sim_met() -> None:
        # TODO: [Sim]:: write sim_met method
        # for each read
        # if include site from met call
        # random seed datetime
        # sim position methylation state (met level as probability)
        ...
