# base_read.py

from abc import abstractmethod
from typing import Dict, List, NoReturn, Tuple, Union


class BaseRead:
    """
    Base class for a bisulfite-sequencing read.

    Attributes
    ----------
    source : str
        The source of the simulation.

    chrom : str
        The chromosome the read is one.

    strand : int
        The strand the read is one. 0 for forward, 1 for reverse.

    start_end : tuple of int
        The start and end positions of the read. [start, end)

    met_calls : dict of {int, int}
        A dictionary of methylation calls. Positional keys and methylation state values.
    """

    def __init__(
        self,
        source: str,
        chrom: str,
        strand: int,
        start_end: Tuple[int, int],
        met_calls: Dict[int, int],
    ) -> None:
        self.source: str = source
        self.chrom: str = chrom
        self.strand: int = strand
        self.read_start: int = start_end[0]
        self.read_end: int = start_end[1]
        self.met_calls: Dict[int, int] = met_calls

    # TODO: [repr]:: add repr method to BaseRead


class BaseReadPair:
    """
    Base class for a pair of paired-end reads on the same strand.

    Attributes
    ----------
    source : str
        The source of the simulation.

    chrom : str
        The chromosome the read pair is on.

    strand : int
        The strand the read pair is on. 0 for forward, 1 for reverse.

    reads : tuple of tuple of int
        A tuple of the read start-end position tuples.

    read_met_calls : tuple of dict of {int, int}
        Methylation states at methylation sites within each of the reads.

    Methods
    -------
    entry()
        Create output file entry/entries for the read pair.
    """

    def __init__(
        self,
        source: str,
        chrom: str,
        strand: int,
        reads: Tuple[Tuple[int, int], Tuple[int, int]],
        read_met_calls: Tuple[Dict[int, int], Dict[int, int]],
    ) -> None:
        self.source: str = source
        self.chrom: str = chrom
        self.strand: int = strand
        self.reads: Tuple[Tuple[int, int], Tuple[int, int]] = reads
        self.read_met_calls: Tuple[Dict[int, int], Dict[int, int]] = read_met_calls

    @abstractmethod
    def entry(self) -> Union[NoReturn, Tuple[str, str], str]:
        """
        Creates the VEF entry/FASTQ entries for the read pair.

        Raises
        ------
        NotImplementedError
            If not implemented in child class.
        """
        raise NotImplementedError

    # TODO: [repr]:: add repr method to BaseReadPair.
