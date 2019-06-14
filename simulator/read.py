# read.py

from typing import Dict, List, NoReturn, Tuple, Union

from Bio.Seq import Seq, MutableSeq


class Read:
    """
    A sequence read.

    Attributes
    ----------
    chrom : str
        The chromosome the read is on.

    strand : int
        The strand the read is on. 0 for forward, 1 for reverse.

    start_end : tuple of int
        The start and end positions of the read. [start, end)

    met_calls : dict of {int, int}
        A dictionary of methylation calls. Positional keys and methylation state values.

    seq : Bio.Seq.MutableSeq
        The mutable sequence object of the read.

    Methods
    -------
    fastq_entry()
        Creates a FASTQ entry for the read.
    """
    def __init__(
        self,
        chrom: str,
        strand: int,
        start_end: Tuple[int],
        met_calls: Dict[int, int],
        seq: MutableSeq,
    ) -> None:
        self.chrom: str = chrom
        self.strand: int = strand
        self.start: int = start_end[0]
        self.end: int = start_end[1]
        self.seq: MutableSeq = seq

    def fastq_entry(self) -> str:
        """
        Creates a FASTQ entry for the read.

        Returns
        -------
        str
            FASTQ entry for the read in this format:
            @Read name
            Sequence
            +
            Quality
        """
        self._update_seq_met()
        return "{}\n{}\n+\n{}\n".format(self._read_name(), self._read_qual(), self.seq)

    def _update_seq_met(self) -> None:
        """Updates read sequence methylation states based on methylation calls."""
        met_states: List[chr] = ["T", "C"]
        for k, v in self.met_calls.iteritems():
            self.seq[k] = met_states[v]

        # Reverse the complemented strand so that it's 5' -> 3'.
        if self.strand:
            self.seq.reverse()

    def _read_name(self) -> str:
        """
        Creates the read name from the chromosome, the strand, and the start and end positions.
        Note that start and end positions are based on forward strand but the sequence will be
        5' -> 3' for reverse strands.

        Returns
        -------
        str
            The read name in the following format:
            @[chromosome], [direction] strand: [start]-[end]
        """
        strands: List[str] = ["forward", "reverse"]
        return "@{}, {} strand: {}-{}".format(
            self.chrom, strands[self.strand], self.start, self.end
        )

    def _read_qual(self) -> str:
        """
        Creates the read quality string (PHRED33). Currently set to 40.
        TODO: [FEATURE]:: Make quality user-specifiable.

        Returns
        -------
        str
            The read quality string.
        """
        return "I" * len(self.seq)


class ReadPair:
    """
    A pair of paired end reads on the same strand.

    Attributes
    ----------
    strand : int
        The strand the reads are on. 0 for forward, 1 for reverse.

    reads : tuple of tuple of int
        A tuple of the read start-end position tuples.

    read_met_calls : tuple of dict of {int, int}
        Methylation states at methylation sites within each of the reads.

    Methods
    -------
    create_read()
        Create a Read object for a specified read in the read pair.
    """

    def __init__(
        self,
        strand: int,
        reads: Tuple[Tuple[int, int], Tuple[int, int]],
        read_met_calls: Tuple[Dict[int, int], Dict[int, int]],
    ) -> None:
        self.strand: int = strand
        self.reads: Tuple[Tuple[int, int], Tuple[int, int]] = reads
        self.read_met_calls: Tuple[Dict[int, int], Dict[int, int]] = read_met_calls

    def create_read(self, chrom: str, read: int, chrom_seq: Seq) -> Union[NoReturn, Read]:
        """
        Create a Read objected from the specified read in the read pair.

        Parameters
        ----------
        chrom : str
            The chromosome the read pair is on.

        read : int
            Which read in the pair to create the Read object from. 0 or 1 only.

        chrom_seq : Bio.Seq.Seq
            The sequence of the chromosome.

        Returns
        -------
        Read
            Read object.

        Throws
        ------
        ValueError
            Read param isn't 0 or 1.
        """
        if read not in [0, 1]:
            raise ValueError("Invalid read in pair specified.")

        curr_read: Tuple[int, int] = self.reads[read]
        read_seq: MutableSeq = chrom_seq[curr_read[0] : curr_read[1]].tomutable()
        if self.strand:
            read_seq.complement()

        return Read(chrom, self.strand, curr_read, self.read_met_calls[read], read_seq)
