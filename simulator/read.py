# read.py

from typing import Dict, List, NoReturn, Optional, Tuple, Union

from Bio.Seq import Seq, MutableSeq

from base import BaseRead, BaseReadPair


class FastqRead(BaseRead):
    """
    A FastQ sequencing read.

    Attributes
    ----------
    seq : Bio.Seq.MutableSeq
        The mutable sequence object of the read.

    Methods
    -------
    fastq_entry()
        Creates a FastQ entry for the read.
    """

    def __init__(
        self,
        source: str,
        chrom: str,
        strand: int,
        start_end: Tuple[int],
        met_calls: Dict[int, int],
        seq: MutableSeq,
    ) -> None:
        super().__init__(
            source=source, chrom=chrom, strand=strand, start_end=start_end, met_calls=met_calls
        )
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
        return "{}\n{}\n+\n{}\n".format(self._read_name(), self.seq, self._read_qual())

    def _update_seq_met(self) -> None:
        """Updates read sequence methylation states based on methylation calls."""
        met_states: List[chr] = ["T", "C"]
        for k, v in self.met_calls.items():
            self.seq[k - self.read_start] = met_states[v]

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
        sorted_met_calls: List[int] = sorted(self.met_calls.keys())
        return "@{}, {} strand: {}-{}, methylation sites: {}".format(
            self.chrom, strands[self.strand], self.read_start, self.read_end, sorted_met_calls
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


class FastqReadPair(BaseReadPair):
    """A pair of FastQ paired-end reads on the same strand."""

    def __init__(
        self,
        source: str,
        chrom: str,
        strand: int,
        reads: Tuple[Tuple[int, int], Tuple[int, int]],
        read_met_calls: Tuple[Dict[int, int], Dict[int, int]],
    ) -> None:
        super().__init__(source, chrom, strand, reads, read_met_calls)
        self.read_objs: Tuple[FastqRead, FastqRead]

    def set_read_objs(self, chrom_seq: Seq) -> None:
        """
        Sets the FastQ read objects for the read pair.

        Parameters
        ----------
        chrom_seq : Bio.Seq.Seq
            The sequence of the chromosome.
        """
        # Set sequences.
        read_seqs: List[MutableSeq] = [
            chrom_seq[self.reads[0][0] : self.reads[0][1]].tomutable(),
            chrom_seq[self.reads[1][0] : self.reads[1][1]].tomutable(),
        ]

        # Reverse complement if reverse strand.
        if self.strand:
            read_seqs[0].complement()
            read_seqs[1].complement()

        # Set read objects.
        self.read_objs = (
            FastqRead(
                source=self.source,
                chrom=self.chrom,
                strand=self.strand,
                start_end=self.reads[0],
                met_calls=self.read_met_calls[0],
                seq=read_seqs[0],
            ),
            FastqRead(
                source=self.source,
                chrom=self.chrom,
                strand=self.strand,
                start_end=self.reads[1],
                met_calls=self.read_met_calls[1],
                seq=read_seqs[1],
            ),
        )

    def entry(self) -> Tuple[str, str]:
        """
        Creates the Fastq entries for the read pair.

        Returns
        -------
        tuple of str
            A tuple of the Fastq read entries.
        """
        return self.read_objs[0].fastq_entry(), self.read_objs[1].fastq_entry()


class VefReadPair(BaseReadPair):
    """A pair of VEF paired-end reads on the same strand."""

    def __init__(
        self,
        source: str,
        chrom: str,
        strand: int,
        reads: Tuple[Tuple[int, int], Tuple[int, int]],
        read_met_calls: Tuple[Dict[int, int], Dict[int, int]],
    ) -> None:
        super().__init__(source, chrom, strand, reads, read_met_calls)
        self.merged_met_calls: Dict[int, int] = {**self.read_met_calls[0], **self.read_met_calls[1]}

    def entry(self) -> List[str]:
        return [
            self._read_pair_name(),
            self._read_pair_alleles(),
            "//",
        ] + self._read_pair_coordinates()

    def _read_pair_name(self) -> str:
        strands: List[str] = ["forward", "reverse"]
        return "@Haplotype_{}_chrom{}_{}_strand:{}_{}".format(
            self.source, self.chrom, strands[self.strand], self.reads[0][0], self.reads[1][0]
        )

    def _read_pair_alleles(self) -> str:
        alleles: str = ""
        for k in sorted(self.merged_met_calls.keys()):
            alleles += "{}={};".format(k, self.merged_met_calls[k])

        return alleles

    def _read_pair_coordinates(self) -> List[str]:
        return [self.reads[0][0], self.reads[0][1], self.reads[1][0], self.reads[1][1]]
