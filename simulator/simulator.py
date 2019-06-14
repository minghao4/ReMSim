# simulator.py

import random
import time
from pathlib import Path
from typing import Dict, Optional, Tuple

from Bio.Seq import Seq
from numpy.random import normal

from .read import ReadPair, Read
import utils.util as util


class Simulator:
    """
    Simulate given number of reads for a given chromosome.

    Attributes
    ----------
    chrom : str
        The chromosome to simulate reads from.

    seq : Bio.Seq.Seq
        Biopython immutable sequence object; the DNA sequence of the current chromosome.

    chrom_met_calls : tuple of dict of {int, int}
        Methylation sites and states for each strand of the current chromosome.

    num_reads : int
        Number of read pairs to simulate.

    read_len : int
        Length of the simulated reads.

    insert_mutupr designated mean insert size.

    insert_sigma : int
        User designated insert size standard deviation.

    output_folder : pathlib.Path
        Output directory path.

    Methods
    -------
    sim_all_reads()
        Simulate all read pairs.
    """
    def __init__(
        self,
        chrom: str,
        seq: Seq,
        chrom_met_calls: Tuple[Dict[int, int]],
        num_reads: int,
        read_len: int,
        insert_mu: int,
        insert_sigma: int,
        output_dir: Path,
    ) -> None:
        # Chromosome information.
        self.chrom = chrom
        self.sequence: str = seq
        self.seq_len: int = len(seq) - 1  # 0-based

        # Methylation calls.
        self.met_calls: Tuple[Dict[int, int]] = chrom_met_calls

        # Read simulation properties.
        self.num_reads: int = num_reads
        self.read_len: int = read_len
        self.insert_mu: int = insert_mu
        self.insert_sigma: int = insert_sigma

        # Output file paths.
        self.fq_1: Path = output_dir.joinpath("sim_reads_1.fq")
        self.fq_2: Path = output_dir.joinpath("sim_reads_2.fq")

    def sim_all_reads(self) -> None:
        """Simulate all read pairs and write to file."""
        simmed_reads: int = 0  # accumulator
        while simmed_reads <= self.num_reads:
            curr_read_pair: ReadPair = self._sim_read()  # sim read pair

            # If pair is discarded, continue loop without updating accumulator.
            if curr_read_pair is None:
                continue

            # Create Read objects for each read in the pair.
            read1: Read = curr_read_pair.create_read(self.chrom, 0, self.sequence)
            read2: Read = curr_read_pair.create_read(self.chrom, 1, self.sequence)

            # Write read entry to output fastq file.
            with self.fq_1.open(mode="a") as f1, self.fq_2.open(mode="a") as f2:
                f1.write(read1.fastq_entry())
                f2.write(read2.fastq_entry())

            simmed_reads += 1  # update accumulator

    def _sim_read(self) -> Optional[ReadPair]:
        """
        Simulate a pair of paired end reads.

        Returns
        -------
        read.PairedRead or None
            A PairedRead if the read coordinates remain within chromosomal boundaires. None
            otherwise.
        """
        random.seed(time.time())  # seed current time

        # Simulate read properties.
        strand: int = random.randint(0, 1)
        read1_start: int = random.randint(0, self.length)
        insert_size: int = round(normal(self.insert_mu, self.insert_sigma, 1)[0])

        # Set direction as per the strand.
        directed_read_len: int = -self.read_len if strand else self.read_len
        directed_insert_size: int = -insert_size if strand else insert_size

        # Early guard against chromosomal boundaries. Discard if bounds are violated.
        read2_end: int = read1_start + directed_insert_size + (2 * directed_read_len)
        if read2_end < 0 or read2_end > (self.seq_len + 1):
            return None

        # Set rest of read values.
        reads: Tuple[Tuple[int, int], Tuple[int, int]] = (
            (read1_start, read1_start + directed_read_len),
            (read2_end - directed_read_len, read2_end),
        )

        # Ensure ascending positional order (for start and end) for reverse strand.
        if strand:
            reads = (util.reverse_tuple_pair(reads[0]), util.reverse_tuple_pair(reads[1]))

        # Extract subset of methylation calls.
        read_met_calls: Dict[int, int] = self._subset_met_calls(strand, reads)

        return ReadPair(strand, reads, read_met_calls)

    def _subset_met_calls(
        self, strand: int, reads: Tuple[Tuple[int, int], Tuple[int, int]]
    ) -> Tuple[Dict[int, int], Dict[int, int]]:
        """
        Discover the subset of methylation calls within a pair of paired end reads.

        Parameters
        ----------
        strand : int
            The strand the read pair is on. 0 for forward, 1 for reverse.

        reads : tuple of tuple of int
            A tuple of the read start-end position tuples.

        Returns
        -------
        tuple of dict of {int, int}
            Dictionaries of position to methylation state of all positions within the paired reads.
        """
        return (
            util.subset_dict(reads[0], self.met_calls[strand]),
            util.subset_dict(reads[1], self.met_calls[strand]),
        )
