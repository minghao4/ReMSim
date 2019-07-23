# simulator.py

from pathlib import Path
from typing import Dict, Optional, Tuple

from Bio.Seq import Seq

from base import BaseSimulator
from .read import FastqRead, FastqReadPair, VefReadPair


class FastqSimulator(BaseSimulator):
    """
    Simulate given number of fastq reads for a given chromosome.

    Attributes
    ----------
    fq_1 : pathlib.Path
        Fastq file containing the first read in each read pair.

    fq_2 : pathlib.Path
        Fastq file containing the second read in each read pair.

    Methods
    -------
    sim_all_reads()
        Simulate all fastq read pairs.
    """

    def __init__(
        self,
        source: str,
        chrom: str,
        chrom_seq: Seq,
        chrom_met_calls: Tuple[Dict[int, int]],
        num_reads: int,
        read_len: int,
        inner_dist_mu: int,
        inner_dist_sigma: int,
        output_dir: Path,
        file_prefix: str,
        window_start: int = 0,
        sim_window: int = 0,
    ) -> None:

        # Output file paths.
        self.fq_1: Path
        self.fq_2: Path

        super().__init__(
            source=source,
            chrom=chrom,
            chrom_seq=chrom_seq,
            chrom_met_calls=chrom_met_calls,
            num_reads=num_reads,
            read_len=read_len,
            inner_dist_mu=inner_dist_mu,
            inner_dist_sigma=inner_dist_sigma,
            output_dir=output_dir,
            file_prefix=file_prefix,
            window_start=window_start,
            sim_window=sim_window,
        )

    def sim_all_reads(self) -> None:
        """Simulate all fastq read pairs and write to file."""
        simmed_reads: int = 0  # accumulator
        while simmed_reads < self.num_reads:
            # TODO: [Logging]:: move to logger.
            print("Number of simmed reads: {}".format(simmed_reads))

            curr_read_pair: Optional[FastqReadPair] = self._sim_read()  # sim read pair

            # If pair is discarded, continue loop without updating accumulator.
            if curr_read_pair is None:
                print("Discard\n")  # TODO: [Logging]:: move to logger
                continue

            # Create Read objects for each read in the pair.
            curr_read_pair.set_read_objs(self.chrom_seq)
            read1: FastqRead = curr_read_pair.read_objs[0]
            read2: FastqRead = curr_read_pair.read_objs[1]

            # TODO: [Logging]:: move to logger.
            print("Accept")
            print("Chrom: {}, strand: {}".format(curr_read_pair.chrom, curr_read_pair.strand))
            print("Read 1: {} - {}".format(read1.read_start, read1.read_end))
            print("Read 1 number of met sites: {}".format(len(read1.met_calls)))
            print("Read 2: {} - {}".format(read2.read_start, read2.read_end))
            print("Read 2 number of met sites: {}\n".format(len(read2.met_calls)))

            # Write read entry to output fastq file.
            with self.fq_1.open(mode="a") as f1, self.fq_2.open(mode="a") as f2:
                fq_entries: Tuple[str, str] = curr_read_pair.entry()
                f1.write(fq_entries[0])
                f2.write(fq_entries[1])

            simmed_reads += 2  # update accumulator

    def _set_output_file_path(self) -> None:
        """Set output file path."""
        self.fq_1 = self.output_dir.joinpath("{}_sim_reads_1.fq".format(self.file_prefix))
        self.fq_2 = self.output_dir.joinpath("{}_sim_reads_2.fq".format(self.file_prefix))

    def _sim_read(self) -> Optional[FastqReadPair]:
        """
        Simulate a pair of paired-end fastq reads.

        Returns
        -------
        read.FastqReadPair or None
            A FastqReadPair if the read coordinates remain within chromosomal boundaries. None
            otherwise.
        """
        read_pair_in: Optional[
            Tuple[
                str, Tuple[Tuple[int, int], Tuple[int, int]], Tuple[Dict[int, int], Dict[int, int]]
            ]
        ] = self._sim_properties()

        if read_pair_in is None:
            return read_pair_in
        elif len(read_pair_in[2][0]) == 0 or len(read_pair_in[2][1]) == 0:
            return None
        else:
            return FastqReadPair(
                self.source, self.chrom, read_pair_in[0], read_pair_in[1], read_pair_in[2]
            )


class VefSimulator(BaseSimulator):
    """
    Simulate a given number of VEF read pairs for a given chromosome.

    Attributes
    ----------
    output_file : pathlib.Path
        Output VEF file.

    Methods
    -------
    sim_all_reads()
        Simulate all VEF read pairs.
    """

    def __init__(
        self,
        source: str,
        chrom: str,
        chrom_seq: Seq,
        chrom_met_calls: Tuple[Dict[int, int]],
        num_reads: int,
        read_len: int,
        inner_dist_mu: int,
        inner_dist_sigma: int,
        output_dir: Path,
        file_prefix: str,
        window_start: int = 0,
        sim_window: int = 0,
    ) -> None:

        # Output file path.
        self.output_file: Path

        super().__init__(
            source=source,
            chrom=chrom,
            chrom_seq=chrom_seq,
            chrom_met_calls=chrom_met_calls,
            num_reads=num_reads,
            read_len=read_len,
            inner_dist_mu=inner_dist_mu,
            inner_dist_sigma=inner_dist_sigma,
            output_dir=output_dir,
            file_prefix=file_prefix,
            window_start=window_start,
            sim_window=sim_window,
        )

    def sim_all_reads(self) -> None:
        """Simulate all VEF read pairs and write to file."""
        simmed_reads: int = 0  # accumulator
        while simmed_reads < self.num_reads:
            # TODO: [Logging]:: move to logger.
            print("Number of simmed reads: {}".format(simmed_reads))

            curr_read_pair: Optional[VefReadPair] = self._sim_read()  # sim read pair

            # If pair is discarded, continue loop without updating accumulator.
            if curr_read_pair is None:
                print("Discard\n")  # TODO: [Logging]:: move to logger
                continue

            # TODO: [Logging]:: move to logger.
            print("Accept")
            print("Chrom: {}, strand: {}".format(curr_read_pair.chrom, curr_read_pair.strand))
            print("Read 1: {} - {}".format(curr_read_pair.reads[0][0], curr_read_pair.reads[0][1]))
            print("Read 1 number of met sites: {}".format(len(curr_read_pair.read_met_calls[0])))
            print("Read 2: {} - {}".format(curr_read_pair.reads[1][0], curr_read_pair.reads[1][1]))
            print("Read 2 number of met sites: {}\n".format(len(curr_read_pair.read_met_calls[1])))

            # Write read entry to output vef file.
            with self.output_file.open(mode="a") as o_f:
                o_f.write(curr_read_pair.entry())

            simmed_reads += 2  # update accumulator...

    def _set_output_file_path(self) -> None:
        """Set the output file path."""
        self.output_file = self.output_dir.joinpath("{}_sim_reads.vef".format(self.file_prefix))

    def _sim_read(self) -> Optional[VefReadPair]:
        """
        Simulate a VEF read pair.

        Returns
        -------
        read.VefReadPair or None
            A VefReadPair if the read coordinates remain within chromosomal boundaires. None
            otherwise.
        """
        read_pair_in: Optional[
            Tuple[
                str, Tuple[Tuple[int, int], Tuple[int, int]], Tuple[Dict[int, int], Dict[int, int]]
            ]
        ] = self._sim_properties()

        if read_pair_in is None:
            return read_pair_in
        elif len(read_pair_in[2][0]) == 0 and len(read_pair_in[2][1]) == 0:
            return None
        else:
            return VefReadPair(
                self.source, self.chrom, read_pair_in[0], read_pair_in[1], read_pair_in[2]
            )
