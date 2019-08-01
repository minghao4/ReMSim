# base_simulator.py

from abc import abstractmethod
import random
import time
from pathlib import Path
from typing import Dict, NoReturn, Optional, Tuple, Union

from Bio.Seq import Seq
from numpy.random import normal

from .base_read import BaseRead, BaseReadPair
import utils.util as util


class BaseSimulator:
    """
    Base simulator class.

    Attributes
    ----------
    source : str
        The source organism.

    chrom : str
        The chromosome to simulate data from.

    chrom_met_calls : tuple of dict of {int, int}
        Methylation sites and states for each strand of the current chromosome.

    num_reads : int
        Number of read pairs to simulate.

    read_len : int
        Length of the simulated reads.

    inner_dist_mu : float
        User designated mean of the inner distance.

    inner_dist_sigma : float
        User designated standard deviation of the inner distance.

    window_start : int
        The start position (0-based) of the window to simulate reads from. If not specified, set to
        beginning of chromosome.

    window_end : int
        The end position (0-based) of the window to simulate reads from. If not specified, set to
        the end of the chromosome.

    file_prefix : str
        Output file prefix.

    Methods
    -------
    sim_all_Reads()
        Simulate all read pairs.
    """

    def __init__(
        self,
        source: str,
        chrom: str,
        chrom_seq: Seq,
        chrom_met_calls: Tuple[Dict[int, int], Dict[int, int]],
        num_reads: int,
        read_len: int,
        inner_dist_mu: int,
        inner_dist_sigma: int,
        output_dir: Path,
        file_prefix: str,
        window_start: int = 0,
        sim_window: int = 0,
    ) -> None:
        # Simulation source information.
        self.chrom: str = chrom
        self.chrom_seq: Seq = chrom_seq
        self.met_calls = chrom_met_calls
        self.window_start: int = window_start

        if sim_window:  # if not 0
            self.window_end: int = self.window_start + sim_window
        else:
            self.window_end = len(chrom_seq) - 1

        self.source: str = source

        # Check if specified window has methylation calls...
        self._check_window_met_calls()

        # Output file paths.
        self.output_dir: Path = output_dir
        self.file_prefix: str = file_prefix
        util.ensure_dir(output_dir)
        self._set_output_file_path()

        # Read simulation properties.
        self.num_reads: int = num_reads
        self.read_len: int = read_len
        self.inner_dist_mu: int = inner_dist_mu  # mean inner distance
        self.inner_dist_sigma: int = inner_dist_sigma  # inner distance standard deviation

    @abstractmethod
    def sim_all_reads(self, *kwargs) -> Optional[NoReturn]:
        """
        Simulate all read pairs and write to file.

        Raises
        ------
        NotImplementedError
            If not implemented in child class.
        """
        raise NotImplementedError

    def _check_window_met_calls(self) -> Optional[NoReturn]:
        """
        """
        if (
            len(
                util.subset_dict(
                    (self.window_start, self.window_end), {**self.met_calls[0], **self.met_calls[1]}
                )
            )
            < 2
        ):
            raise LookupError(
                "Not enough (<2) methylation calls found for chr{}, {}-{} in {}.".format(
                    self.chrom, self.window_start, self.window_end - 1, self.source
                )
            )

    @abstractmethod
    def _set_output_file_path(self) -> Optional[NoReturn]:
        """
        Set output file path.

        Raises
        ------
        NotImplementedError
            If not implemented in child class.
        """
        raise NotImplementedError

    def _sim_properties(
        self
    ) -> Optional[
        Tuple[int, Tuple[Tuple[int, int], Tuple[int, int]], Tuple[Dict[int, int], Dict[int, int]]]
    ]:
        """"""
        random.seed(time.time())  # seed current time

        # Simulate read properties.
        strand: int = random.randint(0, 1)
        read1_start: int = random.randint(self.window_start, self.window_end)
        inner_distance: int = int(round(normal(self.inner_dist_mu, self.inner_dist_sigma, 1)[0]))

        # Set direction as per the strand. 0 = forward, 1 = reverse.
        directed_read_len: int = -self.read_len if strand else self.read_len
        directed_inner_dist: int = -inner_distance if strand else inner_distance

        # Early guard against simulation window boundaries. Discard if bounds are violated.
        read2_end: int = read1_start + directed_inner_dist + (2 * directed_read_len)
        if read2_end < self.window_start or read2_end > (self.window_end + 1):
            return None

        # Set rest of the read values.
        reads: Tuple[Tuple[int, int], Tuple[int, int]] = (
            (read1_start, read1_start + directed_read_len),
            (read2_end - directed_read_len, read2_end),
        )

        # Ensure ascending positional order (for start and end) for reverse strand.
        if strand:
            reads = (util.reverse_tuple_pair(reads[0]), util.reverse_tuple_pair(reads[1]))

        # Extract subset of methylation calls covered by reads.
        read_met_calls: Tuple[Dict[int, int], Dict[int, int]] = self._subset_met_calls(
            strand, reads
        )

        return strand, reads, read_met_calls

    @abstractmethod
    def _sim_read(self) -> Union[NoReturn, Optional[BaseReadPair]]:
        """
        Simulate a pair of paired-end reads.

        Returns
        -------
        base.BaseReadPair or None
            A BaseReadPair if the read coordinates remain within simulation window boundaries. None
            otherwise.
        """
        raise NotImplementedError

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
            util.subset_dict(reads[0], {**self.met_calls[0], **self.met_calls[1]}),
            util.subset_dict(reads[1], {**self.met_calls[0], **self.met_calls[1]}),
        )
