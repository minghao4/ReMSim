# formatter.py

import csv
from pathlib import Path
from typing import List, NoReturn, Optional, Union

from utils.util import ensure_dir


def fmt_met_calls(
    mammal: bool,
    min_cov: int,
    pos_start: int,
    cols: List[int],
    in_fpath: Path,
    out_dpath: Path,
    met_reads_col: Optional[int] = None,
) -> Union[NoReturn, Path]:
    """
    Format methylation calls TSV file.

    Parameters
    ----------
    mammal : bool
        Whether the data is mammalian; false if plant data.

    min_cov : int
        The minimum coverage filter for methylation calls.

    pos_start : int
        The base of the position value. Positions start at either 0 or 1.

    cols : list of int
        An array of column indices to retain from the raw methylation call TSV.
        We need the following columns:
            - Chromosome
            - Position
            - Strand
            - Context (plant only; ideally methylation context, but nucleotide context works as
              well)
            - Coverage (interchangeable with total read count)
            - Methylation level (sometimes there is no methylation level column, in that case set
              as 999 to calculate methylation level from the # of methylated reads divided by the
              total of reads)

    in_fpath : pathlib.Path
        Raw TSV input file path.

    out_dpath : pathlib.Path
        Processed TSV output directory path (intermediate directory).

    met_reads_col: int, optional
        If the methylation level column doesn't exist in the input file (i.e. the methylation level
        column index is 999), indicate the column index of the number of methylated reads.

    Returns
    -------
    out_fpath : pathlib.Path
        Processed TSV output file path.

    Raises
    ------
    ValueError
        If the parameter values are not as expected.
    """

    # Throw error if parameter values are not as expected.
    if pos_start not in [0, 1]:
        raise ValueError("Chromosomal position value must be 0-based or 1-based.")
    elif cols[6] == 999 and isinstance(met_reads_col, int):
        raise ValueError(
            "Methylation level has been indicated to not be present in the input raw methylation "
            + "calls file, but the methylation read count column index has not been properly "
            + "provided."
        )

    # Initialize new header, remove context column if mammalian data.
    new_header: List[str] = ["chrom", "0-pos", "strand", "context", "coverage" "methylation_level"]
    _rm_context_col(mammal, 3, new_header)

    # Set output file path.
    ensure_dir(out_dpath)
    out_fpath: Path = out_dpath.joinpath("processed_" + in_fpath.name)
    with in_fpath.open() as in_f, out_fpath.open(mode="w") as out_f:
        # Reader/writer objects.
        rdr = csv.reader(in_f, delimiter="\t")
        wtr = csv.writer(out_f, delimiter="\t")
        wtr.writerow(new_header)  # write new header
        next(rdr)  # skip old header

        # Read through input file.
        line: List[str]
        for line in rdr:
            coverage: int = int(line[cols[4]])

            if coverage >= min_cov:
                # Set up output row values.
                chrom: str = line[cols[0]]
                pos: int = \
                    int(line[cols[1]]) - pos_start  # minus 1 from pos if input file is 1-based

                strand: str = line[cols[2]]
                context: str = "" if mammal else _cytosine_context(line[cols[3]])

                # Calculate methylation level from read counts if not available. Round to determine
                # methylation state.
                if cols[5] == 999:
                    assert isinstance(met_reads_col, int)  # make mypy linter happy
                    methylation_state: int = _calc_met_state(line, met_reads_col, coverage)
                else:
                    methylation_state = round(float(line[cols[5]]) + 0.00001)  # guard against 0.5.

                # Write new row to output file; remove context column if mammalian data.
                new_row: list = [chrom, pos, strand, context, coverage, methylation_state]
                _rm_context_col(mammal, 3, new_header)
                wtr.writerow(new_row)

    return out_fpath


def _rm_context_col(mammal: bool, context_col_idx: int, row: list) -> None:
    """
    Remove the context column when working with mammalian data.

    Parameters
    ----------
    mammal : bool
        Whether the data is mammalian.

    context_col_idx : int
        The index of the context column.

    row : list
        The row to remove the context column from.
    """
    if mammal:
        del row[context_col_idx]


def _cytosine_context(context: str) -> str:
    """
    Determine the methylation context for plant data.

    Parameters
    ----------
    context : str
        Methylation site context.

    Returns
    -------
    str
        Cytosine methylation context.
    """
    methylation_context: str = "CHH"
    if context[1] == "G":
        methylation_context = "CG"
    elif context[2] == "G":
        methylation_context = "CHG"

    return methylation_context


def _calc_met_state(curr_row: list, met_reads_col: int, coverage: int) -> int:
    """
    Calculate the methylation level from the number of methylated reads and the coverage, then round
    to determine the methylation state.

    Parameters
    ----------
    curr_row : list of str
        Current row being read in the input file.

    met_reads_col : int
        The index of the methylated read count column.

    coverage : int
        The index of the coverage column (or total read count).

    Returns
    -------
    int
        The methylation state.
    """
    return round((int(curr_row[met_reads_col]) / coverage) + 0.00001)  # guard against 0.5
