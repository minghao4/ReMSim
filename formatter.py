# formatter.py

import csv
from pathlib import Path
from typing import List, NoReturn, Optional

"""
TODO: [Docstring]:: add description for formatter.py
"""


def fmt_met_calls(
    mammal: bool,
    min_cov: int,
    pos_start: int,
    cols: List[int],
    in_fpath: Path,
    out_fpath: Path,
    met_reads_col: Optional[int] = None,
) -> Optional[NoReturn]:
    """
    Format methylation calls TSV file.

    Parameters
    ----------
    mammal : bool
        Whether the data is mammalian; false if plant data.

    min_cov : int
        The minimum coverage filter for methylation alls.

    pos_start : int
        The base of the position value. Positions start at 0 or 1.

    cols : list of int
        An array of columns to keep from the raw methylation call TSV.
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

    in_fpath : path
        Raw TSV input file path.

    out_fpath : path
        Processed TSV output file path in intermediate directory.

    met_reads_col: list of int, optional
        If methylation level column doesn't exist in the input file (i.e. methylation level column
        index is 999), indicate the column index of the number of methylated reads.

    Raises
    ------
    ValueError
        If the parameter values are not as expected.
    """

    # Throw error if parameter values are not as expected.
    if pos_start not in [0, 1]:
        raise ValueError("Chromosomal position value must be 0-based or 1-based.")
    elif cols[6] == 999 and met_reads_col is None:
        raise ValueError(
            "Methylation level has been indicated to not be present in the input raw methylation "
            + "calls file but methylation read count column index has not been properly provided."
        )

    # Initialize new header, remove context column if mammalian data.
    new_header: List[str] = ["chrom", "0-pos", "strand", "context", "coverage" "methylation_level"]
    _rm_context_col(mammal, 3, new_header)

    with open(file=in_fpath, mode="r") as in_f, open(file=out_fpath, mode="w") as out_f:
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

                # Calculate methylation level from read counts if not available.
                if cols[5] == 999:
                    methylation_level: float = _calc_met_level(line, met_reads_col, coverage)
                else:
                    methylation_level = float(line[cols[5]])

                # Write new row to output file.
                new_row: list = [chrom, pos, strand, context, coverage, methylation_level]
                _rm_context_col(mammal, 3, new_header)
                wtr.writerow(new_row)


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


def _calc_met_level(curr_row: list, met_reads_col: int, coverage: int) -> float:
    """
    Calculate the methylation level from the number of methylated reads and the coverage.

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
    float
        The methylation level.
    """
    return int(curr_row[met_reads_col]) / coverage
