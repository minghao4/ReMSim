# parsers.py
import csv
import json
from pathlib import Path
from typing import Dict, List, Tuple

from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser


def parse_config(config_fpath: Path) -> dict:
    """
    Load configuration settings from config JSON into a dictionary.

    Parameters
    ----------
    config_fpath : pathlib.Path

    Returns
    -------
    dict
        Dictionary of configuration settings to values.
    """
    with config_fpath.open() as f:
        return json.load(f)


def parse_fasta(fasta_fpath: Path) -> Dict[str, Seq]:
    """
    Load contigs/chromosomes in FASTA file into a dictionary.

    Parameters
    ----------
    fasta_fpath : pathlib.Path

    Returns
    -------
    dict of {str, Bio.Seq.Seq}
        Dictionary of contig/chromosome name to sequence.
    """
    with fasta_fpath.open() as f:
        contig_seq: Dict[str, Seq] = dict(
            (contig, seq) for contig, seq in SimpleFastaParser(f)
        )

    return contig_seq


def parse_met_call(met_calls_fpath: Path) -> Dict[str, Tuple[Dict[int, int], Dict[int, int]]]:
    """
    Load processed methylation calls TSV file into a dictionary.

    Each chromosome is attached to a tuple of dictionaries, one for each strand. These dictionaries
    store information on the position of methylation calls and the methylation state at the position
    in question.

    Parameters
    ----------
    met_calls_fpath : Path

    Returns
    -------
    dict of {str, tuple of dict of {int, int}}
        Dictionary of contig/chromosome name tuple of dictionaries of position to methylation state.
    """
    with met_calls_fpath.open() as f:
        rdr = csv.reader(f, delimiter="\t")
        met_calls: Dict[str, Tuple[Dict[int, int], Dict[int, int]]] = {}

        # Read through file.
        line: List[str]
        for line in rdr:
            chrom: str = line[0]

            # New dict entry if new chromosome.
            if chrom not in met_calls.keys():
                met_calls[chrom] = ({}, {})

            strand: int = 0 if line[2] == "+" else 1
            pos: int = int(line[1])
            met_state: int = int(line[-1])  # might not include context col, met state always last

            # Add methylation state as value to the position key on the specified strand and
            # chromosome.
            met_calls[chrom][strand][pos] = met_state

        return met_calls
