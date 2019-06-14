# parsers.py
import csv
import json
from pathlib import Path
from typing import Dict, Tuple

from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser


def parse_config(config_fpath: Path) -> dict:
    """Load configuration settings from config json into a dictionary."""
    with config_fpath.open() as f:
        return json.load(f)


def parse_fasta(fasta_fpath: Path) -> Dict[str, Seq]:
    """Load contigs/chromosomes in fasta file into a dictionary."""
    with fasta_fpath.open() as f:
        contig_seq: Dict[str, Seq] = dict(
            (contig, seq) for contig, seq in SimpleFastaParser(f)
        )

    return contig_seq


def parse_met_call(met_calls_fpath: Path) -> Dict[str, Tuple[Dict[int, int], Dict[int, int]]]:
    """"""
    with met_calls_fpath.open() as f:
        rdr = csv.reader(f, delimiter="\t")
        met_calls: Dict[str, Tuple[Dict[int, int], Dict[int, int]]] = {}

        for line in rdr:
            chrom: str = line[0]
            if chrom not in met_calls.keys():
                met_calls[chrom] = ({}, {})

            strand: int = 0 if line[2] == "+" else 1
            pos: int = int(line[1])
            met_state: int = int(line[-1])

            met_calls[chrom][strand][pos] = met_state

        return met_calls
