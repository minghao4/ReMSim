# remsim.py

import argparse
from multiprocessing import Process
from pathlib import Path
from typing import Dict, Tuple

from Bio.Seq import Seq

from io_tools import (
    fmt_met_calls as fmt,
    parse_config,
    parse_fasta as fa_parser,
    parse_met_call as met_parser,
)
from simulator import Simulator


def main(config_fpath: Path) -> None:
    config: dict = parse_config(config_fpath)

    # Run formatter with formatter settings.
    fmt_conf: dict = config["formatter"]
    met_calls_fpath: Path = fmt(
        fmt_conf["mammal"],
        fmt_conf["minimum_coverage"],
        fmt_conf["start_position"],
        Path(fmt_conf["raw_methylation_calls"]),
        Path(fmt_conf["intermediate_dir"]),
        fmt_conf["met_reads_col_idx"],
    )

    # Store reference genome fasta info and methylation calls info by chromosome in dictionaries.
    chrom_seq: Dict[str, Seq] = fa_parser(Path(config["reference"]))
    chrom_met_calls: Dict[str, Tuple[Dict[int, int], Dict[int, int]]] = met_parser(met_calls_fpath)

    sim_conf: dict = config["simulator"]
    processes: Dict[str, Process] = {}
    for key in chrom_seq.keys():
        simulator = Simulator(
            chrom=key,
            seq=chrom_seq[key],
            chrom_met_calls=chrom_met_calls[key],
            num_reads=sim_conf["read_count"],
            read_len=sim_conf["read_length"],
            insert_mu=sim_conf["mean_insert_size"],
            insert_sigma=sim_conf["insert_size_standard_deviation"],
            output_dir=sim_conf["output_dir"],
        )
        process: Process = Process(target=simulator.sim_all_reads, args=())
        processes[key] = process

    for process in processes.values():
        process.start()

    for process in processes.values():
        process.join()


if __name__ == "__main__":
    args = argparse.ArgumentParser(description="Reference Methylation Simulator")
    required = args.add_argument_group("required arguments")
    required.add_argument(
        "-c",
        "--config",
        default=None,
        type=str,
        help="Configuration JSON file path.",
        required=True
    )

    main(Path(args.parse_args().config))
