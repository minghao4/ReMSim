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

    print("Begin formatting...")  # TODO: [Logging]:: move to logger

    met_calls_fpath: Path = fmt(
        fmt_conf["mammal"],
        fmt_conf["minimum_coverage"],
        fmt_conf["start_position"],
        fmt_conf["cols_idx"],
        Path(fmt_conf["raw_methylation_calls"]),
        Path(fmt_conf["intermediate_dir"]),
        fmt_conf["met_reads_col_idx"],
    )

    print("Formatting end.\n")  # TODO: [Logging]:: move to logger
    print("Reading reference genome and processed methylation calls...")

    # Store reference genome fasta info and methylation calls info by chromosome in dictionaries.
    chrom_seq: Dict[str, Seq] = fa_parser(Path(config["reference"]))
    chrom_met_calls: Dict[str, Tuple[Dict[int, int], Dict[int, int]]] = met_parser(met_calls_fpath)

    print("Creating simulator objects...")  # TODO: [Logging]:: move to logger

    # Generate new process for each chromosome to simulate separately.
    sim_conf: dict = config["simulator"]  # simulator params
    try:
        seq_start: int = sim_conf["seq_start"]
        sim_window: int = sim_conf["sim_window"]
    except KeyError:
        seq_start = 0
        sim_window = 0

    processes: Dict[str, Tuple[Simulator, Process]] = {}
    for key in chrom_seq.keys():
        simulator = Simulator(
            chrom=key,
            seq=chrom_seq[key],
            chrom_met_calls=chrom_met_calls[key],
            num_reads=sim_conf["read_count"],
            read_len=sim_conf["read_length"],
            insert_mu=sim_conf["mean_insert_size"],
            insert_sigma=sim_conf["insert_size_standard_deviation"],
            output_dir=Path(sim_conf["output_dir"]),
            seq_start=seq_start,
            sim_window=sim_window,
            file_prefix=sim_conf["file_prefix"]
        )
        sim_proc: Process = Process(target=simulator.sim_all_reads, args=())
        processes[key] = (simulator, sim_proc)

    print("Simulation start...")  # TODO: [Logging]:: move to logger

    # Start processes.
    sim_proc: Tuple[Simulator, Process]
    for sim_proc in processes.values():
        sim_proc[1].start()

    # Exit completed processes.
    for sim_proc in processes.values():
        sim_proc[1].join()

    print("Simulation end.\n")  # TODO: [Logging]:: move to logger


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
