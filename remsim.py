# remsim.py

import argparse
from multiprocessing import Process
from pathlib import Path
from typing import Dict, NoReturn, Optional, Tuple

from Bio.Seq import Seq

from base import BaseSimulator
from io_tools import (
    fmt_met_calls as fmt,
    parse_config,
    parse_fasta as fa_parser,
    parse_met_call as met_parser,
)
from simulator import FastqSimulator, VefSimulator


def formatting(fmt_conf: dict) -> None:
    """"""
    fmt(
        fmt_conf["mammal"],
        fmt_conf["minimum_coverage"],
        fmt_conf["start_position"],
        fmt_conf["cols_idx"],
        Path(fmt_conf["raw_methylation_calls"]),
        Path(fmt_conf["intermediate_dir"]),
        fmt_conf["met_reads_col_idx"],
    )


def simulating(
    sim_conf: dict,
    chrom_seqs: Dict[str, Seq],
    chrom_met_calls: Dict[str, Tuple[Dict[int, int], Dict[int, int]]],
) -> None:
    """"""
    try:
        window_start: int = sim_conf["window_start"]
        sim_window: int = sim_conf["sim_window"]
    except KeyError:
        window_start = 0
        sim_window = 0

    sim_input: dict = {
        "source": sim_conf["source"],
        "num_reads": sim_conf["read_count"],
        "read_len": sim_conf["read_length"],
        "inner_dist_mu": sim_conf["mean_inner_distance"],
        "inner_dist_sigma": sim_conf["inner_distance_standard_deviation"],
        "output_dir": Path(sim_conf["output_dir"]),
        "window_start": window_start,
        "sim_window": sim_window,
    }

    processes: Dict[str, Tuple[BaseSimulator, Process]] = {}
    for key in chrom_seqs.keys():
        sim_input["chrom"] = key
        sim_input["chrom_seq"] = chrom_seqs[key]
        sim_input["chrom_met_calls"] = chrom_met_calls[key]

        if sim_conf["mode"] == "fastq":
            simulator: BaseSimulator = FastqSimulator(**sim_input)
        else:
            simulator = VefSimulator(**sim_input)
        sim_proc: Process = Process(target=simulator.sim_all_reads, args=())
        processes[key] = (simulator, sim_proc)

    print("Simulation start...")  # TODO: [Logging]:: move to logger

    # Start processes.
    sim_proc: Tuple[BaseSimulator, Process]
    for sim_proc in processes.values():
        sim_proc[1].start()

    # Exit completed processes.
    for sim_proc in processes.values():
        sim_proc[1].join()


def main(config_fpath: Path, function: str) -> Optional[NoReturn]:
    """
    """
    config: dict = parse_config(config_fpath)
    if function not in ["format", "simulate"]:
        raise ValueError('Invalid function specified, please enter "format" or "simulate"')

    fmt_conf: dict = config["formatter"]
    if function == "format":
        print("Begin formatting...")  # TODO: [Logging]:: move to logger
        formatting(fmt_conf)
        print("Formatting end.\n")  # TODO: [Logging]:: move to logger

    else:
        if config["simulator"]["mode"] not in ["vef", "fastq"]:
            raise ValueError(
                'Invalid simulation mode specified, please use either "vef" or "fastq"'
            )

        out_dpath: Path = Path(fmt_conf["intermediate_dir"])
        in_fpath: Path = Path(fmt_conf["raw_methylation_calls"])
        met_calls_fpath: Path = out_dpath.joinpath("processed_" + in_fpath.name)

        # Store reference genome fasta info and methylation calls info by chromosome in
        # dictionaries.
        print("Reading reference genome and processed methylation calls...")
        chrom_seqs: Dict[str, Seq] = fa_parser(Path(config["reference"]))
        chrom_met_calls: Dict[str, Tuple[Dict[int, int], Dict[int, int]]] = met_parser(
            met_calls_fpath
        )

        # Generate new process for each chromosome to simulate separately.
        sim_conf: dict = config["simulator"]  # simulator params
        print("Creating simulator objects...")  # TODO: [Logging]:: move to logger
        simulating(sim_conf, chrom_seqs, chrom_met_calls)
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
        required=True,
    )
    required.add_argument(
        "-f", "--function", default=None, type=str, help="Function to perform.", required=True
    )

    main(Path(args.parse_args().config), args.parse_args().function)
