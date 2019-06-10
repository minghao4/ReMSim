# main.py

import csv
import sys
from pathlib import Path

import Bio

from formatter import fmt_met_calls as fmt
# from ref_met_sim import read_simulator as sim

# Input parameters
# TODO: [Main]:: write config json parsing logic
config_fpath: Path = Path(sys.argv[1])  # configuration json

# TODO: [Main]:: write rest of main logic
# 1. Format raw methylation call file
# fmt()

# 2. for chrom in ref, create new process to simulate reads
# sim()
