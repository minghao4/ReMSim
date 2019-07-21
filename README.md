# ReMSim
Reference-based bisulfite sequencing read simulator.

## Requirements
* Python >=3.6
* Biopython

## Usage
```
python remsim.py
usage: remsim.py [-h] -c CONFIG

Reference Methylation Simulator

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  -c CONFIG, --config CONFIG
                        Configuration JSON file path.

  -f FUNCTION, --function FUNCTION
                        Function to perform.
```

## TODO
- [ ] Logging module
- [ ] Object reprs
- [ ] More intuitive configuration for simulated sequencing depth/coverage
- [ ] More detailed description in README
