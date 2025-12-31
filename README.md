# qmof-thermo

A toolkit for calculating thermodynamic stability (energy above hull) of Metal-Organic Frameworks (MOFs) using machine learning interatomic potentials (MLIPs). This repository also includes scripts to reproduce figures found in "Predicting the Thermodynamic Limits of Metal–Organic Framework Metastability" (Dallmann et al.).

## Overview

This package provides a streamlined workflow to:
- Load MOF structures (as ASE `Atoms` objects)
- Optionally relax structures using MLIP methods (eSEN, UMA)
- Calculate the energy above hull in a single line

This respository also includes scripts and CSV to reproduce:
- Figures 8, S15, S16, S17, and S18 in the manuscript.

To reproduce the manuscript figures, one can simply clone the repository and follow **Figure Reproducability Installation** directions
denoted below. In order to utilize the energy-above-hull calculation method, one must follow the **Energy-Above-Hull Calculator Installation** directions and pip install the repository.

## Figure Reproducibility Installation
### 1. Clone and Install the Repository
```bash
git clone https://github.com/Quantum-Accelerators/qmof_thermo.git
cd qmof_thermo
```

### 2. Construct Figures
Figures 8, S15, S16, S17, S18 each have scripts located in `qmof_thermo/figures`. Running any of the scripts will produce a corresponding folder `qmof_thermo/figures/figure_#` containing the appropriate plot.
```bash
python figures/figure_<N>.py
```

## Energy-Above-Hull Calculator Installation
Upon installation of the qmof_thermo library, users must setup phase diagram constructions.
To obtain energy_above_hull calculations using DFT values from the qmof_thermo database, users must utilize JSON files found in the qmof-thermo [Figshare](https://doi.org/10.6084/m9.figshare.13147324). All necessary steps are outlined below.

### 1. Clone and Install the Package

```bash
pip install git+https://github.com/Quantum-Accelerators/qmof_thermo.git
```

### 2. Download the QMOF Thermodynamics Database

Download the qmof-thermo database files from [Figshare](https://doi.org/10.6084/m9.figshare.13147324).

Place the following files anywhere within your accessible directory:
- `reference_thermo_structures.json`
- `reference_structures.json`

### 3. Initialize the Phase Diagram

Run the setup script to construct the phase diagram from the reference data. 
Set `structures_path` to the path of the `reference_thermo_structures.json file`, and set 
`thermo_path` to the path of the `reference_thermo.json file`. Set `pd_dir` to an accessible empty directory. `pd_dir` will be accessed everytime `calc.energy_above_hull_from_structure()` is called. 

Here is a sample script one can use to initialize the phase_diagrams. This script will fill `pd_dir` with PhaseDiagram JSON files by chemical space, along with a chemical space to ID mapping JSON file. 

```python
import logging


import qmof_thermo
from qmof_thermo.core import setup_pd

qmof_thermo.set_log_level(logging.INFO)

structures_path = "reference_thermo_structures.json"
thermo_path = "reference_thermo.json"
pd_dir = "phase_diagrams"

setup_pd.setup_phase_diagrams(structures_path, thermo_path, pd_dir)

```

## MLIP Relaxation Setup (Optional)

If you want to perform structure relaxation using MLIPs, additional setup is required. Refer to [FairChem's documentation](https://fair-chem.github.io/) for detailed instructions on using their models.

### Using eSEN

Place the eSEN model checkpoint file in the `models/` directory:

```bash
models/
└── esen_checkpoint.pt
```

### Using UMA

You have two options:

1. **Local checkpoint**: Place the UMA model checkpoint in the `models/` directory.

2. **HuggingFace (recommended)**: Log in with your HuggingFace credentials to download the model automatically:
   ```bash
   huggingface-cli login
   ```
   Enter your HuggingFace token when prompted.

### Relaxation Output

When running `relax.run_calc()`, the following files are generated in `data/relaxations/qmof-XXXXX/`:
- Trajectory file (`.traj`)
- Log file
- Relaxed structure (`.cif`)

## Usage

The following script allows users to relax a CIF file using an MLIP, as well as obtain an energy-above-hull calculation in eV/atom. 
### Relax Structure and Calculate Energy Above Hull

```python
from ase.io import read

import logging
import qmof_thermo
from qmof_thermo.core import calc, relax

#Specify level of logging. Choose between INFO, WARNING, DEBUG
qmof_thermo.set_log_level(logging.INFO)

# Load your structure
atoms = read('data/inputs/qmof-XXXXX.cif')

# Relax the structure and get energy
struct, energy = relax.run_calc(atoms, mpid="qmof-XXXXX")

# Path to directory containing PhaseDiagram JSON files
pd_dir = "phase_diagrams"

# Calculate energy above hull
e_above_hull = calc.energy_above_hull_from_structure(struct, energy, pd_dir)
print(f"Energy above hull: {e_above_hull} eV/atom")
```

<!-- 
## Citation

If you use this package, please cite (tbd)
 -->
