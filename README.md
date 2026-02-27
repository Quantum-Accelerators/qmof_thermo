# qmof-thermo

A toolkit for calculating thermodynamic stability (i.e. formation energy, energy above hull) of metal–organic frameworks (MOFs) using machine-learned interatomic potentials (MLIPs).

Reference: B. Dallmann, A. Saha, A.S. Rosen, "Predicting the Thermodynamic Limits of Metal–Organic Framework Metastability" (2026).

## Overview

This package provides a streamlined workflow to calculate the formation energy and energy above hull for MOFs. In order to utilize the energy-above-hull calculation method, one must follow the **Setup Instructions** directions.

This respository also includes scripts to reproduce key figures in the manuscript. To reproduce the manuscript figures, one can simply clone the repository and follow **Figure Reproducability** directions denoted below.

## Usage

The following script allows users to relax a CIF file using an MLIP, as well as obtain an energy-above-hull calculation in eV/atom.

```python
import logging

from ase.io import read
from qmof_thermo import set_log_level
from qmof_thermo.core import calc, relax

# Specify level of logging. Choose between INFO, WARNING, DEBUG
set_log_level(logging.INFO)

# Load your structure
atoms = read("data/inputs/qmof-XXXXX.cif")

# Relax the structure and get energy
struct, energy = relax.run_calc(atoms, label="qmof-XXXXX")

# Path to directory containing PhaseDiagram JSON files
pd_dir = "phase_diagrams"

# Calculate energy above hull
e_above_hull = calc.energy_above_hull_from_structure(struct, energy, pd_dir)
print(f"Energy above hull: {e_above_hull} eV/atom")
```

## Setup Instructions

To obtain energy_above_hull calculations using DFT values from the qmof_thermo database, users must utilize JSON files found in the qmof-thermo [Figshare](https://doi.org/10.6084/m9.figshare.13147324). All necessary steps are outlined below.

### 1. Install the Package

```bash
pip install git+https://github.com/Quantum-Accelerators/qmof_thermo.git
```

### 2. Download the QMOF-Thermo Database

Download the qmof-thermo database files from [Figshare](https://doi.org/10.6084/m9.figshare.13147324).

Place the following files anywhere within an accessible directory:

- `reference_thermo_structures.json`
- `reference_structures.json`

### 3. Initialize the Phase Diagram

Construct and save the phase diagram from the reference data. An example is provided below.

```python
from qmof_thermo.core import setup_pd

structures_path = "/path/to/reference_thermo_structures.json"
thermo_path = "/path/to/reference_thermo.json"
output_dir = "phase_diagrams" # path to store cached phase diagrams

setup_pd.setup_phase_diagrams(structures_path, thermo_path, output_dir=output_dir)
```

### 4. MLIP Relaxation Setup (Optional)

If you want to perform structure relaxation using MLIPs, additional setup is required. Refer to [FairChem's documentation](https://fair-chem.github.io/) for detailed instructions on using their models.

#### Using eSEN

Place the eSEN model checkpoint file in the `models/` directory:

```bash
models/
└── esen_checkpoint.pt
```

#### Using UMA

You have two options:

1. **Local checkpoint**: Place the UMA model checkpoint in the `models/` directory.

2. **HuggingFace (recommended)**: Log in with your HuggingFace credentials to download the model automatically:
   ```bash
   huggingface-cli login
   ```
   Enter your HuggingFace token when prompted.

## Figure Reproducibility

Scripts to reproduce the figures in the manuscript are also included in this repository and can be run as follows:

### 1. Clone and Install the Repository

```bash
git clone https://github.com/Quantum-Accelerators/qmof_thermo.git
cd qmof_thermo
```

### 2. Construct Figures

```bash
python figures/figure_<N>.py
```
