# qmof-thermo

A toolkit for calculating thermodynamic stability (i.e. formation energy, energy above hull) of metal–organic frameworks (MOFs) using machine-learned interatomic potentials (MLIPs).

Reference: B. Dallmann, A. Saha, A.S. Rosen, "Predicting the Thermodynamic Limits of Metal–Organic Framework Metastability" (2026).

## Overview

This package provides a streamlined workflow to calculate the formation energy and energy above hull for MOFs. In order to utilize the energy-above-hull calculation method, one must follow the **Setup Instructions** directions.

This repository also includes scripts to reproduce key figures in the manuscript. To reproduce the manuscript figures, one can simply clone the repository and follow **Figure Reproducability** directions denoted below.

## Usage

### Phase Diagram Construction

Construct and save the phase diagram from the reference data. An example is provided below.

```python
from qmof_thermo import setup_phase_diagrams

structures_path = "reference_thermo_structures.json"
thermo_path = "reference_thermo.json"
output_dir = "phase_diagrams"

setup_phase_diagrams(structures_path, thermo_path, output_dir=output_dir)
```

### Energy Above Hull Calculation

The following script allows users to relax a CIF file using an MLIP, as well as obtain an energy-above-hull calculation in eV/atom.

```python
from ase.io import read
from qmof_thermo import set_log_level, relax_mof, get_energy_above_hull

# Set logging level
set_log_level("INFO")

# Load your structure
atoms = read("mof.cif")

# Relax the structure and get energy
struct, energy = relax_mof(atoms, model="uma-s-1p1.pt", fmax=0.01, label="mymof")

# Path to directory containing serialized PatchedPhaseDiagram
references_dir = "phase_diagrams"

# Calculate energy above hull
e_above_hull = get_energy_above_hull(struct, energy, references_dir=references_dir)
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

### 3. MLIP Setup

Refer to [FairChem's documentation](https://fair-chem.github.io/) for detailed instructions on using their models.

#### UMA-ODAC

You have two options:

1. **Local checkpoint**: Download the [UMA model checkpoint](https://huggingface.co/facebook/UMA) directly.

2. **HuggingFace**: Log in with your HuggingFace credentials to download the model automatically:

   ```bash
   hf auth login
   ```

   Enter your HuggingFace token when prompted.

#### eSEN-ODAC

**Local checkpoint**: Download the [eSEN-ODAC model checkpoint](https://huggingface.co/facebook/ODAC25) directly.

## Figure Reproducibility

Scripts to reproduce the figures in the manuscript are also included in this repository and can be run as follows:

### 1. Clone and Install the Repository

```bash
git clone https://github.com/Quantum-Accelerators/qmof_thermo.git
cd qmof_thermo
```

### 2. Construct Figures

Run the corresponding Python scripts:

```bash
python figures/figure_<N>.py
```
