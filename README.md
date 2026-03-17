# qmof-thermo

A toolkit for calculating thermodynamic stability (i.e. formation energy, energy above hull) of metal–organic frameworks (MOFs) using machine-learned interatomic potentials (MLIPs).

Note: For the QMOF-Thermo Database, please refer to the data available on [Figshare](https://doi.org/10.6084/m9.figshare.13147324).

Reference: B. Dallmann, A. Saha, A.S. Rosen, "Predicting the Thermodynamic Limits of Metal–Organic Framework Metastability" (2026).

## Overview

This package provides a streamlined workflow to calculate the formation energy and energy above hull for MOFs. In order to utilize the energy-above-hull calculation method, one must follow the **Setup Instructions** directions.

This repository also includes scripts to reproduce key figures in the manuscript. To reproduce the manuscript figures, one can simply clone the repository and follow **Figure Reproducability** directions denoted below.

## Usage

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
energy = relax_mof(atoms, model="uma-s-1p1.pt", fmax=0.05, label="mymof")

# Calculate energy above hull
e_above_hull = get_energy_above_hull(atoms, energy)
print(f"Energy above hull: {e_above_hull} eV/atom")
```

## Setup Instructions

### 1. Install the Package

```bash
pip install git+https://github.com/Quantum-Accelerators/qmof_thermo.git
```

### 2. MLIP Setup

Refer to [FairChem's documentation](https://fair-chem.github.io/) for detailed instructions on using their models.

#### UMA-ODAC

You have two options:

1. _Local checkpoint_: Download the [UMA model checkpoint](https://huggingface.co/facebook/UMA) directly.

2. _HuggingFace_: Log in with your HuggingFace credentials to download the model automatically:

   ```bash
   hf auth login
   ```

   Enter your HuggingFace token when prompted. If using this method, do not include the `.pt` extension in the `model` keyword argument.

#### eSEN-ODAC

_Local checkpoint_: Download the [eSEN-ODAC model checkpoint](https://huggingface.co/facebook/ODAC25) directly.

## Advanced: Phase Diagram Reconstruction

While we provide the phase diagram data with the package, to re-construct it (e.g. to add new phases), you can do as follows:

```python
from qmof_thermo import setup_phase_diagrams

structures_path = "reference_thermo_structures.json"  # from QMOF-Thermo figshare
thermo_path = "reference_thermo.json"  # from QMOF-Thermo figshare
output_dir = "phase_diagrams"  # directory to store patched_phase_diagram.json

setup_phase_diagrams(structures_path, thermo_path, output_dir=output_dir)
```

The resulting `phase_diagrams/patched_phase_diagram.json` can then be passed to the `serialized_phase_diagram` keyword argument of `qmof_thermo.get_energy_above_hull()`.

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
