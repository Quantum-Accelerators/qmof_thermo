# qmof-thermo

A toolkit for calculating thermodynamic stability (energy above hull) of Metal-Organic Frameworks (MOFs) using machine learning interatomic potentials (MLIPs).

## Overview

This package provides a streamlined workflow to:
- Load MOF structures (as ASE `Atoms` objects)
- Optionally relax structures using MLIP methods (eSEN, UMA)
- Calculate the energy above hull in a single line

## Installation

### 1. Clone and Install the Package

```bash
git clone https://github.com/Quantum-Accelerators/qmof_thermo.git
cd qmof_thermo
pip install -e .
```

### 2. Download the QMOF Thermodynamics Database

Download the qmof-thermo database files from [Figshare](https://doi.org/10.6084/m9.figshare.13147324).

Place the following files in `data/external/`:
- `reference_thermo_structures.json`
- `reference_structures.json`

```bash
# Your directory structure should look like:
data/
└── external/
    ├── reference_thermo_structures.json
    └── reference_structures.json
```

### 3. Initialize the Phase Diagram

Run the setup script to construct the phase diagram from the reference data:

```bash
python src/qmof_thermo/core/setup_pd.py
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

### 1. Prepare Your Input

Place your MOF structure as a `.cif` file in `data/inputs/`:

```bash
data/
└── inputs/
    └── qmof-XXXXX.cif
```

### 2. Calculate Energy Above Hull

```python
from ase.io import read
from qmof_thermo import relax, calc

# Load your structure
atoms = read('data/inputs/qmof-XXXXX.cif')

# Relax the structure (optional) and get energy
struct, energy = relax.run_calc(atoms, id="qmof-XXXXX")

# Calculate energy above hull
e_above_hull = calc.energy_above_hull_from_structure(struct, energy)
print(f"Energy above hull: {e_above_hull} eV/atom")
```

## License

[Add your license here]

## Citation

If you use this package, please cite:

[Add citation information here]