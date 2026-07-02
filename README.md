# qmof-thermo

<p align="center">
<img width="443" height="222" alt="qmof_thermo" src="https://github.com/user-attachments/assets/b473f604-79ea-47c2-9bd4-0022d438d3e8" />
</p>

A toolkit for calculating thermodynamic stability (i.e. formation energy, energy above hull) of metal–organic frameworks (MOFs) using machine-learned interatomic potentials (MLIPs).

Note: For the QMOF-Thermo Database, please refer to the data available on [Figshare](https://doi.org/10.6084/m9.figshare.13147324).

Reference: B. Dallmann, A. Saha, A.S. Rosen, "Predicting the Thermodynamic Limits of Metal–Organic Framework Metastability", _J. Am. Chem. Soc._, 19, 19487--19501 (2026). DOI: [https://doi.org/10.1021/jacs.5c20253](https://doi.org/10.1021/jacs.5c20253)

## Overview

This package provides a streamlined workflow to calculate the formation energy and energy above hull for MOFs. In order to utilize the energy-above-hull calculation method, one must follow the **Setup Instructions** directions.

This repository also includes scripts to reproduce key figures in the manuscript. To reproduce the manuscript figures, one can simply clone the repository and follow **Figure Reproducability** directions denoted below.

## Usage

### Energy Above Hull Calculation

If you plan to calculate the energy of your MOF with VASP, first [set up quacc](https://quantum-accelerators.github.io/quacc/install/codes.html#vasp) and then run the following:

```python
from ase.io import read
from qmof_thermo import set_log_level, relax_mof, get_energy_above_hull
from quacc.recipes.vasp.core import static_job

# Set logging level
set_log_level("INFO")

# Load your structure
atoms = read("mof.cif")

# Relax the structure with QMOF settings and get the energy
output = static_job(atoms, preset="QMOFSet")
energy = output["results"]["energy"]

# Calculate energy above hull
e_above_hull = get_energy_above_hull(atoms, energy, energy_type="DFT")
print(f"Energy above hull: {e_above_hull} eV/atom")
```

If you plan to calculate the energy of your MOF with an ODAC-trained MLIP, first [set up the MLIP](https://github.com/Quantum-Accelerators/qmof_thermo/tree/main#2-mlip-setup) and then run the following:

```python
from ase.io import read
from qmof_thermo import set_log_level, relax_mof, get_energy_above_hull

# Set logging level
set_log_level("INFO")

# Load your structure
atoms = read("mof.cif")

# Relax the structure with an ODAC-trained MLIP and get the energy
energy = relax_mof(atoms, model="uma-s-1p1.pt", fmax=0.05, label="mymof")

# Calculate energy above hull
e_above_hull = get_energy_above_hull(atoms, energy, energy_type="ODAC_MLIP")
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

While we provide the phase diagram data with the package, to re-construct it (e.g. to add new phases), you can do so as follows:

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

## FAQ

1. Can I use an MLIP not trained on ODAC to relax my MOF when using the `qmof_thermo` package?

Not without some work. It is critical that the energy for the MOFs and the structures that compose the convex hull are all at the same level of theory, including functional, pseudopotentials, and so on. The DFT-computed energies of both the materials composing the hull and the MOFs in our work are obtained using QMOF settings, which is at the PBE-D3(BJ) level of theory. This is the same functional as UMA-ODAC and eSEN-ODAC, after filtering out elements with incompatible pseudopotentials between QMOF and ODAC. If you would like to use a completely different MLIP, you will need to to ensure internal consistency between all materials composing the the convex hull diagram.

You can use an MLIP that is compatible with the Materials Project (e.g. any [foundation MLIP fine-tuned on OAM](https://matbench-discovery.materialsproject.org/)) to obtain the energy of your MOF, but you would then need to construct the convex hull using [Pymatgen](https://github.com/materialsproject/pymatgen-core) and the Materials Project data available via the [MP API](https://github.com/materialsproject/api). Our package would not be needed in that scenario. Refer to the [Materials Project Thermodynamic Stability documentation](https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability) for details about constructing the convex hull using Materials Project-compatible data.

An MLIP that is not compatible with ODAC or the Materials Project would require recomputing the energies of all the materials composing the convex hull as well as the MOF(s) of interest, which is not recommended.
