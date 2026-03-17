from qmof_thermo.core import calc, setup_pd
from pymatgen.core import Structure
import logging

logging.basicConfig(level=logging.DEBUG)

# Paths
structures_path = "data/external/reference_thermo_structures.json"
thermo_path = "data/external/reference_thermo.json"
pd_dir = "data/references"

# Step 1: Build the PatchedPhaseDiagram
setup_pd.setup_phase_diagrams(structures_path, thermo_path, pd_dir)

# Step 2: Load structure and compute energy above hull
struct = Structure.from_file("mof.cif")
energy = 0.1

e_above_hull = calc.energy_above_hull_from_structure(struct, energy, pd_dir)
print(f"Energy above hull: {e_above_hull:.6f} eV/atom")
