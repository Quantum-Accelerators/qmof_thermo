import logging

from ase.io import read

import qmof_thermo
from qmof_thermo.core import calc, relax, setup_pd

qmof_thermo.set_log_level(logging.INFO)

structures_path = "reference_thermo_structures.json"
thermo_path = "reference_thermo.json"
pd_dir = "phase_diagrams"

setup_pd.setup_phase_diagrams(structures_path, thermo_path, pd_dir)

atoms = read("data/inputs/qmof-XXXXXX.cif")
struct, energy = relax.run_calc(
    atoms,
    id="qmof-XXXXXX",
)
print(calc.energy_above_hull_from_structure(struct, energy, pd_dir))