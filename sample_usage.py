import logging

from ase.io import read

import qmof_thermo
from qmof_thermo.core import calc, relax

qmof_thermo.set_log_level(logging.INFO)

atoms = read("data/inputs/qmof-XXXXXX.cif")
struct, energy = relax.run_calc(
    atoms,
    id="qmof-XXXXXX",
)
print(calc.energy_above_hull_from_structure(struct, energy))
