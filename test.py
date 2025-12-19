from qmof_thermo.core import calc, relax
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import read

atoms = read('data/inputs/qmof-bda2f7d.cif')
struct, energy = relax.run_calc(atoms, id="qmof-bda2f7d",)
print(calc.energy_above_hull_from_structure(struct, energy))
