from qmof_thermo.core import calc, relax
from ase.io import read
import pytest

@pytest.fixture(scope="module")
def relaxed_structure():
    atoms = read('tests/qmof-bda2f7d.cif')
    return relax.run_calc(atoms, id="qmof-bda2f7d", out_dir="tests/relaxations_test/")


def test_energy_above_hull(relaxed_structure):
    struct, energy = relaxed_structure
    e_above_hull = calc.energy_above_hull_from_structure(struct, energy)
    assert e_above_hull == pytest.approx(0.19212944022861933)

