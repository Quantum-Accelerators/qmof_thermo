from pathlib import Path

import pytest
from ase.io import read
from pymatgen.core import Structure

from qmof_thermo.core import calc, relax

out_dir = Path("tests/relaxations/")
out_dir.mkdir(exist_ok=True)


@pytest.fixture()
def relaxed_structure():
    atoms = read("tests/qmof-bda2f7d_relaxed.cif")
    return Structure.from_ase_atoms(atoms)


@pytest.fixture()
def unrelaxed_atoms():
    return read("tests/qmof-bda2f7d.cif")


def test_relax(unrelaxed_atoms, tmp_path):
    struct, energy = relax.run_calc(
        unrelaxed_atoms, id="qmof-bda2f7d", out_dir=out_dir, device="cpu"
    )
    assert struct is not None
    assert energy == pytest.approx(-1191.972702)


def test_energy_above_hull(relaxed_structure):
    energy = -1191.972702
    e_above_hull = calc.energy_above_hull_from_structure(relaxed_structure, energy)
    assert e_above_hull == pytest.approx(0.19212944022861933)
