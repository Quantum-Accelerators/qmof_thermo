from __future__ import annotations

import json
from pathlib import Path

import pytest
from ase.io import read
from pymatgen.core import Structure

from qmof_thermo import get_energy_above_hull, relax_mof, setup_phase_diagrams

FILE_DIR = Path(__file__).parent
TEST_DATA_DIR = FILE_DIR / "test_data"


@pytest.fixture(scope="module")
def out_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("relaxations")


@pytest.fixture(scope="module")
def pd_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("phase_diagrams")


@pytest.fixture
def relaxed_structure():
    atoms = read(TEST_DATA_DIR / "qmof-bda2f7d_relaxed.cif")
    return Structure.from_ase_atoms(atoms)  # type: ignore


@pytest.fixture
def unrelaxed_atoms():
    return read(TEST_DATA_DIR / "qmof-bda2f7d.cif")


def test_make_phase_diagram(pd_dir):
    structures_path = TEST_DATA_DIR / "test_reference_thermo_structures.json"
    thermo_path = TEST_DATA_DIR / "test_reference_thermo.json"
    setup_phase_diagrams(structures_path, thermo_path, output_dir=pd_dir)

    with open(pd_dir / "('C', 'H', 'N', 'O', 'Zn')_phase_diagram.json") as f:
        new_chem_space = json.load(f)

    with open(
        TEST_DATA_DIR
        / "test_phase_diagrams"
        / "('C', 'H', 'N', 'O', 'Zn')_phase_diagram.json"
    ) as f:
        test_chem_space = json.load(f)

    with open(pd_dir / "chemical_space_to_mpids.json") as f:
        new_chem_map = json.load(f)

    with open(
        TEST_DATA_DIR / "test_phase_diagrams" / "chemical_space_to_mpids.json"
    ) as f:
        test_chem_map = json.load(f)

    assert new_chem_space == test_chem_space
    assert new_chem_map == test_chem_map


def test_relax(unrelaxed_atoms, out_dir):
    struct, energy = relax_mof(unrelaxed_atoms, label="qmof-bda2f7d", out_dir=out_dir)
    assert struct.volume == pytest.approx(5284.412604266308)
    assert energy == pytest.approx(-1191.972703923097)


def test_energy_above_hull(relaxed_structure, pd_dir):
    energy = -1191.972703923097
    e_above_hull = get_energy_above_hull(relaxed_structure, energy, pd_dir)
    assert e_above_hull == pytest.approx(0.1921294352092806)
