from __future__ import annotations

from pathlib import Path

import pytest
from ase.io import read
from monty.serialization import loadfn
from pymatgen.core import Structure

from qmof_thermo import get_energy_above_hull, relax_mof, setup_phase_diagrams
from qmof_thermo.phase_diagram import _DEFAULT_PD_FILENAME

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
    pd_path = pd_dir / _DEFAULT_PD_FILENAME
    assert pd_path.is_file()

    ppd = loadfn(pd_path)
    assert len(ppd.all_entries) > 0
    assert len(ppd.elements) > 0


def test_relax(unrelaxed_atoms, out_dir):
    struct, energy = relax_mof(
        unrelaxed_atoms, label="qmof-bda2f7d", fmax=0.03, out_dir=out_dir
    )
    assert struct.volume == pytest.approx(5284.412604266308)
    assert energy == pytest.approx(-1191.972703923097)


def test_energy_above_hull(relaxed_structure, pd_dir):
    energy = -1191.972703923097
    e_above_hull = get_energy_above_hull(relaxed_structure, energy, pd_dir)
    assert e_above_hull == pytest.approx(0.1921294352092806)
