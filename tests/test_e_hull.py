from __future__ import annotations

from pathlib import Path

import pytest
from ase import Atoms
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
    atoms = unrelaxed_atoms.copy()
    energy = relax_mof(atoms, label="qmof-bda2f7d", fmax=0.03, out_dir=out_dir)
    assert atoms.get_volume() != unrelaxed_atoms.get_volume()
    assert atoms.get_volume() == pytest.approx(5284.412604266308)
    assert energy == pytest.approx(-1191.972703923097)


def test_energy_above_hull_default(relaxed_structure):
    energy = -1191.972703923097
    thermo = get_energy_above_hull(relaxed_structure, energy)
    assert thermo["energy_above_hull"] == pytest.approx(0.1921294352092806)
    assert isinstance(thermo["formation_energy"], float)
    assert isinstance(thermo["decomposition_products"], dict)


def test_energy_above_hull(relaxed_structure, pd_dir):
    energy = -1191.972703923097
    e_above_hull = get_energy_above_hull(
        relaxed_structure,
        energy,
        serialized_phase_diagram=pd_dir / _DEFAULT_PD_FILENAME,
    )["energy_above_hull"]
    assert e_above_hull == pytest.approx(0.1921294352092806)


def test_warning_for_ooqmof_chem_space(caplog):
    atoms = Atoms("He", positions=[[0, 0, 0]], cell=[3, 3, 3], pbc=True)
    energy = -1.0

    with caplog.at_level("WARNING"), pytest.raises(
        ValueError, match="Unable to get decomposition"
    ):
        get_energy_above_hull(atoms, energy)

    assert "elements that are not present in the QMOF chemical space" in caplog.text
