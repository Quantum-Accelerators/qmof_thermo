import json
from pathlib import Path

import pytest
from ase.io import read
from pymatgen.core import Structure

from qmof_thermo.core import calc, relax, setup_pd

FILE_DIR = Path(__file__).parent

out_dir = FILE_DIR / "relaxations"
out_dir.mkdir(exist_ok=True)

pd_dir = FILE_DIR / "phase_diagrams"
pd_dir.mkdir(exist_ok=True)


@pytest.fixture()
def relaxed_structure():
    atoms = read(FILE_DIR / "qmof-bda2f7d_relaxed.cif")
    return Structure.from_ase_atoms(atoms)


@pytest.fixture()
def unrelaxed_atoms():
    return read(FILE_DIR / "qmof-bda2f7d.cif")


def test_make_phase_diagram():
    structures_path = FILE_DIR / "test_reference_thermo_structures.json"
    thermo_path = FILE_DIR / "test_reference_thermo.json"
    pd_dir = FILE_DIR / "phase_diagrams"
    setup_pd.setup_phase_diagrams(structures_path, thermo_path, pd_dir)

    new_chem_space = json.load(
        open(pd_dir / "('C', 'H', 'N', 'O', 'Zn')_phase_diagram.json")
    )
    test_chem_space = json.load(
        open(
            FILE_DIR
            / "test_phase_diagrams"
            / "('C', 'H', 'N', 'O', 'Zn')_phase_diagram.json"
        )
    )

    new_chem_map = json.load(open(pd_dir / "chemical_space_to_mpids.json"))
    test_chem_map = json.load(
        open(FILE_DIR / "test_phase_diagrams" / "chemical_space_to_mpids.json")
    )

    assert new_chem_space == test_chem_space
    assert new_chem_map == test_chem_map


def test_relax(unrelaxed_atoms):
    struct, energy = relax.run_calc(
        unrelaxed_atoms, id="qmof-bda2f7d", out_dir=out_dir, device="cpu"
    )
    assert struct
    assert energy == pytest.approx(-1191.972702)


def test_energy_above_hull(relaxed_structure):
    energy = -1191.972702
    e_above_hull = calc.energy_above_hull_from_structure(
        relaxed_structure, energy, pd_dir
    )
    assert e_above_hull == pytest.approx(0.19212944022861933)
