"""
Module for calculating energy above hull.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from ase import Atoms
from monty.serialization import loadfn
from pymatgen.analysis.phase_diagram import PatchedPhaseDiagram, PDEntry
from pymatgen.io.ase import AseAtomsAdaptor

from qmof_thermo.phase_diagram import _DEFAULT_PD_FILENAME

if TYPE_CHECKING:
    from pymatgen.core import Structure

_DEFAULT_PD_JSON = Path(__file__).parent.resolve() / _DEFAULT_PD_FILENAME


def _load_patched_phase_diagram(pd_dir: Path | str) -> PatchedPhaseDiagram:
    """
    Load a precomputed PatchedPhaseDiagram from disk.

    Parameters
    ----------
    pd_dir
        Directory containing the ``patched_phase_diagram.json`` file.

    Returns
    -------
    PatchedPhaseDiagram
        The loaded PatchedPhaseDiagram object.

    Raises
    ------
    FileNotFoundError
        If the PatchedPhaseDiagram JSON file cannot be found.
    """
    pd_dir = Path(pd_dir)
    pd_path = pd_dir / _DEFAULT_PD_FILENAME

    if not pd_path.is_file():
        raise FileNotFoundError(
            f"PatchedPhaseDiagram not found at {pd_path}. "
            "Run qmof_thermo.setup_phase_diagrams() to build it first."
        )

    return loadfn(pd_path)


def get_energy_above_hull(
    struct: Structure | Atoms,
    energy: float,
    references_dir: Path | str = _DEFAULT_PD_JSON,
) -> float:
    """
    Calculate the energy above hull for a structure with a given total energy.

    Parameters
    ----------
    struct
        Input structure as either a pymatgen Structure or ASE Atoms object.
        If an Atoms object is provided, it will be converted to a Structure.
    energy
        Total relaxed energy of the structure in eV.
    references_dir
        Path to the directory containing the precomputed PatchedPhaseDiagram.
        Default is ``"data/references"``.

    Returns
    -------
    float
        Energy above the convex hull in eV/atom.
    """
    if isinstance(struct, Atoms):
        struct = AseAtomsAdaptor.get_structure(struct)

    ppd = _load_patched_phase_diagram(references_dir)

    entry = PDEntry(struct.composition, energy)
    result = ppd.get_decomp_and_e_above_hull(entry)

    if result[1] is None:
        msg = (
            f"Could not compute energy above hull for composition "
            f"{struct.composition.reduced_formula}."
        )
        raise ValueError(msg)

    return float(result[1])
