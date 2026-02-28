"""
Module for calculating energy above hull.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING, Any

from ase import Atoms
from monty.serialization import loadfn
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.io.ase import AseAtomsAdaptor

if TYPE_CHECKING:
    from pymatgen.core import Structure


def get_energy_above_hull(
    struct: Structure | Atoms,
    energy: float,
    references_dir: Path | str = Path("data/references"),
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
        Path to the directory containing precomputed phase diagram
        references. Default is ``"data/references"``.

    Returns
    -------
    float
        Energy above the convex hull in eV/atom.
    """
    if isinstance(struct, Atoms):
        struct = AseAtomsAdaptor.get_structure(struct)

    space = _chemical_space_from_structure(struct)
    pd_obj = _load_phase_diagram_for_space(space, references_dir)

    entry = PDEntry(struct.composition, energy)
    _, ehull = pd_obj.get_decomp_and_e_above_hull(entry)

    return float(ehull)


def _chemical_space_from_structure(struct: Structure) -> tuple[str, ...]:
    """
    Extract the chemical space from a structure as a sorted tuple of element symbols.

    Parameters
    ----------
    struct
        Pymatgen Structure object from which to extract the chemical space.

    Returns
    -------
    tuple[str, ...]
        Sorted tuple of element symbols present in the structure's composition.
    """
    return tuple(sorted(el.symbol for el in struct.composition.elements))


def _load_phase_diagram_for_space(
    space: tuple[str, ...],
    pd_dir: Path | str,
    mapping_filename: str = "chemical_space_to_mpids.json",
) -> PhaseDiagram:
    """
    Load a precomputed PhaseDiagram for an exact chemical space.

    Parameters
    ----------
    space
        Sorted tuple of element symbols defining the chemical space,
        e.g., ``('Ba', 'O', 'V')``.
    pd_dir
        Directory containing the precomputed phase diagram JSON files
        and the mapping file.
    mapping_filename
        Name of the JSON file mapping chemical spaces to material IDs.
        Default is ``"chemical_space_to_mpids.json"``.

    Returns
    -------
    PhaseDiagram
        Pymatgen PhaseDiagram object for the specified chemical space.

    Raises
    ------
    FileNotFoundError
        If the mapping file or the phase diagram JSON file for the specified
        space cannot be found.
    ValueError
        If the specified chemical space is not present in the mapping file.
    """

    pd_dir = Path(pd_dir)
    mapping_path = pd_dir / mapping_filename

    if not mapping_path.is_file():
        raise FileNotFoundError(
            f"Could not find mapping file at {mapping_path}. "
            "Run `qmof_thermo.phase_diagram.setup_phase_diagrams` to build the reference phase diagrams."
        )

    with mapping_path.open() as f:
        space_mapping: dict[str, Any] = json.load(f)

    key = str(space)  # e.g. "('Ba', 'O', 'V')"
    if key not in space_mapping:
        raise ValueError(
            f"No phase diagram found for chemical space {space}. "
            f"Known spaces: {len(space_mapping)}"
        )

    pd_filename = f"{key}_phase_diagram.json"
    pd_path = pd_dir / pd_filename

    if not pd_path.is_file():
        raise FileNotFoundError(
            f"PhaseDiagram JSON for space {space} not found at {pd_path}. "
            "Make sure `qmof_thermo.phase_diagram.setup_phase_diagrams` finished successfully."
        )

    pd_obj: PhaseDiagram = loadfn(pd_path)
    return pd_obj
