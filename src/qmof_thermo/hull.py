"""
Module for calculating energy above hull.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import TYPE_CHECKING

from ase import Atoms
from monty.serialization import loadfn
from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.io.ase import AseAtomsAdaptor

from qmof_thermo.phase_diagram import QMOF_COMPATIBLE_ELEMENTS, _DEFAULT_PD_FILENAME

if TYPE_CHECKING:
    from pymatgen.core import Structure

_DEFAULT_PD_JSON = Path(__file__).parent.resolve() / _DEFAULT_PD_FILENAME


def get_energy_above_hull(
    struct: Structure | Atoms,
    energy: float,
    serialized_phase_diagram: Path | str = _DEFAULT_PD_JSON,
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

    mol_elements = {str(el) for el in struct.composition.elements}
    incompatible = mol_elements - QMOF_COMPATIBLE_ELEMENTS
    if incompatible:
        warnings.warn(
            f"Structure contains elements whose UMA-ODAC default "
            f"pseudopotentials do not match QMOF: {sorted(incompatible)}. "
            f"Energy predictions will not be comparable to QMOF DFT references.",
            stacklevel=2,
        )

    ppd = loadfn(serialized_phase_diagram)

    entry = PDEntry(struct.composition, energy)
    result = ppd.get_decomp_and_e_above_hull(entry)

    if result[1] is None:
        msg = (
            f"Could not compute energy above hull for composition "
            f"{struct.composition.reduced_formula}."
        )
        raise ValueError(msg)

    return float(result[1])
