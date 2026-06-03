"""
Module for calculating energy above hull.
"""

from __future__ import annotations

from logging import getLogger
from pathlib import Path
from typing import TYPE_CHECKING

from ase import Atoms
from monty.serialization import loadfn
from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.io.ase import AseAtomsAdaptor

from qmof_thermo.phase_diagram import _DEFAULT_PD_FILENAME, QMOF_COMPATIBLE_ELEMENTS, QMOF_ELEMENTS

LOGGER = getLogger(__name__)

if TYPE_CHECKING:
    from pymatgen.core import Structure

_DEFAULT_PD_JSON = Path(__file__).parent.resolve() / _DEFAULT_PD_FILENAME


def get_energy_above_hull(
    struct: Structure | Atoms,
    energy: float,
    serialized_phase_diagram: Path | str = _DEFAULT_PD_JSON,
    return_dict: bool = False,
) -> float | dict[str, float | dict[str, float]]:
    """
    Calculate the energy above hull for a structure with a given total energy.

    Parameters
    ----------
    struct
        Input structure as either a pymatgen Structure or ASE Atoms object.
        If an Atoms object is provided, it will be converted to a Structure.
    energy
        Total relaxed energy of the structure in eV.
    serialized_phase_diagram
        Path to the directory containing the precomputed PatchedPhaseDiagram.
        Default is the phase diagram bundled with the package.
    return_dict
        If True, return energy above hull, formation energy, and decomposition
        products in a dictionary. If False, return only energy above hull.

    Returns
    -------
    float or dict
        Energy above the convex hull in eV/atom, or a dictionary of
        thermodynamic results.
    """

    if isinstance(struct, Atoms):
        struct = AseAtomsAdaptor.get_structure(struct)

    mol_elements = {str(el) for el in struct.composition.elements}
    incompatible = mol_elements - QMOF_COMPATIBLE_ELEMENTS

    if incompatible:
        LOGGER.warning(
            "Structure contains elements whose UMA-ODAC "
            f"pseudopotentials do not match QMOF's: {sorted(incompatible)}. "
            "Relaxations ran via UMA-ODAC will not be comparable to QMOF DFT references."
        )

    out_of_chemical_space = mol_elements - QMOF_ELEMENTS
    if out_of_chemical_space:
        LOGGER.warning(
            "Structure contains elements not present in the QMOF chemical space: "
            f"{sorted(out_of_chemical_space)}. Relaxations possibly are not be comparable to QMOF references."
        )

    ppd = loadfn(serialized_phase_diagram)

    entry = PDEntry(struct.composition, energy)
    decomp, e_above_hull = ppd.get_decomp_and_e_above_hull(entry)

    if e_above_hull is None:
        msg = (
            f"Could not compute energy above hull for composition "
            f"{struct.composition.reduced_formula}."
        )
        raise ValueError(msg)

    if return_dict:
        return {
            "energy_above_hull": float(e_above_hull),
            "formation_energy": float(ppd.get_form_energy_per_atom(entry)),
            "decomposition_products": {
                decomp_entry.composition.reduced_formula: float(amount)
                for decomp_entry, amount in decomp.items()
            },
        }

    return float(e_above_hull)
