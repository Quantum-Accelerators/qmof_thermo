"""
Module for calculating energy above hull.
"""

from __future__ import annotations

from logging import getLogger
from pathlib import Path
from typing import TYPE_CHECKING, Literal, TypedDict

from ase import Atoms
from monty.serialization import loadfn
from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.io.ase import AseAtomsAdaptor

from qmof_thermo import _DEFAULT_PD_JSON, QMOF_ELEMENTS, QMOF_ODAC_COMPATIBLE_ELEMENTS

LOGGER = getLogger(__name__)

if TYPE_CHECKING:
    from pymatgen.core import Structure


class HullOutput(TypedDict):
    energy_above_hull: float
    formation_energy: float
    decomposition_products: dict[str, float]


def get_energy_above_hull(
    struct: Structure | Atoms,
    energy: float,
    *,
    energy_type: Literal["DFT", "ODAC_MLIP"],
    serialized_phase_diagram: Path | str = _DEFAULT_PD_JSON,
) -> HullOutput:
    """
    Calculate the energy above hull for a structure with a given total energy.

    Parameters
    ----------
    struct
        Input structure as either a pymatgen Structure or ASE Atoms object.
        If an Atoms object is provided, it will be converted to a Structure.
    energy
        Total relaxed energy of the structure in eV.
    energy_type
        How the energy data was sourced, either from "DFT" or "ODAC_MLIP"
    serialized_phase_diagram
        Path to the directory containing the precomputed PatchedPhaseDiagram.
        Default is the phase diagram bundled with the package.

    Returns
    -------
    dict
        A dictionary containing energy above the convex hull, formation energy, and decomposition products.
    """

    energy_type = energy_type.lower()
    if isinstance(struct, Atoms):
        struct = AseAtomsAdaptor.get_structure(struct)

    mol_elements = {str(el) for el in struct.composition.elements}

    out_of_chemical_space = mol_elements - QMOF_ELEMENTS
    if out_of_chemical_space:
        LOGGER.warning(
            "Your structure contains elements that are not present in the QMOF chemical space: "
            f"{sorted(out_of_chemical_space)}. The convex hull phase diagram will be incomplete."
        )

    if energy_type == "odac_mlip" and (
        incompatible := mol_elements - QMOF_ODAC_COMPATIBLE_ELEMENTS
    ):
        LOGGER.warning(
            "Structure contains elements whose ODAC "
            f"pseudopotentials do not match QMOF's: {sorted(incompatible)}. "
            "If `energy` is obtained from an ODAC MLIP, the results will likely not be compatible with "
            "the QMOF DFT reference data."
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

    return {
        "energy_above_hull": float(e_above_hull),
        "formation_energy": float(ppd.get_form_energy_per_atom(entry)),
        "decomposition_products": {
            decomp_entry.composition.reduced_formula: float(amount)
            for decomp_entry, amount in decomp.items()
        },
    }
