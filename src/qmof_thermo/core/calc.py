# src/qmof_thermo/core/e_above_hull.py

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Tuple

from monty.serialization import loadfn
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.core import Structure


def _chemical_space_from_structure(struct: Structure) -> Tuple[str, ...]:
    """Sorted tuple of elements,"""
    return tuple(sorted(el.symbol for el in struct.composition.elements))


def _load_phase_diagram_for_space(
    space: Tuple[str, ...],
    references_dir: Path | str = "data/references",
    mapping_filename: str = "chemical_space_to_mpids.json",
) -> PhaseDiagram:
    """
    Load precomputed PhaseDiagram for exact chemical space.

    `setup_pd.py` creates:
      - data/references/chemical_space_to_mpids.json
      - data/references/"('Ba', 'O', 'V')_phase_diagram.json"
    """
    references_dir = Path(references_dir)
    mapping_path = references_dir / mapping_filename

    if not mapping_path.is_file():
        raise FileNotFoundError(
            f"Could not find mapping file at {mapping_path}. "
            "Run setup_pd.py to build the reference phase diagrams."
        )

    with mapping_path.open() as f:
        space_mapping: Dict[str, Any] = json.load(f)

    key = str(space)  # e.g. "('Ba', 'O', 'V')"
    if key not in space_mapping:
        raise ValueError(
            f"No phase diagram found for chemical space {space}. "
            f"Known spaces: {len(space_mapping)}"
        )

    pd_filename = f"{key}_phase_diagram.json"
    pd_path = references_dir / pd_filename

    if not pd_path.is_file():
        raise FileNotFoundError(
            f"PhaseDiagram JSON for space {space} not found at {pd_path}. "
            "Make sure setup_pd.py finished successfully."
        )

    pd_obj: PhaseDiagram = loadfn(pd_path)
    return pd_obj


def energy_above_hull_from_structure(
    struct: Structure,
    energy: float,
    references_dir: Path | str = "data/references",
) -> float:
    """
    Energy above hull (eV/atom) for a pymatgen Structure with a given total energy.

    struct -> input structure
    energy -> total relaxed energy of such structure
    references_dir -> filled by setup_pd.py
    """

    space = _chemical_space_from_structure(struct)
    pd_obj = _load_phase_diagram_for_space(space, references_dir)

    entry = PDEntry(struct.composition, energy)
    _, ehull = pd_obj.get_decomp_and_e_above_hull(entry)

    return float(ehull)
