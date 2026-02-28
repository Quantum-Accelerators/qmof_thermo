from __future__ import annotations

import json
from dataclasses import dataclass
from logging import getLogger
from pathlib import Path

import pandas as pd
from monty.serialization import dumpfn
from pymatgen.analysis.phase_diagram import PatchedPhaseDiagram, PDEntry
from pymatgen.core import Structure

LOGGER = getLogger(__name__)

DEFAULT_PD_FILENAME = "patched_phase_diagram.json"


@dataclass
class HullEntry:
    """
    Container for a single reference hull entry.

    Attributes
    ----------
    mpid
        Materials Project ID for this entry.
    structure
        Pymatgen Structure object for this material.
    energy
        Total energy in eV.
    elements
        Frozenset of element symbols present in the structure.
    """

    mpid: str
    structure: Structure
    energy: float  # total energy (eV)
    elements: frozenset[str]  # e.g. set ({"Ba", "O", "V"})


def chemical_space_from_structure(struct: Structure) -> set[str]:
    """
    Extract the chemical space from a structure as a frozenset of element symbols.

    Parameters
    ----------
    struct
        Pymatgen Structure object.

    Returns
    -------
    frozenset[str]
        Frozenset of element symbols present in the structure's composition.
    """
    return frozenset(str(el.symbol) for el in struct.composition.elements)


def _load_hull_entries(
    structures_path: Path,
    thermo_path: Path,
    mpid_key: str = "mpid",
    energy_key: str = "energy_total",
    ehull_key: str = "energy_above_hull",
) -> list[HullEntry]:
    """
    Load all hull entries (energy_above_hull == 0) with structures and energies.

    Reads structure and thermodynamic data from separate JSON files, filters
    for materials on the convex hull, and returns a list of HullEntry objects
    containing matched data.

    Parameters
    ----------
    structures_path
        Path to a JSON file containing structure records. Each record should
        have an ID field (matching ``mpid_key``) and a ``"structure"`` field
        containing a Structure object.
    thermo_path
        Path to a JSON file containing thermodynamic data.
        Must include columns of ID, total energy, and energy above hull.
    mpid_key
        Column/key for the material ID in both data sources.
    energy_key
        Column name for total energy (eV) in thermo data.
    ehull_key
        Column name for energy above hull (eV) in thermo data.

    Returns
    -------
    list[HullEntry]
        List of HullEntry objects for all valid materials
        with ``energy_above_hull = 0``.

    Raises
    ------
    KeyError
        If required columns (``mpid_key``, ``energy_key``, or ``ehull_key``)
        not found in the thermo JSON.
    """

    LOGGER.info(f"Loading structures from: {structures_path}")
    with structures_path.open() as f:
        struct_records = json.load(f)
    LOGGER.info(f"Loaded {len(struct_records)} structure records.")

    LOGGER.info(f"Loading thermo data from: {thermo_path}")
    df = pd.read_json(thermo_path)

    if ehull_key not in df.columns:
        raise KeyError(f"Column '{ehull_key}' not found in thermo JSON.")
    if energy_key not in df.columns:
        raise KeyError(f"Column '{energy_key}' not found in thermo JSON.")
    if mpid_key not in df.columns:
        raise KeyError(f"Column '{mpid_key}' not found in thermo JSON.")

    # Only hull entries
    hull_df = df[df[ehull_key] == 0].copy()
    hull_mpids = hull_df[mpid_key].tolist()
    LOGGER.info(f"Found {len(hull_mpids)} hull MPIDs with {ehull_key} == 0.")
    LOGGER.info(f"Using {len(hull_mpids)} MPIDs as reference hull entries.")

    # Lookup from mpid -> energy_total
    hull_energy_lookup: dict[str, float] = dict(
        zip(hull_df[mpid_key], hull_df[energy_key], strict=True)
    )

    # Build mpid -> structure mapping
    struct_lookup: dict[str, Structure] = {}
    missing_struct_count = 0

    for rec in struct_records:
        if mpid_key not in rec:
            continue
        mpid = rec[mpid_key]
        if mpid not in hull_energy_lookup:
            # not on hull; we don't need this entry
            continue

        if "structure" not in rec:
            missing_struct_count += 1
            if missing_struct_count <= 10:
                LOGGER.warning(f"Warning: structure missing for {mpid}, skipping.")
            elif missing_struct_count == 11:
                LOGGER.warning("Further missing structure warnings suppressed...")
            continue

        struct_obj = rec["structure"]
        if isinstance(struct_obj, dict):
            struct = Structure.from_dict(struct_obj)
        elif isinstance(struct_obj, Structure):
            struct = struct_obj
        else:
            LOGGER.warning(
                f"Warning: structure for {mpid} is not a dict or Structure "
                f"(type={type(struct_obj)}), skipping."
            )
            continue

        struct_lookup[mpid] = struct

    LOGGER.info(
        f"Structures available for {len(struct_lookup)} "
        f"of {len(hull_mpids)} hull MPIDs."
    )

    # Assemble final HullEntry list
    all_entries: list[HullEntry] = []
    used_count = 0
    for mpid in hull_mpids:
        if mpid not in struct_lookup:
            continue

        struct = struct_lookup[mpid]
        energy = float(hull_energy_lookup[mpid])
        elements = chemical_space_from_structure(struct)

        all_entries.append(HullEntry(mpid, struct, energy, elements))
        used_count += 1

    LOGGER.info(f"\nTotal hull entries with both energy and structure: {used_count}\n")

    return all_entries


def setup_phase_diagrams(
    structures_path: str | Path,
    thermo_path: str | Path,
    output_dir: str | Path,
    id_key: str = "mpid",
    energy_key: str = "energy_total",
    ehull_key: str = "energy_above_hull",
) -> None:
    """
    Load reference hull data and construct a PatchedPhaseDiagram.

    Builds a single PatchedPhaseDiagram from all stable compounds
    (energy_above_hull = 0) and saves it to disk. The PatchedPhaseDiagram
    internally partitions entries by chemical space for efficient
    energy-above-hull queries.

    Parameters
    ----------
    structures_path : str | Path
        Path to a JSON file containing structure records. Each record should
        have an ID field and a "structure" field (pymatgen Structure).
    thermo_path : str | Path
        Path to a JSON file containing thermo data.
        Must include columns for ID, total energy, and energy above hull.
    output_dir : str | Path
        Directory where the PatchedPhaseDiagram JSON file will be saved.
        Created if it does not exist.
    id_key : str, default "mpid"
        Column/key name for the material ID in both data sources.
    energy_key : str, default "energy_total"
        Column name for total energy (eV) in the thermo data.
    ehull_key : str, default "energy_above_hull"
        Column name for energy above hull (eV) in the thermo data.

    Returns
    -------
    None
        Outputs ``patched_phase_diagram.json`` to ``output_dir``.
    """
    structures_path = Path(structures_path)
    thermo_path = Path(thermo_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    hull_entries = _load_hull_entries(
        structures_path, thermo_path, id_key, energy_key, ehull_key
    )

    pd_entries = [PDEntry(e.structure.composition, e.energy) for e in hull_entries]

    LOGGER.info(f"Building PatchedPhaseDiagram from {len(pd_entries)} entries...")
    ppd = PatchedPhaseDiagram(pd_entries)
    n_elements = len(ppd.elements) if ppd.elements else 0
    LOGGER.info(
        f"PatchedPhaseDiagram built with {n_elements} elements "
        f"and {len(ppd)} chemical sub-spaces."
    )

    pd_path = output_dir / DEFAULT_PD_FILENAME
    dumpfn(ppd, pd_path)
    LOGGER.info(f"Saved PatchedPhaseDiagram to: {pd_path}")
