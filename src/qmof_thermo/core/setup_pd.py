from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import pandas as pd
from pymatgen.core import Structure
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from monty.serialization import dumpfn
from logging import getLogger

LOGGER = getLogger(__name__)


@dataclass
class HullEntry:
    """Container for a single reference hull entry."""

    mpid: str
    structure: Structure
    energy: float  # total energy (eV)
    elements: frozenset[str]  # e.g. set ({"Ba", "O", "V"})


def chemical_space_from_structure(struct: Structure) -> set[str]:
    """Return a set of element symbols in the structure."""
    return frozenset(str(el.symbol) for el in struct.composition.elements)


def load_hull_entries(
    structures_path: Path,
    thermo_path: Path,
    mpid_key: str = "mpid",
    energy_key: str = "energy_total",
    ehull_key: str = "energy_above_hull",
) -> list[HullEntry]:
    """
    Load all hull entries (energy_above_hull == 0) with structures and energies.
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
    hull_energy_lookup: Dict[str, float] = dict(
        zip(hull_df[mpid_key], hull_df[energy_key])
    )

    # Build mpid -> structure mapping
    struct_lookup: Dict[str, Structure] = {}
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
    all_entries: List[HullEntry] = []
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


def build_phase_diagrams_by_space(
    entries: List[HullEntry], output_dir: Path
) -> None:
    """
    For each unique chemical space S (set of elements), build a PhaseDiagram
    using *all* hull entries whose element set is a subset of S.

    Also writes a chemical_space_to_mpids.json mapping.
    """

    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. All chemical spaces present (by element set)
    all_spaces = {e.elements for e in entries}

    # Keep only multi-element spaces for now
    spaces = {s for s in all_spaces if len(s) > 1}

    LOGGER.info(f"Total unique chemical spaces (multi-element): {len(spaces)}")

    # 3. Build PhaseDiagram per chemical space
    chemical_space_to_mpids: Dict[str, List[str]] = {}

    all_entries = entries

    # start sorting through different chemical spaces
    for i, space in enumerate(
        sorted(spaces, key=lambda s: (len(s), sorted(s))), start=1
    ):
        space_tuple = tuple(sorted(space))
        LOGGER.info(
            f"Building PhaseDiagram for chemical space {space_tuple} "
            f"({i}/{len(spaces)})..."
        )

        entries_for_space: list[PDEntry] = []
        mpids_for_space: list[str] = []

        # add entries/mpids to lists
        for e in all_entries:
            if e.elements.issubset(space):
                entries_for_space.append(PDEntry(e.structure.composition, e.energy))
                mpids_for_space.append(e.mpid)

        if len(entries_for_space) < len(space):
            # Not enough entries to possibly have all elemental references
            LOGGER.info(
                f"  Skipping {space_tuple}: only {len(entries_for_space)} entries, "
                f"need at least {len(space)} for elemental refs."
            )
            continue

        # Check we have pure element references for each element in the space
        elemental_refs_present = set()
        for pd_entry in entries_for_space:
            comp = pd_entry.composition
            if len(comp.elements) == 1:
                elemental_refs_present.add(str(comp.elements[0].symbol))

        missing_elements = [el for el in space if el not in elemental_refs_present]
        if missing_elements:
            LOGGER.debug(
                f"  Skipping {space_tuple}: missing elemental references for "
                f"{missing_elements}"
            )
            continue

        # Try to build PhaseDiagram
        try:
            pd = PhaseDiagram(entries_for_space)
        except Exception as exc:
            LOGGER.info(f"  Error constructing PhaseDiagram for {space_tuple}: {exc}")
            continue

        # Save PD as JSON
        filename = f"{space_tuple}_phase_diagram.json"
        pd_path = output_dir / filename
        dumpfn(pd, pd_path)
        LOGGER.info(f"Saved PhaseDiagram to: {pd_path}")

        # Store mapping for later use
        chemical_space_to_mpids[str(space_tuple)] = sorted(set(mpids_for_space))

    # Save chemical space → MPIDs mapping
    mapping_path = output_dir / "chemical_space_to_mpids.json"
    with mapping_path.open("w") as f:
        json.dump(chemical_space_to_mpids, f, indent=2)
    LOGGER.info(f"\nSaved chemical space to MPIDs mapping to: {mapping_path}")