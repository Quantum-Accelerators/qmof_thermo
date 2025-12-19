from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
from monty.json import MontyEncoder
from pymatgen.core import Structure
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from monty.serialization import dumpfn


@dataclass
class HullEntry:
    """Container for a single reference hull entry."""
    mpid: str
    structure: Structure
    energy: float           # total energy (eV)
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
) -> List[HullEntry]:
    """
    Load all hull entries (energy_above_hull == 0) with structures and energies.
    """

    print(f"Loading structures from: {structures_path}")
    with structures_path.open() as f:
        struct_records = json.load(f)
    print(f"Loaded {len(struct_records)} structure records.")

    print(f"Loading thermo data from: {thermo_path}")
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
    print(f"Found {len(hull_mpids)} hull MPIDs with {ehull_key} == 0.")
    print(f"Using {len(hull_mpids)} MPIDs as reference hull entries.")

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
                print(f"Warning: structure missing for {mpid}, skipping.")
            elif missing_struct_count == 11:
                print("Further missing structure warnings suppressed...")
            continue

        struct_obj = rec["structure"]
        if isinstance(struct_obj, dict):
            struct = Structure.from_dict(struct_obj)
        elif isinstance(struct_obj, Structure):
            struct = struct_obj
        else:
            print(
                f"Warning: structure for {mpid} is not a dict or Structure "
                f"(type={type(struct_obj)}), skipping."
            )
            continue

        struct_lookup[mpid] = struct

    print(
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

    print(
        f"\nTotal hull entries with both energy and structure: {used_count}\n"
    )

    return all_entries


def build_phase_diagrams_by_space(
    entries: List[HullEntry],
    output_dir: Path,
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

    print(f"Total unique chemical spaces (multi-element): {len(spaces)}")

    mpid_to_entry: Dict[str, HullEntry] = {e.mpid: e for e in entries}

    # 3. Build PhaseDiagram per chemical space
    chemical_space_to_mpids: Dict[str, List[str]] = {}

    all_entries = entries

    # start sorting through different chemical spaces 
    for i, space in enumerate(sorted(spaces, key=lambda s: (len(s), sorted(s))), start=1):
        space_tuple = tuple(sorted(space))  # for pretty printing / filenames
        print(
            f"Building PhaseDiagram for chemical space {space_tuple} "
            f"({i}/{len(spaces)})..."
        )

        # set up dicts, entries in chem space + mpid in chem space
        entries_for_space: List[PDEntry] = []
        mpids_for_space: List[str] = []

        # add entries/mpids to lists
        for e in all_entries:
            if e.elements.issubset(space):
                entries_for_space.append(PDEntry(e.structure.composition, e.energy))
                mpids_for_space.append(e.mpid)

        if len(entries_for_space) < len(space):
            # Not enough entries to possibly have all elemental references
            print(
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
            print(
                f"  Skipping {space_tuple}: missing elemental references for "
                f"{missing_elements}"
            )
            continue

        # Try to build PhaseDiagram
        try:
            pd = PhaseDiagram(entries_for_space)
        except Exception as exc:
            print(
                f"  Error constructing PhaseDiagram for {space_tuple}: {exc}"
            )
            continue

        # Save PD as JSON
        filename = f"{space_tuple}_phase_diagram.json"
        pd_path = output_dir / filename
        dumpfn(pd, pd_path)
        print(f"  Saved PhaseDiagram to: {pd_path}")

        # Store mapping for later use
        chemical_space_to_mpids[str(space_tuple)] = sorted(set(mpids_for_space))

    # Save chemical space → MPIDs mapping
    mapping_path = output_dir / "chemical_space_to_mpids.json"
    with mapping_path.open("w") as f:
        json.dump(chemical_space_to_mpids, f, indent=2)
    print(
        f"\nSaved chemical space → MPIDs mapping to: {mapping_path}"
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build reference PhaseDiagram objects for each chemical space "
            "from reference_thermo_structures.json and reference_thermo.json."
        )
    )
    parser.add_argument(
        "--structures_json",
        type=str,
        default="data/external/reference_thermo_structures.json",
        help="Path to data/external/reference_thermo_structures.json",
    )
    parser.add_argument(
        "--thermo_json",
        type=str,
        default="data/external/reference_thermo.json",
        help="Path to data/external/reference_thermo.json",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="data/references",
        help="Directory to write phase diagram JSONs and mapping file.",
    )
    parser.add_argument(
        "--mpid-key",
        type=str,
        default="mpid",
        help="Key name for MPID in both JSON files (default: 'mpid').",
    )
    parser.add_argument(
        "--energy-key",
        type=str,
        default="energy_total",
        help="Column name for total energy in thermo JSON (default: 'energy_total').",
    )
    parser.add_argument(
        "--ehull-key",
        type=str,
        default="energy_above_hull",
        help="Column name for energy above hull in thermo JSON "
             "(default: 'energy_above_hull').",
    )

    args = parser.parse_args()

    structures_path = Path(args.structures_json).resolve()
    thermo_path = Path(args.thermo_json).resolve()
    output_dir = Path(args.output_dir).resolve()

    hull_entries = load_hull_entries(
        structures_path=structures_path,
        thermo_path=thermo_path,
        mpid_key=args.mpid_key,
        energy_key=args.energy_key,
        ehull_key=args.ehull_key,
    )

    build_phase_diagrams_by_space(
        entries=hull_entries,
        output_dir=output_dir,
    )


if __name__ == "__main__":
    main()


"""
i need you to help make a readme for my github repository.

the main purpose is for one to easily be able to pass in a Atoms strucutre of a MOF, optionally run a relaxation via an MLIP method, and then to be able to get it's energy above hull value in one line.

There are many steps needed to set up. the first one is to install this package as an editable package. it's located here: https://github.com/Quantum-Accelerators/qmof_thermo

After this, I need them to load in the qmof-thermo database files. These can be accessed on https://doi.org/10.6084/m9.figshare.13147324). One needs to download the qmof-thermo database and place the reference_thermo_structures.json and reference_structures.json in the the data/external/ folder. 
"""