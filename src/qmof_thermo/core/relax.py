from __future__ import annotations

import json
from logging import getLogger
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.optimize import BFGS
from fairchem.core import FAIRChemCalculator
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

if TYPE_CHECKING:
    from ase import Atoms
    from pymatgen.core import Structure

LOGGER = getLogger(__name__)


def run_calc(
    atoms: Atoms,
    label: None | str = None,
    model: str | Path = "uma-s-1p1",
    uma_task_name: str | None = "odac",
    fmax: float = 0.03,
    max_steps: int = 10000,
    device: str = "cuda",
    out_dir: Path | str = Path("data/relaxations"),
) -> tuple[Structure, float]:
    """
    Relax an ASE Atoms structure using a FAIRChem MLIP calculator.

    Performs a full relaxation of both atomic positions and cell parameters
    using the BFGS optimizer with a FrechetCellFilter. Outputs are saved to
    ``<out_dir>/<label>/`` including optimization log, trajectory, relaxed CIF file,
    and a JSON summary of results.

    Parameters
    ----------
    atoms
        ASE Atoms object to be relaxed.
    label
        Unique identifier for the relaxation job. Used to name output files.
    model
        Model name or path to checkpoint file. For UMA models, pass the model
        name (e.g., ``"uma-s-1p1"``) to load from Hugging Face, or provide a
        local checkpoint path. For eSEN models, pass the local checkpoint
        file path.
    uma_task_name
        Task name for UMA models. Set to "odac" for UMA models, or
        ``None`` for non-UMA models like eSEN.
    fmax
        Convergence criteria, set maximum force on atoms in in eV/Å. Optimization
        stops when maximum force on any atom falls below this value.
    max_steps
        Maximum number of optimization steps allowed.
    device
        Device to run calculation on, e.g., "cpu" or "cuda".
    out_dir
        Base directory for output files. Subdirectory ``<label>``
        created and stores all relaxation specific outputs.

    Returns
    -------
    tuple[Structure, float]
        - The final relaxed structure as a pymatgen Structure object.
        - The final relaxed total energy in eV.

    Notes
    -----
    Following files are written to ``<out_dir>/<label>/``:
        - ``opt.log``: Optimization log
        - ``opt.traj``: Trajectory file with all optimization steps
        - ``<label>.cif``: Final relaxed structure in CIF format
        - ``results.json``: Summary including energy, volume, forces, and steps
    """
    if label is None:
        label = "output"

    atoms.calc = FAIRChemCalculator.from_model_checkpoint(
        name_or_path=model,
        task_name=uma_task_name if "uma" in str(model) else None,
        device=device,
    )

    filter_atoms = FrechetCellFilter(atoms)

    out_dir = Path(out_dir) / label
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "opt.log"
    traj_path = out_dir / "opt.traj"

    opt = BFGS(filter_atoms, trajectory=traj_path, logfile=log_path)  # type: ignore
    opt.run(fmax=fmax, steps=max_steps)  # runs the optimization until max|F| <= fmax

    final_forces = atoms.get_forces()
    final_fmax = float(np.max(np.linalg.norm(final_forces, axis=1)))
    nsteps = opt.get_number_of_steps()
    final_volume = atoms.get_volume()
    final_energy = atoms.get_potential_energy()
    LOGGER.info(
        f"Energy: {final_energy}, Volume: {final_volume}, fmax: {final_fmax}, steps: {nsteps}"
    )

    final_atoms = read(traj_path, index=-1)  # last frame
    final_struct = AseAtomsAdaptor.get_structure(final_atoms)  # type: ignore
    cif_path = out_dir / f"{label}.cif"
    write(cif_path, final_atoms)
    LOGGER.info(f"Final relaxed structure written to: {cif_path}")

    summary = {
        "id": label,
        "model": model,
        "fmax_target": fmax,
        "max_steps": max_steps,
        "nsteps": nsteps,
        "final_energy": final_energy,
        "final_volume": final_volume,
        "final_fmax": final_fmax,
    }
    summary_path = out_dir / "results.json"
    with summary_path.open("w") as f:
        json.dump(summary, f, indent=2)
    LOGGER.info(f"Summary written to: {summary_path}")

    return final_struct, final_energy
