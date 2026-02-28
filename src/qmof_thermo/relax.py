"""
Module for relaxations.
"""

from __future__ import annotations

from logging import getLogger
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import torch
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.optimize import BFGS
from fairchem.core import FAIRChemCalculator
from fairchem.core.units.mlip_unit.api.inference import UMATask
from monty.serialization import dumpfn
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

if TYPE_CHECKING:
    from typing import Literal

    from ase import Atoms
    from ase.optimize.optimize import Optimizer
    from pymatgen.core import Structure

LOGGER = getLogger(__name__)


def relax_mof(
    atoms: Atoms,
    label: str = "output",
    model: str | Path = "uma-s-1p1",
    uma_task_name: UMATask | None = UMATask.ODAC,
    fmax: float = 0.01,
    max_steps: int = 10000,
    optimizer: type[Optimizer] = BFGS,
    device: Literal["cpu", "cuda"] | None = None,
    out_dir: Path | str = Path("data/relaxations"),
) -> tuple[Atoms, float]:
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
    tuple[Atoms, float]
        - The final relaxed structure as an ASE Atoms object.
        - The final relaxed total energy in eV.

    Notes
    -----
    Following files are written to ``<out_dir>/<label>/``:
        - ``opt.log``: Optimization log
        - ``opt.traj``: Trajectory file with all optimization steps
        - ``<label>.cif``: Final relaxed structure in CIF format
        - ``results.json``: Summary including energy, volume, forces, and steps
    """
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    atoms.calc = FAIRChemCalculator.from_model_checkpoint(
        name_or_path=model,
        task_name=uma_task_name if "uma" in str(model) else None,
        device=device,
    )

    filter_atoms = FrechetCellFilter(atoms)

    out_dir = (Path(out_dir) / label).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    traj_path = out_dir / "opt.traj"

    opt = optimizer(
        filter_atoms,  # type: ignore
        trajectory=traj_path,
    )
    opt.run(fmax=fmax, steps=max_steps)

    final_forces = atoms.get_forces()
    final_fmax = float(np.max(np.linalg.norm(final_forces, axis=1)))
    nsteps = opt.get_number_of_steps()
    final_volume = atoms.get_volume()
    final_energy = atoms.get_potential_energy()
    LOGGER.info(
        f"Energy: {final_energy}, Volume: {final_volume}, fmax: {final_fmax}, steps: {nsteps}"
    )

    final_atoms = read(traj_path, index=-1)
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
    dumpfn(summary, summary_path)
    LOGGER.info(f"Summary written to: {summary_path}")

    return final_atoms, final_energy
