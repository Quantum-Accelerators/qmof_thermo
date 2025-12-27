import json
from logging import getLogger
from pathlib import Path

from ase import Atoms
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.optimize import BFGS
from fairchem.core import FAIRChemCalculator, pretrained_mlip
from fairchem.core.units.mlip_unit import load_predict_unit
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core import Structure
import numpy as np

LOGGER = getLogger(__name__)


def run_calc(
    atoms: Atoms,
    id: str,
    model: str | Path = "uma-s-1p1",
    uma_task_name: str | None = "odac",
    fmax: float = 0.03,
    max_steps: int = 10000,
    device: str = "cuda",
    out_dir: Path | str = Path("data/relaxations"),
) -> tuple[Structure, float]:
    """
    Relax an ASE `Atoms` structure (positions + cell) with a FAIRChem MLIP and save outputs to
    `data/relaxations/<id>/` (log, traj, CIF, and a results.json summary).

    Model usage:
    - eSEN: pass the local checkpoint file path via `model_path`.
        - task_name
    - UMA: either pass a local checkpoint file, or load from Hugging Face by passing the UMA
    model name (e.g. "uma-s-1p1") to `pretrained_mlip.get_predict_unit(...)`.

    log_level - specify level of logging, either INFO or DEBUG

    Returns: (final_struct: pymatgen Structure, final_energy: float)
    """
    
    atoms.calc = FAIRChemCalculator.from_model_checkpoint(name_or_path=model, task_name=uma_task_name)

    filter_atoms = FrechetCellFilter(atoms)

    out_dir = out_dir / id
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "opt.log"
    traj_path = out_dir / "opt.traj"

    opt = BFGS(
        filter_atoms, trajectory=traj_path, logfile=log_path
    )  # selects the optimizer to use for the geometry optimization
    opt.run(fmax=fmax, steps=max_steps)  # runs the optimization until max|F| <= fmax

    final_forces = atoms.get_forces()
    final_fmax = np.max(np.linalg.norm(final_forces, axis=1)) 
    nsteps = opt.get_number_of_steps()
    final_volume = atoms.get_volume()
    final_energy = atoms.get_potential_energy()
    LOGGER.info(
        f"Energy: {final_energy}, Volume: {final_volume}, fmax: {final_fmax}, steps: {nsteps}"
    )

    try:
        final_atoms = read(traj_path, index=-1)  # last frame
        final_struct = AseAtomsAdaptor.get_structure(final_atoms)
        cif_path = out_dir / f"{id}.cif"
        write(cif_path, final_atoms)
        # optional: print path so you know where it went
        LOGGER.info(f"Final relaxed structure written to: {cif_path}")
    except Exception as e:
        LOGGER.debug(f"Warning: failed to write CIF from trajectory: {e}")

    summary = {
        "id": id,
        "model_ckpt": model_path,
        "fmax_target": fmax,
        "max_steps": max_steps,
        "nsteps": int(nsteps),
        "final_energy": float(final_energy),
        "final_volume": float(final_volume),
        "final_fmax": float(final_fmax),
    }
    summary_path = out_dir / "results.json"
    with summary_path.open("w") as f:
        json.dump(summary, f, indent=2)
    LOGGER.info(f"Summary written to: {summary_path}")

    return final_struct, final_energy
