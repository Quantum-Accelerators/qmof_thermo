from qmof_thermo import setup_pd
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)

structures_path = Path("data/external/reference_thermo_structures.json")
thermo_path = Path("data/external/reference_thermo.json")
output_dir = Path("data/references")

hull_entries = setup_pd.load_hull_entries(
    structures_path=structures_path,
    thermo_path=thermo_path,
    mpid_key="mpid",
    energy_key="energy_total",
    ehull_key="energy_above_hull"
)

setup_pd.build_phase_diagrams_by_space(
    entries=hull_entries,
    output_dir=output_dir,
)