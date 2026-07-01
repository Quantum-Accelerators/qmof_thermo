from __future__ import annotations

import logging

from qmof_thermo.hull import get_energy_above_hull
from qmof_thermo.phase_diagram import setup_phase_diagrams
from qmof_thermo.relax import relax_mof

__all__ = [
    "get_energy_above_hull",
    "relax_mof",
    "set_log_level",
    "setup_phase_diagrams",
]

logger = logging.getLogger(__name__)


def set_log_level(level=logging.INFO):
    logger.setLevel(level)
