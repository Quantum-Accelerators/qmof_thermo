from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

def set_log_level(level=logging.INFO):
    logger.setLevel(level)
