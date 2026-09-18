"""Rest- and observed-frame wavelength bands for continuum / ratio summaries.

Shared constants align with plot limits in analyze_ehvo (3800–6000 Å obs) and the
1450 Å anchor used in convert_fits. Import REGIONS where you need the same masks.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import NamedTuple

import numpy as np

LYA_REST = 1215.67
LYA_EXCL = 15.0
WAVE_OBS_LO = 3800.0
WAVE_OBS_HI = 6000.0

MaskFn = Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]


class ContinuumRegion(NamedTuple):
    region_id: str
    label: str
    mask_fn: MaskFn


def _obs_window(wr: np.ndarray, wo: np.ndarray, good: np.ndarray) -> np.ndarray:
    return (wo >= WAVE_OBS_LO) & (wo <= WAVE_OBS_HI) & good


REGIONS: tuple[ContinuumRegion, ...] = (
    ContinuumRegion("full", "Full (obs 3800–6000 Å)", _obs_window),
    ContinuumRegion(
        "anchor_1450",
        "Continuum anchor (1450 Å rest)",
        lambda wr, wo, g: (wr > 1445) & (wr < 1455) & g,
    ),
    ContinuumRegion(
        "red_wing",
        "Red wing (1250–1350 Å rest)",
        lambda wr, wo, g: (wr > 1250) & (wr < 1350) & g,
    ),
    ContinuumRegion(
        "nv",
        "N V (1235–1245 Å rest)",
        lambda wr, wo, g: (wr > 1235) & (wr < 1245) & g,
    ),
    ContinuumRegion(
        "siiv",
        "Si IV (1395–1405 Å rest)",
        lambda wr, wo, g: (wr > 1395) & (wr < 1405) & g,
    ),
    ContinuumRegion(
        "civ",
        "C IV (1540–1560 Å rest)",
        lambda wr, wo, g: (wr > 1540) & (wr < 1560) & g,
    ),
    ContinuumRegion(
        "excl_lya",
        "Excl. Lyα ±15 Å (obs 3800–6000)",
        lambda wr, wo, g: _obs_window(wr, wo, g)
        & ~((wr > LYA_REST - LYA_EXCL) & (wr < LYA_REST + LYA_EXCL)),
    ),
)
