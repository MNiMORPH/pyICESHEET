"""Build a basal shear-stress field from a substrate/class map.

Rather than hand-drawing shear-stress polygons, the basal shear stress can be
*derived* from mapped geology: a classified substrate raster (e.g. the sediment /
grain-size / bedrock datasets of Gowan et al. 2019, ESSD 11, 375) plus a per-class
yield-stress table. Deformable sediment classes take a low shear stress (fast flow
/ ice streams); hard-bedrock classes take a high one.

This is a composable upstream operator: ``class map + per-class table -> tau``.
The per-class table is a natural handle for calibration (see
:mod:`pyicesheet.calibrate`).
"""

from __future__ import annotations

import numpy as np

__all__ = ["tau_from_classes"]


def tau_from_classes(class_ids, class_tau, nodata=np.nan):
    """Map an integer class raster to a shear-stress raster.

    Parameters
    ----------
    class_ids : 2-D int array
        Substrate class per cell, values ``1..K`` (``0`` or negative = no data).
    class_tau : mapping or 1-D array
        Shear stress (Pa) for each class. If an array, class ``k`` uses
        ``class_tau[k - 1]``; if a mapping, ``class_tau[k]``.
    nodata : float
        Value for cells whose class is missing.

    Returns
    -------
    tau : 2-D float array
    """
    ci = np.asarray(class_ids)
    if hasattr(class_tau, "keys"):                       # mapping
        table = np.full(max(class_tau) + 1, np.nan)
        for k, v in class_tau.items():
            table[k] = v
        valid = ci >= 1
        tau = np.where(valid, table[np.clip(ci, 0, len(table) - 1)], nodata)
    else:
        table = np.asarray(class_tau, dtype=float)
        valid = (ci >= 1) & (ci <= table.size)
        tau = np.where(valid, table[np.clip(ci, 1, table.size) - 1], nodata)
    return tau
