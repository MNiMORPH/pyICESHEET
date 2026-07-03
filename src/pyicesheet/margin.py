"""The ice margin: the boundary condition for the reconstruction.

The mapped margin is the outermost (near-zero-thickness) contour of the ice
sheet. :class:`IceMargin` wraps a closed polygon, resamples its boundary to an
even point spacing, and computes the inward-pointing normal at each point — the
initial along-flow direction from which flowlines are integrated inward.
"""

from __future__ import annotations

import numpy as np
from shapely.geometry import Polygon

from ._geometry import resample_ring, inward_normals

__all__ = ["IceMargin"]


class IceMargin:
    """A resampled, oriented ice-margin loop.

    Attributes
    ----------
    x, y : ndarray
        Resampled boundary point coordinates (map units).
    inward_x, inward_y : ndarray
        Components of the inward-pointing unit normal at each point (the initial
        along-flow direction for the flowline seeded there).
    direction : ndarray
        Azimuth of the inward normal, ``atan2(inward_y, inward_x)`` (radians),
        for convenience.
    polygon : shapely.geometry.Polygon
        The original polygon (used for inside/outside tests).
    spacing : float
        Target point spacing used for resampling (map units).
    """

    def __init__(self, x, y, inward_x, inward_y, polygon, spacing):
        self.x = np.asarray(x, dtype=float)
        self.y = np.asarray(y, dtype=float)
        self.inward_x = np.asarray(inward_x, dtype=float)
        self.inward_y = np.asarray(inward_y, dtype=float)
        self.direction = np.arctan2(self.inward_y, self.inward_x)
        self.polygon = polygon
        self.spacing = float(spacing)

    def __len__(self):
        return self.x.size

    @classmethod
    def from_polygon(cls, polygon: Polygon, spacing: float):
        """Resample a polygon's exterior to ``spacing`` and orient inward.

        Parameters
        ----------
        polygon : shapely.geometry.Polygon
            Closed margin polygon (map units). Its exterior ring is used.
        spacing : float
            Target spacing between resampled boundary points (map units).
        """
        if not isinstance(polygon, Polygon):
            raise TypeError("IceMargin.from_polygon expects a shapely Polygon")
        x, y = resample_ring(polygon, spacing)
        inx, iny = inward_normals(x, y, polygon, spacing)
        return cls(x, y, inx, iny, polygon, spacing)

    def seed_points(self):
        """Yield ``(x, y, direction)`` tuples for each margin point."""
        for xi, yi, di in zip(self.x, self.y, self.direction):
            yield float(xi), float(yi), float(di)
