"""The reconstructed ice surface: the output of a solve.

:class:`IceSurface` holds the point cloud of reconstructed surface samples
``(x, y, elevation, thickness, bed)`` accumulated over all contours, and gridding
helpers to turn it into rasters.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import griddata

__all__ = ["IceSurface"]


@dataclass
class IceSurface:
    """Reconstructed surface as a scattered point cloud plus gridding helpers.

    Attributes
    ----------
    x, y : ndarray
        Sample coordinates (map units).
    elevation : ndarray
        Reconstructed ice-surface elevation at each sample (m).
    thickness : ndarray
        Ice thickness ``elevation - bed`` (m).
    bed : ndarray
        Bed elevation at each sample (m).
    """

    x: np.ndarray
    y: np.ndarray
    elevation: np.ndarray
    thickness: np.ndarray
    bed: np.ndarray

    def __len__(self):
        return self.x.size

    @property
    def bounds(self):
        """``(xmin, ymin, xmax, ymax)`` of the sample cloud."""
        return (float(self.x.min()), float(self.y.min()),
                float(self.x.max()), float(self.y.max()))

    def to_grid(self, resolution, bounds=None, value="thickness", method="linear"):
        """Grid a chosen field onto a regular raster.

        Parameters
        ----------
        resolution : float
            Grid cell size (map units).
        bounds : (xmin, ymin, xmax, ymax), optional
            Extent to grid; defaults to the sample bounds.
        value : {"thickness", "elevation", "bed"}
            Which field to grid.
        method : str
            ``scipy.interpolate.griddata`` method.

        Returns
        -------
        grid : 2-D ndarray
            Rows north-to-south (top-down), NaN outside the sample hull.
        x_coords, y_coords : 1-D ndarray
            Cell-centre coordinates (``y_coords`` descending).
        """
        if bounds is None:
            bounds = self.bounds
        xmin, ymin, xmax, ymax = bounds
        xc = np.arange(xmin + resolution / 2, xmax, resolution)
        yc = np.arange(ymin + resolution / 2, ymax, resolution)
        XX, YY = np.meshgrid(xc, yc)
        values = getattr(self, value)
        grid = griddata((self.x, self.y), values, (XX, YY), method=method)
        # Return top-down (north first) to match raster conventions.
        return grid[::-1, :], xc, yc[::-1]

    def to_geodataframe(self, crs=None):
        """Return the sample cloud as a GeoDataFrame (requires geopandas)."""
        import geopandas as gpd
        from shapely.geometry import Point

        geom = [Point(xi, yi) for xi, yi in zip(self.x, self.y)]
        return gpd.GeoDataFrame(
            {
                "elevation": self.elevation,
                "thickness": self.thickness,
                "bed": self.bed,
            },
            geometry=geom,
            crs=crs,
        )
