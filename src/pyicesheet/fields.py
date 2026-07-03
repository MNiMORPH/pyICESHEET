"""Smoothly-sampleable raster fields (bed topography, basal shear stress).

The flowline integrator differentiates the bed and shear-stress fields, so it
needs *smooth* values and gradients at arbitrary points, not just at grid nodes.
The original Fortran ICESHEET achieved this with a hand-rolled bicubic
interpolation over a GMT binary grid. pyICESHEET uses
:class:`scipy.interpolate.RectBivariateSpline`, a bicubic (by default) spline on
a regular grid that provides analytic partial derivatives — cleaner, and
continuous through the second derivative.

Field *smoothing/regularization is deliberately not done here*: it is a separate
upstream step (e.g. GRASS ``r.resamp.bspline`` / ``r.mapcalc``). A
:class:`RasterField` assumes it is handed an already-smooth field to interpolate.

Coordinate convention
---------------------
A field is defined on a regular grid with strictly increasing 1-D coordinates
``x`` (columns) and ``y`` (rows), and a 2-D ``values`` array of shape
``(ny, nx)`` such that ``values[j, i]`` is the field at ``(x[i], y[j])``. Rasters
that are stored top-down (y decreasing) are flipped on construction by
:meth:`RasterField.from_arrays`.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import RectBivariateSpline

__all__ = ["RasterField"]


class RasterField:
    """A smooth scalar field on a regular grid, with analytic gradients.

    Parameters
    ----------
    x, y : 1-D array_like
        Strictly increasing grid coordinates (map units, e.g. metres).
    values : 2-D array_like, shape (ny, nx)
        Field values; ``values[j, i]`` is the value at ``(x[i], y[j])``.
    degree : int, optional
        Spline degree in each direction (default 3, bicubic).
    """

    def __init__(self, x, y, values, degree: int = 3):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        values = np.asarray(values, dtype=float)
        if x.ndim != 1 or y.ndim != 1:
            raise ValueError("x and y must be 1-D")
        if values.shape != (y.size, x.size):
            raise ValueError(
                f"values shape {values.shape} != (ny, nx) = {(y.size, x.size)}"
            )
        if not (np.all(np.diff(x) > 0) and np.all(np.diff(y) > 0)):
            raise ValueError("x and y must be strictly increasing")

        self.x = x
        self.y = y
        self.values = values
        self._xmin, self._xmax = float(x[0]), float(x[-1])
        self._ymin, self._ymax = float(y[0]), float(y[-1])
        # RectBivariateSpline wants z[i, j] = f(x[i], y[j]) -> transpose.
        self._spline = RectBivariateSpline(x, y, values.T, kx=degree, ky=degree)

    @classmethod
    def from_arrays(cls, x, y, values, degree: int = 3):
        """Build from arrays, flipping a top-down (y-decreasing) raster.

        Convenience for raster data whose rows run north-to-south (``y``
        decreasing), as GDAL/GeoTIFF typically store them.
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        values = np.asarray(values, dtype=float)
        if y[0] > y[-1]:
            y = y[::-1]
            values = values[::-1, :]
        if x[0] > x[-1]:
            x = x[::-1]
            values = values[:, ::-1]
        return cls(x, y, values, degree=degree)

    # -- sampling -------------------------------------------------------- #

    def value(self, x, y):
        """Field value(s) at point(s) ``(x, y)``.

        Scalar inputs give a 0-d array (``float()``-able); array inputs give an
        array of matching shape.
        """
        return np.squeeze(self._spline.ev(x, y))[()]

    __call__ = value

    def gradient(self, x, y):
        """Cartesian gradient ``(df/dx, df/dy)`` at point(s) ``(x, y)``."""
        df_dx = np.squeeze(self._spline.ev(x, y, dx=1, dy=0))[()]
        df_dy = np.squeeze(self._spline.ev(x, y, dx=0, dy=1))[()]
        return df_dx, df_dy

    def directional_derivative(self, x, y, ux, uy):
        """Derivative along unit vector ``(ux, uy)``: ``grad f . (ux, uy)``.

        No check that ``(ux, uy)`` is normalized — callers pass unit vectors.
        """
        df_dx, df_dy = self.gradient(x, y)
        return df_dx * ux + df_dy * uy

    # -- bounds ---------------------------------------------------------- #

    @property
    def bounds(self):
        """``(xmin, ymin, xmax, ymax)`` of the grid extent."""
        return (self._xmin, self._ymin, self._xmax, self._ymax)

    def contains(self, x, y):
        """Boolean: whether point(s) lie within the grid extent."""
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        return (
            (x >= self._xmin) & (x <= self._xmax)
            & (y >= self._ymin) & (y <= self._ymax)
        )
