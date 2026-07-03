"""The reconstruction driver: seed the margin, march contours inward.

:class:`IceSheetModel` couples the fields, the flowline integrator, and the
contour manager. :meth:`IceSheetModel.solve` seeds a starting contour on the ice
margin, then marches contours inward one elevation interval at a time — recursing
into each lobe produced when a contour splits — accumulating the reconstructed
surface samples into an :class:`~pyicesheet.surface.IceSurface`.
"""

from __future__ import annotations

import numpy as np
from shapely.geometry import Polygon, MultiPolygon

from . import physics
from .config import ModelConfig
from .contour import Contour, ContourManager
from .fields import RasterField
from .flowline import FlowlineIntegrator
from ._geometry import inward_normals, resample_ring
from .surface import IceSurface

__all__ = ["IceSheetModel"]


class IceSheetModel:
    """Equilibrium perfectly-plastic ice-sheet reconstruction.

    Parameters
    ----------
    bed : RasterField
        Bed topography (m).
    tau : RasterField
        Basal shear stress (Pa).
    margin : shapely Polygon or MultiPolygon
        The mapped ice margin (map units, same CRS as the fields).
    config : ModelConfig, optional
    """

    def __init__(self, bed: RasterField, tau: RasterField, margin,
                 config: ModelConfig | None = None):
        self.bed = bed
        self.tau = tau
        self.margin = margin
        self.config = config or ModelConfig()
        self.constants = self.config.constants

        self.integrator = FlowlineIntegrator(
            bed, tau, self.constants,
            min_thickness=self.config.min_thickness,
            rtol=self.config.rtol, atol=self.config.atol,
        )
        self.manager = ContourManager(
            self.integrator, self.config.spacing, self.config.elevation_interval,
            min_area=self.config.min_area, constants=self.constants,
            require_inside=self.config.require_inside,
        )

    # -- public ---------------------------------------------------------- #

    def solve(self) -> IceSurface:
        """Run the reconstruction and return the ice surface."""
        polygons = self._margin_polygons()
        xs, ys, Es = [], [], []
        for poly in polygons:
            self._march_polygon(poly, xs, ys, Es)

        x = np.array(xs)
        y = np.array(ys)
        E = np.array(Es)
        bed_vals = np.asarray(self.bed.value(x, y), dtype=float)
        thickness = E - bed_vals
        return IceSurface(x=x, y=y, elevation=E, thickness=thickness, bed=bed_vals)

    # -- internals ------------------------------------------------------- #

    def _margin_polygons(self):
        if isinstance(self.margin, MultiPolygon):
            return list(self.margin.geoms)
        if isinstance(self.margin, Polygon):
            return [self.margin]
        raise TypeError("margin must be a shapely Polygon or MultiPolygon")

    def _seed_contour(self, polygon: Polygon) -> Contour:
        """Seed a starting contour on the margin of one polygon.

        Each margin point is placed at the bed elevation plus the nominal
        margin thickness; below sea level, a near-flotation surface is used.
        """
        x, y = resample_ring(polygon, self.config.spacing)
        inx, iny = inward_normals(x, y, polygon, self.config.spacing)
        Bv = np.asarray(self.bed.value(x, y), dtype=float)
        S = physics.marine_flotation_surface(Bv, self.constants)
        E = S + self.config.min_thickness
        level = float(np.min(E))
        return Contour(x=x, y=y, E=E, inward_x=inx, inward_y=iny,
                       polygon=polygon, level=level)

    def _march_polygon(self, polygon: Polygon, xs, ys, Es):
        """March one margin polygon inward, accumulating surface samples."""
        interval = self.config.elevation_interval
        max_elev = self.config.max_elevation

        seed = self._seed_contour(polygon)
        _accumulate(seed, xs, ys, Es)

        step0 = int(np.floor(np.min(seed.E) / interval)) + 1
        stack = [(seed, step0)]
        while stack:
            contour, step = stack.pop()
            target = step * interval
            if target > max_elev:
                continue
            children = self.manager.advance(contour, target)
            for child in children:
                _accumulate(child, xs, ys, Es)
                stack.append((child, step + 1))


def _accumulate(contour: Contour, xs, ys, Es):
    xs.extend(np.asarray(contour.x, dtype=float).tolist())
    ys.extend(np.asarray(contour.y, dtype=float).tolist())
    Es.extend(np.asarray(contour.E, dtype=float).tolist())
