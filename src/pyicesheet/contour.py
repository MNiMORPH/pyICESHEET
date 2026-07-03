"""Contours and the front-advancing step (the crowding/divide handling).

A :class:`Contour` is one closed loop of the reconstruction: ordered boundary
points, each with a surface elevation ``E``, plus inward normals giving the local
along-flow direction. :class:`ContourManager` advances a contour inward by one
elevation interval and returns the resulting contour(s).

Advancing is where converging flowlines, ice divides, and dome/saddle splitting
are handled. Following the design (``docs/design-note-02``), the topology is
delegated to GEOS: every point is integrated up its flowline to the target
elevation, the advanced points form a ring, and :func:`clean_polygon`
(``make_valid``) resolves any self-intersections and splits pinched lobes into
separate polygons automatically. Each output polygon is resampled to even
spacing and becomes a new contour to recurse on.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree
from shapely.geometry import Polygon, Point
from shapely.prepared import prep

from ._geometry import resample_ring, inward_normals, clean_polygon
from .constants import DEFAULT_CONSTANTS, PhysicalConstants
from .flowline import FlowlineIntegrator

__all__ = ["Contour", "ContourManager"]


@dataclass
class Contour:
    """A closed contour loop with per-point surface elevation and normals."""

    x: np.ndarray
    y: np.ndarray
    E: np.ndarray
    inward_x: np.ndarray
    inward_y: np.ndarray
    polygon: Polygon
    level: float          # nominal elevation level of this contour (m)

    @property
    def direction(self):
        """Along-flow azimuths (radians) = inward-normal directions."""
        return np.arctan2(self.inward_y, self.inward_x)

    def __len__(self):
        return self.x.size

    @classmethod
    def build(cls, polygon, src_x, src_y, src_E, level, spacing):
        """Resample ``polygon`` to ``spacing``; assign E from source points.

        ``src_*`` are the pre-cleaning advanced points carrying known elevations;
        each resampled boundary point takes the elevation of its nearest source
        point (points that GEOS inserted at intersections thereby take ~the
        target level, which is where those intersections physically sit).
        """
        x, y = resample_ring(polygon, spacing)
        inx, iny = inward_normals(x, y, polygon, spacing)
        tree = cKDTree(np.column_stack([src_x, src_y]))
        _, idx = tree.query(np.column_stack([x, y]))
        E = np.asarray(src_E, dtype=float)[idx]
        return cls(x, y, E, inx, iny, polygon, float(level))


class ContourManager:
    """Advance contours inward, delegating divide/dome topology to GEOS.

    Parameters
    ----------
    integrator : FlowlineIntegrator
    spacing : float
        Target along-contour point spacing (m).
    elevation_interval : float
        Elevation increment between contours (m).
    min_area : float, optional
        Discard advanced polygons smaller than this (m^2); default
        ``4 * spacing**2`` (a few points' worth).
    constants : PhysicalConstants, optional
    require_inside : bool, optional
        If True (default), reject advanced points that fall outside the parent
        contour polygon (a flowline should march inward, not escape).
    survivor_rule : {"geos"}, optional
        How converging-flowline crossings are resolved. Only ``"geos"`` (GEOS
        make_valid) is implemented; the original ICESHEET distance-rule
        ("shorter flowline wins") is a planned configurable alternative.
    """

    def __init__(self, integrator: FlowlineIntegrator, spacing: float,
                 elevation_interval: float, min_area: float | None = None,
                 constants: PhysicalConstants = DEFAULT_CONSTANTS,
                 require_inside: bool = True, survivor_rule: str = "geos"):
        self.integrator = integrator
        self.spacing = float(spacing)
        self.elevation_interval = float(elevation_interval)
        self.min_area = (4.0 * spacing ** 2) if min_area is None else float(min_area)
        self.constants = constants
        self.require_inside = require_inside
        if survivor_rule != "geos":
            raise NotImplementedError(
                f"survivor_rule={survivor_rule!r}: only 'geos' is implemented; "
                "the ICESHEET distance-rule is a planned alternative"
            )
        self.survivor_rule = survivor_rule

    def advance(self, contour: Contour, target: float):
        """Advance ``contour`` to elevation ``target``; return new contour(s).

        Points already at/above ``target`` are kept in place; the rest are
        integrated up their flowlines. The advanced ring is cleaned with GEOS,
        which resolves crossings and splits pinched lobes; each resulting polygon
        is resampled into a new :class:`Contour`.
        """
        px, py, pE = self._advance_points(contour, target)
        if len(px) < 3:
            return []
        polys = clean_polygon(list(zip(px, py)), min_area=self.min_area)
        children = []
        for poly in polys:
            if poly.exterior is None or poly.area < self.min_area:
                continue
            children.append(
                Contour.build(poly, px, py, pE, target, self.spacing)
            )
        return children

    # -- internals ------------------------------------------------------- #

    def _advance_points(self, contour: Contour, target: float):
        """Advance each contour point to ``target`` (or keep it if already up)."""
        prepared = prep(contour.polygon) if self.require_inside else None
        px, py, pE = [], [], []
        direction = contour.direction
        for i in range(len(contour)):
            Ei = float(contour.E[i])
            if Ei >= target - 1e-9:
                px.append(float(contour.x[i]))
                py.append(float(contour.y[i]))
                pE.append(Ei)
                continue
            step = self.integrator.step_to_elevation(
                float(contour.x[i]), float(contour.y[i]), float(direction[i]),
                Ei, target, q0=0.0,
            )
            if not step.reached_target:
                continue  # flowline stalled (too thin / singularity): drop point
            if prepared is not None and not prepared.contains(Point(step.x, step.y)):
                continue  # escaped the parent contour: drop
            px.append(step.x)
            py.append(step.y)
            pE.append(target)
        return px, py, pE
