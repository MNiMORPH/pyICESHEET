"""Shared geometry helpers: boundary resampling, inward normals, cleaning.

These back both :mod:`pyicesheet.margin` (the initial margin) and
:mod:`pyicesheet.contour` (the advancing contours). Keeping them in one place
avoids duplicating the resample/normal logic.
"""

from __future__ import annotations

import numpy as np
from shapely.geometry import Polygon, MultiPolygon, GeometryCollection, Point
from shapely import make_valid

__all__ = [
    "resample_ring",
    "inward_normals",
    "polygon_components",
    "clean_polygon",
]


def resample_ring(polygon: Polygon, spacing: float):
    """Even arc-length resample of a polygon's exterior ring.

    Returns ``(x, y)`` arrays of boundary points spaced ~``spacing`` apart, with
    the duplicate closing vertex dropped.
    """
    ring = polygon.exterior
    length = ring.length
    n = max(3, int(round(length / spacing)))
    dists = np.linspace(0.0, length, n, endpoint=False)
    pts = [ring.interpolate(d) for d in dists]
    x = np.array([p.x for p in pts])
    y = np.array([p.y for p in pts])
    return x, y


def inward_normals(x, y, polygon: Polygon, spacing: float):
    """Unit inward normals at each boundary point of a closed loop.

    The normal is perpendicular to the local tangent (from neighbouring points,
    wrapping around the loop), with sign chosen so a small step lands inside
    ``polygon``.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = x.size
    inx = np.empty(n)
    iny = np.empty(n)
    eps = max(1e-6, 0.01 * spacing)
    for i in range(n):
        tx = x[(i + 1) % n] - x[(i - 1) % n]
        ty = y[(i + 1) % n] - y[(i - 1) % n]
        norm = np.hypot(tx, ty)
        if norm == 0.0:
            tx, ty, norm = 1.0, 0.0, 1.0
        tx, ty = tx / norm, ty / norm
        nx, ny = -ty, tx                      # left normal of the tangent
        if not polygon.contains(Point(x[i] + eps * nx, y[i] + eps * ny)):
            nx, ny = -nx, -ny
        inx[i] = nx
        iny[i] = ny
    return inx, iny


def polygon_components(geom, min_area: float = 0.0):
    """Flatten a geometry into its constituent polygons above ``min_area``.

    ``make_valid`` on a self-intersecting ring can return a Polygon,
    MultiPolygon, or GeometryCollection (with stray lines/points); this pulls out
    just the polygonal pieces of significant area.
    """
    polys = []
    if isinstance(geom, Polygon):
        candidates = [geom]
    elif isinstance(geom, (MultiPolygon, GeometryCollection)):
        candidates = list(geom.geoms)
    else:
        candidates = []
    for g in candidates:
        if isinstance(g, (MultiPolygon, GeometryCollection)):
            polys.extend(polygon_components(g, min_area))
        elif isinstance(g, Polygon) and g.area >= min_area:
            polys.append(g)
    return polys


def clean_polygon(points_xy, min_area: float = 0.0):
    """Build a valid (multi)polygon from a possibly self-intersecting ring.

    This is the geometric heart of the crowding/divide handling: a ring of
    inward-advanced flowline endpoints may self-intersect (converging flowlines)
    or pinch (domes/saddles). ``make_valid`` resolves the self-intersections and
    splits pinched lobes into separate polygons automatically — replacing the
    original ICESHEET's hand-rolled motorcycle-graph crossover elimination and
    polygon-splitting bookkeeping.

    Returns a list of :class:`shapely.geometry.Polygon`.
    """
    if len(points_xy) < 3:
        return []
    poly = Polygon(points_xy)
    valid = make_valid(poly)
    return polygon_components(valid, min_area=min_area)
