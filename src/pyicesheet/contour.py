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

from ._geometry import (
    resample_ring, inward_normals, clean_polygon, polygon_components,
)
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

    Physically-based thinning and stopping (near the summit)
    -------------------------------------------------------
    As flowlines converge toward an ice divide, the surface flattens: gaining one
    contour interval requires a horizontal distance of the Nye length
    ``L = interval * H / Hf`` (thickness ``H``, stress length ``Hf = tau/(rho g)``).
    This single length scale drives two resolution-independent behaviours:

    * **Stopping.** A contour whose equivalent radius ``sqrt(area/pi)`` is smaller
      than ``climb_factor * L`` cannot climb another interval — it is the summit —
      so it is not advanced further. This replaces the original resolution-scaled
      ``4*spacing**2`` area cutoff, which made the summit height depend on the
      input spacing.
    * **Thinning.** Where ``L`` is large (flat interior) the surface is smooth over
      that scale, so the contour is resampled at ``spacing_growth * L`` (clamped
      between the base ``spacing`` and ``spacing_cap_factor * spacing``) — far
      fewer flowlines are integrated in the flat interior, full resolution is kept
      at the steep margins.

    Parameters
    ----------
    integrator : FlowlineIntegrator
    spacing : float
        Base (finest) along-contour point spacing (m), used near steep margins.
    elevation_interval : float
        Elevation increment between contours (m).
    min_area : float, optional
        Numerical floor: discard degenerate polygons below this area (m^2);
        default ``spacing**2``. The real stopping is the physical Nye criterion.
    constants : PhysicalConstants, optional
    require_inside : bool, optional
        Reject advanced points that fall outside the parent contour polygon.
    survivor_rule : {"geos", "distance"}, optional
        How converging-flowline crossings are resolved (GEOS make_valid, or the
        ICESHEET motorcycle-graph distance-rule).
    climb_factor : float, optional
        Stop advancing a contour when ``sqrt(area/pi) < climb_factor * L``.
    spacing_growth : float, optional
        Interior spacing target as a fraction of the Nye length ``L``.
    spacing_cap_factor : float, optional
        Maximum interior spacing as a multiple of the base ``spacing``.
    clip_convergence : bool, optional
        Clip the advanced front so it cannot extend past where converging
        flowlines meet (the divide / medial axis). ``make_valid`` alone takes the
        *outer envelope* of a front whose flowlines have crossed, over-enclosing
        area beyond the divide (a spurious summit plateau). When True (default)
        each advanced polygon is intersected with the parent contour eroded
        inward by the local advance distance — geometrically "advance the whole
        parent inward by one step" — which collapses to nothing exactly where the
        opposing fronts have crossed, removing the over-enclosure. See
        ``docs/design-note-04``.
    """

    def __init__(self, integrator: FlowlineIntegrator, spacing: float,
                 elevation_interval: float, min_area: float | None = None,
                 constants: PhysicalConstants = DEFAULT_CONSTANTS,
                 require_inside: bool = True, survivor_rule: str = "geos",
                 climb_factor: float = 1.0, spacing_growth: float = 0.5,
                 spacing_cap_factor: float = 8.0, clip_convergence: bool = True):
        self.integrator = integrator
        self.spacing = float(spacing)
        self.elevation_interval = float(elevation_interval)
        self.min_area = (spacing ** 2) if min_area is None else float(min_area)
        self.constants = constants
        self.require_inside = require_inside
        if survivor_rule not in ("geos", "distance"):
            raise ValueError(
                f"survivor_rule={survivor_rule!r}: expected 'geos' or 'distance'"
            )
        self.survivor_rule = survivor_rule
        self.climb_factor = float(climb_factor)
        self.spacing_growth = float(spacing_growth)
        self.spacing_cap_factor = float(spacing_cap_factor)
        self.clip_convergence = bool(clip_convergence)

    # -- physically-based length scale, stopping, thinning --------------- #

    def _nye_length(self, poly, level):
        """Local Nye length ``L = interval * H / Hf`` and thickness at ``poly``.

        Sampled at an interior point of the contour. Returns ``(H, L)``.
        """
        p = poly.representative_point()
        B = float(self.integrator.bed.value(p.x, p.y))
        H = float(level) - B
        tau = float(self.integrator.tau.value(p.x, p.y))
        hf = tau / self.constants.rho_g
        L = (self.elevation_interval * H / hf) if hf > 0 else float("inf")
        return H, L

    def can_advance(self, poly, level):
        """Whether a contour at ``level`` can physically climb another interval."""
        H, L = self._nye_length(poly, level)
        if H <= self.integrator.min_thickness:
            return False
        radius = (poly.area / np.pi) ** 0.5
        return radius >= self.climb_factor * L

    def _effective_spacing(self, poly, level):
        """Interior-adaptive point spacing: coarsen toward the Nye scale."""
        _H, L = self._nye_length(poly, level)
        eff = self.spacing_growth * L
        eff = max(self.spacing, eff)
        return min(eff, self.spacing_cap_factor * self.spacing)

    def advance(self, contour: Contour, target: float):
        """Advance ``contour`` to elevation ``target``; return new contour(s).

        Points already at/above ``target`` are kept in place; the rest are
        integrated up their flowlines. Converging flowlines are then resolved by
        the chosen ``survivor_rule`` (GEOS ``make_valid``, or the ICESHEET
        motorcycle-graph distance-rule), and pinched lobes are split; each
        resulting polygon is resampled into a new :class:`Contour`.
        """
        ox, oy, nx, ny, pE, advanced = self._advance_points(contour, target)
        if len(nx) < 3:
            return []

        if self.survivor_rule == "distance":
            keep = _prune_crossing_flowlines(ox, oy, nx, ny, advanced)
            nx, ny, pE = nx[keep], ny[keep], pE[keep]
            if len(nx) < 3:
                return []

        # GEOS make_valid still runs: it splits pinched lobes (domes/saddles) and
        # cleans any residual self-intersection among the survivors.
        polys = clean_polygon(list(zip(nx, ny)), min_area=self.min_area)
        if self.clip_convergence:
            polys = self._clip_to_advance(polys, contour.polygon,
                                          ox, oy, nx, ny, advanced)
        children = []
        for poly in polys:
            if poly.exterior is None or poly.area < self.min_area:
                continue
            eff = self._effective_spacing(poly, target)
            children.append(Contour.build(poly, nx, ny, pE, target, eff))
        return children

    def _clip_to_advance(self, polys, parent, ox, oy, nx, ny, advanced):
        """Clip advanced polygons to the parent eroded by the advance distance.

        ``make_valid`` on the advanced ring keeps the *outer envelope* where
        converging flowlines have crossed, over-enclosing area past the divide
        (a spurious plateau near ice divides). Eroding the parent contour inward
        by the local advance distance ``d`` (``parent.buffer(-d)``) is exactly
        "advance the whole parent inward by one step": it reproduces the true
        front where fronts have *not* crossed, and collapses to nothing exactly
        where opposing fronts have crossed — so intersecting the ``make_valid``
        result with it removes only the over-enclosed fold-over.

        ``d`` is the median of the per-point flowline advance distances this
        step (robust to the odd long flowline; on a flat bed with constant tau
        the advance is uniform along a contour, so the choice is immaterial
        there and matters only where the fields vary).
        """
        if advanced.any():
            d = float(np.median(np.hypot(nx - ox, ny - oy)[advanced]))
        else:
            d = 0.0
        if d <= 0.0:
            return polys
        eroded = parent.buffer(-d)
        if eroded.is_empty:
            return []
        clipped = []
        for poly in polys:
            clipped.extend(
                polygon_components(poly.intersection(eroded), self.min_area)
            )
        return clipped

    # -- internals ------------------------------------------------------- #

    def _advance_points(self, contour: Contour, target: float):
        """Advance each contour point to ``target`` (or keep it if already up).

        Returns parallel arrays ``(old_x, old_y, new_x, new_y, E, advanced)`` for
        the surviving points: ``advanced`` flags points that were integrated up a
        flowline (vs. kept in place because already at/above ``target``). The
        old→new pairing is the flowline segment used by the distance-rule.
        """
        prepared = prep(contour.polygon) if self.require_inside else None
        ox, oy, nx, ny, pE, adv = [], [], [], [], [], []
        direction = contour.direction
        for i in range(len(contour)):
            xi, yi, Ei = float(contour.x[i]), float(contour.y[i]), float(contour.E[i])
            if Ei >= target - 1e-9:
                ox.append(xi); oy.append(yi); nx.append(xi); ny.append(yi)
                pE.append(Ei); adv.append(False)
                continue
            step = self.integrator.step_to_elevation(
                xi, yi, float(direction[i]), Ei, target, q0=0.0,
            )
            if not step.reached_target:
                continue  # flowline stalled (too thin / singularity): drop point
            if prepared is not None and not prepared.contains(Point(step.x, step.y)):
                continue  # escaped the parent contour: drop
            ox.append(xi); oy.append(yi); nx.append(step.x); ny.append(step.y)
            pE.append(target); adv.append(True)
        return (np.array(ox), np.array(oy), np.array(nx), np.array(ny),
                np.array(pE), np.array(adv, dtype=bool))


def _prune_crossing_flowlines(ox, oy, nx, ny, advanced):
    """Motorcycle-graph pruning: eliminate the shorter of each crossing pair.

    Reproduces the original ICESHEET rule: where two flowline segments (old
    contour point -> advanced point) cross, the line with the *shorter* distance
    from its start to the crossover is eliminated, processing crossings in order
    of increasing distance so a line pruned by one crossing auto-resolves its
    others. Only advanced points can be pruned; kept points always survive.

    Returns a boolean survivor mask aligned to the input arrays.
    """
    from shapely.geometry import LineString
    from shapely.strtree import STRtree

    n = len(nx)
    keep = np.ones(n, dtype=bool)
    idx_adv = np.where(advanced)[0]
    if len(idx_adv) < 2:
        return keep

    segs = [LineString([(ox[i], oy[i]), (nx[i], ny[i])]) for i in idx_adv]
    tree = STRtree(segs)

    crossings = []          # (distance_to_crossover, global_line_index, crossing_id)
    other_of = {}           # crossing_id -> (line_a_global, line_b_global)
    cid = 0
    for a_local, i in enumerate(idx_adv):
        for b_local in tree.query(segs[a_local]):
            if b_local <= a_local:
                continue
            if not segs[a_local].crosses(segs[b_local]):
                continue
            pt = segs[a_local].intersection(segs[b_local])
            if pt.geom_type != "Point":
                continue
            j = idx_adv[b_local]
            di = float(np.hypot(pt.x - ox[i], pt.y - oy[i]))
            dj = float(np.hypot(pt.x - ox[j], pt.y - oy[j]))
            other_of[cid] = (i, j)
            crossings.append((di, i, cid))
            crossings.append((dj, j, cid))
            cid += 1

    if not crossings:
        return keep

    crossings.sort(key=lambda e: e[0])
    eliminated = set()
    resolved = set()
    for _d, line, c in crossings:
        if c in resolved:
            continue
        if line in eliminated:
            resolved.add(c)
            continue
        a, b = other_of[c]
        other = b if line == a else a
        if other in eliminated:
            resolved.add(c)          # the other line already gone; keep this one
            continue
        eliminated.add(line)         # shorter line (smallest distance) is pruned
        resolved.add(c)

    for i in eliminated:
        keep[i] = False
    return keep
