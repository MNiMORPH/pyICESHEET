"""Model configuration."""

from __future__ import annotations

from dataclasses import dataclass, field

from .constants import DEFAULT_CONSTANTS, PhysicalConstants

__all__ = ["ModelConfig"]


@dataclass
class ModelConfig:
    """Parameters controlling a pyICESHEET reconstruction.

    Attributes
    ----------
    spacing : float
        Target along-contour point spacing (m); the flowline-seed density. The
        original ICESHEET ``minimum_spacing``.
    elevation_interval : float
        Elevation increment between successive contours (m).
    min_thickness : float
        Nominal ice thickness assigned at the margin (m).
    min_elevation : float
        Bed elevation (m) below which the margin is treated as marine and seeded
        at a near-flotation surface. Default 0 (sea level).
    max_elevation : float
        Stop marching once contours exceed this surface elevation (m).
    require_inside : bool
        Reject advanced flowline points that escape the parent contour.
    survivor_rule : {"geos", "distance"}
        How converging flowlines are resolved near ice divides. "geos" (default)
        lets GEOS make_valid resolve crossings (smoother, closer to observed);
        "distance" applies the original ICESHEET motorcycle-graph pruning
        (reproduces Gowan's more aggressive interior pruning). See
        docs/design-note-03.
    rtol, atol : float
        Flowline integrator tolerances.
    min_area : float or None
        Numerical floor on polygon area (m^2); None -> spacing^2. The real
        summit stopping is the physical Nye criterion below.
    climb_factor : float
        Physical stopping: a contour stops advancing when its equivalent radius
        drops below climb_factor * (Nye length L = interval*H/Hf). See
        ContourManager and docs/design-note-04.
    spacing_growth : float
        Interior point spacing target as a fraction of the Nye length L.
    spacing_cap_factor : float
        Maximum interior spacing as a multiple of the base spacing.
    clip_convergence : bool
        Clip the advancing front so it cannot extend past where converging
        flowlines meet (the ice divide / medial axis). GEOS make_valid alone
        takes the outer envelope of a crossed front, over-enclosing area beyond
        the divide (a spurious summit plateau); this trims it. See
        ContourManager and docs/design-note-04.
    constants : PhysicalConstants
    """

    spacing: float = 5000.0
    elevation_interval: float = 100.0
    min_thickness: float = 1.0
    min_elevation: float = 0.0
    max_elevation: float = 5000.0
    require_inside: bool = True
    survivor_rule: str = "geos"
    rtol: float = 1e-8
    atol: float = 1e-6
    min_area: float | None = None
    climb_factor: float = 1.0
    spacing_growth: float = 0.5
    spacing_cap_factor: float = 8.0
    clip_convergence: bool = True
    constants: PhysicalConstants = field(default=DEFAULT_CONSTANTS)
