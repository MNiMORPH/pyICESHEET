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
    rtol, atol : float
        Flowline integrator tolerances.
    min_area : float or None
        Discard advanced polygons below this area (m^2); None -> 4*spacing^2.
    constants : PhysicalConstants
    """

    spacing: float = 5000.0
    elevation_interval: float = 100.0
    min_thickness: float = 1.0
    min_elevation: float = 0.0
    max_elevation: float = 5000.0
    require_inside: bool = True
    rtol: float = 1e-8
    atol: float = 1e-6
    min_area: float | None = None
    constants: PhysicalConstants = field(default=DEFAULT_CONSTANTS)
