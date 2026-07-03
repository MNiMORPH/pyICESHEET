"""Plan B: the Eulerian grid (eikonal) solver.

Two correctness checks: the analytic Nye cap (a point divide), and a line-divide
case (an ellipse) where the Lagrangian flowline solver over-encloses but the grid
solver — being single-valued — reproduces the distance-transform ground truth.
"""

import numpy as np
import pytest
from shapely.geometry import Point
from shapely.affinity import scale

from pyicesheet.grid import solve_surface_grid
from pyicesheet.constants import DEFAULT_CONSTANTS as C


def _flat_fields(x, y, tau=1.0e5):
    XX, _ = np.meshgrid(x, y)
    return np.zeros_like(XX), np.full_like(XX, tau)


def test_grid_reproduces_nye_cap():
    tau, R, h = 1.0e5, 2.0e5, 5000.0
    x = np.arange(-2.6e5, 2.6e5, h)
    y = np.arange(-2.6e5, 2.6e5, h)
    bed, tauA = _flat_fields(x, y, tau)
    margin = Point(0, 0).buffer(R, quad_segs=256)
    S = solve_surface_grid(x, y, bed, tauA, margin).surface

    XX, YY = np.meshgrid(x, y)
    d = R - np.hypot(XX, YY)                      # distance from margin
    analytic = np.sqrt(np.maximum(2 * tau * d / C.rho_g, 0.0))
    m = np.isfinite(S) & (d > 1e4) & (d < R - 1e4)
    resid = S[m] - analytic[m]

    assert np.sqrt((resid ** 2).mean()) < 40.0          # ~0.7% of a 2 km cap
    assert np.nanmax(S) == pytest.approx(np.sqrt(2 * tau * R / C.rho_g), abs=100.0)


def test_grid_handles_line_divide_without_over_enclosure():
    """Ellipse -> a line divide. The grid solver matches the distance-transform
    ground truth (no over-enclosure)."""
    from rasterio.features import rasterize
    from rasterio.transform import from_origin
    from scipy.ndimage import distance_transform_edt

    tau, h = 1.0e5, 5000.0
    x = np.arange(-5e5, 5e5, h)
    y = np.arange(-3e5, 3e5, h)
    bed, tauA = _flat_fields(x, y, tau)
    margin = scale(Point(0, 0).buffer(1.0, quad_segs=256), 3.5e5, 1.5e5)
    S = solve_surface_grid(x, y, bed, tauA, margin).surface

    tr = from_origin(x.min() - h / 2, y.max() + h / 2, h, h)
    inside = rasterize([(margin, 1)], out_shape=bed.shape, transform=tr,
                       fill=0, dtype="uint8").astype(bool)[::-1, :]
    d = distance_transform_edt(inside) * h
    gt = np.where(inside, np.sqrt(2 * tau * d / C.rho_g), np.nan)

    m = np.isfinite(S) & np.isfinite(gt) & (gt > 200)
    resid = S[m] - gt[m]
    assert np.sqrt((resid ** 2).mean()) < 40.0
    assert np.nanmax(S) == pytest.approx(np.nanmax(gt), abs=60.0)

    # Hypsometry: the grid must NOT enclose much more area than ground truth.
    cell = h * h
    for th in (1000.0, 1500.0):
        a_grid = np.nansum(S > th) * cell
        a_gt = np.nansum(gt > th) * cell
        assert a_grid <= 1.1 * a_gt                     # no over-enclosure


def test_grid_surface_above_bed_and_masked_outside():
    x = np.arange(-2.6e5, 2.6e5, 5000.0)
    y = np.arange(-2.6e5, 2.6e5, 5000.0)
    bed, tauA = _flat_fields(x, y)
    res = solve_surface_grid(x, y, bed, tauA, Point(0, 0).buffer(2.0e5, quad_segs=128))
    ice = np.isfinite(res.surface)
    assert np.all(res.surface[ice] >= res.bed[ice])     # surface at/above bed
    assert np.any(~ice)                                 # cells outside are masked
    assert np.nanmax(res.thickness) > 1500.0
