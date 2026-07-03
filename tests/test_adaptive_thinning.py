"""Physically-based summit stopping and interior thinning (design-note-04).

The Nye length ``L = interval * H / Hf`` gives a resolution-independent stopping
criterion (a contour smaller than ``L`` cannot climb another interval) and an
interior-adaptive point spacing (coarsen toward ``L`` where the ice is flat).
The key property is that the reconstructed summit no longer depends on the input
spacing.
"""

import numpy as np
import pytest
from shapely.geometry import Point

from pyicesheet import RasterField, IceSheetModel, ModelConfig
from pyicesheet.contour import ContourManager
from pyicesheet.flowline import FlowlineIntegrator


def _cap(spacing):
    half, n = 4.0e5, 41
    x = np.linspace(-half, half, n)
    y = np.linspace(-half, half, n)
    XX, _ = np.meshgrid(x, y)
    bed = RasterField(x, y, np.zeros_like(XX))
    tau = RasterField(x, y, np.full_like(XX, 1.0e5))
    margin = Point(0.0, 0.0).buffer(2.0e5, quad_segs=64)
    cfg = ModelConfig(spacing=spacing, elevation_interval=100.0, max_elevation=5000.0)
    return IceSheetModel(bed, tau, margin, cfg).solve()


def test_summit_is_resolution_independent():
    """The cap apex is the same at coarse and fine spacing (the whole point)."""
    apex_20 = _cap(20000.0).elevation.max()
    apex_10 = _cap(10000.0).elevation.max()
    assert abs(apex_20 - apex_10) <= 50.0


def _manager(spacing=5000.0, interval=100.0):
    half, n = 4.0e5, 41
    x = np.linspace(-half, half, n)
    y = np.linspace(-half, half, n)
    XX, _ = np.meshgrid(x, y)
    bed = RasterField(x, y, np.zeros_like(XX))
    tau = RasterField(x, y, np.full_like(XX, 1.0e5))
    return ContourManager(FlowlineIntegrator(bed, tau), spacing, interval)


def test_can_advance_stops_small_contour():
    mgr = _manager()
    level = 2000.0                      # H=2000, Hf~11.1 -> L = 100*2000/11.1 ~ 18 km
    big = Point(0.0, 0.0).buffer(40_000.0)    # radius 40 km > L: can climb
    small = Point(0.0, 0.0).buffer(8_000.0)   # radius 8 km  < L: cannot
    assert mgr.can_advance(big, level)
    assert not mgr.can_advance(small, level)


def test_effective_spacing_coarsens_where_ice_is_thick():
    mgr = _manager(spacing=5000.0)
    thin_level = 200.0                  # H small -> L small -> spacing stays at base
    thick_level = 2500.0                # H large -> L large -> spacing coarsens
    poly = Point(0.0, 0.0).buffer(100_000.0)
    eff_thin = mgr._effective_spacing(poly, thin_level)
    eff_thick = mgr._effective_spacing(poly, thick_level)
    assert eff_thin == pytest.approx(5000.0)      # clamped to base
    assert eff_thick > 5000.0                     # coarsened toward the Nye length
