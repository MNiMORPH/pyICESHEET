"""M3: the full 2-D solver on a radial dome, checked against the Nye profile.

A circular margin over a flat bed with uniform shear stress must reconstruct a
dome whose radial surface profile follows the analytic Nye parabola. This
exercises the whole chain (seed -> flowline advance -> GEOS cleaning -> resample
-> recurse) and ties it back to the M1 analytic anchor.
"""

import numpy as np
import pytest
from shapely.geometry import Point

from pyicesheet import RasterField, IceSheetModel, ModelConfig, physics
from pyicesheet.constants import DEFAULT_CONSTANTS as C


def _flat_uniform_model(radius=2.0e5, tau=1.0e5, spacing=10_000.0, interval=200.0):
    half = 2.0 * radius
    n = 41
    x = np.linspace(-half, half, n)
    y = np.linspace(-half, half, n)
    XX, _ = np.meshgrid(x, y)
    bed = RasterField(x, y, np.zeros_like(XX))
    tau_f = RasterField(x, y, np.full_like(XX, tau))
    margin = Point(0.0, 0.0).buffer(radius, quad_segs=64)
    cfg = ModelConfig(spacing=spacing, elevation_interval=interval, max_elevation=5000.0)
    return IceSheetModel(bed, tau_f, margin, cfg)


def test_radial_dome_matches_nye_profile():
    radius, tau = 2.0e5, 1.0e5
    surf = _flat_uniform_model(radius=radius, tau=tau).solve()

    r = np.hypot(surf.x, surf.y)
    dist_in = radius - r
    interior = dist_in > 2.0e4          # avoid margin toe and the singular apex
    Enye = physics.nye_thickness(dist_in[interior], tau, C, thickness0=1.0)
    resid = surf.elevation[interior] - Enye

    # Unbiased, and small relative to a ~2 km dome.
    assert abs(resid.mean()) < 5.0
    assert resid.std() < 40.0
    assert np.abs(resid).max() < 200.0


def test_radial_dome_apex_near_analytic():
    radius, tau = 2.0e5, 1.0e5
    surf = _flat_uniform_model(radius=radius, tau=tau).solve()
    apex_analytic = float(physics.nye_thickness(radius, tau, C, thickness0=1.0))
    # Highest reconstructed sample within ~one interval of the analytic apex.
    assert surf.elevation.max() == pytest.approx(apex_analytic, abs=250.0)


def test_surface_thickness_and_bounds():
    surf = _flat_uniform_model().solve()
    # Flat bed at zero -> thickness == elevation.
    np.testing.assert_allclose(surf.thickness, surf.elevation, atol=1e-6)
    xmin, ymin, xmax, ymax = surf.bounds
    assert xmin >= -2.05e5 and xmax <= 2.05e5


def test_to_grid_produces_raster():
    surf = _flat_uniform_model().solve()
    grid, xc, yc = surf.to_grid(resolution=10_000.0, value="thickness")
    assert grid.ndim == 2
    assert np.nanmax(grid) > 1500.0        # dome interior present
    assert yc[0] > yc[-1]                   # top-down
