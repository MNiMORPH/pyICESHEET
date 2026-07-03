"""The configurable divide-handling rule (GEOS vs motorcycle-graph distance-rule).

The distance-rule reproduces the original ICESHEET crossover pruning; GEOS
make_valid is the default. Both must produce valid reconstructions; on a clean
radial cap the distance-rule prunes slightly more, giving a marginally lower apex.
See docs/design-note-03.
"""

import numpy as np
import pytest
from shapely.geometry import Point

from pyicesheet import RasterField, IceSheetModel, ModelConfig
from pyicesheet.contour import ContourManager
from pyicesheet.flowline import FlowlineIntegrator


def _radial(rule):
    half, n = 4.0e5, 41
    x = np.linspace(-half, half, n)
    y = np.linspace(-half, half, n)
    XX, _ = np.meshgrid(x, y)
    bed = RasterField(x, y, np.zeros_like(XX))
    tau = RasterField(x, y, np.full_like(XX, 1.0e5))
    margin = Point(0.0, 0.0).buffer(2.0e5, quad_segs=48)
    cfg = ModelConfig(spacing=15000.0, elevation_interval=200.0,
                      max_elevation=5000.0, survivor_rule=rule)
    return IceSheetModel(bed, tau, margin, cfg).solve()


def test_distance_rule_runs_and_is_valid():
    surf = _radial("distance")
    assert len(surf) > 100
    assert 1800.0 < surf.elevation.max() < 2300.0   # sane cap apex


def test_distance_rule_prunes_at_least_as_much_as_geos():
    """On a clean radial cap the distance-rule apex is <= the GEOS apex."""
    apex_geos = _radial("geos").elevation.max()
    apex_dist = _radial("distance").elevation.max()
    assert apex_dist <= apex_geos + 1e-6


def test_invalid_survivor_rule_raises():
    integ = FlowlineIntegrator(
        RasterField(np.linspace(0, 40, 5), np.linspace(0, 40, 5), np.zeros((5, 5))),
        RasterField(np.linspace(0, 40, 5), np.linspace(0, 40, 5), np.ones((5, 5))),
    )
    with pytest.raises(ValueError):
        ContourManager(integ, spacing=5000.0, elevation_interval=100.0,
                       survivor_rule="bogus")
