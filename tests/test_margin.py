"""IceMargin: resampling and inward normals."""

import numpy as np
from shapely.geometry import Point

from pyicesheet.margin import IceMargin


def _circle(radius=100_000.0, n=400, cx=0.0, cy=0.0):
    t = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    return Point(cx, cy).buffer(radius, quad_segs=n // 4)


def test_resampling_spacing():
    poly = _circle(radius=100_000.0)
    spacing = 5_000.0
    m = IceMargin.from_polygon(poly, spacing)
    # Point count roughly perimeter / spacing.
    perim = poly.exterior.length
    assert abs(len(m) - round(perim / spacing)) <= 1
    # Consecutive spacings are close to target.
    dx = np.diff(np.r_[m.x, m.x[0]])
    dy = np.diff(np.r_[m.y, m.y[0]])
    seg = np.hypot(dx, dy)
    assert np.allclose(seg, spacing, rtol=0.05)


def test_inward_normals_point_to_centre():
    cx, cy, R = 12_000.0, -3_000.0, 80_000.0
    poly = _circle(radius=R, cx=cx, cy=cy)
    m = IceMargin.from_polygon(poly, 4_000.0)
    # For a circle, the inward normal at each point should point toward centre.
    to_centre_x = cx - m.x
    to_centre_y = cy - m.y
    norm = np.hypot(to_centre_x, to_centre_y)
    dot = (m.inward_x * to_centre_x + m.inward_y * to_centre_y) / norm
    # Cosine of angle between inward normal and radial-inward direction ~ 1.
    assert np.all(dot > 0.99)


def test_inward_normals_are_unit():
    poly = _circle(radius=50_000.0)
    m = IceMargin.from_polygon(poly, 3_000.0)
    assert np.allclose(np.hypot(m.inward_x, m.inward_y), 1.0, atol=1e-9)
