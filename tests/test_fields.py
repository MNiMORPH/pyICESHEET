"""RasterField: smooth sampling and analytic gradients."""

import numpy as np
import pytest

from pyicesheet.fields import RasterField


def _grid(f, n=25, lo=-100.0, hi=100.0):
    x = np.linspace(lo, hi, n)
    y = np.linspace(lo, hi, n)
    XX, YY = np.meshgrid(x, y)          # shape (ny, nx); XX[j,i]=x[i]
    return x, y, f(XX, YY)


def test_value_linear_field():
    x, y, v = _grid(lambda X, Y: 2.0 * X + 3.0 * Y + 5.0)
    fld = RasterField(x, y, v)
    for px, py in [(0.0, 0.0), (10.0, -20.0), (55.5, 33.3)]:
        assert float(fld.value(px, py)) == pytest.approx(2 * px + 3 * py + 5, abs=1e-6)


def test_gradient_linear_field():
    x, y, v = _grid(lambda X, Y: 2.0 * X + 3.0 * Y + 5.0)
    fld = RasterField(x, y, v)
    gx, gy = fld.gradient(12.0, -7.0)
    assert float(gx) == pytest.approx(2.0, abs=1e-6)
    assert float(gy) == pytest.approx(3.0, abs=1e-6)


def test_gradient_quadratic_field():
    # f = x^2 + y^2 -> grad = (2x, 2y)
    x, y, v = _grid(lambda X, Y: X ** 2 + Y ** 2)
    fld = RasterField(x, y, v)
    gx, gy = fld.gradient(30.0, -40.0)
    assert float(gx) == pytest.approx(60.0, rel=1e-3)
    assert float(gy) == pytest.approx(-80.0, rel=1e-3)


def test_directional_derivative():
    x, y, v = _grid(lambda X, Y: 2.0 * X + 3.0 * Y + 5.0)
    fld = RasterField(x, y, v)
    # along (1,0) -> 2; along (0,1) -> 3; along (0,-1) -> -3
    assert float(fld.directional_derivative(0.0, 0.0, 1.0, 0.0)) == pytest.approx(2.0, abs=1e-6)
    assert float(fld.directional_derivative(0.0, 0.0, 0.0, -1.0)) == pytest.approx(-3.0, abs=1e-6)


def test_from_arrays_flips_topdown_raster():
    x = np.linspace(-100, 100, 20)
    y_topdown = np.linspace(100, -100, 20)   # decreasing (GDAL style)
    XX, YY = np.meshgrid(x, y_topdown)
    v = 2.0 * XX + 3.0 * YY
    fld = RasterField.from_arrays(x, y_topdown, v)
    assert float(fld.value(10.0, 20.0)) == pytest.approx(2 * 10 + 3 * 20, abs=1e-6)


def test_bounds_and_contains():
    x, y, v = _grid(lambda X, Y: X + Y, lo=-100, hi=100)
    fld = RasterField(x, y, v)
    assert fld.bounds == (-100.0, -100.0, 100.0, 100.0)
    assert bool(fld.contains(0.0, 0.0))
    assert not bool(fld.contains(200.0, 0.0))
