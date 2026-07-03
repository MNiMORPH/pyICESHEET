"""tau-from-substrate builder and the calibration layer."""

import numpy as np
import pytest
from shapely.geometry import Point

from pyicesheet.grid import solve_surface_grid
from pyicesheet.tau import tau_from_classes
from pyicesheet.calibrate import SurfaceMisfit, calibrate, calibrate_scalar


def test_tau_from_classes_array_and_dict():
    cls = np.array([[1, 2, 0], [3, 2, 1]])
    tau = tau_from_classes(cls, [100e3, 80e3, 50e3])
    assert tau[0, 0] == 100e3 and tau[0, 1] == 80e3 and tau[1, 0] == 50e3
    assert np.isnan(tau[0, 2])                      # class 0 -> nodata
    tau2 = tau_from_classes(cls, {1: 100e3, 2: 80e3, 3: 50e3})
    np.testing.assert_array_equal(np.nan_to_num(tau), np.nan_to_num(tau2))


def _cap(h=10000.0, R=1.2e5):
    x = np.arange(-1.6e5, 1.6e5, h)
    y = x.copy()
    XX, YY = np.meshgrid(x, y)
    bed = np.zeros_like(XX)
    return x, y, bed, XX, YY, Point(0, 0).buffer(R, quad_segs=128)


def test_calibrate_scalar_recovers_a_known_multiplier():
    x, y, bed, XX, _, margin = _cap()
    base = np.full_like(XX, 1.0e5)
    forward = lambda a: solve_surface_grid(x, y, bed, a * base, margin).surface
    truth = forward(0.8)                            # synthetic "observation"
    obj = SurfaceMisfit(truth, mask=np.isfinite(truth))
    res = calibrate_scalar(forward, obj, bounds=(0.4, 1.4))
    assert res.x == pytest.approx(0.8, abs=0.03)
    assert res.fun < 2.0                            # near-zero misfit at the optimum


def test_calibrate_per_class_recovers_a_two_class_field():
    x, y, bed, XX, YY, margin = _cap()
    cls = np.where(np.hypot(XX, YY) < 6.0e4, 1, 2)  # inner / outer class
    forward = lambda m: solve_surface_grid(
        x, y, bed, tau_from_classes(cls, [m[0] * 1e5, m[1] * 1e5]), margin).surface
    truth = forward(np.array([0.7, 1.1]))
    obj = SurfaceMisfit(truth, mask=np.isfinite(truth))
    res = calibrate(forward, [1.0, 1.0], obj, bounds=[(0.4, 1.5)] * 2,
                    method="L-BFGS-B", eps=0.02, maxiter=40)
    assert res.x[0] == pytest.approx(0.7, abs=0.08)
    assert res.x[1] == pytest.approx(1.1, abs=0.08)
