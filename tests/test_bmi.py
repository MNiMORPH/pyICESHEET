"""BMI adapter: equilibrium single-shot lifecycle on the radial case."""

import numpy as np
import pytest
from shapely.geometry import Point

from pyicesheet import RasterField
from pyicesheet.bmi import IceSheetBMI


def _config():
    half, n = 4.0e5, 41
    x = np.linspace(-half, half, n)
    y = np.linspace(-half, half, n)
    XX, _ = np.meshgrid(x, y)
    bed = RasterField(x, y, np.zeros_like(XX))
    tau = RasterField(x, y, np.full_like(XX, 1.0e5))
    margin = Point(0.0, 0.0).buffer(2.0e5, quad_segs=48)
    return {
        "bed": bed, "tau": tau, "margin": margin,
        "spacing": 15000.0, "elevation_interval": 200.0,
        "max_elevation": 5000.0, "output_resolution": 20000.0,
    }


def test_bmi_lifecycle_and_time_is_degenerate():
    bmi = IceSheetBMI()
    bmi.initialize(_config())
    # Equilibrium: time functions are degenerate.
    assert bmi.get_start_time() == bmi.get_end_time() == bmi.get_current_time() == 0.0
    bmi.update()
    assert "elevation" in bmi.get_component_name() or "ice sheet" in bmi.get_component_name().lower()
    bmi.finalize()


def test_bmi_get_value_returns_thickness_grid():
    bmi = IceSheetBMI()
    bmi.initialize(_config())
    bmi.update()
    assert set(bmi.get_output_var_names()) == {
        "land_ice_surface__elevation", "land_ice__thickness"
    }
    thk = bmi.get_value("land_ice__thickness")
    assert thk.ndim == 1
    assert np.nanmax(thk) > 1500.0            # dome interior present
    assert bmi.get_var_units("land_ice__thickness") == "m"
    shape = bmi.get_grid_shape()
    assert len(shape) == 2 and shape[0] * shape[1] == thk.size


def test_bmi_read_before_update_raises():
    bmi = IceSheetBMI()
    bmi.initialize(_config())
    with pytest.raises(RuntimeError):
        bmi.get_value("land_ice__thickness")
