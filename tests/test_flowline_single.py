"""FlowlineIntegrator: single-flowline behaviour through real fields.

Ties M2 back to the M1 analytic anchor: a flowline through a uniform (flat-bed,
constant-tau) field must reproduce the Nye distance, travel straight, and keep
q = 0. A cross-flow bed gradient must make q and the lateral offset develop with
the correct sign.
"""

import numpy as np
import pytest

from pyicesheet import physics
from pyicesheet.constants import DEFAULT_CONSTANTS as C
from pyicesheet.fields import RasterField
from pyicesheet.flowline import FlowlineIntegrator


def _uniform_fields(tau_value=1.0e5, bed_slope_x=0.0, bed_slope_y=0.0,
                    n=41, half=6.0e5):
    x = np.linspace(-half, half, n)
    y = np.linspace(-half, half, n)
    XX, YY = np.meshgrid(x, y)
    bed = bed_slope_x * XX + bed_slope_y * YY   # flat if both slopes 0
    tau = np.full_like(XX, tau_value)
    return RasterField(x, y, bed), RasterField(x, y, tau)


def _lateral_offset(step, x0, y0, theta0):
    """Signed offset of the endpoint along the cross-flow axis n=(-sin,cos)."""
    nx, ny = -np.sin(theta0), np.cos(theta0)
    return (step.x - x0) * nx + (step.y - y0) * ny


def test_radial_flowline_reproduces_nye_distance():
    tau = 1.0e5
    bed, tau_f = _uniform_fields(tau_value=tau)
    integ = FlowlineIntegrator(bed, tau_f, C)

    x0, y0 = 4.0e5, 0.0
    theta0 = np.pi                       # inward normal points toward -x
    E0, target = 1.0, 1000.0
    step = integ.step_to_elevation(x0, y0, theta0, E0, target)

    assert step.reached_target and step.status == "reached"
    assert step.elevation == pytest.approx(target, abs=1e-3)

    # Flat bed: thickness == elevation, so arc length == Nye distance.
    expected = physics.nye_distance(E0, target, tau, C)
    assert step.arc_length == pytest.approx(expected, rel=1e-4)

    # Straight, radial: endpoint at x0 - arc_length, no drift, q ~ 0.
    assert step.x == pytest.approx(x0 - step.arc_length, abs=1.0)
    assert step.y == pytest.approx(0.0, abs=1e-3)
    assert abs(step.q) < 1e-9
    assert step.direction == pytest.approx(theta0, abs=1e-9)


def test_uniform_field_direction_unchanged_over_multiple_steps():
    """Marching a uniform field inward should keep a straight radial flowline."""
    tau = 1.0e5
    bed, tau_f = _uniform_fields(tau_value=tau)
    integ = FlowlineIntegrator(bed, tau_f, C)

    x, y, theta, E = 4.0e5, 0.0, np.pi, 1.0
    for target in (200.0, 500.0, 900.0, 1400.0):
        step = integ.step_to_elevation(x, y, theta, E, target)
        assert step.reached_target
        x, y, theta, E = step.x, step.y, step.direction, step.elevation
        assert abs(step.q) < 1e-8
        assert step.y == pytest.approx(0.0, abs=1e-2)

    # Total inward distance matches Nye from 1 -> 1400 m.
    assert (4.0e5 - x) == pytest.approx(physics.nye_distance(1.0, 1400.0, tau, C), rel=1e-4)


def test_cross_flow_bed_gradient_develops_q():
    """A cross-flow bed gradient develops q, lateral drift, and turns the flowline.

    q and the lateral offset must share a sign (dy/dx = q/p, p > 0). We do not
    assert an absolute turn direction here; antisymmetry and the flow-aligned
    control (below) pin down the gradient-coupling sign instead.
    """
    tau = 1.0e5
    x0, y0, theta0 = 4.0e5, 0.0, np.pi
    bed, tau_f = _uniform_fields(tau_value=tau, bed_slope_y=2.0e-3)
    step = FlowlineIntegrator(bed, tau_f, C).step_to_elevation(x0, y0, theta0, 5.0, 1200.0)

    assert step.reached_target
    assert abs(step.q) > 1e-4
    lat = _lateral_offset(step, x0, y0, theta0)
    assert abs(lat) > 1.0
    assert np.sign(lat) == np.sign(step.q)          # dy/dx = q/p, p>0
    assert abs(step.direction - np.pi) > 1e-4       # flowline turned


def test_flow_aligned_bed_gradient_keeps_q_zero():
    """A bed gradient ALONG the flow axis has zero cross-flow component -> q~0.

    This specifically checks that only the component of grad(B) along n (not t)
    drives q, i.e. the cross-flow projection is correct.
    """
    tau = 1.0e5
    # theta0 = pi -> flow axis t = (-1, 0); a bed sloping in x is flow-aligned.
    # Seed at x=0 where bed = slope_x * 0 = 0 (so thickness starts positive); the
    # flowline runs into -x, where the bed drops.
    bed, tau_f = _uniform_fields(tau_value=tau, bed_slope_x=2.0e-3)
    step = FlowlineIntegrator(bed, tau_f, C).step_to_elevation(0.0, 0.0, np.pi, 5.0, 1200.0)
    assert step.reached_target
    assert abs(step.q) < 1e-8
    assert abs(_lateral_offset(step, 0.0, 0.0, np.pi)) < 1e-2


def test_cross_flow_gradient_sign_is_antisymmetric():
    """Mirroring the cross-flow bed gradient flips the sign of q and the drift."""
    tau = 1.0e5
    x0, y0, theta0 = 4.0e5, 0.0, np.pi
    s_plus = FlowlineIntegrator(*_uniform_fields(tau_value=tau, bed_slope_y=+2.0e-3), C) \
        .step_to_elevation(x0, y0, theta0, 5.0, 1200.0)
    s_minus = FlowlineIntegrator(*_uniform_fields(tau_value=tau, bed_slope_y=-2.0e-3), C) \
        .step_to_elevation(x0, y0, theta0, 5.0, 1200.0)
    assert np.sign(s_plus.q) == -np.sign(s_minus.q)
    assert s_plus.q == pytest.approx(-s_minus.q, rel=1e-3)


def test_flat_bed_control_keeps_q_zero():
    bed, tau_f = _uniform_fields(tau_value=1.0e5)
    step = FlowlineIntegrator(bed, tau_f, C).step_to_elevation(4.0e5, 0.0, np.pi, 5.0, 1200.0)
    assert abs(step.q) < 1e-9
