"""M1 validation: perfectly-plastic physics against the analytic Nye solution.

The key test integrates the flowline ODE right-hand side
(:func:`pyicesheet.physics.slope_derivatives`) numerically for the constant-shear-
stress, flat-bed case and confirms it reproduces the closed-form Nye parabola
``H(x) = sqrt(H0**2 + 2 tau x / (rho_i g))``. This is the M1 "radial solver":
integrating the physics along a straight inward flowline. It exercises the actual
ODE form, not the closed form against itself.
"""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from pyicesheet import physics
from pyicesheet.constants import DEFAULT_CONSTANTS as C


# --------------------------------------------------------------------------- #
# Analytic relations                                                          #
# --------------------------------------------------------------------------- #

def test_nye_thickness_closed_form():
    """nye_thickness matches H = sqrt(H0^2 + 2 tau x / (rho g))."""
    tau = 1.0e5  # 100 kPa
    x = np.linspace(0.0, 2.0e5, 50)
    H0 = 1.0
    expected = np.sqrt(H0 ** 2 + 2.0 * tau * x / C.rho_g)
    got = physics.nye_thickness(x, tau, C, thickness0=H0)
    np.testing.assert_allclose(got, expected, rtol=0, atol=1e-9)


def test_nye_thickness_at_margin_is_H0():
    assert physics.nye_thickness(0.0, 1.5e5, C, thickness0=3.0) == pytest.approx(3.0)


def test_nye_distance_inverts_nye_thickness():
    """nye_distance is the exact inverse of nye_thickness."""
    tau = 8.0e4
    H0, H1 = 2.0, 1200.0
    x = physics.nye_distance(H0, H1, tau, C)
    assert physics.nye_thickness(x, tau, C, thickness0=H0) == pytest.approx(H1, rel=1e-12)


def test_stress_length_and_slope_magnitude():
    tau = 1.0e5
    hf = physics.stress_length(tau, C)
    assert hf == pytest.approx(tau / (C.rho_ice * C.g))
    # |grad S| = Hf / H
    H = 500.0
    assert physics.slope_magnitude(tau, H, C) == pytest.approx(hf / H)


# --------------------------------------------------------------------------- #
# The M1 radial solve: integrate the flowline RHS, compare to Nye             #
# --------------------------------------------------------------------------- #

def _integrate_straight_flowline(tau, x_max, H0, constants, n_eval=200):
    """Integrate (y, E, q) inward along a straight flowline; flat bed, const tau.

    Bed B = 0 so thickness H = E. With q = 0 and dHf/dy = 0 the flowline stays
    straight (dy/dx = 0, dq/dx = 0) and only E grows. Returns (x, E).
    """
    hf = float(physics.stress_length(tau, constants))

    def rhs(x, state):
        y, E, q = state
        H = E  # flat bed at zero
        dy, dE, dq = physics.slope_derivatives(hf, H, q, dB_dy=0.0, dHf_dy=0.0)
        return [dy, dE, dq]

    x_eval = np.linspace(0.0, x_max, n_eval)
    sol = solve_ivp(
        rhs, (0.0, x_max), [0.0, H0, 0.0],
        t_eval=x_eval, rtol=1e-10, atol=1e-12, method="RK45",
    )
    assert sol.success, sol.message
    return sol.t, sol.y


def test_flowline_integration_reproduces_nye_parabola():
    """Numerically integrating slope_derivatives reproduces the Nye parabola."""
    tau = 1.0e5
    x_max = 4.0e5   # 400 km
    H0 = 1.0
    x, state = _integrate_straight_flowline(tau, x_max, H0, C)
    y, E, q = state

    expected = physics.nye_thickness(x, tau, C, thickness0=H0)
    # Surface elevation (= thickness here) matches the analytic parabola.
    np.testing.assert_allclose(E, expected, rtol=1e-4, atol=1e-2)

    # A straight flowline stays straight: no lateral drift, no cross-flow slope.
    assert np.max(np.abs(y)) < 1e-6
    assert np.max(np.abs(q)) < 1e-12


def test_flowline_matches_nye_at_specific_distance():
    """Spot-check a physically meaningful number: thickness at 300 km."""
    tau = 1.0e5
    x, state = _integrate_straight_flowline(tau, 3.0e5, 1.0, C)
    y, E, q = state
    # analytic: sqrt(1 + 2*1e5*3e5/8992.7) ~ 2582 m
    analytic = physics.nye_thickness(3.0e5, tau, C, thickness0=1.0)
    assert E[-1] == pytest.approx(analytic, rel=1e-4)
    assert 2500.0 < E[-1] < 2650.0  # sane ice thickness


# --------------------------------------------------------------------------- #
# Singularity guard and marine margin                                         #
# --------------------------------------------------------------------------- #

def test_slope_derivatives_raises_at_turning_singularity():
    """q >= Hf/H is the flowline-turning singularity; must raise, not NaN."""
    hf, H = 100.0, 500.0
    q_crit = hf / H  # = |grad S|; p = 0 here
    with pytest.raises(ValueError):
        physics.slope_derivatives(hf, H, q_crit * 1.001)


def test_slope_derivatives_constant_tau_keeps_q_zero_rhs():
    """With q=0, flat bed, const tau: dy=0, dq=0, dE=|grad S|."""
    hf, H = 100.0, 500.0
    dy, dE, dq = physics.slope_derivatives(hf, H, 0.0, dB_dy=0.0, dHf_dy=0.0)
    assert dy == pytest.approx(0.0)
    assert dq == pytest.approx(0.0)
    assert dE == pytest.approx(hf / H)  # = |grad S|


def test_marine_flotation_surface():
    # Below sea level -> near-flotation surface B*(1 - rho_w/rho_i) > 0.
    B = -500.0
    S = float(physics.marine_flotation_surface(B, C))
    assert S == pytest.approx(B * (1.0 - C.rho_water / C.rho_ice))
    assert S > 0.0
    # Above sea level -> unchanged.
    assert float(physics.marine_flotation_surface(250.0, C)) == pytest.approx(250.0)
