"""Flowline integration: step one flowline inward to the next elevation contour.

Given a start point on a contour, the local inward (along-flow) direction, and the
current surface elevation, :class:`FlowlineIntegrator` integrates the perfectly-
plastic flowline equations (:func:`pyicesheet.physics.slope_derivatives`) inward
until the surface reaches a target elevation, sampling the bed and shear-stress
fields along the way.

Local frame
-----------
The integration runs along a *fixed* along-flow axis ``t = (cos theta, sin theta)``
set at the start of the step; ``n = (-sin theta, cos theta)`` is the cross-flow
axis. The independent variable is arc length ``x_l`` along ``t``; the state is
``(y_l, E, q)`` — lateral offset along ``n``, surface elevation, and cross-flow
slope. The map position is

    X = x0 + x_l * cos(theta) + y_l * (-sin(theta))
    Y = y0 + x_l * sin(theta) + y_l * ( cos(theta)).

After the step, the along-flow direction is updated to follow steepest ascent:
``theta_new = theta + atan2(q_end, p_end)`` (as in the original ICESHEET).

Where the Fortran hand-rolled an adaptive Runge–Kutta with step doubling, this
uses :func:`scipy.integrate.solve_ivp` (adaptive RK45) with a terminal event at
the target elevation — the same numerics, better tested.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from . import physics
from .constants import DEFAULT_CONSTANTS, PhysicalConstants
from .fields import RasterField

__all__ = ["FlowlineIntegrator", "FlowlineStep"]


@dataclass
class FlowlineStep:
    """Result of integrating one flowline to a target elevation.

    Attributes
    ----------
    x, y : float
        Map coordinates of the endpoint.
    elevation : float
        Surface elevation at the endpoint (== target if ``reached_target``).
    q : float
        Cross-flow slope at the endpoint.
    direction : float
        Updated along-flow azimuth (radians) for the next step.
    arc_length : float
        Along-axis arc length ``x_l`` travelled.
    reached_target : bool
        Whether the target elevation was reached (vs. stalling on a stop
        condition such as the ice thinning to nothing).
    status : str
        One of ``"reached"``, ``"too_thin"``, ``"too_long"``, ``"singularity"``.
    """

    x: float
    y: float
    elevation: float
    q: float
    direction: float
    arc_length: float
    reached_target: bool
    status: str


class FlowlineIntegrator:
    """Integrate flowlines through given bed and shear-stress fields.

    Parameters
    ----------
    bed : RasterField
        Bed topography ``B`` (m).
    tau : RasterField
        Basal shear stress ``tau_b`` (Pa).
    constants : PhysicalConstants, optional
    max_arc_length : float, optional
        Safety cap on a single step's arc length (m); default 5e6 (5000 km).
    min_thickness : float, optional
        Stop the step if ice thins to this (m); default 0.1.
    rtol, atol : float, optional
        Integrator tolerances.
    """

    def __init__(self, bed: RasterField, tau: RasterField,
                 constants: PhysicalConstants = DEFAULT_CONSTANTS,
                 max_arc_length: float = 5.0e6, min_thickness: float = 0.1,
                 rtol: float = 1e-8, atol: float = 1e-6):
        self.bed = bed
        self.tau = tau
        self.constants = constants
        self.max_arc_length = float(max_arc_length)
        self.min_thickness = float(min_thickness)
        self.rtol = rtol
        self.atol = atol

    def step_to_elevation(self, x0, y0, theta0, elevation0, target_elevation,
                          q0=0.0) -> FlowlineStep:
        """Integrate from ``(x0, y0)`` at ``elevation0`` up to ``target_elevation``.

        Returns a :class:`FlowlineStep`. ``target_elevation`` must exceed
        ``elevation0``.
        """
        if target_elevation <= elevation0:
            raise ValueError("target_elevation must exceed elevation0")

        rho_g = self.constants.rho_g
        cos_t, sin_t = np.cos(theta0), np.sin(theta0)
        # cross-flow (left-normal) unit vector
        nx, ny = -sin_t, cos_t

        def to_map(x_l, y_l):
            X = x0 + x_l * cos_t + y_l * nx
            Y = y0 + x_l * sin_t + y_l * ny
            return X, Y

        def rhs(x_l, state):
            y_l, E, q = state
            X, Y = to_map(x_l, y_l)
            B = float(self.bed.value(X, Y))
            H = E - B
            if H <= 0.0:
                # Degenerate; let the thin-ice event catch it. Return ~flat.
                H = self.min_thickness
            tau_val = float(self.tau.value(X, Y))
            hf = tau_val / rho_g
            dB_dy = float(self.bed.directional_derivative(X, Y, nx, ny))
            dtau_dy = float(self.tau.directional_derivative(X, Y, nx, ny))
            dHf_dy = dtau_dy / rho_g
            try:
                dy, dE, dq = physics.slope_derivatives(hf, H, q, dB_dy, dHf_dy)
            except ValueError:
                # Flowline-turning singularity: signal via a huge dq so the
                # solver shrinks the step; the singularity event will fire.
                return [0.0, 0.0, np.sign(q) * 1e12 if q else 1e12]
            return [dy, dE, dq]

        def reach_target(x_l, state):
            return state[1] - target_elevation
        reach_target.terminal = True
        reach_target.direction = 1.0

        def too_thin(x_l, state):
            X, Y = to_map(x_l, state[0])
            return (state[1] - float(self.bed.value(X, Y))) - self.min_thickness
        too_thin.terminal = True
        too_thin.direction = -1.0

        def singular(x_l, state):
            y_l, E, q = state
            X, Y = to_map(x_l, y_l)
            H = E - float(self.bed.value(X, Y))
            H = max(H, self.min_thickness)
            hf = float(self.tau.value(X, Y)) / rho_g
            # p^2 = (Hf/H)^2 - q^2; stop just before it hits zero.
            return (hf / H) ** 2 - q * q - 1e-9
        singular.terminal = True
        singular.direction = -1.0

        sol = solve_ivp(
            rhs, (0.0, self.max_arc_length), [0.0, elevation0, q0],
            events=(reach_target, too_thin, singular),
            rtol=self.rtol, atol=self.atol, method="RK45", dense_output=False,
            max_step=self.max_arc_length,
        )

        # Endpoint from the terminating event (or the final state).
        status, reached = self._classify(sol)
        y_l_end, E_end, q_end = self._endpoint_state(sol)
        x_l_end = self._endpoint_t(sol)
        X_end, Y_end = to_map(x_l_end, y_l_end)

        # Updated along-flow direction: rotate toward steepest ascent.
        B_end = float(self.bed.value(X_end, Y_end))
        H_end = max(E_end - B_end, self.min_thickness)
        hf_end = float(self.tau.value(X_end, Y_end)) / rho_g
        ratio2 = (hf_end / H_end) ** 2
        if q_end * q_end < ratio2:
            p_end = np.sqrt(ratio2 - q_end * q_end)
            theta_new = theta0 + np.arctan2(q_end, p_end)
        else:
            theta_new = theta0

        return FlowlineStep(
            x=X_end, y=Y_end, elevation=E_end, q=q_end,
            direction=float(theta_new), arc_length=float(x_l_end),
            reached_target=reached, status=status,
        )

    # -- helpers --------------------------------------------------------- #

    @staticmethod
    def _classify(sol):
        # events order: (reach_target, too_thin, singular)
        if len(sol.t_events[0]) > 0:
            return "reached", True
        if len(sol.t_events[1]) > 0:
            return "too_thin", False
        if len(sol.t_events[2]) > 0:
            return "singularity", False
        return "too_long", False

    @staticmethod
    def _endpoint_state(sol):
        for ye in sol.y_events:
            if len(ye) > 0:
                return ye[-1]
        return sol.y[:, -1]

    @staticmethod
    def _endpoint_t(sol):
        for te in sol.t_events:
            if len(te) > 0:
                return te[-1]
        return sol.t[-1]
