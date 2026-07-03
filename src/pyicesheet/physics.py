"""Perfectly-plastic ice-sheet physics.

This module holds the *physics only* — a small set of pure functions with no
knowledge of grids, geometry, or I/O. Everything the reconstruction needs about
ice mechanics lives here, which makes it directly testable against analytic
solutions (see ``tests/test_physics_nye.py``).

Governing assumption (perfect plasticity / shallow-ice limit): basal shear stress
equals driving stress,

    tau_b = rho_i * g * H * |grad S|,

with ``H = S - B`` the ice thickness (surface ``S`` minus bed ``B``). Hence the
surface-slope magnitude is ``|grad S| = tau_b / (rho_i g H)``. It is convenient to
carry the *stress length* ``Hf = tau_b / (rho_i g)`` (a length), so that
``|grad S| = Hf / H``.

Flowline frame
--------------
Following Gowan's ICESHEET (after Reeh 1982 and Fisher et al. 1985), the surface
is built by integrating inward along flowlines. In a local frame fixed to the
inward-normal at a point, two surface-slope components are tracked:

* ``p`` — along-flow slope,
* ``q`` — cross-flow slope,

with ``|grad S|**2 = p**2 + q**2``, so ``p = sqrt((Hf/H)**2 - q**2)``. Integrating
along the local along-flow axis ``x_l``, the state ``(y, E, q)`` (lateral offset,
surface elevation, cross-flow slope) evolves as

    dy/dx_l = q / p
    dE/dx_l = Hf**2 / (H**2 * p)                                  (= |grad S|**2 / p)
    dq/dx_l = Hf**2 * (dB/dy - q) / (H**3 * p)
              + Hf * (dHf/dy) / (p * H**2)                        (Fisher eq. A8)

The second ``dq`` term (the shear-stress gradient) is what distinguishes the
ICESHEET reconstruction from a constant-stress one; drop it (``dHf/dy = 0``) and
the equations reduce to the classic constant-shear-stress flowline.
"""

from __future__ import annotations

import numpy as np

from .constants import DEFAULT_CONSTANTS, PhysicalConstants

__all__ = [
    "stress_length",
    "slope_magnitude",
    "along_flow_slope",
    "slope_derivatives",
    "nye_thickness",
    "nye_distance",
    "marine_flotation_surface",
]


# --------------------------------------------------------------------------- #
# Fundamental stress / slope relations                                        #
# --------------------------------------------------------------------------- #

def stress_length(tau, constants: PhysicalConstants = DEFAULT_CONSTANTS):
    """Stress length ``Hf = tau / (rho_i g)`` (metres).

    ``Hf`` has units of length and equals the ice thickness at which a unit
    surface slope would balance the basal shear stress ``tau`` (Pa).
    """
    return np.asarray(tau, dtype=float) / constants.rho_g


def slope_magnitude(tau, thickness, constants: PhysicalConstants = DEFAULT_CONSTANTS):
    """Surface-slope magnitude ``|grad S| = tau / (rho_i g H)`` (dimensionless).

    Parameters
    ----------
    tau : array_like
        Basal shear stress (Pa).
    thickness : array_like
        Ice thickness ``H = S - B`` (m); must be positive.
    """
    thickness = np.asarray(thickness, dtype=float)
    return stress_length(tau, constants) / thickness


def along_flow_slope(hf, thickness, q):
    """Along-flow slope ``p = sqrt((Hf/H)**2 - q**2)``.

    ``p`` is the along-flow component of the surface gradient in the local
    flowline frame. It becomes zero (a turning/tangent condition) when the
    cross-flow slope ``q`` reaches the full slope magnitude ``Hf/H``; callers
    must guard against ``q**2 >= (Hf/H)**2``.
    """
    hf = np.asarray(hf, dtype=float)
    thickness = np.asarray(thickness, dtype=float)
    q = np.asarray(q, dtype=float)
    return np.sqrt((hf / thickness) ** 2 - q ** 2)


def slope_derivatives(hf, thickness, q, dB_dy=0.0, dHf_dy=0.0):
    """Flowline ODE right-hand side.

    Returns the along-``x_l`` derivatives ``(dy_dx, dE_dx, dq_dx)`` of the state
    ``(y, E, q)``, given the local ice mechanics. This is the pure-physics heart
    of the reconstruction; :mod:`pyicesheet.flowline` wraps it with field
    sampling for ``B``, ``Hf``, ``dB/dy`` and ``dHf/dy``.

    Parameters
    ----------
    hf : float
        Stress length ``Hf = tau / (rho_i g)`` (m) at the point.
    thickness : float
        Ice thickness ``H = E - B`` (m); must be positive.
    q : float
        Cross-flow surface slope (dimensionless).
    dB_dy : float, optional
        Cross-flow bed slope ``dB/dy`` (dimensionless).
    dHf_dy : float, optional
        Cross-flow gradient of the stress length ``dHf/dy`` (dimensionless);
        equals ``(dtau/dy) / (rho_i g)``. Zero for constant shear stress.

    Returns
    -------
    (dy_dx, dE_dx, dq_dx) : tuple of float

    Raises
    ------
    ValueError
        If ``q**2 >= (Hf/H)**2`` (the flowline-turning singularity), where the
        along-flow slope ``p`` is undefined.
    """
    ratio2 = (hf / thickness) ** 2
    if q * q >= ratio2:
        raise ValueError(
            "flowline-turning singularity: q**2 >= (Hf/H)**2 "
            f"(q={q!r}, Hf/H={np.sqrt(ratio2)!r})"
        )
    p = np.sqrt(ratio2 - q * q)

    dy_dx = q / p
    dE_dx = hf ** 2 / (thickness ** 2 * p)
    dq_dx = (
        hf ** 2 * (dB_dy - q) / (thickness ** 3 * p)
        + hf * dHf_dy / (p * thickness ** 2)
    )
    return dy_dx, dE_dx, dq_dx


# --------------------------------------------------------------------------- #
# Nye (constant shear stress, flat bed) analytic solution                     #
# --------------------------------------------------------------------------- #

def nye_thickness(x, tau, constants: PhysicalConstants = DEFAULT_CONSTANTS,
                  thickness0=0.0):
    """Analytic Nye thickness profile for uniform ``tau`` and a flat bed.

    Integrating the flowline equations with ``q = 0`` and ``dB/dy = 0`` gives the
    parabolic profile

        H(x)**2 = H0**2 + (2 tau / (rho_i g)) * x,

    i.e. ``H = sqrt(H0**2 + 2 tau x / (rho_i g))``, where ``x`` is distance inward
    from the margin and ``H0`` the (nominal) margin thickness. See Paterson,
    *The Physics of Glaciers*, 3rd ed., ch. 11, eq. 4.

    Parameters
    ----------
    x : array_like
        Distance inward from the margin (m), ``x >= 0``.
    tau : float
        Basal shear stress (Pa).
    thickness0 : float, optional
        Margin thickness ``H0`` (m); default 0.
    """
    x = np.asarray(x, dtype=float)
    return np.sqrt(thickness0 ** 2 + 2.0 * float(tau) * x / constants.rho_g)


def nye_distance(thickness0, thickness1, tau,
                 constants: PhysicalConstants = DEFAULT_CONSTANTS):
    """Inward distance to grow from ``thickness0`` to ``thickness1`` under Nye.

    Inverse of :func:`nye_thickness`:

        x = (H1**2 - H0**2) * rho_i g / (2 tau).

    Used to estimate along-contour point spacing. This includes ``g`` — a
    deliberate correction over the original Fortran, which drops ``g`` from the
    equivalent resampling estimate (see ``docs/design-note-01`` section 6).
    """
    return (
        (float(thickness1) ** 2 - float(thickness0) ** 2)
        * constants.rho_g
        / (2.0 * float(tau))
    )


# --------------------------------------------------------------------------- #
# Marine margin                                                               #
# --------------------------------------------------------------------------- #

def marine_flotation_surface(bed, constants: PhysicalConstants = DEFAULT_CONSTANTS):
    """Surface elevation of grounded ice near flotation at a marine margin.

    Where the bed is below sea level, ICESHEET seeds the margin ice at a
    thickness in near-hydrostatic balance with the water column rather than at
    zero. For bed elevation ``B < 0`` this gives surface

        S = B * (1 - rho_w / rho_i),

    a crude grounded-but-near-flotation approximation (a true ice shelf is not
    modelled). For ``B >= 0`` the surface is just ``B``.
    """
    bed = np.asarray(bed, dtype=float)
    flot = bed * (1.0 - constants.rho_water / constants.rho_ice)
    return np.where(bed < 0.0, flot, bed)
