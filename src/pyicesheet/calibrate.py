"""Calibrate a reconstruction against data.

The reconstruction is a *forward operator* — parameters -> ice surface — and this
module fits the parameters to observations. Because the grid solver is fast (a
full Greenland solve is ~0.1 s), a genuine multi-parameter inversion is cheap
(hundreds of solves in a minute).

The design is deliberately open:

* the caller supplies a ``forward(params) -> surface`` closure (it builds the tau
  field from the parameters — e.g. a global scale, or a per-substrate-class table
  via :func:`pyicesheet.tau.tau_from_classes` — and runs the solver);
* the *objective* is a callable ``surface -> misfit`` that can be composed from
  several data constraints. :class:`SurfaceMisfit` (RMS vs an observed surface) is
  provided; relative-sea-level and drainage-routing misfits are intended to plug
  in as additional weighted terms via :func:`composite`.

This keeps the solver, the tau builder, and the data constraints decoupled.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["SurfaceMisfit", "composite", "calibrate", "calibrate_scalar"]


@dataclass
class SurfaceMisfit:
    """RMS misfit between a reconstructed and an observed surface (m)."""

    observed: np.ndarray
    mask: np.ndarray | None = None          # where to compare (e.g. observed ice)
    weight: float = 1.0

    def __call__(self, surface):
        m = np.isfinite(surface)
        if self.mask is not None:
            m = m & self.mask
        r = surface[m] - self.observed[m]
        return self.weight * float(np.sqrt((r ** 2).mean()))


def composite(*objectives):
    """Combine objectives into a single weighted misfit ``surface -> float``.

    Each objective carries its own weight; the total is their sum. Add
    relative-sea-level or drainage-routing misfit terms here as they are built.
    """
    def total(surface):
        return sum(obj(surface) for obj in objectives)
    return total


def calibrate(forward, x0, objective, bounds=None, method="L-BFGS-B", **options):
    """Fit a parameter vector by minimizing ``objective(forward(params))``.

    Parameters
    ----------
    forward : callable
        ``forward(params) -> surface`` (builds tau, runs the solver).
    x0 : array_like
        Initial parameters.
    objective : callable
        ``surface -> misfit``; see :class:`SurfaceMisfit` / :func:`composite`.
    bounds, method, **options
        Passed to :func:`scipy.optimize.minimize`.

    Returns
    -------
    scipy.optimize.OptimizeResult
    """
    from scipy.optimize import minimize

    def f(x):
        return objective(forward(np.asarray(x, dtype=float)))

    return minimize(f, np.asarray(x0, dtype=float), method=method,
                    bounds=bounds, options=options)


def calibrate_scalar(forward, objective, bounds=(0.5, 1.5)):
    """Fit a single scalar parameter (e.g. a global shear-stress multiplier)."""
    from scipy.optimize import minimize_scalar

    return minimize_scalar(lambda a: objective(forward(float(a))),
                           bounds=bounds, method="bounded")
