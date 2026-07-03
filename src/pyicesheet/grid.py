"""Plan B: an Eulerian (grid) solver for the perfectly-plastic surface.

Instead of tracing flowlines inward (the Lagrangian approach in
:mod:`pyicesheet.solver`), this solves the *same* stress balance directly on the
grid as a stationary Hamilton–Jacobi / eikonal equation:

    |grad S| = tau_b / (rho_i g (S - B)),

with the surface fixed at the margin (``S = B + min_thickness``, or a marine
flotation value) and marched inward (``S`` increasing). Because ``S`` is a single-
valued field on the grid, ice **divides form automatically** where fronts from
different directions meet — there is no flowline convergence to track and no
over-enclosure to clip. It is also fast (fast-sweeping/iterative, ~O(N)) and
GRASS-native (raster in, raster out).

The one subtlety is that the local speed ``F = tau/(rho_i g (S - B))`` depends on
the unknown ``S`` (through the ice thickness ``H = S - B``). Each cell update
therefore solves the *coupled* local equation, not a fixed-speed eikonal:

* 1-D (one upwind neighbour ``a``): closed form
  ``U = [(a + B) + sqrt((a - B)^2 + 4 tau h / (rho g))] / 2``;
* 2-D (both neighbours ``a <= b`` upwind): a few Newton iterations on
  ``((U-a)^2 + (U-b)^2)(U-B)^2 = (tau h / (rho g))^2``.

For a flat bed and uniform ``tau`` this reduces to ``|grad(S^2)| = 2 tau / (rho g)``,
whose solution is the Nye parabola ``S = sqrt(2 tau d / (rho g))`` with ``d`` the
distance from the margin — the analytic check in the tests.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import physics
from .constants import DEFAULT_CONSTANTS, PhysicalConstants

__all__ = ["GridSurface", "solve_surface_grid"]


@dataclass
class GridSurface:
    """Result of a grid solve."""

    x: np.ndarray            # 1-D cell-centre coordinates (ascending)
    y: np.ndarray            # 1-D cell-centre coordinates (ascending)
    surface: np.ndarray      # 2-D ice-surface elevation, NaN outside the ice
    bed: np.ndarray          # 2-D bed elevation
    iterations: int

    @property
    def thickness(self):
        H = self.surface - self.bed
        return np.where(np.isfinite(self.surface), H, np.nan)


def _coupled_update(a, b, B, tau, h, rho_g, newton_iter=4):
    """Godunov eikonal update with speed F = tau/(rho_g (U - B)).

    ``a <= b`` are the smaller upwind neighbour values per axis (may be +inf).
    Returns the updated value ``U`` (+inf where no finite neighbour exists).
    """
    C = tau * h / rho_g                       # (tau h / rho g); F h = C/(U-B)
    # 1-D update (uses only the smaller neighbour a): (U-a)(U-B) = C
    disc1 = (a - B) ** 2 + 4.0 * C
    U1 = 0.5 * ((a + B) + np.sqrt(np.maximum(disc1, 0.0)))

    # Where the second neighbour is also upwind (b < U1), solve the 2-D equation.
    two_d = np.isfinite(b) & (b < U1)
    U = U1.copy()
    if np.any(two_d):
        aa, bb, BB, CC = a[two_d], b[two_d], B[two_d], C[two_d]
        u = U1[two_d].copy()
        for _ in range(newton_iter):
            s = (u - aa) ** 2 + (u - bb) ** 2
            d = u - BB
            g = s * d * d - CC * CC
            gp = (2.0 * (u - aa) + 2.0 * (u - bb)) * d * d + s * 2.0 * d
            step = np.where(gp != 0.0, g / gp, 0.0)
            u = u - step
        U[two_d] = u

    U = np.where(np.isfinite(a), U, np.inf)
    return U


def _solve(bed, tau, inside, boundary, boundary_S, h,
           constants=DEFAULT_CONSTANTS, min_thickness=1.0,
           tol=1e-2, max_iter=20000):
    """Fast-iterative solve of the state-dependent eikonal on a regular grid."""
    rho_g = constants.rho_g
    ny, nx = bed.shape
    S = np.full((ny, nx), np.inf)
    S[boundary] = boundary_S[boundary]

    update = inside & ~boundary
    INF = np.inf
    it = 0
    for it in range(1, max_iter + 1):
        up = np.pad(S, ((1, 0), (0, 0)), constant_values=INF)[:-1, :]
        down = np.pad(S, ((0, 1), (0, 0)), constant_values=INF)[1:, :]
        left = np.pad(S, ((0, 0), (1, 0)), constant_values=INF)[:, :-1]
        right = np.pad(S, ((0, 0), (0, 1)), constant_values=INF)[:, 1:]
        Uy = np.minimum(up, down)
        Ux = np.minimum(left, right)
        a = np.minimum(Ux, Uy)
        b = np.maximum(Ux, Uy)

        U = _coupled_update(a, b, bed, tau, h, rho_g)
        prev = S[update]
        cand = np.minimum(prev, U[update])
        # Convergence: max change over reached cells; a cell that just became
        # finite (was +inf) counts as not-yet-converged.
        both = np.isfinite(prev) & np.isfinite(cand)
        newly = np.isfinite(cand) & ~np.isfinite(prev)
        diff = np.zeros_like(cand)
        diff[both] = np.abs(cand[both] - prev[both])
        diff[newly] = np.inf
        change = diff.max() if diff.size else 0.0
        S_new = S.copy()
        S_new[update] = cand
        S = S_new
        if change < tol:
            break

    # Enforce the physical floor (S >= B + min_thickness) and mask outside.
    S = np.where(inside, np.maximum(S, bed + min_thickness), np.nan)
    S[~inside] = np.nan
    return S, it


def solve_surface_grid(x, y, bed, tau, margin, config=None,
                       constants: PhysicalConstants = DEFAULT_CONSTANTS):
    """Reconstruct the ice surface on a grid from bed, shear stress, and margin.

    Parameters
    ----------
    x, y : 1-D arrays
        Cell-centre coordinates (ascending), square cells (dx == dy).
    bed, tau : 2-D arrays, shape (len(y), len(x))
        Bed elevation (m) and basal shear stress (Pa), with ``bed[j, i]`` at
        ``(x[i], y[j])``.
    margin : shapely Polygon or MultiPolygon
        Ice margin; cells inside are solved, cells outside are masked.
    config : ModelConfig, optional
        Uses ``min_thickness`` and ``min_elevation`` (marine flotation).

    Returns
    -------
    GridSurface
    """
    from rasterio.features import rasterize
    from rasterio.transform import from_origin

    from .config import ModelConfig
    cfg = config or ModelConfig()

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    bed = np.asarray(bed, dtype=float)
    tau = np.asarray(tau, dtype=float)
    dx = float(abs(x[1] - x[0]))
    dy = float(abs(y[1] - y[0]))
    if not np.isclose(dx, dy, rtol=1e-3):
        raise ValueError(f"grid must have square cells (dx={dx}, dy={dy})")

    # Rasterize the margin. transform is top-left origin, so build a top-down grid.
    y_desc = y[::-1]
    bed_td = bed[::-1, :]
    tau_td = tau[::-1, :]
    west = x.min() - dx / 2
    north = y.max() + dy / 2
    transform = from_origin(west, north, dx, dy)
    inside = rasterize([(margin, 1)], out_shape=bed.shape, transform=transform,
                       fill=0, dtype="uint8").astype(bool)

    # Near-margin boundary band. Fixing S on the inner ring of cells introduces a
    # half-cell error (that ring is already ~dx inside the margin, where the true
    # surface is well above the margin value). Instead, seed the first ring with
    # the local plastic toe using each cell's true sub-cell distance to the margin
    # d: S ~= B + sqrt(2 tau d / (rho g)) (valid locally, bed ~flat over a cell).
    from scipy.ndimage import distance_transform_edt
    d_cells = distance_transform_edt(inside)          # distance to margin, in cells
    band = inside & (d_cells < 1.5)
    d_m = d_cells * dx
    S_base = physics.marine_flotation_surface(bed_td, constants)
    toe = S_base + np.sqrt(2.0 * tau_td * d_m / constants.rho_g)
    boundary = band
    boundary_S = np.where(band, np.maximum(toe, S_base + cfg.min_thickness), np.nan)

    S_td, it = _solve(bed_td, tau_td, inside, boundary, boundary_S, dx,
                      constants=constants, min_thickness=cfg.min_thickness)

    # Return in ascending-y orientation.
    return GridSurface(x=x, y=y, surface=S_td[::-1, :], bed=bed,
                       iterations=it)
