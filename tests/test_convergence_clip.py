"""Divide over-enclosure: the convergence clip (docs/design-note-04).

Where flowlines converge on an ice divide whose medial axis is a *line* (not a
single point), GEOS ``make_valid`` takes the outer envelope of the crossed
front and encloses area *past* the divide — a spurious summit plateau of ice
above the physically maximum possible elevation. The clip
(``ContourManager.clip_convergence``) trims the advanced front to the parent
eroded by the local advance distance, removing that over-enclosure.

These use a flat bed and constant shear stress, for which the exact steady
surface is the eikonal ``S(x) = sqrt(2 tau / (rho_i g) * d(x))`` with ``d`` the
distance inward from the margin. The maximum possible elevation is therefore
fixed by the largest distance-to-margin, and *no* ice may exist above it.
"""

import numpy as np
from shapely.geometry import LineString, Point

from pyicesheet import RasterField, IceSheetModel, ModelConfig, physics
from pyicesheet.constants import DEFAULT_CONSTANTS as C

TAU = 1.0e5


def _flat_uniform_fields(hx, hy, n=121):
    x = np.linspace(-hx, hx, n)
    y = np.linspace(-hy, hy, n)
    XX, _ = np.meshgrid(x, y)
    bed = RasterField(x, y, np.zeros_like(XX))
    tau_f = RasterField(x, y, np.full_like(XX, TAU))
    return bed, tau_f


def _stadium(half_len, half_width):
    """A capsule: its medial axis is the central segment (a line divide)."""
    return LineString([(-half_len, 0.0), (half_len, 0.0)]).buffer(
        half_width, quad_segs=64)


def _solve_stadium(clip, half_len=3.0e5, half_width=1.2e5):
    margin = _stadium(half_len, half_width)
    bed, tau_f = _flat_uniform_fields(5.0e5, 3.0e5)
    cfg = ModelConfig(spacing=8000.0, elevation_interval=300.0,
                      max_elevation=4000.0, clip_convergence=clip)
    surf = IceSheetModel(bed, tau_f, margin, cfg).solve()
    # Exact eikonal cap: the deepest interior point is half_width from the margin.
    cap = float(physics.nye_thickness(half_width, TAU, C, thickness0=1.0))
    return surf, cap


def test_clip_removes_over_enclosure_at_a_line_divide():
    """With the clip on, no ice is built above the physical eikonal cap."""
    surf, cap = _solve_stadium(clip=True)
    # The reconstruction may undershoot the cap by up to one interval (it stops
    # at the last full contour), but must never exceed it: ice above the maximum
    # possible elevation is unphysical over-enclosure.
    assert surf.elevation.max() <= cap + 1.0, (
        f"clip on: summit {surf.elevation.max():.0f} exceeds physical cap "
        f"{cap:.0f} -- over-enclosure not removed")


def test_unclipped_front_over_encloses_past_the_divide():
    """Regression guard: without the clip the front over-builds past the cap.

    This documents the bug the clip fixes; if this ever stops holding, the
    over-enclosure has changed character and the clip test above needs review.
    """
    surf, cap = _solve_stadium(clip=False)
    assert surf.elevation.max() > cap, (
        f"clip off: summit {surf.elevation.max():.0f} did NOT exceed cap "
        f"{cap:.0f} -- expected the make_valid over-enclosure")


def test_clip_leaves_the_radial_cap_unchanged():
    """A circular margin converges to a *point*; the clip must not touch it."""
    radius = 2.0e5
    half = 2.0 * radius
    x = np.linspace(-half, half, 41)
    y = np.linspace(-half, half, 41)
    XX, _ = np.meshgrid(x, y)
    bed = RasterField(x, y, np.zeros_like(XX))
    tau_f = RasterField(x, y, np.full_like(XX, TAU))
    margin = Point(0.0, 0.0).buffer(radius, quad_segs=64)

    def solve(clip):
        cfg = ModelConfig(spacing=10_000.0, elevation_interval=200.0,
                          max_elevation=5000.0, clip_convergence=clip)
        return IceSheetModel(bed, tau_f, margin, cfg).solve()

    on, off = solve(True), solve(False)
    assert on.x.size == off.x.size
    np.testing.assert_allclose(np.sort(on.elevation), np.sort(off.elevation))
