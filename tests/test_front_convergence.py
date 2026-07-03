"""M3: crowding / divide formation in isolation (a dumbbell margin).

Where flowlines converge, ice divides form and pinched contours split into
separate lobes. This is the part of the algorithm we are most wary of, so it is
tested on a controlled two-lobe (dumbbell) margin *before* the full Greenland
case. The GEOS ``make_valid`` step is expected to split a pinching contour into
two polygons automatically.
"""

import numpy as np
from shapely.geometry import Point

from pyicesheet import RasterField, IceSheetModel, ModelConfig
from pyicesheet.contour import Contour, ContourManager
from pyicesheet.flowline import FlowlineIntegrator
from pyicesheet._geometry import resample_ring, inward_normals


def _flat_uniform_fields(half=6.0e5, n=61, tau=1.0e5):
    x = np.linspace(-half, half, n)
    y = np.linspace(-half, half, n)
    XX, _ = np.meshgrid(x, y)
    bed = RasterField(x, y, np.zeros_like(XX))
    tau_f = RasterField(x, y, np.full_like(XX, tau))
    return bed, tau_f


def _dumbbell(radius=1.0e5, half_sep=9.9e4):
    """Two overlapping circles -> a single peanut polygon with a narrow neck."""
    left = Point(-half_sep, 0.0).buffer(radius, quad_segs=64)
    right = Point(half_sep, 0.0).buffer(radius, quad_segs=64)
    return left.union(right)


def test_manager_splits_a_pinching_contour():
    """Advancing a narrow-necked peanut contour eventually yields two lobes."""
    bed, tau_f = _flat_uniform_fields()
    integ = FlowlineIntegrator(bed, tau_f)
    mgr = ContourManager(integ, spacing=8000.0, elevation_interval=200.0)

    poly = _dumbbell(radius=1.0e5, half_sep=9.9e4)
    x, y = resample_ring(poly, mgr.spacing)
    inx, iny = inward_normals(x, y, poly, mgr.spacing)
    E0 = 1000.0
    contour = Contour(x=x, y=y, E=np.full(x.size, E0),
                      inward_x=inx, inward_y=iny, polygon=poly, level=E0)

    contours = [contour]
    n_lobes = 1
    for step in range(1, 8):
        target = E0 + step * 200.0
        nxt = []
        for c in contours:
            nxt.extend(mgr.advance(c, target))
        if not nxt:
            break
        contours = nxt
        n_lobes = max(n_lobes, len(contours))
        if n_lobes >= 2:
            break
    assert n_lobes >= 2, "the pinching neck never split into separate lobes"


def test_dumbbell_reconstructs_two_domes():
    """End-to-end: a dumbbell margin gives two spatially-separated ice domes."""
    bed, tau_f = _flat_uniform_fields()
    margin = _dumbbell(radius=1.2e5, half_sep=1.15e5)
    cfg = ModelConfig(spacing=10_000.0, elevation_interval=200.0, max_elevation=5000.0)
    surf = IceSheetModel(bed, tau_f, margin, cfg).solve()

    # Look at the upper part of the reconstruction; it should occur over both
    # lobe centres (x<0 and x>0) and be absent over the neck (x~0).
    hi = surf.elevation > 0.7 * surf.elevation.max()
    xhi = surf.x[hi]
    assert np.any(xhi < -5.0e4), "no high ice over the left lobe"
    assert np.any(xhi > 5.0e4), "no high ice over the right lobe"
    # The two highs are separated: few high samples straddle the neck.
    near_neck = np.sum(np.abs(xhi) < 2.0e4)
    assert near_neck < 0.15 * xhi.size
