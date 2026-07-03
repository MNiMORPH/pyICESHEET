"""Calibrate the basal shear stress against the observed Greenland surface.

Demonstrates the calibration layer (:mod:`pyicesheet.calibrate`) and the
tau-from-substrate builder (:mod:`pyicesheet.tau`) on the fast grid solver:

  * a single global multiplier (alpha) on Gowan's tau, then
  * a per-substrate-class multiplier vector, using Gowan's discrete shear-stress
    polygons as a proxy substrate classification (his values are geology-derived).

The per-class fit attacks the spatial pattern the scalar cannot. With the real
Gowan et al. (2019, ESSD 11, 375) geology rasters the class map would come from
sediment/bedrock type instead of the polygons; the machinery is identical.

Usage:
    python calibrate_tau.py --fields greenland_fields.npz
"""

from __future__ import annotations

import argparse

import numpy as np
import geopandas as gpd
from scipy.ndimage import gaussian_filter, distance_transform_edt
from rasterio.features import rasterize
from rasterio.transform import from_origin

from pyicesheet.grid import solve_surface_grid
from pyicesheet.tau import tau_from_classes
from pyicesheet.calibrate import SurfaceMisfit, calibrate, calibrate_scalar

SS_DIR = "/home/awickert/models/icesheet/Greenland_final/shear_stress/qgis"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fields", default="greenland_fields.npz")
    ap.add_argument("--tau-sigma", type=float, default=1.5)
    args = ap.parse_args()

    d = np.load(args.fields)
    x = d["xb"].astype(float)
    y = d["yb"][::-1].astype(float)
    bed = d["bed"][::-1, :]
    obs = d["surf_obs"][::-1, :]
    margin = gpd.read_file(f"{SS_DIR}/outline5.shp").geometry.iloc[0]

    # Substrate class map (proxy): rasterize Gowan's discrete shear-stress polygons.
    ss = gpd.read_file(f"{SS_DIR}/shear_stress.shp")
    classes = np.array(sorted(set(ss["shear_stre"].astype(float))))
    h = abs(x[1] - x[0])
    tr = from_origin(x.min() - h / 2, y.max() + h / 2, h, h)
    shapes = [(g, int(np.where(classes == float(v))[0][0]) + 1)
              for g, v in zip(ss.geometry, ss["shear_stre"])]
    cls = rasterize(shapes, out_shape=bed.shape, transform=tr, fill=0,
                    dtype="int32")[::-1, :]
    idx = distance_transform_edt(cls == 0, return_distances=False, return_indices=True)
    cls = np.clip(np.where(cls == 0, cls[tuple(idx)], cls), 1, classes.size)

    mask = obs > 1.0
    obj = SurfaceMisfit(obs, mask=mask)

    def forward_scalar(alpha):
        tau = gaussian_filter(alpha * classes[cls - 1], args.tau_sigma)
        return solve_surface_grid(x, y, bed, tau, margin).surface

    def forward_perclass(m):
        tau = gaussian_filter((classes * m)[cls - 1], args.tau_sigma)
        return solve_surface_grid(x, y, bed, tau, margin).surface

    print(f"Gowan (alpha=1):      RMS {obj(forward_scalar(1.0)):.0f} m")
    rs = calibrate_scalar(forward_scalar, obj, bounds=(0.5, 1.3))
    print(f"scalar optimum:       alpha={rs.x:.3f}  RMS {rs.fun:.0f} m")
    n = classes.size
    rp = calibrate(forward_perclass, np.full(n, rs.x), obj,
                   bounds=[(0.4, 1.6)] * n, method="L-BFGS-B",
                   eps=0.02, maxiter=40)
    print(f"per-class optimum:    RMS {rp.fun:.0f} m  ({rp.nfev} solves)")
    print("  class tau (kPa) -> multiplier:")
    for c, mult in sorted(zip(classes, rp.x), key=lambda t: t[1]):
        print(f"    {int(c/1e3):>3} kPa  x{mult:.2f}")


if __name__ == "__main__":
    main()
