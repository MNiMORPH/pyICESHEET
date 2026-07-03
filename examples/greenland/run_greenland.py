"""Reconstruct the Greenland ice sheet with pyICESHEET.

Inputs (see README.md for provenance):
  * bed topography   : BedMachineGreenland-v6.nc  (NSIDC; EPSG:3413, 150 m)
  * ice margin       : outline5.shp               (from the ICESHEET repo)
  * basal shear stress: shear_stress.shp          (polygons, Pa; ICESHEET repo)

All three are in the same NSIDC polar-stereographic CRS (EPSG:3413), so no
reprojection is needed — coordinates are metres throughout.

The shear-stress polygons are rasterized and then *smoothed* (a Gaussian filter):
smoothing is a deliberate upstream operator, kept out of the solver, because the
flowline integrator differentiates the field.

Usage:
    python run_greenland.py [--bed PATH] [--resolution-factor N] [--spacing M]
"""

from __future__ import annotations

import argparse
import time

import numpy as np
from scipy.ndimage import gaussian_filter

from pyicesheet import RasterField, IceSheetModel, ModelConfig
from pyicesheet.io.raster import read_netcdf_downsampled, fill_invalid
from pyicesheet.io.vector import read_polygon, read_polygons, rasterize_polygons

DEFAULT_BED = "/home/awickert/Downloads/BedMachineGreenland-v6.nc"
SS_DIR = "/home/awickert/models/icesheet/Greenland_final/shear_stress/qgis"


def build_fields(bed_path, factor, tau_sigma=1.5, bed_sigma=1.0):
    """Load bed + shear stress onto a common grid; smooth (upstream operator)."""
    xb, yb, bed = read_netcdf_downsampled(bed_path, "bed", factor=factor)
    bed = fill_invalid(bed, method="nearest")

    ss = read_polygons(f"{SS_DIR}/shear_stress.shp")
    tau = rasterize_polygons(ss.geometry, ss["shear_stre"], xb, yb, fill=np.nan)
    tau = fill_invalid(tau, method="nearest")

    tau = gaussian_filter(tau, sigma=tau_sigma)
    bed = gaussian_filter(bed, sigma=bed_sigma)

    bed_f = RasterField.from_arrays(xb, yb, bed)
    tau_f = RasterField.from_arrays(xb, yb, tau)
    return bed_f, tau_f, (xb, yb)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", default=DEFAULT_BED)
    ap.add_argument("--resolution-factor", type=int, default=33,
                    help="downsample factor from the 150 m bed (33 ~ 5 km)")
    ap.add_argument("--spacing", type=float, default=15000.0)
    ap.add_argument("--interval", type=float, default=200.0)
    ap.add_argument("--out", default="greenland_recon.npz")
    args = ap.parse_args()

    bed_f, tau_f, _ = build_fields(args.bed, args.resolution_factor)
    margin = read_polygon(f"{SS_DIR}/outline5.shp")

    cfg = ModelConfig(spacing=args.spacing, elevation_interval=args.interval,
                      max_elevation=4000.0, min_thickness=1.0)
    model = IceSheetModel(bed_f, tau_f, margin, cfg)

    t0 = time.time()
    surf = model.solve(progress=True)
    print(f"solve: {time.time() - t0:.1f} s, {len(surf)} samples, "
          f"max surface {surf.elevation.max():.0f} m")

    np.savez(args.out, x=surf.x, y=surf.y, elevation=surf.elevation,
             thickness=surf.thickness, bed=surf.bed)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
