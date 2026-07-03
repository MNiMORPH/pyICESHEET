"""Greenland head-to-head: original Fortran ICESHEET vs pyICESHEET.

Runs the original Fortran ICESHEET on the SAME bed + shear-stress fields and the
same margin as a pyICESHEET reconstruction, then compares the two surfaces. This
is the method-to-method validation ("basically reproduces Gowan"), not bit-for-bit
(see docs/design-note-02 for the deliberate physics corrections).

The Fortran binary-grid format is written by _fortran_io.write_grid, VALIDATED on
the analytic circular cap (fortran_cap_test.py).

Usage:
    python make_fortran_reference.py --fields greenland_fields.npz \\
        --recon greenland_recon.npz --spacing 40000 --interval 400
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import numpy as np

from _fortran_io import write_grid, write_grid_params, parse_contours


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fields", default="greenland_fields.npz",
                    help="npz with xb, yb (north-first), bed, tau (row0=north)")
    ap.add_argument("--recon", default="greenland_recon.npz",
                    help="pyICESHEET reconstruction to compare against (optional)")
    ap.add_argument("--margin",
                    default="/home/awickert/models/icesheet/Greenland_final/"
                            "shear_stress/qgis/outline5.shp")
    ap.add_argument("--icesheet", default="/home/awickert/models/icesheet/icesheet")
    ap.add_argument("--spacing", type=float, default=40000.0)
    ap.add_argument("--interval", type=float, default=400.0)
    ap.add_argument("--simplify", type=float, default=15000.0)
    ap.add_argument("--workdir", default="fortran_ref")
    args = ap.parse_args()

    import geopandas as gpd

    workdir = Path(args.workdir)
    workdir.mkdir(exist_ok=True)

    d = np.load(args.fields)
    xb, yb = d["xb"].astype(float), d["yb"].astype(float)
    bed, tau = d["bed"], d["tau"]          # (ny, nx), row 0 = north (yb descending)
    assert yb[0] > yb[-1], "expected yb descending (north first)"
    spacing = int(round(abs(xb[1] - xb[0])))
    xmin, xmax = int(round(xb[0])), int(round(xb[-1]))
    ymin, ymax = int(round(yb[-1])), int(round(yb[0]))

    write_grid(workdir / "bed.bin", bed)
    write_grid(workdir / "ss.bin", tau)
    write_grid_params(workdir / "elev_parameters.txt", "bed.bin",
                      xmin, xmax, ymin, ymax, spacing)
    write_grid_params(workdir / "ss_parameters.txt", "ss.bin",
                      xmin, xmax, ymin, ymax, spacing)

    margin = gpd.read_file(args.margin).geometry.iloc[0].simplify(args.simplify)
    coords = np.asarray(margin.exterior.coords)[:, :2]
    np.savetxt(workdir / "outline.xyz", coords, fmt="%.3f")

    (workdir / "params.txt").write_text(
        f"outline.xyz\nelev_parameters.txt\nss_parameters.txt\n"
        f"{args.interval}\n{args.spacing}\n"
    )

    print("Running the Fortran ICESHEET on Greenland...")
    r = subprocess.run([args.icesheet], cwd=workdir, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:])
        raise SystemExit(f"Fortran icesheet failed (rc={r.returncode})")

    fx, fy, fE = parse_contours(workdir / "contours.txt")
    print(f"Fortran: {len(fE)} points, max surface {fE.max():.0f} m")

    if args.recon and Path(args.recon).exists():
        from scipy.interpolate import griddata
        rec = np.load(args.recon)
        # Compare at Fortran points: interpolate pyICESHEET surface there.
        py_at = griddata((rec["x"], rec["y"]), rec["elevation"], (fx, fy),
                         method="linear")
        good = np.isfinite(py_at)
        resid = fE[good] - py_at[good]
        print(f"pyICESHEET: {len(rec['x'])} points, max surface "
              f"{rec['elevation'].max():.0f} m")
        print(f"Fortran - pyICESHEET (n={good.sum()} common-support points):")
        print(f"  mean {resid.mean():+.0f} m   RMS {np.sqrt((resid**2).mean()):.0f} m"
              f"   median|.| {np.median(np.abs(resid)):.0f} m")


if __name__ == "__main__":
    main()
