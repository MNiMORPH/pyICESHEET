"""Generate a reference reconstruction from the ORIGINAL Fortran ICESHEET.

Status: scaffold for a supervised run. The intended validation is method-to-
method — run the original Fortran ICESHEET and pyICESHEET on the *same* bed,
margin, and shear-stress fields and compare the surfaces within a tolerance.

This script writes the inputs the Fortran expects (GMT native-binary grids via
``gmt xyz2grd ... =bf``, the ``*_parameters.txt`` files, an ``outline.xyz``
margin, and ``params.txt``) and invokes the compiled ``icesheet`` binary.

CAVEATS to verify before trusting the output (why this is not yet wired into the
test suite): the Fortran's grid I/O has subtleties that must be checked against a
known case first —
  * gridline vs pixel registration of the ``=bf`` grid;
  * row order / byte order of GMT native binary vs the Fortran's top-left-origin
    record indexing (a mismatch silently flips the grid);
  * the Fortran under-allocates its in-memory grid store by one row/column, so
    the domain must extend a few cells beyond the margin (buffer) or edge
    flowlines read out of bounds.
Validate on the flat-bed / constant-tau circular cap (analytic Nye answer known)
before running Greenland, so any format error is caught where the answer is
known.

Usage (from a working directory):
    python make_fortran_reference.py --fields greenland_fields.npz \\
        --icesheet /path/to/icesheet --spacing 40000 --interval 250
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import numpy as np


def write_bf_grid(name, x, y, values, workdir):
    """Write a GMT native-float (=bf) grid from a regular field via xyz2grd."""
    xmin, xmax = int(round(x.min())), int(round(x.max()))
    ymin, ymax = int(round(y.min())), int(round(y.max()))
    res = int(round(abs(x[1] - x[0])))

    XX, YY = np.meshgrid(x, y)
    xyz = np.column_stack([XX.ravel(), YY.ravel(), np.asarray(values).ravel()])
    xyz_path = workdir / f"{name}.xyz"
    np.savetxt(xyz_path, xyz, fmt="%.6g")

    bin_path = workdir / f"{name}.bin"
    subprocess.run(
        ["gmt", "xyz2grd", str(xyz_path), f"-G{bin_path}=bf",
         f"-I{res}", f"-R{xmin}/{xmax}/{ymin}/{ymax}"],
        check=True,
    )
    # elev/ss parameter file: filename, xmin, xmax, ymin, ymax, spacing.
    params = workdir / f"{name}_parameters.txt"
    params.write_text(f"{bin_path.name}\n{xmin}\n{xmax}\n{ymin}\n{ymax}\n{res}\n")
    return params, (xmin, xmax, ymin, ymax, res)


def write_margin_xyz(margin, path):
    """Write the margin exterior as a plain 'x y' (metres) file."""
    coords = np.asarray(margin.exterior.coords)
    np.savetxt(path, coords[:, :2], fmt="%.3f")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fields", default="greenland_fields.npz",
                    help="npz with xb, yb, bed, tau (integer coords, metres)")
    ap.add_argument("--margin",
                    default="/home/awickert/models/icesheet/Greenland_final/"
                            "shear_stress/qgis/outline5.shp")
    ap.add_argument("--icesheet", default="/home/awickert/models/icesheet/icesheet")
    ap.add_argument("--spacing", type=float, default=40000.0)
    ap.add_argument("--interval", type=float, default=250.0)
    ap.add_argument("--workdir", default="fortran_ref")
    args = ap.parse_args()

    import geopandas as gpd

    workdir = Path(args.workdir)
    workdir.mkdir(exist_ok=True)

    d = np.load(args.fields)
    xb, yb = d["xb"].astype(float), d["yb"].astype(float)
    # y must ascend for xyz2grd; reorder bed/tau to match.
    order = np.argsort(yb)
    yb = yb[order]
    bed = d["bed"][order, :]
    tau = d["tau"][order, :]

    elev_par, _ = write_bf_grid("elev", xb, yb, bed, workdir)
    ss_par, _ = write_bf_grid("ss", xb, yb, tau, workdir)

    margin = gpd.read_file(args.margin).geometry.iloc[0]
    margin_path = workdir / "outline.xyz"
    write_margin_xyz(margin, margin_path)

    (workdir / "params.txt").write_text(
        f"{margin_path.name}\n{elev_par.name}\n{ss_par.name}\n"
        f"{args.interval}\n{args.spacing}\n"
    )

    print(f"Wrote Fortran inputs to {workdir}/. Running icesheet...")
    subprocess.run([args.icesheet], cwd=workdir, check=True)
    print(f"Done. See {workdir}/contours.txt")


if __name__ == "__main__":
    main()
