"""Validate the Fortran ICESHEET binary-grid format on the analytic circular cap.

Before trusting a Fortran-vs-pyICESHEET comparison on Greenland (where there is no
analytic answer), we validate the whole Fortran pipeline — the binary grid I/O,
the parameter files, the margin file — on a case whose answer is known exactly:
a flat-bed, constant-shear-stress circular ice cap, whose surface is the Nye
parabola ``H(r) = sqrt(2 tau (R - r) / (rho g))``.

The Fortran reads a headered native binary (it skips 896 bytes, then reads
4-byte floats via direct access, indexed from the top-left / north-west corner,
row-major). GMT 6 writes an 892-byte header for ``=bf`` — a 4-byte mismatch — so
we write the raw binary ourselves in exactly the layout the Fortran expects,
rather than depend on GMT's (version-dependent) header.

Usage:
    python fortran_cap_test.py --icesheet /path/to/icesheet [--workdir cap_ref]
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import numpy as np

from pyicesheet import RasterField, IceSheetModel, ModelConfig, physics
from pyicesheet.constants import DEFAULT_CONSTANTS as C
from _fortran_io import write_grid, write_grid_params, parse_contours

# Case definition (integer coords/spacing, metres).
RADIUS = 300_000.0
TAU = 1.0e5
XMIN, XMAX = -400_000, 400_000
YMIN, YMAX = -400_000, 400_000
SPACING = 5_000
INTERVAL = 100.0
POINT_SPACING = 20_000.0


def build_grids(workdir):
    ncols = (XMAX - XMIN) // SPACING + 1
    nrows = (YMAX - YMIN) // SPACING + 1
    # Row 0 = north (y=YMAX); col 0 = west (x=XMIN).
    x = XMIN + np.arange(ncols) * SPACING
    y = YMAX - np.arange(nrows) * SPACING          # descending (north first)
    bed = np.zeros((nrows, ncols), dtype=float)
    tau = np.full((nrows, ncols), TAU, dtype=float)

    write_grid(workdir / "bed.bin", bed)
    write_grid(workdir / "ss.bin", tau)
    write_grid_params(workdir / "elev_parameters.txt", "bed.bin",
                      XMIN, XMAX, YMIN, YMAX, SPACING)
    write_grid_params(workdir / "ss_parameters.txt", "ss.bin",
                      XMIN, XMAX, YMIN, YMAX, SPACING)
    return x, y, bed, tau


def write_margin(workdir):
    t = np.linspace(0, 2 * np.pi, 400, endpoint=True)
    xy = np.column_stack([RADIUS * np.cos(t), RADIUS * np.sin(t)])
    np.savetxt(workdir / "outline.xyz", xy, fmt="%.3f")


def nye(r):
    return physics.nye_thickness(np.maximum(RADIUS - r, 0.0), TAU, C, thickness0=1.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--icesheet", default="/home/awickert/models/icesheet/icesheet")
    ap.add_argument("--workdir", default="cap_ref")
    args = ap.parse_args()

    workdir = Path(args.workdir)
    workdir.mkdir(exist_ok=True)

    x, y, bed, tau = build_grids(workdir)
    write_margin(workdir)
    (workdir / "params.txt").write_text(
        f"outline.xyz\nelev_parameters.txt\nss_parameters.txt\n{INTERVAL}\n{POINT_SPACING}\n"
    )

    print("Running Fortran ICESHEET on the circular cap...")
    subprocess.run([args.icesheet], cwd=workdir, check=True,
                   stdout=subprocess.DEVNULL)

    fx, fy, fE = parse_contours(workdir / "contours.txt")
    fr = np.hypot(fx, fy)
    fnye = nye(fr)
    interior = fr < RADIUS - 2 * SPACING
    fresid = fE[interior] - fnye[interior]
    print(f"\nFORTRAN: {len(fE)} points, apex {fE.max():.0f} m (analytic apex {nye(0):.0f} m)")
    print(f"  vs analytic Nye: mean {fresid.mean():+.1f}  RMS {np.sqrt((fresid**2).mean()):.1f} m")
    # Radial symmetry / angular spread: within each radial bin, detrend by Nye(r)
    # first (so we measure angular variation, not the parabola's steepness across
    # the bin). A flipped/shifted/garbled grid would break symmetry -> large spread.
    bins = np.linspace(0, RADIUS - 2 * SPACING, 25)
    idx = np.digitize(fr, bins)
    detrended = fE - fnye
    spread = np.nanmax([detrended[idx == b].std()
                        for b in range(1, len(bins)) if np.sum(idx == b) > 3])
    print(f"  angular spread (Nye-detrended, max within-bin std): {spread:.1f} m")

    # pyICESHEET on the same case.
    bedF = RasterField.from_arrays(x, y, bed)
    tauF = RasterField.from_arrays(x, y, tau)
    from shapely.geometry import Point
    margin = Point(0, 0).buffer(RADIUS, quad_segs=100)
    cfg = ModelConfig(spacing=POINT_SPACING, elevation_interval=INTERVAL,
                      max_elevation=3000.0, rtol=1e-6, atol=1e-3)
    surf = IceSheetModel(bedF, tauF, margin, cfg).solve()
    pr = np.hypot(surf.x, surf.y)
    pint = pr < RADIUS - 2 * SPACING
    presid = surf.elevation[pint] - nye(pr[pint])
    print(f"\npyICESHEET: {len(surf)} points, apex {surf.elevation.max():.0f} m")
    print(f"  vs analytic Nye: mean {presid.mean():+.1f}  RMS {np.sqrt((presid**2).mean()):.1f} m")

    ok = np.sqrt((fresid ** 2).mean()) < 30 and spread < 30
    print("\nVERDICT: Fortran binary format is",
          "VALIDATED — the cap follows the Nye parabola and is radially symmetric"
          if ok else "SUSPECT — check byte order / row order / header offset")


if __name__ == "__main__":
    main()
