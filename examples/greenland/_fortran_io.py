"""I/O for the original Fortran ICESHEET, in the exact format it expects.

VALIDATED on the analytic circular cap (see fortran_cap_test.py): the Fortran
reads a native binary grid by skipping a 896-byte header, then reading 4-byte
floats via direct access, indexed from the top-left (north-west) corner,
row-major. We write that layout directly (rather than depend on GMT's ``=bf``,
whose header is 892 bytes in GMT 6 — a 4-byte mismatch).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

HEADER_BYTES = 896


def write_grid(path, values_north_first):
    """Write a Fortran-ICESHEET native grid.

    ``values_north_first`` is a 2-D array with row 0 = north (y = ymax), col 0 =
    west (x = xmin): 896 zero header bytes + native float32, row-major.
    """
    with open(path, "wb") as fh:
        fh.write(b"\x00" * HEADER_BYTES)
        fh.write(np.ascontiguousarray(values_north_first, dtype="<f4").tobytes())


def write_grid_params(path, binfile, xmin, xmax, ymin, ymax, spacing):
    """Write an elevation/shear-stress parameter file."""
    Path(path).write_text(
        f"{binfile}\n{int(xmin)}\n{int(xmax)}\n{int(ymin)}\n{int(ymax)}\n{int(spacing)}\n"
    )


def parse_contours(path):
    """Parse a Fortran ``contours.txt`` into ``(x, y, surface_elevation)`` arrays."""
    xs, ys, Es = [], [], []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        p = line.split()
        xs.append(float(p[0])); ys.append(float(p[1])); Es.append(float(p[2]))
    return np.array(xs), np.array(ys), np.array(Es)
