"""Sanity-check a pyICESHEET Greenland reconstruction vs the observed surface.

A perfectly-plastic reconstruction is NOT expected to match the modern observed
ice surface exactly (real ice is not perfectly plastic, and the shear-stress
field is an assumption). This is a *sanity* reference: the reconstruction should
be Greenland-scale ice with a broadly comparable surface. The rigorous validation
is the method-to-method comparison against the original Fortran ICESHEET
(see make_fortran_reference.py).

Usage:
    python compare_observed.py --recon greenland_recon.npz \\
        --bed /path/to/BedMachineGreenland-v6.nc [--figure out.png]
"""

from __future__ import annotations

import argparse

import numpy as np

from pyicesheet import RasterField
from pyicesheet.io.raster import read_netcdf_downsampled, fill_invalid


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recon", default="greenland_recon.npz")
    ap.add_argument("--bed", default="/home/awickert/Downloads/BedMachineGreenland-v6.nc")
    ap.add_argument("--factor", type=int, default=33)
    ap.add_argument("--figure", default=None)
    args = ap.parse_args()

    r = np.load(args.recon)
    rx, ry, rE = r["x"], r["y"], r["elevation"]

    xs, ys, obs = read_netcdf_downsampled(args.bed, "surface", factor=args.factor)
    obs_f = RasterField.from_arrays(xs, ys, fill_invalid(obs, "nearest"))
    obs_at = np.asarray(obs_f.value(rx, ry), dtype=float)

    valid = obs_at > 1.0
    resid = rE[valid] - obs_at[valid]
    print(f"reconstruction: {len(rE)} samples, max surface {rE.max():.0f} m")
    print(f"observed max surface: {obs_at.max():.0f} m")
    print(f"residual (recon - observed) over ice, n={valid.sum()}:")
    print(f"  mean {resid.mean():+.0f} m   RMS {np.sqrt((resid**2).mean()):.0f} m")
    print(f"  median|.| {np.median(np.abs(resid)):.0f} m   "
          f"p90|.| {np.percentile(np.abs(resid), 90):.0f} m")

    if args.figure:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 6))
        sc = ax[0].scatter(rx, ry, c=rE, s=4, cmap="viridis")
        ax[0].set_title("pyICESHEET reconstructed surface [m]")
        ax[0].set_aspect("equal"); plt.colorbar(sc, ax=ax[0], shrink=0.6)
        sc2 = ax[1].scatter(rx[valid], ry[valid], c=resid, s=4, cmap="RdBu_r",
                            vmin=-800, vmax=800)
        ax[1].set_title("recon - observed [m]")
        ax[1].set_aspect("equal"); plt.colorbar(sc2, ax=ax[1], shrink=0.6)
        fig.tight_layout()
        fig.savefig(args.figure, dpi=120)
        print(f"wrote {args.figure}")


if __name__ == "__main__":
    main()
