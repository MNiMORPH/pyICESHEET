#!/usr/bin/env python3

##############################################################################
# MODULE:    r.icesheet
# AUTHOR(S): pyICESHEET developers; method after Evan J. Gowan (ICESHEET)
# PURPOSE:   Reconstruct an equilibrium perfectly-plastic ice-sheet surface
#            from a basal shear-stress field and a mapped ice margin.
# COPYRIGHT: (C) 2026 pyICESHEET developers. GPL-3.0-or-later.
##############################################################################

# %module
# % description: Reconstructs an equilibrium perfectly-plastic ice-sheet surface from a basal shear-stress raster and a mapped ice margin (pyICESHEET).
# % keyword: raster
# % keyword: glaciology
# % keyword: ice sheet
# % keyword: reconstruction
# %end
# %option G_OPT_R_INPUT
# % key: bed
# % description: Bed topography raster [m]
# %end
# %option G_OPT_R_INPUT
# % key: shear_stress
# % description: Basal shear stress raster [Pa] (should be pre-smoothed, e.g. with r.resamp.bspline)
# %end
# %option G_OPT_V_INPUT
# % key: margin
# % description: Ice margin (area vector)
# %end
# %option G_OPT_R_OUTPUT
# % key: surface
# % description: Output ice-surface elevation raster [m]
# %end
# %option G_OPT_R_OUTPUT
# % key: thickness
# % required: no
# % description: Output ice-thickness raster [m]
# %end
# %option
# % key: spacing
# % type: double
# % required: no
# % answer: 5000
# % description: Along-contour point spacing / flowline seed spacing [m]
# %end
# %option
# % key: interval
# % type: double
# % required: no
# % answer: 100
# % description: Elevation interval between contours [m]
# %end
# %option
# % key: max_elevation
# % type: double
# % required: no
# % answer: 5000
# % description: Maximum surface elevation to march to [m]
# %end
# %option
# % key: min_thickness
# % type: double
# % required: no
# % answer: 1
# % description: Nominal ice thickness at the margin [m]
# %end

import sys

import numpy as np
import grass.script as gs


def _region_coords(region):
    """Cell-centre coordinate arrays for the current region (x asc, y desc)."""
    cols, rows = int(region["cols"]), int(region["rows"])
    w, n = float(region["w"]), float(region["n"])
    ewres, nsres = float(region["ewres"]), float(region["nsres"])
    x = w + ewres / 2.0 + np.arange(cols) * ewres
    y = n - nsres / 2.0 - np.arange(rows) * nsres
    return x, y


def _read_margin_polygon(margin):
    """Return a shapely (Multi)Polygon for the margin area vector."""
    from shapely import wkt
    from shapely.geometry import MultiPolygon, Polygon
    from shapely.ops import unary_union

    txt = gs.read_command("v.out.ascii", input=margin, format="wkt",
                          type="area", quiet=True)
    geoms = []
    for line in txt.splitlines():
        line = line.strip()
        if line.upper().startswith(("POLYGON", "MULTIPOLYGON")):
            g = wkt.loads(line)
            geoms.append(g)
    if not geoms:
        gs.fatal(f"no area geometry found in vector <{margin}>")
    merged = unary_union(geoms)
    if isinstance(merged, Polygon):
        return merged
    if isinstance(merged, MultiPolygon):
        return merged
    # GeometryCollection: keep polygons only
    polys = [g for g in getattr(merged, "geoms", []) if isinstance(g, Polygon)]
    return MultiPolygon(polys) if len(polys) > 1 else polys[0]


def main():
    options, flags = gs.parser()

    try:
        import grass.script.array as garray
        from scipy.interpolate import griddata
        from pyicesheet import RasterField, IceSheetModel, ModelConfig
    except ImportError as exc:
        gs.fatal(
            "pyICESHEET (and its dependencies) must be installed in the Python "
            f"used by GRASS. Import failed: {exc}"
        )

    region = gs.region()
    x, y = _region_coords(region)

    # Read input rasters as NumPy arrays aligned to the region.
    bed_arr = np.asarray(garray.array(mapname=options["bed"]), dtype=float)
    tau_arr = np.asarray(garray.array(mapname=options["shear_stress"]), dtype=float)
    if np.isnan(bed_arr).any():
        gs.warning("bed raster contains NULLs; filling with the median")
        bed_arr = np.where(np.isnan(bed_arr), np.nanmedian(bed_arr), bed_arr)
    if np.isnan(tau_arr).any():
        gs.warning("shear-stress raster contains NULLs; filling with the median")
        tau_arr = np.where(np.isnan(tau_arr), np.nanmedian(tau_arr), tau_arr)

    bed = RasterField.from_arrays(x, y, bed_arr)
    tau = RasterField.from_arrays(x, y, tau_arr)
    margin = _read_margin_polygon(options["margin"])

    cfg = ModelConfig(
        spacing=float(options["spacing"]),
        elevation_interval=float(options["interval"]),
        max_elevation=float(options["max_elevation"]),
        min_thickness=float(options["min_thickness"]),
    )
    gs.message("Running pyICESHEET reconstruction...")
    surf = IceSheetModel(bed, tau, margin, cfg).solve(
        progress=lambda n, s, t: gs.percent(min(int(t), int(cfg.max_elevation)),
                                            int(cfg.max_elevation), 5)
    )
    gs.message(f"Reconstruction: {len(surf)} samples, "
               f"max surface {surf.elevation.max():.0f} m")

    # Grid the sample cloud onto the region and write outputs.
    XX, YY = np.meshgrid(x, y)
    def _write(mapname, values):
        grid = griddata((surf.x, surf.y), values, (XX, YY), method="linear")
        out = garray.array()
        out[...] = grid.astype(np.float64)
        out.write(mapname=mapname, overwrite=gs.overwrite())

    _write(options["surface"], surf.elevation)
    if options["thickness"]:
        _write(options["thickness"], surf.thickness)

    gs.raster_history(options["surface"], overwrite=True)
    gs.message("Done.")


if __name__ == "__main__":
    sys.exit(main())
