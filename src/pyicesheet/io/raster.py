"""Raster input/output.

Readers for gridded fields (bed topography, and optionally a pre-rasterized
shear-stress field), plus a GeoTIFF writer for results. Heavy geospatial
dependencies (``netCDF4``, ``rasterio``) are imported lazily so the numerical
core does not require them.
"""

from __future__ import annotations

import numpy as np

from ..fields import RasterField

__all__ = [
    "read_netcdf_downsampled",
    "fill_invalid",
    "grid_transform",
    "write_geotiff",
]


def read_netcdf_downsampled(path, var, factor=1, x_name="x", y_name="y",
                            fill_below=-9990.0):
    """Read a NetCDF variable, strided by ``factor`` for downsampling.

    Returns ``(x, y, values)`` with 1-D coordinate arrays and a 2-D array shaped
    ``(len(y), len(x))``. Values ``<= fill_below`` are set to NaN (BedMachine
    uses -9999 as a fill). Reading a coarse stride keeps a multi-GB file
    tractable and is itself a mild anti-alias smoothing.
    """
    import netCDF4

    ds = netCDF4.Dataset(path)
    try:
        x = np.asarray(ds[x_name][::factor], dtype=float)
        y = np.asarray(ds[y_name][::factor], dtype=float)
        v = np.asarray(ds[var][::factor, ::factor], dtype=float)
    finally:
        ds.close()
    if fill_below is not None:
        v = np.where(v <= fill_below, np.nan, v)
    return x, y, v


def fill_invalid(values, method="median"):
    """Fill NaNs in a 2-D array.

    ``method="median"`` fills with the array median (cheap; fine for fill cells
    that lie outside the ice). ``method="nearest"`` fills each NaN from its
    nearest valid neighbour (better near the ice boundary).
    """
    values = np.array(values, dtype=float)
    nan = np.isnan(values)
    if not nan.any():
        return values
    if method == "median":
        values[nan] = np.nanmedian(values)
    elif method == "nearest":
        from scipy.ndimage import distance_transform_edt
        idx = distance_transform_edt(nan, return_distances=False,
                                     return_indices=True)
        values = values[tuple(idx)]
    else:
        raise ValueError(f"unknown fill method {method!r}")
    return values


def grid_transform(x, y):
    """Affine transform (top-left origin) for a regular grid given 1-D coords."""
    from rasterio.transform import from_origin
    resx = abs(float(x[1] - x[0]))
    resy = abs(float(y[1] - y[0]))
    west = min(float(x[0]), float(x[-1]))
    north = max(float(y[0]), float(y[-1]))
    return from_origin(west - resx / 2, north + resy / 2, resx, resy), resx, resy


def write_geotiff(path, grid, x, y, crs=None, nodata=np.nan):
    """Write a 2-D array to a GeoTIFF. ``grid`` must be top-down (north first)."""
    import rasterio
    transform, _, _ = grid_transform(x, y)
    grid = np.asarray(grid, dtype="float32")
    with rasterio.open(
        path, "w", driver="GTiff", height=grid.shape[0], width=grid.shape[1],
        count=1, dtype="float32", crs=crs, transform=transform, nodata=nodata,
    ) as dst:
        dst.write(grid, 1)


def field_from_netcdf(path, var, factor=1, fill="median", **kw):
    """Convenience: read a NetCDF variable and return a :class:`RasterField`."""
    x, y, v = read_netcdf_downsampled(path, var, factor=factor, **kw)
    v = fill_invalid(v, method=fill)
    return RasterField.from_arrays(x, y, v)
