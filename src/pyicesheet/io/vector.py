"""Vector input: ice margins and shear-stress polygons.

Ease of use with vector data is a design goal: the margin is naturally a polygon,
and basal shear stress is naturally a set of polygons each carrying a value.
Vector shear stress is rasterized onto the model grid (a fixed resolution) here,
mirroring what a GRASS ``v.to.rast`` step would do. ``geopandas``/``rasterio``
are imported lazily.
"""

from __future__ import annotations

import numpy as np

from .raster import grid_transform

__all__ = ["read_polygon", "read_polygons", "rasterize_polygons"]


def read_polygon(path, index=0):
    """Read a single polygon geometry from a vector file (e.g. the margin)."""
    import geopandas as gpd
    gdf = gpd.read_file(path)
    return gdf.geometry.iloc[index]


def read_polygons(path):
    """Read a vector file, returning the GeoDataFrame (attributes + geometry)."""
    import geopandas as gpd
    return gpd.read_file(path)


def rasterize_polygons(geometries, values, x, y, fill=np.nan, dtype="float64"):
    """Burn polygon ``values`` onto the regular grid defined by ``x``, ``y``.

    Parameters
    ----------
    geometries : iterable of shapely geometries
    values : iterable of float
        Value to burn for each geometry (e.g. shear stress in Pa).
    x, y : 1-D arrays
        Grid coordinates (regular spacing).
    fill : float
        Value for cells not covered by any polygon.

    Returns
    -------
    grid : 2-D ndarray, shape ``(len(y), len(x))``, top-down (north first).
    """
    from rasterio.features import rasterize

    transform, _, _ = grid_transform(x, y)
    ny, nx = len(y), len(x)
    shapes = [(g, float(v)) for g, v in zip(geometries, values)]
    grid = rasterize(shapes, out_shape=(ny, nx), transform=transform,
                     fill=fill, dtype=dtype)
    return grid
