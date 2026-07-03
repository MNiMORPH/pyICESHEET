"""Input/output adapters for pyICESHEET.

Raster and vector readers/writers live here so the numerical core stays free of
geospatial dependencies. GIS libraries (rasterio, geopandas) are imported lazily
inside the adapter modules, so importing :mod:`pyicesheet` does not require them.
"""
