"""A CSDMS-style Basic Model Interface (BMI) adapter for pyICESHEET.

pyICESHEET is an *equilibrium* model: there is no time evolution. The BMI is
therefore an "outside-looking-in" wrapper (kept out of the numerical core) with
single-shot semantics: :meth:`initialize` builds the model, :meth:`update` runs
the *entire* reconstruction to its steady state, and the time functions
degenerate (start = end = current = 0). This is a common way to expose a steady
solver through BMI.

This is a BMI-*shaped* adapter (the lifecycle and the ``get_value`` grid access
that couplers need). Full ``bmipy.Bmi`` conformance — every ``get_grid_*`` /
``get_var_*`` method — is a later, mechanical step.

``initialize`` takes a configuration mapping (or a YAML file path) with:

    bed:     path to a bed raster (NetCDF var via {path, var, factor}) or an
             in-memory pyicesheet.RasterField
    tau:     likewise for basal shear stress
    margin:  path to a polygon vector file, or a shapely geometry
    output_resolution: grid cell size for get_value rasters (m)
    ... plus any ModelConfig fields (spacing, elevation_interval, ...)
"""

from __future__ import annotations

import numpy as np

from .config import ModelConfig
from .fields import RasterField
from .solver import IceSheetModel

__all__ = ["IceSheetBMI"]

_OUTPUT_VARS = ("land_ice_surface__elevation", "land_ice__thickness")
_VAR_TO_ATTR = {
    "land_ice_surface__elevation": "elevation",
    "land_ice__thickness": "thickness",
}
_VAR_UNITS = {
    "land_ice_surface__elevation": "m",
    "land_ice__thickness": "m",
}


class IceSheetBMI:
    """Equilibrium (single-shot) BMI wrapper around :class:`IceSheetModel`."""

    def __init__(self):
        self._model = None
        self._surface = None
        self._grids = {}
        self._grid_x = None
        self._grid_y = None
        self._resolution = None
        self._done = False

    # -- lifecycle ------------------------------------------------------- #

    def initialize(self, config):
        """Build the model from a config mapping or YAML file path."""
        cfg = self._load_config(config)
        bed = self._as_field(cfg["bed"])
        tau = self._as_field(cfg["tau"])
        margin = self._as_margin(cfg["margin"])
        self._resolution = float(cfg.get("output_resolution", cfg.get("spacing", 5000.0)))

        model_cfg = ModelConfig(**{
            k: cfg[k] for k in (
                "spacing", "elevation_interval", "min_thickness", "min_elevation",
                "max_elevation", "require_inside", "rtol", "atol",
            ) if k in cfg
        })
        self._model = IceSheetModel(bed, tau, margin, model_cfg)
        self._done = False

    def update(self):
        """Run the full equilibrium reconstruction (there is only one step)."""
        if self._model is None:
            raise RuntimeError("initialize() must be called before update()")
        self._surface = self._model.solve()
        self._grids = {}
        self._done = True

    def update_until(self, then):
        """Equilibrium model: any target time resolves to the single solve."""
        if not self._done:
            self.update()

    def finalize(self):
        self._model = None
        self._surface = None
        self._grids = {}

    # -- identity / vars ------------------------------------------------- #

    def get_component_name(self):
        return "pyICESHEET (equilibrium perfectly-plastic ice sheet)"

    def get_input_var_names(self):
        return ()

    def get_output_var_names(self):
        return _OUTPUT_VARS

    def get_var_units(self, name):
        return _VAR_UNITS[name]

    def get_var_type(self, name):
        return "float64"

    def get_var_grid(self, name):
        return 0

    # -- time (degenerate: equilibrium) ---------------------------------- #

    def get_start_time(self):
        return 0.0

    def get_end_time(self):
        return 0.0

    def get_current_time(self):
        return 0.0

    def get_time_step(self):
        return 0.0

    def get_time_units(self):
        return "1"

    # -- values / grid --------------------------------------------------- #

    def get_value(self, name):
        """Return the named output as a flattened gridded array."""
        return self._grid_for(name).ravel()

    def get_value_ptr(self, name):
        return self._grid_for(name)

    def get_grid_shape(self, grid=0):
        self._ensure_grid("land_ice__thickness")
        return self._grids["land_ice__thickness"].shape

    def get_grid_spacing(self, grid=0):
        return (self._resolution, self._resolution)

    def get_grid_type(self, grid=0):
        return "uniform_rectilinear"

    def get_grid_rank(self, grid=0):
        return 2

    # -- internals ------------------------------------------------------- #

    def _grid_for(self, name):
        self._ensure_grid(name)
        return self._grids[name]

    def _ensure_grid(self, name):
        if self._surface is None:
            raise RuntimeError("update() must be called before reading values")
        if name not in self._grids:
            attr = _VAR_TO_ATTR[name]
            grid, xc, yc = self._surface.to_grid(self._resolution, value=attr)
            self._grids[name] = np.asarray(grid, dtype=float)
            self._grid_x, self._grid_y = xc, yc
        return self._grids[name]

    @staticmethod
    def _load_config(config):
        if isinstance(config, dict):
            return config
        import yaml  # optional; only needed for file-based config
        with open(config) as fh:
            return yaml.safe_load(fh)

    @staticmethod
    def _as_field(spec):
        if isinstance(spec, RasterField):
            return spec
        from .io.raster import field_from_netcdf
        if isinstance(spec, dict):
            return field_from_netcdf(spec["path"], spec["var"],
                                     factor=spec.get("factor", 1))
        raise TypeError("field spec must be a RasterField or {path, var, factor}")

    @staticmethod
    def _as_margin(spec):
        from shapely.geometry.base import BaseGeometry
        if isinstance(spec, BaseGeometry):
            return spec
        from .io.vector import read_polygon
        return read_polygon(spec)
