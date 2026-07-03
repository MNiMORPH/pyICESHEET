"""pyICESHEET: equilibrium perfectly-plastic ice-sheet surface reconstruction.

A Pythonic, GIS-enabled port of Evan J. Gowan's ICESHEET. Given a mapped ice
margin, a basal shear-stress field, and bed topography, pyICESHEET reconstructs
the steady-state ice-surface elevation from the perfectly-plastic stress balance.

See ``docs/`` for the design notes describing the physics and architecture.
"""

from . import constants
from . import physics
from .config import ModelConfig
from .fields import RasterField
from .margin import IceMargin
from .flowline import FlowlineIntegrator
from .contour import Contour, ContourManager
from .surface import IceSurface
from .solver import IceSheetModel
from .grid import GridSurface, solve_surface_grid

__all__ = [
    "constants",
    "physics",
    "ModelConfig",
    "RasterField",
    "IceMargin",
    "FlowlineIntegrator",
    "Contour",
    "ContourManager",
    "IceSurface",
    "IceSheetModel",
    "GridSurface",
    "solve_surface_grid",
    "__version__",
]

__version__ = "0.0.1.dev0"
