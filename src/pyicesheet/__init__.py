"""pyICESHEET: equilibrium perfectly-plastic ice-sheet surface reconstruction.

A Pythonic, GIS-enabled port of Evan J. Gowan's ICESHEET. Given a mapped ice
margin, a basal shear-stress field, and bed topography, pyICESHEET reconstructs
the steady-state ice-surface elevation from the perfectly-plastic stress balance.

See ``docs/`` for the design notes describing the physics and architecture.
"""

from . import constants
from . import physics

__all__ = ["constants", "physics", "__version__"]

__version__ = "0.0.1.dev0"
