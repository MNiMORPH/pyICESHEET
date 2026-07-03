"""Physical constants for pyICESHEET.

These are defined once, here, and imported everywhere. This is a deliberate
improvement over the original Fortran ICESHEET, which defined the ice density in
``global_parameters.f90`` as 917 kg m^-3 but then hardcoded 920 kg m^-3 inside
its Runge-Kutta routine. pyICESHEET uses a single, documented value.

All values are SI. Densities in kg m^-3, acceleration in m s^-2.
"""

from __future__ import annotations

from dataclasses import dataclass

#: Density of glacier ice (kg m^-3). Matches the primary value in the original
#: ICESHEET ``global_parameters.f90``; see module docstring re: the 917/920 fix.
RHO_ICE: float = 917.0

#: Density of seawater (kg m^-3), used for the marine (below-sea-level) margin
#: flotation approximation.
RHO_WATER: float = 1025.0

#: Standard gravitational acceleration (m s^-2).
G: float = 9.80665


@dataclass(frozen=True)
class PhysicalConstants:
    """Bundle of physical constants, so a run can override them explicitly.

    Defaults reproduce :data:`RHO_ICE`, :data:`RHO_WATER`, and :data:`G`. Pass a
    customized instance through the model configuration to vary them without
    mutating module globals.
    """

    rho_ice: float = RHO_ICE
    rho_water: float = RHO_WATER
    g: float = G

    @property
    def rho_g(self) -> float:
        """Convenience product ``rho_ice * g`` (Pa m^-1), the driving-stress scale."""
        return self.rho_ice * self.g


#: Module-level default instance.
DEFAULT_CONSTANTS = PhysicalConstants()
