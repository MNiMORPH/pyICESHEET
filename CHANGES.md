# Changelog

All notable changes to pyICESHEET are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/), and the
project aims to follow [Semantic Versioning](https://semver.org/).

## Authorship and provenance

**pyICESHEET is a Python port and modernization of ICESHEET, whose author and the
originator of the reconstruction method and algorithm is Evan J. Gowan.** This
project reproduces and improves upon his work; it is not new science. Please cite
the original ICESHEET papers:

- Gowan, E.J., Tregoning, P., Purcell, A., Lea, J., Fransner, O.J., Noormets, R.,
  Dowdeswell, J.A. (2016). ICESHEET 1.0: a program to produce paleo-ice sheet
  reconstructions with minimal assumptions. *Geoscientific Model Development*,
  9(5), 1673–1682. https://doi.org/10.5194/gmd-9-1673-2016
- Gowan, E.J., et al. (2021). A new global ice sheet reconstruction for the past
  80000 years. *Nature Communications*, 12, 1199.
  https://doi.org/10.1038/s41467-021-21469-w

Original ICESHEET (Fortran) is GPL-3.0; pyICESHEET inherits GPL-3.0.

## [Unreleased]

### Added
- Initial repository scaffold, design notes (physics/numerics/I/O of the Fortran;
  the pyICESHEET architecture and flowline-crowding plan).

### Notes on intended improvements over the Fortran
These are deliberate changes from the original, documented as they land:
- A single, documented ice density (ρ_i = 917 kg m⁻³) throughout — the Fortran
  mixes 917 and a hardcoded 920 inside its RK routine.
- Correct gravitational acceleration in the contour-resampling (Nye) distance
  estimate — the Fortran drops `g` there.
- Flowline-crowding, divide formation, and dome/saddle splitting handled via
  GEOS geometry operations (Shapely) instead of the original hand-rolled
  motorcycle-graph + polygon-splitting bookkeeping. The physics (RK integration
  of the Fisher slope equations) is preserved faithfully.
