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
- Initial repository scaffold and design notes (physics/numerics/I/O of the
  Fortran; the pyICESHEET architecture and flowline-crowding plan).
- `physics`: perfectly-plastic mechanics — stress/slope relations, the flowline
  ODE right-hand side (Fisher eq. A8), the Nye analytic profile/distance, and the
  marine-margin flotation seed. Validated against the analytic Nye parabola.
- `fields.RasterField`: smooth (bicubic-spline) sampling with analytic gradients.
- `margin.IceMargin`: boundary resampling and inward normals.
- `flowline.FlowlineIntegrator`: adaptive-RK (scipy) integration of a flowline to
  the next elevation contour, with the turning-singularity handled by an event.
- `contour.ContourManager`: the inward-advance step, with converging-flowline
  crossings and dome/saddle splitting delegated to GEOS `make_valid` (replacing
  the Fortran motorcycle-graph + polygon-split machinery).
- `solver.IceSheetModel`: seed the margin, march contours inward (recursing into
  split lobes), accumulate the surface.
- `surface.IceSurface`: sample cloud + gridding / GeoDataFrame export.
- `io`: NetCDF/GeoTIFF raster adapters and vector rasterization.
- `bmi.IceSheetBMI`: equilibrium (single-shot) CSDMS-style adapter.
- GRASS addon `r.icesheet`.
- Greenland example (bed = BedMachine v6; bundled margin + shear-stress polygons).

### Known limitations
- **Performance.** The solve integrates one adaptive-RK flowline per contour point
  per elevation interval; a full Greenland reconstruction currently takes several
  minutes. Loosening integrator tolerances, vectorizing field sampling, and
  simplifying the (very detailed) margin are the obvious next optimizations.
- Summit maxima differ between the two codes on Greenland (Fortran ~2800 m vs
  pyICESHEET ~3600 m, which reached its elevation cap); under investigation —
  likely divide placement (GEOS vs the original motorcycle-graph) and/or
  integrator tolerance.

### Validation
- Analytic Nye parabola (unit tests), a radial dome, and a dumbbell divide split.
- **Fortran ICESHEET method-to-method comparison.** The Fortran binary-grid format
  is reproduced and verified on the analytic circular cap (Fortran matches Nye to
  RMS 2.7 m, radially symmetric to 0.1 m). On Greenland, on identical fields,
  Fortran − pyICESHEET is mean +65 m / RMS 259 m / median 152 m — agreement to
  ~6–10 % on a 2–3 km ice sheet.
- Greenland reconstruction vs the observed BedMachine surface: mean +165 m,
  RMS 373 m (the positive bias is the expected perfectly-plastic signature).

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
