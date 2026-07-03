# pyICESHEET

A Pythonic, GIS-enabled port of **ICESHEET** — a program to reconstruct the
equilibrium surface of a perfectly-plastic ice sheet from a **basal shear-stress
field** and a **mapped ice margin**, with no assumptions about climate or mass
balance.

pyICESHEET always produces a **steady-state (equilibrium)** solution: given the
ice margin, the basal shear stress, and the bed topography, the ice-surface
elevation is fully determined by the stress balance
`τ_b = ρ_i g H |∇S|` (driving stress equals basal shear stress).

> **Status: pre-alpha, under active development.** APIs will change.

## Provenance and credit

**ICESHEET, and the reconstruction method and algorithm reproduced here, are the
work of Evan J. Gowan.** pyICESHEET is a port and modernization — not new
science. If you use it, please cite the original ICESHEET papers:

- Gowan, E.J., Tregoning, P., Purcell, A., Lea, J., Fransner, O.J., Noormets, R.,
  Dowdeswell, J.A. (2016). ICESHEET 1.0: a program to produce paleo-ice sheet
  reconstructions with minimal assumptions. *Geoscientific Model Development*,
  9(5), 1673–1682. https://doi.org/10.5194/gmd-9-1673-2016
- Gowan, E.J., et al. (2021). A new global ice sheet reconstruction for the past
  80000 years. *Nature Communications*, 12, 1199.
  https://doi.org/10.1038/s41467-021-21469-w

## Why a port?

The original ICESHEET is Fortran. A reading of the source (see
[`docs/design-note-01`](docs/design-note-01-understanding-icesheet.md)) shows that
only a small fraction of the code is physics; the bulk is computational geometry —
detecting where inward-marching flowlines cross and where contours pinch off ice
domes. Modern geometry libraries (Shapely/GEOS) and a clean object model replace
that machinery, and native raster/vector I/O makes the model usable directly
inside a GIS. The intended endpoint is a smooth [GRASS GIS](https://grass.osgeo.org/)
integration.

Design and rationale:
[`docs/design-note-02`](docs/design-note-02-architecture.md).

## How it works (one paragraph)

The mapped margin is the zero-thickness contour. Flowlines are integrated inward
from the margin, raising the surface according to the perfectly-plastic slope
relation until the next elevation contour is reached; the process recurses inward
until the ice sheet closes. Where flowlines converge, ice **divides** form; where
contours pinch, **domes and saddles** separate. pyICESHEET keeps Gowan's flowline
physics (the Fisher slope equations, integrated with adaptive Runge–Kutta) and
delegates the divide/dome topology to GEOS.

## Inputs

- **Ice margin** — a closed polygon (vector).
- **Basal shear stress** τ_b — vector polygons (each with a value) or a raster.
- **Bed topography** — a raster.

Vector inputs are rasterized to a fixed model resolution on ingest. Smoothing of
the shear-stress and bed fields is a **separate upstream step**, by design — the
solver assumes it is handed smooth, sampleable fields.

## Installation (development)

```bash
git clone https://github.com/MNiMORPH/pyICESHEET.git
cd pyICESHEET
pip install -e ".[gis,bmi,test]"
```

## License

GPL-3.0-only, inherited from the original ICESHEET.
