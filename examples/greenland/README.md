# Greenland example

Reconstruct the Greenland ice sheet from a mapped margin and a basal shear-stress
field, and compare against Evan Gowan's original Fortran ICESHEET on the same
inputs.

## Data provenance

| Input | File | Source | CRS |
|---|---|---|---|
| Bed topography | `BedMachineGreenland-v6.nc` (variable `bed`) | NSIDC IDBMG4 v6 (Morlighem et al.); **not redistributed** — download separately | EPSG:3413, 150 m |
| Ice margin | `outline5.shp` | ICESHEET repo (`Greenland_final/shear_stress/qgis/`) | EPSG:3413 |
| Basal shear stress | `shear_stress.shp` (attribute `shear_stre`, Pa) | ICESHEET repo, same folder | EPSG:3413 |

All three share the NSIDC polar-stereographic CRS (EPSG:3413), so the model runs
in projected metres with no reprojection. The bed DEM is large (2.8 GB) and is
**not** included in this repository; set `--bed` to your local copy.

The observed modern ice `surface` and `thickness` in the BedMachine file are used
only as a *secondary sanity reference* — a perfectly-plastic reconstruction is not
expected to match the real surface exactly.

## Running

```bash
# pyICESHEET reconstruction
python run_greenland.py --bed /path/to/BedMachineGreenland-v6.nc \
    --resolution-factor 33 --spacing 15000 --interval 200

# Reference from the original Fortran ICESHEET on identical fields, then compare
python make_fortran_reference.py --bed /path/to/BedMachineGreenland-v6.nc
python compare_to_fortran.py
```

## Preprocessing notes

- **Shear-stress smoothing is an explicit upstream step.** The polygons are
  rasterized (piecewise-constant, sharp edges) and then Gaussian-smoothed before
  the solve, because the flowline integrator differentiates the field. This
  mirrors the original ICESHEET's reliance on bicubic interpolation of a coarse
  grid, and keeps smoothing out of the solver (one operator, one task).
- **Resolution.** `--resolution-factor 33` downsamples the 150 m bed to ~5 km.
  Coarser factors run faster; finer ones resolve more of the fjord margins.

## Validation

The intended validation is *method-to-method*: run the original Fortran ICESHEET
and pyICESHEET on the **same** bed, margin, and shear-stress fields, and compare
the reconstructed surfaces within a tolerance (not bit-for-bit — see
`docs/design-note-02` §1 and §6 for the deliberate physics corrections that make
exact agreement neither expected nor desired).
