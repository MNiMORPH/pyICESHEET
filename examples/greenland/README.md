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
    --resolution-factor 33 --spacing 40000 --interval 400

# Validate the Fortran binary-grid format on the analytic circular cap first
python fortran_cap_test.py --icesheet /path/to/icesheet

# Then the Greenland head-to-head: Fortran vs pyICESHEET on identical fields
python make_fortran_reference.py --fields greenland_fields.npz \
    --recon greenland_recon.npz --icesheet /path/to/icesheet
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

The validation is *method-to-method*: the original Fortran ICESHEET and
pyICESHEET run on the **same** bed, margin, and shear-stress fields, and their
surfaces are compared within a tolerance (not bit-for-bit — see
`docs/design-note-02` §1 and §6 for the deliberate physics corrections).

1. **Binary-grid format, on the analytic circular cap** (`fortran_cap_test.py`).
   The Fortran cap follows the Nye parabola to **RMS 2.7 m** with **0.1 m**
   angular spread (radially symmetric), confirming the grid I/O is correct.
   pyICESHEET matches the same analytic answer to RMS 11.5 m.
2. **Greenland head-to-head** (`make_fortran_reference.py`). On identical fields,
   Fortran − pyICESHEET is **mean +65 m, RMS 259 m, median 152 m** — the two
   independent implementations agree to ~6–10 % on a 2–3 km ice sheet.

Open discrepancy: the summit maxima differ (Fortran ~2800 m vs pyICESHEET ~3600 m,
which reached its elevation cap) — under investigation; likely divide-placement
(GEOS vs the original motorcycle-graph) and/or integrator tolerance.

The `_fortran_io.py` grid writer reproduces the Fortran's native format directly
(896-byte header + native float32, north-first row-major), which is more reliable
than GMT's `=bf` (an 892-byte header in GMT 6 — a 4-byte mismatch).
