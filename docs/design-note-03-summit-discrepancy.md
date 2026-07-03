# Design Note 03 — the Greenland summit discrepancy

**Status:** investigation result (2026-07-03)
**Question:** on identical Greenland fields, the Fortran ICESHEET reaches a summit
of ~2800 m while pyICESHEET reaches ~3600 m (capped). Where does the difference
come from? A specific hypothesis was raised: that it relates to how the basal
shear stress is balanced through the flowline integral, and thus how the gridded
stress translates into accumulated surface elevation.

## Result in one line

The stress→elevation **integral is not the cause** (ruled out directly); the codes
agree on the ice geometry up to ~2000 m and then diverge because the Fortran's
**interior contours stop advancing** (their area cliffs at 2000–2500 m) while
pyICESHEET's march smoothly on to the divide. This is the divide-placement /
flowline-pruning behaviour flagged in Design Note 02 §2.

## Evidence

**1. The integral is tolerance- and interpolation-invariant (hypothesis ruled out).**
Integrating a single straight radial flowline through the real Greenland fields:

| τ sampling / tolerance | E @ 300 km | E @ 700 km |
|---|---|---|
| bicubic spline, rtol 1e-4 | 4172 | 5067 |
| bicubic spline, rtol 1e-9 | 4172 | 5067 |
| bilinear, rtol 1e-9 | 4170 | 5066 |

Tightening the integrator tolerance by five orders of magnitude changes the result
by 0 m; switching the τ interpolation from bicubic spline to bilinear changes it by
~1 m (mean τ along the path differs by 0.01 %). So neither the integration nor the
gridding-to-elevation translation is responsible.

**2. Hypsometry: the codes agree up to ~2000 m, then the Fortran cliffs.**
Ice area (10³ km²) above each elevation, vs the observed BedMachine surface:

| elevation | observed | pyICESHEET | Fortran |
|---|---|---|---|
| > 2000 m | 1049 | 1241 | 1229 |
| > 2500 m | 582 | 878 | **95** |
| > 3000 m | 109 | 412 | **0** |
| max | 3233 | 3600 | 2800 |

Up to 2000 m the two reconstructions are nearly identical. Above ~2500 m the
Fortran's contour area collapses (a cliff), while pyICESHEET's shrinks gradually
toward the divide. See `examples/greenland/summit_discrepancy.png`.

**3. Against reality, pyICESHEET is the closer of the two — but both miss.**
The Fortran badly *under*-builds the interior (95k vs 582k km² above 2500 m; summit
2800 vs 3233 m). pyICESHEET *over*-builds (878k; 3600 m, capped) but is much closer
to observed. Reality lies between them.

## Interpretation (leading hypothesis, not yet proven)

The divergence is in the **interior contour geometry near the ice divide**, not the
per-flowline physics. As flowlines converge toward the divide, the original
ICESHEET prunes them aggressively — the motorcycle-graph crossover rule (shorter
flowline wins) and the `check_polygon2` rejection of over-extended contours — which
shrinks the contour rapidly and halts it (the cliff). pyICESHEET instead resolves
converging fronts with GEOS `make_valid`, which preserves more contour area, so its
contours march farther inward and build a higher summit.

This is exactly the GEOS-vs-motorcycle-graph divide-placement difference flagged in
Design Note 02 §2 as a validation item.

**To confirm it:** implement the configurable distance-rule ("shorter flowline
wins") survivor selection in `ContourManager` (currently
`survivor_rule="geos"` only) and check whether it collapses pyICESHEET's summit
toward the Fortran's. If it does, the mechanism is confirmed and becomes a user
choice: GEOS (smoother, closer to observed) vs distance-rule (reproduces Gowan).

## What this is NOT

- Not the integrator tolerance (invariant to 1e-9).
- Not the τ interpolation / gridding (spline ≈ bilinear to 0.01 %).
- Not the ρ_i 917-vs-920 constant (a ~0.3 % effect, ≪ the 800 m seen).
- Not a uniform bias — the medians agree (1200 m both); only the deep interior
  (long integration paths past the divide region) diverges.
