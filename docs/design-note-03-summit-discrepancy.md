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

## Interpretation — confirmed by experiment

The divergence is in the **interior contour geometry near the ice divide**, not the
per-flowline physics. Two candidate Fortran pruning mechanisms were tested:

**(a) The motorcycle-graph crossover rule ("shorter flowline wins") — RULED OUT.**
Implemented as `ContourManager(survivor_rule="distance")` and run on Greenland: it
made essentially no difference (area above 2500 m: 877k km² with GEOS vs 883k with
the distance-rule; both summits capped at 3600 m). On a clean radial cap it does
lower the apex slightly (2000 vs 2200 m), but it is not what caps the Fortran on
Greenland.

**(b) The `check_polygon2` rejection of over-extended contours — CONFIRMED.**
The Fortran writes discarded contours to `contours-rejected.txt`. On this Greenland
run it rejects a growing number of contours with elevation — 20 at 2000 m, 23 at
2400 m, **35 at 2800 m, 8 at 3200 m**. The rejection peaks exactly where the
Fortran's accepted surface caps (2800 m): the Fortran *tries* to build to 3200 m
but rejects those interior contours. pyICESHEET has no such rejection (only a
minimum-area filter), so it keeps building to 3600 m.

So the cliff is the Fortran's **polygon rejection**, not the divide-placement
crossover rule. (Andy's original stress-integral hypothesis and my own
motorcycle-graph hypothesis were both wrong; the evidence points to rejection.)

## Correction — the coarse comparison was resolution-confounded

**The analysis above was run at deliberately coarse settings (40 km spacing,
400 m interval, 15 km-simplified margin) for speed, and that confounded the
conclusion.** Two facts, established afterward:

1. The Greenland shear-stress field is **tuned so that Gowan's ICESHEET reproduces
   the observed surface** (per E. Gowan / A. Wickert). So the field is calibrated
   *against the Fortran*, and the Fortran run — done correctly — should match
   reality, not fall short of it.
2. Re-running the **Fortran at fine settings (5 km spacing, 100 m interval)**
   gives a summit of **3300 m vs observed 3232 m** — it *does* reproduce reality.
   The coarse run's 2800 m was a resolution artifact (a simplified margin and
   400 m steps produce more folded contours, so `check_polygon2` rejects far more
   of the interior). The rejection *cliff* is real but its severity is
   resolution-dependent.

Hypsometry, area (10³ km²) above elevation:

| elevation | observed | Fortran fine (5 km) | Fortran coarse (40 km) | pyICESHEET coarse |
|---|---|---|---|---|
| > 2500 m | 581 | 714 | 96 | 878 |
| > 3000 m | 110 | 221 | 0 | 411 |
| max | 3232 | 3300 | 2800 | 3600 |

So the earlier claim that "the Fortran under-builds reality and pyICESHEET is
closer" is **retracted**: at proper resolution the Fortran ~matches reality. Both
codes' coarse results merely bracketed reality by artifact.

Two things remain genuinely open:

- The fine Fortran still *over*-builds the interior somewhat (714k vs 581k km²
  above 2500 m). This is likely because our reproduction of the τ field
  (rasterize + Gaussian smooth) does not exactly match Gowan's original grid
  pipeline (`nearest_int`/`reduce_dem`) that the values were tuned with — i.e. we
  are not running on *precisely* the tuned field.
- **Does pyICESHEET converge to reality at fine resolution?** Untested: pyICESHEET
  is too slow to run at 5 km spacing today. This makes solver performance the
  critical path for the reality goal, not the contour-rejection algorithm.

## Note on the model physics

The governing balance `τ_b = ρ_i g H |∇S|` is a **force balance** — the left side
is the gravitational driving stress — and is independent of the ice deformation
mechanism (internal creep, basal sliding, till). It is *not* a perfectly-plastic
rheology; ICESHEET with a spatially variable τ_b is a **prescribed-basal-shear-
stress** model. Its one approximation is that all driving stress is balanced by
basal drag (the shallow-ice / local balance, neglecting longitudinal and lateral
stress-gradient terms), which is least accurate at divides — but a τ_b tuned to
reality absorbs that. Earlier references in this project to a "perfectly-plastic
vs. real ice" gap were mis-stated on this point.

The `survivor_rule="distance"` option remains available as a faithful reproduction
of Gowan's crossover handling; GEOS is the default.

## What this is NOT

- Not the integrator tolerance (invariant to 1e-9).
- Not the τ interpolation / gridding (spline ≈ bilinear to 0.01 %).
- Not the ρ_i 917-vs-920 constant (a ~0.3 % effect, ≪ the 800 m seen).
- Not a uniform bias — the medians agree (1200 m both); only the deep interior
  (long integration paths past the divide region) diverges.
