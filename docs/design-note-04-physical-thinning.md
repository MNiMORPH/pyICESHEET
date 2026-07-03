# Design Note 04 — physically-based thinning and stopping near the summit

**Status:** implemented (2026-07-03)
**Motivation:** two coupled problems, both traced to the contour machinery, not the
physics: (1) pyICESHEET's summit height depended on the input spacing (finer →
higher) and over-built vs reality; (2) it was too slow to run at fine resolution.
Both come from the same place — how the solver handles flowlines converging in the
flat interior — so both are fixed by one physical idea.

## The length scale

As flowlines converge toward an ice divide the surface flattens, and gaining one
contour interval requires a horizontal distance of the **Nye length**

    L = interval · H / Hf,     Hf = τ_b / (ρ_i g),

where `H` is the local ice thickness. `L` is large where the ice is thick and flat
(the interior), small where it is thin and steep (the margins). This single,
locally-computed scale drives both fixes.

## Stopping (fixes over-building and resolution dependence)

The original stopping was a **resolution-scaled area cutoff**, `min_area =
4·spacing²`: finer spacing → smaller cutoff → contours survive smaller and march
farther inward → higher summit. That is why the summit was not resolution-
independent (see Design Note 03).

Replace it with a physical test: a contour at elevation `E` enclosing area `A` has
equivalent radius `r = √(A/π)`. If `r < climb_factor · L`, the footprint is too
small to gain another interval — it *is* the summit — so it is not advanced
further. This is resolution-independent (it depends on `interval`, `H`, `τ`, not on
`spacing`).

`ContourManager.can_advance()` implements it; the solver accumulates every contour
but only keeps marching those that can climb. `min_area` is now only a numerical
floor against degenerate polygons.

## Thinning (fixes performance)

Where `L` is large the surface is smooth over that scale, so the contour does not
need fine point spacing there. Each contour is resampled at

    spacing_eff = clamp(spacing_growth · L, spacing, spacing_cap_factor · spacing),

so the flat interior is sampled coarsely (few flowlines to integrate) while the
steep margins keep the full base `spacing`. `ContourManager._effective_spacing()`.

## Why this is better than the original ICESHEET

Gowan's ICESHEET thins with a fixed `distance_factor · minimum_spacing` and stops
via the geometric `check_polygon2` rejection of over-extended contours. Here the
thinning *and* stopping length scale is derived from the **local stress balance**,
so it adapts to the shear stress and thickness automatically, and the stopping is
tied to a physical statement ("can this footprint support another interval?")
rather than a geometric heuristic.

## Results

Radial cap (flat bed, constant τ; analytic Nye apex 2109 m):

| spacing | old apex | new apex |
|---|---|---|
| 20 km | 2000 | 2000 |
| 10 km | (higher) | 2000 |
| 5 km | (higher, slow) | 2000 |

The new apex is **identical across resolutions** (it undershoots the analytic apex
by ~one interval — the cost of stopping when a full interval can no longer be
climbed; reduce `climb_factor` to march closer to the tip).

Greenland (same tuned τ field):

| spacing | summit | area > 2500 m | solve time |
|---|---|---|---|
| 40 km | 3750 m | 885k km² | 110 s |
| 20 km | 3600 m | 845k km² | 346 s |
| observed | 3232 m | 581k km² | — |

The summit now **decreases toward reality** with finer resolution (3750 → 3600),
where before it diverged upward, and the 20 km run is ~2–3× faster than before.

## What remains

The Greenland reconstruction still over-builds vs observed (3600 vs 3232;
845k vs 581k km² above 2500 m). This residual is **localized** — a spurious high
dome in east Greenland over the coastal mountains, where thin ice over high bed
gives a small `L`, so the stopping does not bite. That is a **divide-placement /
convergence** issue in that region (how the fronts from the fjorded east coast
merge), separate from the resolution-independence fixed here. It is the next thread
for the reality goal, along with pushing the resolution finer (now tractable).
