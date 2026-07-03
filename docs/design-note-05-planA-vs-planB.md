# Design Note 05 — Plan A (flowline) vs Plan B (grid)

**Status:** Plan B prototype validated (2026-07-03)

Two ways to solve the same perfectly-plastic stress balance
`τ_b = ρ_i g H |∇S|`:

- **Plan A — Lagrangian flowlines** (`pyicesheet.solver`, the faithful port of
  Gowan's ICESHEET). Trace flowlines inward from the margin, advancing contours.
  Its structural weakness is *front over-enclosure at convergences*: where
  flowlines from different directions meet, GEOS `make_valid` takes the outer
  envelope, enclosing more area than the true divide (Design Notes 03–04). It is
  also slow (minutes for Greenland).
- **Plan B — Eulerian grid** (`pyicesheet.grid`). Solve `|∇S| = τ/(ρ_i g (S−B))`
  directly on the raster as a state-dependent eikonal, fixed at the margin and
  marched inward. `S` is single-valued, so **divides form by construction** — no
  convergence to track, no over-enclosure to clip.

## Head-to-head on Greenland (identical bed, τ, margin)

| | summit | area > 2500 km² | area > 3000 km² | RMS vs observed | time |
|---|---|---|---|---|---|
| observed (BedMachine) | ~3.2–3.5 km | 582k | 109k | — | — |
| **Plan B (grid)** | 3424 m | **712k** | **213k** | **177 m** | **1.9 s** |
| Plan A (flowline) | 3600 m | 821k | 411k | 373 m | ~10 min |
| original Fortran (fine) | 3300 m | 714k | 221k | — | ~min |

Plan B's hypsometry matches the original Fortran (712k vs 714k above 2500 m) —
**the over-enclosure is gone** (Plan A enclosed 821k) — it is **~2× closer to
observed** (RMS 177 vs 373 m) and **~300× faster** (1.9 s vs ~10 min). Figure:
`examples/greenland/planB_grid_greenland.png`.

## Why Plan B works (validated before Greenland)

- Analytic Nye cap (point divide): RMS 15 m (~0.7 % of a 2 km cap).
- Ellipse (line divide — exactly where Plan A over-encloses): matches the
  distance-transform ground truth, hypsometry within ~4 %, no over-enclosure.

Numerics: fast-iterative eikonal; each cell solves the *coupled* local update
(the speed `F = τ/(ρg(S−B))` depends on the unknown `S`) — closed-form 1-D
Godunov plus a Newton step for 2-D cells. The near-margin boundary is seeded with
the plastic toe `S ≈ B + √(2τd/ρg)` from each cell's sub-cell distance to the
margin, removing a half-cell bias.

## Trade-offs and status

- Plan B is a **different discretization**, not the faithful flowline port: it is
  first-order (a slight, resolution-shrinking summit undershoot ~2 %), and it uses
  the local force balance on the grid. The residual +120 m mean over observed is
  the shared τ/model over-build (the fine Fortran shows the same).
- Plan A remains the reference for reproducing Gowan exactly, and the front-clip
  (Task #8) is still worth having for it.
- **Recommendation:** Plan B is the better path for the reality goal and the GRASS
  target — correct at divides by construction, ~300× faster, raster-native, and a
  fraction of the code. Next: a variable-bed validation, a second-order eikonal
  update to remove the undershoot, and wiring Plan B into the GRASS addon.
