# Design Note 06 — pyICESHEET as one module in a reconstruction system

**Status:** plan / roadmap (2026-07-03). The forward solver, the tau-from-substrate
builder, and the scalar + per-class calibration layer are built; the rest is
scoped here as the direction of travel.

pyICESHEET's job is narrow and fast: an **equilibrium, perfectly-plastic forward
operator** `(bed, tau, margin) -> ice surface`, in ~0.1 s for Greenland (grid
solver, Design Note 05). That speed is what makes it useful as the *inner* engine
of a larger, data-constrained ice-sheet-reconstruction workflow — a modernized,
GRASS-native cousin of the ICESHEET → PaleoMIST program. This note lays out that
larger picture and the module boundaries.

## Three layers

```
   observations                          parameters
 (surface, RSL, drainage, dated margins)  (alpha; tau-per-substrate; regional ...)
        |                                        |
        v                                        v
  +-----------------------------+     +-----------------------------+
  | CALIBRATION / INVERSE       |<--->| CONDITIONING (upstream ops) |
  |  pluggable objectives       |     |  smoothing                  |
  |  + parameters, scipy/ensembl|     |  tau-from-substrate geology |
  +-----------------------------+     +-----------------------------+
        |                                        |
        |            +-------------------+       |
        +----------->|  FORWARD SOLVER   |<------+
                     |  grid eikonal, fast|
                     |  (bed,tau,margin)->S|
                     +-------------------+
```

**1. Forward solver** (built, `pyicesheet.grid`). Pure, fast, side-effect-free.
Everything else treats it as a black box called many times.

**2. Conditioning operators** (upstream, composable — "one thing, one task").
These *build the inputs* the solver assumes are clean:
- field smoothing (built, kept external by design);
- **tau from substrate geology** (built, `pyicesheet.tau.tau_from_classes`): a
  classified substrate raster + a per-class yield-stress table → tau. The classes
  come from mapped geology (Gowan et al. 2019, ESSD 11, 375: sediment cover,
  grain size, bedrock type) — deformable sediment → low tau → ice streams, bedrock
  → high tau. This replaces hand-drawn shear-stress polygons with a physically
  grounded, few-parameter field.

**3. Calibration / inverse layer** (scalar + vector built, `pyicesheet.calibrate`).
The caller gives a `forward(params) -> surface` closure (build tau, solve) and an
*objective* (`surface -> misfit`). Objectives compose (weighted sum) so several
data constraints fit simultaneously. Fitting is `scipy.optimize` today; ensembles
/ MCMC later (the fast solver makes thousands of evaluations cheap — a 13-class
Greenland fit was 756 solves in 63 s).

## Data constraints (objectives)

| constraint | status | notes |
|---|---|---|
| modern surface (RMS) | built (`SurfaceMisfit`) | 177 → 137 (scalar) → 125 m (13-class) on Greenland |
| relative sea level / GIA | external module | ice load → Earth deformation → RSL; the PaleoMIST loop; a substantial coupled model (SELEN-class), not in-package |
| drainage routing | near (via GRASS) | proglacial spillways, ice-dammed lakes, meltwater routes constrain margin/surface; `r.watershed` / `r.lake` / `r.fill.dir` |
| dated ice-margin chronology | future | drives a *time series* of equilibrium reconstructions |

## Parameters

- **global scalar** `alpha` (built): `tau -> alpha*tau`, `H ∝ √alpha`. Corrects the
  multiplicative bias; cannot fix spatial pattern.
- **per-substrate-class table** (built): one yield stress per geology class. Attacks
  the *pattern* (puts low tau where sediment/ice-streams are). Physically
  constrained by the geology map, so far fewer effective DOF than free per-cell tau.
- **regional / marine-vs-land multipliers, effective-pressure dependence** (future).

## Module boundaries — what pyICESHEET owns vs. couples to

- **In package:** the forward solver, the conditioning operators (smoothing, tau
  builder), the calibration layer with a surface objective, the GRASS addon.
- **Coupled, external:** the GIA/RSL model (its own repo/tool), the drainage
  routing (GRASS modules), the geology datasets (data, downloaded), and any
  time-stepping driver that sequences equilibrium states through a deglaciation.

pyICESHEET stays a clean, fast, testable equilibrium engine; the "reconstruction
system" is the *composition* of it with those external pieces, orchestrated in
GRASS / a Python driver.

## Roadmap

1. **Now (done):** fast grid solver; tau-from-classes; scalar + per-class
   calibration against the modern surface.
2. **Next:** ingest the real Gowan et al. 2019 geology rasters (not the proxy
   polygon classes) and calibrate the per-class table; a second-order eikonal
   update for cleaner verification; wire the grid solver into `r.icesheet`.
3. **Then:** drainage-routing objective via GRASS; a `forward` that builds tau from
   geology inside the loop.
4. **Later:** GIA/RSL coupling (the paleo inverse), time-dependent reconstruction
   from dated margins, and ensemble / Bayesian calibration over the parameter set.
