# Design Note 02 — pyICESHEET architecture and the flowline-crowding plan

**Status:** agreed design (working prototype name: *pyICESHEET*)
**Companion:** builds on
[`design-note-01-understanding-icesheet.md`](design-note-01-understanding-icesheet.md),
which specifies the physics, numerics, and I/O of the original Fortran ICESHEET.
This note records *what we are building and why*, with the decisions made
2026-07-03.

pyICESHEET is a Pythonic, object-oriented port of Evan Gowan's ICESHEET, intended
to (a) reproduce the analytical-style plastic ice-sheet reconstruction, (b) be
cleaner and easier to reason about than the Fortran, and (c) drop smoothly into
GRASS GIS. It is always an **equilibrium** (steady-state) solution — there is no
time stepping.

---

## 1. Decisions (with rationale)

| Decision | Choice | Why |
|---|---|---|
| **Solver method** | Faithful **Lagrangian flowline** port (not an Eulerian eikonal reformulation) | Reproduces Evan's method and the exact Fisher cross-flow physics; the eikonal was attractive but is a *different* method. Kept as a possible future backend, not the prototype. |
| **Contour/topology handling** | **Contour-as-front + GEOS** (see §2) | Collapses ~2500 lines of hand-rolled crossover/split bookkeeping into a few robust geometry-library calls. |
| **Field smoothing** | **Separate upstream operator**, not in the solver | One-thing-one-task. The solver *assumes* it is handed a smooth, sampleable τ (and bed) field. In GRASS this is `r.resamp.bspline`/`r.mapcalc`; in the library a standalone smoothing step. |
| **Distribution** | **Standalone pip-installable library + thin GRASS addon** (`r.icesheet`) | Physics testable/reusable without GRASS; GRASS is one front-end. Matches "CSDMS outside-looking-in." |
| **CSDMS/BMI** | **External adapter** (`bmi.py`), core stays plain Python | BMI wraps the model from outside; `update()` runs the *entire* equilibrium solve; time functions degenerate (documented). |
| **First milestone** | **1-D / radial Nye first**, validated against the analytic parabola | Cheap executable anchor that pins the physics before any 2-D/topology machinery. |
| **License** | **GPL-3.0** | ICESHEET is GPL-3.0; a faithful derived port must remain GPL-3.0. Compatible with GRASS (GPL-2-or-later). |
| **Authorship** | Evan Gowan credited as originator of the method/algorithm | We are porting and improving, not inventing. Recorded in `CHANGES.md` and `README`. |
| **Validation target** | Greenland — *basically the same as Gowan, not bit-for-bit* | A documented tolerance, not exact reproduction (§4). |

---

## 2. The flowline-crowding plan (the heart of the port)

Nearly all of ICESHEET's bulk exists to manage flowlines converging as the front
marches inward. Three things happen:

1. **Flowlines cross** — physically an **ice divide / ridge**, i.e. a *shock*
   where the PDE's characteristics collide. To be located, not merely suppressed.
2. **The contour self-intersects and pinches off lobes** — domes and saddles
   separating.
3. **Points crowd** without crossing — accuracy and speed degrade.

The Fortran handles these with a motorcycle-graph crossover rule ("when two
flowlines cross, kill the one with the shorter run to the crossing"), explicit
self-intersection detection + polygon splitting + index bookkeeping, and
one-directional oversampling.

**The clean reframing.** ICESHEET is already *contour-centric*: it advances the
whole front one elevation at a time, and flowlines are per-point local
integrations. Leaning into that, GEOS owns the topology:

- **Crossing + lobe-splitting collapse into one GEOS step.** After advancing
  every point to the next contour, build the ring and run
  `make_valid` / `unary_union` (`buffer(0)`). Self-intersections resolve; separate
  lobes fall out as a `MultiPolygon` automatically — **domes pinch off for free**.
  Recurse on each resulting polygon.
- **Crowding is arc-length redistribution.** Reparameterize each valid contour to
  uniform `minimum_spacing` every step (densify where sparse, **decimate where
  crowded** — the Fortran only densifies).

**Why this is principled, not a shortcut.** The motorcycle graph is literally an
algorithm for straight skeletons / medial axes, and the medial axis of the margin
under the τ-weighted metric *is* the locus of ice divides. So the shock
interpretation is preserved; GEOS simply produces the front geometry.

**Honesty flag / validation item.** GEOS-union places the divide at the
*geometric overlap boundary*; the motorcycle rule places it by
*shorter-path-wins*. These agree for symmetric convergence but can differ when
two flowlines reach a crossing having travelled unequal distances. For
"reproduce Greenland, not bit-for-bit" this is very likely fine — but it is a
**validation item, not an assumed equivalence**.

**Decision:** GEOS-union survivor selection is the **default**; the explicit
distance-rule is retained as a **configurable strategy** if Greenland validation
demands it. Either way, **the physics — RK integration of the Fisher ODEs — stays
faithful to Evan**; the elegance is confined to the contour-management layer.

---

## 3. Architecture

The crowding concern is isolated behind a single component, so the rest of the
solver never sees the geometry mess:

```python
class ContourManager:              # a.k.a. FrontTracker  (contour.py)
    def advance(self, front: Contour, fields) -> list[Contour]:
        # 1. integrate each point up its flowline  (-> FlowlineIntegrator)
        # 2. GEOS make_valid / unary_union         -> resolve crossings, split lobes
        # 3. arc-length redistribute               -> control crowding
        # returns one or more contours to recurse on
```

Repository layout:

```
pyICESHEET/
  pyproject.toml                 # src layout; import name: pyicesheet
  LICENSE                        # GPL-3.0
  CHANGES.md                     # authorship (Evan = originator), port + improvements
  README.md
  docs/
    design-note-01-understanding-icesheet.md
    design-note-02-architecture.md            # this file
  src/pyicesheet/
    constants.py     # densities (one chosen rho_i, documented), g
    config.py        # ModelConfig dataclass + from_yaml
    physics.py       # driving stress, p/q slope relations, Fisher q', Nye profile, flotation — pure, tested
    fields.py        # RasterField: sample + gradients  (assumes pre-smoothed input)
    margin.py        # IceMargin: polygon(s), resample, inward normals
    flowline.py      # FlowlineIntegrator: adaptive RK, one flowline -> next contour
    contour.py       # ContourManager: the crowding seam (GEOS topology + redistribution)
    solver.py        # IceSheetModel: inward march; .solve() -> IceSurface
    surface.py       # IceSurface: contour cloud -> raster (elev/thickness) + contour vectors
    io/
      raster.py      # rasterio <-> numpy + affine
      vector.py      # geopandas; rasterize vector tau at fixed resolution
    bmi.py           # IceSheetBMI (bmipy) — external adapter; update() = full solve
  grass/
    r.icesheet/      # thin addon: GRASS raster/vector <-> numpy/shapely -> IceSheetModel
  tests/
    test_physics_nye.py            # M1: analytic parabola
    test_flotation.py
    test_flowline_single.py        # M2
    test_front_convergence.py      # crowding/divide in isolation (dumbbell margin)
    test_greenland_regression.py   # M3: tolerance-based vs Gowan
  examples/greenland/
```

**Vector-data ergonomics.** Margin and τ can be supplied as vector data; `io/vector.py`
rasterizes τ to a fixed model resolution on ingest (the GRASS analog is `v.to.rast`).
The solver works on gridded fields internally regardless of input type.

---

## 4. Milestones and validation

- **M1 — physics + Nye.** `physics.py` and the simplest solver on a flat-bed,
  constant-τ case, checked against the analytic Nye parabola `H² = (2τ/ρg)·x`.
  Pins the physics before any topology.
- **M2 — single flowline.** `FlowlineIntegrator` on real fields; one flowline /
  1-D march, adaptive RK matching the Fortran's behaviour.
- **M3 — full 2-D.** `ContourManager` + `solver.py` on Greenland. Precede the
  full run with `test_front_convergence.py`: exercise the crowding/divide logic
  **in isolation** on a controlled two-lobe (dumbbell) margin, so the one part we
  are wary of is anchored before the full case.

Validation is **tolerance-based** against Gowan's Greenland reconstruction, not
bit-for-bit. The physics discrepancies noted in Design Note 01 §6 (ρ_i 917 vs 920;
missing `g` in a resampling estimate; stale "km" comment) are resolved
deliberately in the port and documented, which is one reason exact reproduction
is not the target.

---

## 5. Authorship and license

- **Evan J. Gowan** is the author of ICESHEET and the originator of the method and
  algorithm reproduced here. pyICESHEET is a **port and improvement**, not new
  science. This is stated in `README.md` and recorded in `CHANGES.md`.
- Cite the ICESHEET papers (Gowan et al. 2016, GMD; Gowan et al. 2021, Nat.
  Commun.) as in the original `readme`.
- **License: GPL-3.0**, inherited from ICESHEET.

---

*End of Design Note 02.*
