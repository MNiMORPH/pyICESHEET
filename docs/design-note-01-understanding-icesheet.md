# Design Note 01 — Understanding ICESHEET, toward a GIS-enabled Python port

**Status:** draft / working document
**Purpose:** capture a faithful reading of Evan Gowan's ICESHEET (Fortran) — its
physics, numerics, algorithm, and I/O contract — as the reference specification
for an eventual Python implementation runnable inside GRASS GIS. This note
documents *what the existing code does*; it does not yet propose the port's
design.

All line references are to the current tree
(`find_flowline_fisher_adaptive_4.f90` unless otherwise noted).

---

## 1. What the model computes

ICESHEET reconstructs a **steady-state, perfectly-plastic ice-sheet surface**
from three spatial inputs:

1. a **mapped ice margin** — a closed polygon (projected, planar metres);
2. a **basal shear-stress field** τ_b(x, y) — originally vector polygons, each
   carrying a τ value, rasterized before the solve;
3. a **bed topography** grid B(x, y).

There is **no climate, mass balance, or time evolution.** Given the margin and
the stress/bed fields, the surface is fully determined by a stress balance. The
approach follows Reeh (1982) and Fisher et al. (1985); the code comments cite
"Fisher's equation A8" for the governing cross-flow slope equation.

Output is the ice **surface elevation** (and derived thickness) as a set of
nested elevation contours, later gridded.

---

## 2. Physical model

### 2.1 Governing assumption

Basal shear stress equals driving stress everywhere (the perfectly-plastic /
shallow-ice limit):

```
τ_b = ρ_i g H |∇S|
```

with H = S − B the ice thickness (surface S minus bed B). Hence the surface
slope magnitude is

```
|∇S| = τ_b / (ρ_i g H).
```

The code carries `Hf = τ_b / (ρ_i g)` (units of length; `grids.f90:270`,
`RK4:1690`), so `|∇S| = Hf / H`.

### 2.2 One-dimensional limit (Nye parabola)

With a flat bed and uniform τ, integrating `dS/dx = τ/(ρ_i g H)` gives the Nye
profile `H² = (2 τ / ρ_i g) x`. This exact relation is used for:

- the initial boundary-consistency check between adjacent margin points
  (`:160-167`), and
- the Nye distance estimates that set point spacing during resampling
  (`:1037-1038`, `:1090-1091`).

### 2.3 Marine (below-sea-level) margins

Where the bed at the margin is below sea level, the margin ice is set to a
flotation-balanced thickness rather than zero:

```
S_margin = B · (1 − ρ_w / ρ_i)       (icesheet init, :111-118)
```

i.e. a grounded-but-near-flotation approximation (a real ice shelf is not
modelled). A subsequent pass (`check_boundary`, :137-176) walks from the lowest
margin point upward and raises any neighbour whose implied slope would exceed
the Nye limit for the smaller of the two local shear stresses.

### 2.4 Constants

| Quantity | Value | Location |
|---|---|---|
| ρ_i | 917 kg m⁻³ | `global_parameters.f90:35` |
| ρ_i | **920 kg m⁻³** | `RK4:1681` (local override — inconsistent, see §6) |
| ρ_w | 1025 kg m⁻³ | `global_parameters.f90:36` |
| g | 9.80665 m s⁻² | both |
| max surface elevation | 5000 m | `global_parameters.f90:42` |
| min margin thickness | 1 m | `global_parameters.f90:44` |

---

## 3. Numerical method — flowline coordinate system

The 2-D surface is built by tracing **flowlines inward from the margin** and
integrating the surface up each one. This is the conceptual heart of the model.

### 3.1 Local frame

At each contour point a local frame is fixed to the current **inward-normal
direction** (an azimuth). The integration variable `x_l` advances along this
fixed axis; `y_l` is the lateral (cross-axis) offset. Two slope components are
tracked:

- `p` — along-axis surface slope,
- `q` — cross-axis surface slope,

with `|∇S|² = p² + q²`, so `p = √(|∇S|² − q²) = √((Hf/H)² − q²)`.

Because the true flowline follows steepest descent (the gradient), it curves
away from the fixed `x_l` axis; `y_l` records that drift and the local direction
is re-rotated by `atan2(q, p)` after each contour advance (`:419-429`).

### 3.2 The ODE system (integrated in `x_l`)

State vector `(y, E, q)` where `E` is surface elevation. From `y_prime`,
`E_prime`, `q_prime_fisher` (`:1964-2025`), with `H = E − B` and
`p = √((Hf/H)² − q²)`:

```
dy/dx_l = q / p
dE/dx_l = (Hf/H)² / p                                   = |∇S|² / p
dq/dx_l = Hf² (dB/dy − q) / (H³ p)  +  Hf (dHf/dy) / (p H²)   [Fisher A8]
```

Sanity check: when `q = 0`, `p = |∇S|`, `dy/dx_l = 0`, and
`dE/dx_l = |∇S|` — pure along-axis ascent, as expected. When `q ≠ 0` the
elevation gain `dE/dx_l = (p² + q²)/p = p + q²/p` correctly accounts for the
lateral drift `dy = (q/p) dx_l` via `dE = ∇S · displacement`.

The `q'` equation is what makes the reconstruction respond to **gradients** in
bed (`dB/dy`) and shear stress (`dHf/dy`): the second term is the ICESHEET
contribution that distinguishes it from a constant-τ reconstruction
(`q_prime` vs `q_prime_fisher`, :1988-2025).

Bed and stress derivatives (`dB/dy`, `dHf/dy`) and the field values themselves
come from **bicubic interpolation** of the input grids (`grids.f90`), using
Skidmore (1989) finite differences on the interpolant. Smooth derivatives are
essential — the integrator differentiates τ and B, so a non-smooth field breaks
it.

### 3.3 Integration control

`step_flow_line` (`:1336-1654`) drives an **adaptive 4th/5th-order Runge–Kutta**
with step doubling (`RK4:1661-1960`): full step vs. two half steps, error
estimated on `(y, E, q)` against fixed tolerances, step size grown/shrunk by the
usual safety-factor powers (`:1595-1620`). Special handling:

- the integrand has a square-root singularity as `q → |∇S|` (the flowline
  turning tangent); `q_prime_fisher` returns a sentinel `-999999` there
  (`:2021-2023`), and `RK4` solves a local quadratic to back off the step
  (`:1737-1774`, and twice more);
- when a step still fails below `stop_value`, the flowline is **rotated 90°**
  (`rotation_amount1 = −sign(q)·π/2`, :1461-1475) and restarted — a recovery
  for near-tangential flow.

Stopping: integrate until `E` reaches the next elevation contour
(`current_elevation = step × elevation_interval`).

---

## 4. Algorithm — contour-advancing recursion

`calculate_polygon5` (`:212-1330`) is recursive: each time a contour pinches into
separate lobes, the routine calls itself on each lobe. One invocation = advance
the whole current contour by one `elevation_interval`.

Per elevation step:

1. **Advance** every non-skipped contour point up its flowline to
   `current_elevation` (`boundary_search`, :356-433). Points already above the
   target, or that step outside the previous polygon, are held/flagged.
2. **Eliminate crossing flowlines** — a *motorcycle-graph* rule
   (`:459-581`): among all pairwise flowline crossings, repeatedly remove the
   flowline with the **shorter** run to the crossover point. Prevents the new
   contour from self-tangling where ice converges.
3. **Detect self-intersections** of the advanced contour
   (`contour_crossover`, :2163-2241) and **split** it into sub-polygons; insert
   the crossover points; assign each point a polygon index (:722-897).
4. **Reject** folded/invalid polygons (`check_polygon2`, :2506-2571): shrink each
   point slightly along its inward normal; if too many land outside the polygon,
   the normal directions are inconsistent and the polygon is discarded (written
   to `contours-rejected.txt`).
5. **Oversample** accepted contours (`oversample_loop`, :1027-1194): insert
   points so along-contour spacing stays below `distance_factor × minimum_spacing`
   *and* the Nye-estimated next-step distance stays bounded. Recompute inward
   normals for the resampled contour (:1203-1259).
6. **Recurse** to `current_step + 1`. Terminate when elevation exceeds
   `max_elevation` (5000 m) or a polygon drops below 3 points.

Most of the ~2650 lines are this geometric bookkeeping (crossover detection,
polygon splitting, dynamic array growth). The physics proper is ~60 lines
(§3.2–3.3). **This split is the single most important fact for the port:** the
hard, voluminous part is computational geometry that mature libraries
(Shapely/GEOS, NumPy) already provide.

---

## 5. I/O contract

### 5.1 `params.txt` (5 lines; read in `icesheet.f90:50-59`)

```
<ice margin file>
<elevation parameter file>
<shear stress parameter file>
<elevation_interval>     # metres between output contours
<minimum_spacing>        # metres; margin resample + flowline seed spacing
```

### 5.2 Elevation / shear-stress parameter files (`grids.f90:87-153, 173-243`)

```
<binary grid filename>
<xmin>          # integer, projected metres
<xmax>
<ymin>
<ymax>
<grid_spacing>  # integer metres
```

### 5.3 Grid binary format

GMT-style binary: **896-byte header** (`header_offset = 896/4` records,
`grids.f90:63`) followed by **4-byte floats**, row-major, ordered from the
**top-left** (record index uses `(ymax − y)`, `:343-345`). Values read via
bicubic interpolation; the whole grid is slurped into memory when
`store_dem = .true.` (default).

### 5.4 Ice margin file

Plain `x  y` (projected metres), one vertex per line, closed polygon. Read and
resampled to `minimum_spacing` in `read_icefile.f90`.

### 5.5 Outputs

- `contours.txt` — GMT multisegment (`> -Z <elev> <polygon#>` headers). Per
  point: `x, y, surface_elev, thickness, bed, skip_flag, shear_stress`
  (`:994-997`).
- `contours-rejected.txt` — rejected polygons, same format.
- `oversample_points.txt` — diagnostic.

Gridding of the contour point cloud into a raster surface is done **downstream**
by GMT scripts, not by ICESHEET itself.

### 5.6 Coordinates

Everything in the solver is **projected planar metres**; the core is Cartesian.
`r_earth` exists in `global_parameters` but is unused by the solver. Projection
choice lives entirely in the preprocessing (margin digitizing, grid creation).

---

## 6. Noted discrepancies (record, do not fix here)

1. **ρ_i inconsistency:** 917 kg m⁻³ globally vs. **920** hardcoded inside `RK4`
   (`:1681`). The two are used in the same solve. The port should choose one
   value deliberately and document it.
2. **Missing g in resampling Nye estimate:** the oversample distance
   (`:1037-1038`, `:1090-1091`) uses `rho_ice / (2·ss)` where the Nye scaling is
   `rho_ice·g / (2·τ)`. This affects only point *spacing*, not the surface
   physics — but it is a genuine unit discrepancy relative to the boundary check
   at `:162`, which does include `g`.
3. **Stale comment:** `icesheet.f90:57` says `minimum_spacing` is "in km"; the
   code uses it directly against metre coordinates. It is metres.

These are flagged for the port's test/validation stage, not as fixes to the
published Fortran.

---

## 7. Implications for a GRASS/Python port (preliminary)

Not a design yet — just what the reading implies. Open questions marked **[Q]**.

- **Clean separation to exploit:** physics (§3.2–3.3, ~60 lines) is small and
  well-defined; the bulk is geometry (§4) that Shapely/GEOS + NumPy replace
  wholesale. The port is mostly *re-expressing the contour-advance loop* on top
  of a geometry library, not re-deriving stress balance.
- **Field access:** ICESHEET's bespoke bicubic reader over GMT binaries becomes
  raster sampling with smooth interpolation. In GRASS this is `r.what` /
  library raster access, or a resampled NumPy array with a spline interpolant.
  The smoothness requirement (§3.2) is a real constraint on how τ and B are
  interpolated. **[Q]** how smooth are the operational τ rasters, given they
  come from hard-edged polygon partitions? The Fortran relies on bicubic
  smoothing of a blocky field — the port must reproduce or improve on that.
- **Inputs as native GIS:** margin = vector polygon; τ = vector (attribute) or
  raster; bed = raster. This matches GRASS data types directly and removes the
  custom binary format entirely.
- **Crossover / motorcycle-graph and polygon splitting:** the delicate part.
  **[Q]** reimplement the same greedy shorter-flowline-wins rule, or lean on a
  contour/offset-curve formulation? Fidelity vs. robustness tradeoff to decide
  explicitly.
- **Validation:** the port should reproduce the Fortran on the bundled
  Greenland case (and a flat-bed uniform-τ case against the analytic Nye
  parabola) before any refactor of the algorithm. The §6 discrepancies mean
  bit-for-bit agreement is not the target; agreement within a documented
  tolerance is.

---

*End of Design Note 01.*
