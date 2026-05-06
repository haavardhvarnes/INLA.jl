# INLASPDERasters.jl

Raster glue for `INLASPDE.jl`: extract covariates at mesh vertices
from a `Rasters.Raster`, and project an SPDE posterior back onto a
raster grid. Lives in its own sub-package because Rasters transitively
pulls GDAL_jll / Proj_jll (~hundreds of MB) — too heavy for a weakdep
that most users will never trigger.

## What's here

- **`extract_at_mesh(raster, mesh; mesh_crs = nothing)`** — barycentric
  / nearest-neighbour sampling of a raster at the mesh vertex
  coordinates. When `mesh_crs` is supplied it is checked against the
  raster's CRS at the API boundary (ADR-041); mismatches raise rather
  than silently mis-locating points.
- **`predict_raster(values, mesh, template)`** — vertex-vector
  primitive: project a length-`num_vertices(mesh)` field onto
  `template`'s cell centres via `INLASPDE.MeshProjector`. Cells
  outside the mesh carry `missingval` (default `NaN`).
- **`predict_raster(model, res, template; component, quantity, level, mesh_crs)`** —
  Gaussian-approximation overload: pulls the SPDE component's
  per-vertex `:mean` / `:sd` / `:lower` / `:upper` summary out of
  `random_effects(model, res; level)` and forwards it through the
  primitive (ADR-040). `component` accepts an `Int` index, a
  component-name string, or `Type{<:SPDE2}` for auto-resolution.
- **`predict_raster(rng, model, res, template; component, quantity, n_samples, …)`** —
  sample-based overload: draws `n_samples` joint posterior samples,
  projects them all in a single sparse-dense GEMM, and reduces
  per cell. `quantity` accepts `:mean`, a `Real ∈ [0, 1]` (empirical
  quantile), or `Exceedance(c)` (per-cell `P(u > c | y)`).
- **`quantile_rasters(mean_v, sd_v, mesh, template; z, …)`** — vertex-
  vector helper that builds a NamedTuple `(; mean, sd, lower, upper)`
  of `Raster`s from a mean / sd pair under a Gaussian assumption.

## Why a sub-package, not a weakdep

Rasters' transitive closure (GDAL_jll, Proj_jll, NetCDF_jll) is large
and license-fragile across platforms. Gating the cost behind an
explicit `Pkg.add("INLASPDERasters")` makes the load-time and install-
size impact visible to users, and means CI for `INLASPDE.jl` does not
have to include the Rasters matrix.

## API

```@autodocs
Modules = [INLASPDERasters]
```
