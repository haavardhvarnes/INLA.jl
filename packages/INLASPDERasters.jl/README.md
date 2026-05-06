# INLASPDERasters.jl

Raster ↔ SPDE glue: extract covariate values from `Rasters.Raster`
sources onto SPDE mesh vertices, and return posterior-field predictions
as `Rasters.Raster` objects.

Companion package to [`INLASPDE.jl`](../INLASPDE.jl/).

## Status

`v0.4.0`. The package ships a working raster-SPDE bridge and is
exercised end-to-end by the
[Meuse SPDE vignette](../../docs/src/vignettes/meuse-spde.md):

- [`extract_at_mesh`](src/extract.jl) — sample a `Raster` at mesh
  vertex coordinates.
- [`predict_raster`](src/predict.jl) — Gaussian-approximation and
  sample-based SPDE → raster projection.
- [`quantile_rasters`](src/predict.jl) — joint mean / sd / lower /
  upper credible-interval rasters from per-vertex summaries.
- [`Exceedance`](src/exceedance.jl) — quantity wrapper for tail
  probabilities under the sample-based path.

CRS mismatches between mesh and raster are surfaced at the API
boundary (see `mesh_crs` keyword, ADR-041).

## API at a glance

```julia
using INLASPDE, INLASPDERasters, Rasters, Random

# 1. Extract a covariate raster onto mesh vertices.
elev_raster = Raster("elevation.tif")
elev_at_vertices = extract_at_mesh(elev_raster, mesh)

# 2. Project a fitted SPDE component back onto a target grid.
template = Raster(zeros(nx, ny), (X(xs), Y(ys)))
r_mean  = predict_raster(model, res, template;
    component = SPDE2, quantity = :mean)
r_lower = predict_raster(model, res, template;
    component = SPDE2, quantity = :lower)

# 3. Sample-based path for exceedance probabilities P(η > c | y).
rng = Xoshiro(123)
r_exc = predict_raster(rng, model, res, template;
    component = SPDE2, quantity = Exceedance(0.5), n_samples = 1000)
```

`predict_raster(model, res, template)` requires the `SPDE2` component
to retain its mesh — construct it via `SPDE2(mesh::INLAMesh; …)`
rather than the back-compat `SPDE2(points, triangles; …)` form
(ADR-036).

## Why a separate package, not a weakdep of INLASPDE

- `Rasters.jl` pulls `GDAL_jll`, `Proj_jll`, `NetCDF_jll`, and other
  binary artifacts totalling hundreds of MB. A weakdep triggered by
  `using Rasters` would inflate the install / precompile cost for
  everyone who happens to load Rasters in the same session, regardless
  of whether they care about SPDE.
- Raster-specific covariate extraction belongs in one dedicated place,
  not sprinkled across an extension.

## See also

- [`INLASPDE.jl`](../INLASPDE.jl/) — the SPDE framework itself.
- [`Rasters.jl`](https://github.com/rafaqz/Rasters.jl) — the raster
  abstraction.
