# Phase O — Raster bridge maturity (`INLASPDERasters.jl`)

## Context

Phase O closes the user-facing geospatial story per
[`plans/replan-2026-04-28.md:434-451`](replan-2026-04-28.md). Phase N
landed at commit `cbfae1e` (`f((east, north, time), KroneckerComponent(...))`
3-tuple `@lgm` syntax); Phase M shipped SPDE expansion at v0.2.0 with
the Cameletti space-time oracle green; Phase G's
`AbstractObservationMapping` seam has been load-bearing since v0.1 and
already accepts the per-obs `MeshProjector` as an `as_matrix(...)`
producer. The replan estimate is 4–6 weeks; the audit below collapses
that to ~2 weeks because most of the surface area is already shipped
and what remains is glue, fixtures, and vignette polish.

The audit also surfaces a re-scoping issue the replan papers over:
two of the three Phase O scope items are essentially complete in
`packages/INLASPDERasters.jl/` v0.2.0. What is genuinely missing is
(a) a `(model, res, …)`-shaped user-facing entry point that wraps the
already-shipped vertex-vector primitives, (b) exceedance-probability
rasters, (c) the Meuse vignette actually rendering predicted rasters,
(d) the long-deferred Meuse R-INLA `predict.inla` oracle fixture.
Phase O is therefore much closer to "publish what we already built"
than to "build a raster prediction surface."

## State of the world (read-only audit)

What ships today, against the replan's Phase O surface:

| Replan scope item | Today | Gap |
|---|---|---|
| **Item 1 — `predict_raster(model, res, raster_template)`** | Vertex-vector form `predict_raster(values::AbstractVector, mesh, template; outside, missingval)` ships at [`packages/INLASPDERasters.jl/src/predict.jl:35`](../packages/INLASPDERasters.jl/src/predict.jl). Builds a barycentric `MeshProjector(mesh, cell_centres)` with `outside=:zero` masking; reproduces linear fields exactly (regression). `quantile_rasters(mean, sd, mesh, template; z, outside, missingval)` ships at line 127, returning `(mean, sd, lower, upper)` rasters via `mean ± z·sd` linear projection. | No `(model, res, ::Raster)` overload; the user must today extract the SPDE component slice from `random_effects(model, res)["..."]` themselves and pass `mean` / `sd` vectors in. The docs page [`docs/src/packages/inlaspderasters.md:15`](../docs/src/packages/inlaspderasters.md) advertises `predict_raster(fit, mesh, template; quantity = :mean)` — that signature does not exist. **Exceedance probability `P(u > c | y)` is not implemented at all.** |
| **Item 2 — bilinear + nearest covariate extraction at mesh nodes** | `extract_at_mesh(raster, mesh; method = :bilinear|:nearest, outside, missingval)` ships at [`packages/INLASPDERasters.jl/src/extract.jl:45`](../packages/INLASPDERasters.jl/src/extract.jl). Both methods land. Ascending and descending coord axes handled, outside-extent policy enforced, regression tests confirm bilinear reproduces affine fields to machine precision. | The replan calls this "currently partial." **The audit disagrees: it is complete except for CRS handling.** No `crs.jl` file exists despite the package plan's `src/crs.jl` slot; the docstring says "the caller must pre-project." `INLAMesh` carries no CRS, so equality assertion has nowhere to read from. The CRS test fixture promised in `CLAUDE.md` (`test/regression/test_extract_crs.jl`) does not exist. |
| **Item 3 — composability with the SPDE projector under `AbstractObservationMapping`** | `MeshProjector` at [`packages/INLASPDE.jl/src/projector.jl:30`](../packages/INLASPDE.jl/src/projector.jl) materialises `A::SparseMatrixCSC`; the LGM constructor wraps any `AbstractMatrix` in `LinearProjector`, so fit-time composability already works (every Phase M oracle uses it). `posterior_predictive(rng, res, model, mapping; n_samples)` at [`packages/LatentGaussianModels.jl/src/inference/diagnostics.jl:146`](../packages/LatentGaussianModels.jl/src/inference/diagnostics.jl) accepts any `AbstractObservationMapping`, including a `LinearProjector` wrapping a mesh→raster projection matrix. | **The audit shows item 3 is structurally complete.** What's missing is a public entry point that *uses* this composability: nothing in `INLASPDERasters.jl` currently calls `posterior_predictive(...)`, so sample-based exceedance probabilities have no ergonomic surface. The "composability under `AbstractObservationMapping`" framing in the replan is ambiguous; the actionable item is "expose a sample-based prediction path through the existing seam," not "build a new seam." |
| **Acceptance — Meuse vignette publishes posterior-mean and exceedance-probability rasters end-to-end** | [`docs/src/vignettes/meuse-spde.md`](../docs/src/vignettes/meuse-spde.md) currently fits the Gaussian + Intercept + dist + SPDE model and prints `fixed_effects(model, res)` plus `(range, σ, τ_noise)`. Last-section text caveats that `@lgm` cannot yet express Meuse — but PR-7c landed coordinate-tuple syntax, so this caveat is now stale. | **No raster output is rendered.** No template-raster construction, no `predict_raster` call, no exceedance map, no `@lgm` form. README at [`packages/INLASPDERasters.jl/README.md:11-13`](../packages/INLASPDERasters.jl/README.md) still says "scaffolding only — not yet implemented" despite shipping working extraction and projection. |

Three cross-cutting realities the audit surfaces:

1. **Replan item #2 is mis-scoped.** "Bilinear + nearest covariate
   extraction (currently partial)" reads as "bilinear missing." It
   isn't — bilinear is the default and is regression-tested to
   machine precision against affine fields. What is partial is
   **CRS handling**, which the package plan tracks separately at
   [`packages/INLASPDERasters.jl/plans/plan.md:45`](../packages/INLASPDERasters.jl/plans/plan.md)
   ("deferred until `INLAMesh` carries CRS metadata").

2. **Replan item #3 is largely vacuous as written.** The `AbstractObservationMapping`
   seam (ADR-017) was designed to host any sparse projector including
   the mesh-barycentric `MeshProjector`. It does, today. The seam-
   level composability story is shipped. The actionable shape is
   "expose Phase K-style `posterior_predictive` plumbing as a
   raster-side primitive," which is a thin wrapper rather than a new
   abstraction.

3. **The Meuse R-INLA `predict.inla` oracle was deferred at M2 close**
   ([`packages/INLASPDERasters.jl/plans/plan.md:58-62`](../packages/INLASPDERasters.jl/plans/plan.md)).
   Phase M's INLASPDE-side oracle (`meuse_spde.jld2`) covers the
   posterior fit; there is no vertex-mean/sd raster fixture for the
   raster-side projection. This is the one numerical gap that
   genuinely deserves new fixture work.

## Pre-Phase-O housekeeping (audit findings)

Before starting Phase O implementation, two pieces of Phase N close
slipped and should be reconciled:

- **ADRs 036–039 not in `plans/decisions.md`.** The Phase N PR-7 plan
  ([`plans/phase-n-pr7.md:24-50`](phase-n-pr7.md)) called for
  ADR-036 (`SPDE2.mesh` retention), ADR-037 (`@lgm` tuple-coord
  syntax), ADR-038 (3-tuple Kronecker route), and ADR-039 (weakdep
  extension boundary). The commits referenced these ADR numbers in
  their messages, but `plans/decisions.md`'s last entry remains
  ADR-035. Either backfill ADRs 036–039 from the PR-7 plan body
  before Phase O starts (preferred) or fold into the Phase O close
  with an explanatory CHANGELOG note.
- **No `v0.2.1` / `v0.2.2` tag.** PR-7's plan targeted `v0.2.2` at
  PR-7c close; the most recent tag is `v0.2.0`. Either tag `v0.2.2`
  retroactively at `cbfae1e` (preferred — clean phase boundaries) or
  fold the PR-7 close into the Phase O close at `v0.3.0`.

The recommended sequencing is: open Phase O with a small housekeeping
PR (PR-0 below) that backfills the ADRs and applies the missing tag,
then move into PR-1.

## Design calls (ADRs)

Phase O surfaces three ADR candidates for the body of `decisions.md`.
Recommended numbers continue from ADR-035; Phase N PR-7's ADRs 036–039
land in the housekeeping PR-0 first, so Phase O's entries are 040–042.

### ADR-040 candidate — `predict_raster(model, res, …)` user-facing surface: posterior mean + sd via `random_effects`, exceedance via `posterior_sample`

The user-facing entry point lifts the existing vertex-vector
`predict_raster` and `quantile_rasters` to `(model, res, …)` shape.
Two semantic axes:

1. **Where do the vertex statistics come from?** Two paths:
   - (a) **`random_effects(model, res)` Gaussian-approximation slice.**
     `random_effects` returns per-component mean / sd / lower / upper
     vectors, where the SPDE component slice is exactly `length(mesh.points)`
     long. Cheap, no sampling. This is the path that drives the
     existing `quantile_rasters` interval rasters.
   - (b) **`posterior_sample(rng, res, model)` joint draws.** Slice
     the SPDE block out of `x_samples` (using `model.latent_ranges[i]`),
     compute per-vertex empirical quantiles or `P(u > c)` exceedance
     across draws, project the resulting vertex vector through the
     same barycentric `MeshProjector`. Sample-based; correctly
     captures non-Gaussian posterior tails. Costs `n_samples` Cholesky
     back-solves.

2. **Quantile semantics.** The current `quantile_rasters` projects
   `mean ± z·sd` linearly through `P`. That is sharp at vertices and
   linear inside triangles — a "vertex Gaussian" view, not a
   pixel-level Gaussian. For posterior **exceedance** `P(u(s) > c | y)`
   the linear-projection trick fails (probabilities don't compose
   linearly under barycentric averaging). The honest path is sampling
   per (a), then `mean(samples .> c, dims = 2)`.

**Recommendation.** Both paths land:

- `predict_raster(model, res, template; component, quantity = :mean,
  outside, missingval)` — Gaussian-approximation summary. `quantity ∈
  (:mean, :sd, :lower, :upper)`. `component` is the SPDE component to
  slice (`Int` index, or a name `String` matching `random_effects`
  keys). Wraps existing primitives; ~30 LOC.
- `predict_raster(rng, model, res, template; component, quantity, n_samples,
  outside, missingval)` — sample-based. `quantity` accepts numeric
  quantiles in `[0, 1]`, the symbol `:mean`, or an `Exceedance(c)`
  wrapper for `P(u > c | y)`. Calls `posterior_sample(rng, res, model,
  n_samples)`, slices the SPDE block, projects each draw through the
  per-cell `MeshProjector` (or applies the reduction at vertices then
  projects — see "performance" below).

The `Exceedance(c)` wrapper is a tiny `struct Exceedance{T}; c::T end`
exported from `INLASPDERasters` so the dispatch is type-stable and
the API self-documents.

**Performance note.** Build the `MeshProjector` once per
`predict_raster` call and reuse `P.A` across draws —
`P.A * x_samples[spde_range, :]` is one sparse-dense GEMM and gives a
`n_cells × n_samples` matrix from which any reduction (mean, sd,
quantile, exceedance) is a column-wise reduction. For 10⁶-cell rasters
× 10³ draws this is the only tolerable shape; element-wise per-draw
projection is 1000× slower.

This is a contract change visible only to `INLASPDERasters.jl`;
nothing upstream changes.

### ADR-041 candidate — CRS policy: `predict_raster` and `extract_at_mesh` reject mismatched CRS at the API boundary; reproject is opt-in via a `mesh_crs` keyword

The `INLAMesh` struct does not carry CRS metadata
([`packages/INLASPDE.jl/src/mesh/inla_mesh.jl:19-25`](../packages/INLASPDE.jl/src/mesh/inla_mesh.jl)).
`Rasters.Raster` does, via `Rasters.crs(raster)`. Three shapes:

(a) **Reject mismatch silently — caller pre-projects.** The current
    behavior. Cheap, but silent CRS mismatches are the single most
    common geo foot-gun. A WGS84 raster (degrees) and a UTM mesh
    (metres) cannot be flagged.

(b) **Add `mesh_crs::Union{Nothing, CRS} = nothing` keyword.** When
    supplied, assert `mesh_crs == Rasters.crs(raster)` at the API
    boundary; `nothing` keeps current "trust the caller" behavior.
    Errors point at the mismatch concretely. No `INLAMesh` change.
    `CoordRefSystems.jl` is already in `INLASPDE`'s deps
    (`packages/INLASPDE.jl/src/INLASPDE.jl:31`); the type is
    available for keyword-arg use without new deps.

(c) **Promote `INLAMesh` to carry `crs::Union{Nothing, CRS}`.** Field
    addition is a contract change (PR-7a-style audit risk); every
    Phase M oracle fixture needs to round-trip the new field.
    Heavy.

**Recommendation: (b).** The `mesh_crs` keyword adds zero cost when
omitted (back-compat) and the cheapest possible safety net when
supplied. Reprojection itself is **out of scope** for Phase O —
reprojecting the *raster* requires `Proj_jll` which is already in the
INLASPDERasters transitive closure, but reprojecting the *mesh* would
require changing `mesh.points` in place, which mutates SPDE FEM
matrices that have already been assembled. Defer raster reprojection
to v0.4 / Phase O+1 with a `reproject = true` opt-in and an explicit
pre-condition that the caller reassemble FEM matrices on the
reprojected mesh.

Promote (c) — `INLAMesh.crs` field — when CRS becomes load-bearing
elsewhere (sphere SPDE, geodesic distance priors, multi-CRS joint
likelihoods). All deferred to v0.4+.

### ADR-042 candidate — INLASPDERasters takes a soft dep on LatentGaussianModels via the `(model, res, …)` overload

The current `INLASPDERasters/Project.toml` already lists
`LatentGaussianModels = "0.2"` — this dep was added at package
inception even though no `src/` file uses LGM types today. Phase O's
`predict_raster(model, res, …)` overloads turn that latent dep into
an active one: the function dispatches on `LatentGaussianModel` and
reads `res.x_mean`, `res.x_var`, `model.latent_ranges`,
`random_effects(model, res)`, and `posterior_sample(rng, res, model)`.

The dep was foreseen and accepted at scaffolding time
([`plans/dependencies.md:101`](dependencies.md)). No new ADR
is strictly required — this is an "ADR for the file" only because the
dep moves from "declared but unused" to "load-bearing for the public
API." If the writer prefers, this can land as an inline note in
`CHANGELOG.md` rather than a numbered ADR. Recommend the inline note
unless the LGM-shape coupling reveals a layering concern during PR-1
implementation.

## PR sequence

Five PRs (plus an optional housekeeping PR-0), with PR-5 as a stretch
tail. Total ~2 weeks at solo developer pace. The replan's 4–6 week
estimate assumed building extraction and projection from scratch; the
audit shows those are shipped, which compresses the surface to "wrap,
fixture, document, release."

### PR-0 — Phase N close housekeeping (ADR backfill + tag)

Reconcile the slipped pieces of Phase N close:

**Files:**
- `plans/decisions.md` — append ADR-036 (`SPDE2.mesh` retention),
  ADR-037 (`@lgm` tuple-coord syntax), ADR-038 (3-tuple Kronecker
  route), ADR-039 (weakdep extension boundary). Body text taken
  from `plans/phase-n-pr7.md` lines 24–50 plus the corresponding PR
  commit messages.
- `CHANGELOG.md` — `v0.2.2` entry for Phase N close.

**Tag:** `v0.2.2` retroactively on `cbfae1e` (PR-7c commit).

**Tests:** none — pure documentation / metadata.

### PR-1 — `predict_raster(model, res, …)` Gaussian-approximation overload (ADR-040 first half)

User-facing entry point lifting the vertex-vector primitive.

**Files:**
- `packages/INLASPDERasters.jl/src/predict.jl` — append two new
  methods on `predict_raster`:
  - `predict_raster(model::LatentGaussianModel, res::INLAResult,
    template::Raster; component, quantity = :mean, level = 0.95,
    outside = :missing, missingval = NaN, mesh_crs = nothing)`. Uses
    `random_effects(model, res; level)` to fetch mean / sd / lower /
    upper, slices the named or indexed SPDE component, locates the
    component's mesh via `model.components[i].mesh::INLAMesh`
    (Phase N PR-7a's mesh field), forwards to the existing vertex
    primitive.
  - Same method with `quantity = :sd` / `:lower` / `:upper` switching
    on which `random_effects` slot to project.
- `packages/INLASPDERasters.jl/src/INLASPDERasters.jl` — update
  `using` block to import `LatentGaussianModel`, `INLAResult`,
  `random_effects`, `LatentGaussianModels`. Re-export `Exceedance`
  symbol (lands in PR-2).
- `packages/INLASPDERasters.jl/src/component_resolution.jl` (new) —
  `_resolve_spde_component(model, component_id) -> (i, mesh)` helper
  that accepts `Int` (1-based component index), `String` (name match
  via `_component_name`), or `Type{<:SPDE2}` (auto-locate the
  unique SPDE2 component, error if 0 or ≥2 found). Returns
  the component index and the retained `INLAMesh`. Validates the
  component is `SPDE2`; rejects `SPDE1D`, `SPDE2NonStationary`
  (the latter has its own mesh but the projector geometry differs;
  defer to PR-5 if demand surfaces).

**Tests** (`packages/INLASPDERasters.jl/test/regression/test_predict_model.jl`):
- Build a small Meuse-shape model (Gaussian + Intercept + dist + SPDE
  on a 50-vertex mesh, fit on synthetic data with `inla(model, y)`),
  call `predict_raster(model, res, template; component = "spde")`.
  Compare with the existing vertex-vector form
  `predict_raster(re["spde"].mean, mesh, template)` for exact
  agreement.
- `quantity = :sd` / `:lower` / `:upper` agree similarly.
- `component` resolution: pass `Int`, name-`String`, and `SPDE2`
  type; pass an out-of-range index, an unknown name, an ambiguous
  type when two SPDE2 components exist — each returns a
  user-readable error pointing at the component list.
- `mesh_crs` keyword: omitted is back-compat; supplied with a
  matching `EPSG:31370` (Meuse Lambert72) passes; supplied with
  `EPSG:4326` against an unprojected `template` raises
  `ArgumentError` naming both sides of the mismatch.

**Tests** (`packages/INLASPDERasters.jl/test/regression/test_extract_crs.jl`,
new — fills the gap promised in `CLAUDE.md:34`):
- `extract_at_mesh(raster, mesh; mesh_crs = ...)` rejects WGS84
  raster against UTM `mesh_crs` with a mismatch error; accepts
  matched.
- Omitting `mesh_crs` preserves the current "trust the caller"
  behavior (no regression on the M1 tests).

**Compat bumps:** `INLASPDERasters.jl` 0.2.0 → 0.3.0 (additive
contract change: new public method shapes; new `mesh_crs` kwarg).
The dep on `LatentGaussianModels = "0.2"` was already declared.

### PR-2 — Sample-based `predict_raster` for exceedance + true posterior quantiles (ADR-040 second half)

The Phase O acceptance criterion's exceedance-probability rasters
land here.

**Files:**
- `packages/INLASPDERasters.jl/src/predict.jl` — append a third
  method:
  - `predict_raster(rng::AbstractRNG, model::LatentGaussianModel,
    res::INLAResult, template::Raster; component, quantity,
    n_samples = 1000, outside = :missing, missingval = NaN,
    mesh_crs = nothing)`. Calls `posterior_sample(rng, res, model;
    n_samples)`, slices `x_samples[latent_ranges[i], :]`, projects
    `n_samples × n_v` block through `P.A` in one sparse GEMM
    (`P.A * X = n_cells × n_samples`), reduces column-wise:
    - `quantity == :mean` → `mean(η_samples, dims = 2)`.
    - `quantity isa Real && 0 ≤ quantity ≤ 1` → `quantile(η_samples,
      quantity)` per cell (call `Statistics.quantile` row-wise).
    - `quantity isa Exceedance` → `mean(η_samples .> quantity.c, dims
      = 2)`.
- `packages/INLASPDERasters.jl/src/exceedance.jl` (new) —
  `struct Exceedance{T <: Real}; c::T end`, exported.

**Tests** (`packages/INLASPDERasters.jl/test/regression/test_predict_sample.jl`):
- Synthetic Meuse-shape fit; for `quantity = :mean` with
  `n_samples = 5000`, sample-based mean agrees with
  Gaussian-approximation mean to MC tolerance (3 sd of the sample
  mean).
- `quantity = 0.5` (median) agrees with `:mean` for symmetric
  Gaussian posterior to MC tolerance.
- `quantity = Exceedance(c)` for a constant-vertex `c = -∞`
  returns 1.0 in every cell; `c = +∞` returns 0.0.
- For a known threshold and a constant-mean / non-trivial-sd vertex
  field, exceedance equals `1 - Φ((c - μ) / σ)` to MC tolerance.
- RNG reproducibility: `rng_seeded` produces bitwise-identical
  rasters across two calls.

**Compat bumps:** `INLASPDERasters.jl` 0.3.0 → 0.4.0 (new public
type `Exceedance`; new method signature with `rng` first arg).

### PR-3 — Meuse vignette renders posterior-mean + exceedance-probability rasters end-to-end (acceptance gate)

The Phase O acceptance close.

**Files:**
- `docs/src/vignettes/meuse-spde.md` — append three new sections
  after "Comparing to R-INLA":
  - **"Posterior surface" `@example`.** Construct a template raster
    over the Meuse mesh extent at 50m resolution, call
    `predict_raster(model, res, template; component = "SPDE2_3",
    quantity = :mean)`. Render with `Makie` (existing
    `INLASPDEMakieExt`) showing the posterior-mean log-zinc surface.
  - **"Posterior credible interval" `@example`.** Two more
    `predict_raster` calls for `:lower` / `:upper` and a Makie
    panel comparing the three surfaces.
  - **"Exceedance probability" `@example`.** Sample-based
    `predict_raster(rng, model, res, template; component = "SPDE2_3",
    quantity = Exceedance(log(500.0)), n_samples = 2000)` rendering
    `P(zinc > 500 ppm | y)` as a heatmap. Threshold `500 ppm` is the
    standard Meuse environmental-risk cutoff used in the gstat
    literature.
- `docs/src/vignettes/meuse-spde.md:126-140` — replace the "What
  about `@lgm`?" stub paragraphs with a live `@example` block
  showing `@lgm log_zinc ~ 1 + dist + f((east, north), spde) data=meuse_df`
  matching the explicit constructor — Phase N PR-7c made this
  expressible. Re-uses the same `predict_raster` calls verbatim.
- `packages/INLASPDERasters.jl/README.md` — replace the
  "scaffolding only" wording (lines 9-21) with a status block
  pointing at the docs vignette and the live API surface.
- `docs/src/packages/inlaspderasters.md:11-20` — sync the
  function-signature sketches with the actual shipped surface
  (`predict_raster(model, res, template; component, quantity, …)`,
  not the stale `predict_raster(fit, mesh, template; quantity = :mean)`).

**Tests** (no new src tests — this is documentation-shaped). Doctest
coverage:
- The vignette is built under `Documenter.makedocs(strict = true)` in
  CI; the `@example` blocks failing breaks doc build, which is the
  acceptance gate. No separate "rendered raster" assertion.

### PR-4 — Meuse R-INLA `predict.inla` oracle fixture (deferred from M2 close)

Long-deferred numerical-correctness gate. The Phase M Meuse oracle
fixture (`packages/INLASPDE.jl/test/oracle/fixtures/meuse_spde.jld2`)
currently carries `summary_fixed`, `summary_hyperpar`, `mlik`, and
`A_field` but **not** R-INLA's posterior-mean raster from
`inla.mesh.project` or `predict.inla`.

**Files:**
- `scripts/generate-fixtures/meuse_spde_predict.R` (new) — extends
  the existing `meuse_spde.R` fixture script with:
  - A regular-grid prediction surface over the Meuse hull at 50m
    resolution (matching what the vignette will render).
  - `inla.mesh.project(mesh, loc = grid_centres)` — R-INLA's mesh→pixel
    projector (the load-bearing reference for `predict_raster`).
  - Pull `result$summary.random$field$mean` and `$sd` for the
    SPDE component.
  - Save grid coordinates, R-INLA's projected mean and sd at each
    grid cell, plus the vertex-level mean and sd vectors, to
    `meuse_spde_predict.jld2`.
- `scripts/generate-fixtures/_jld2_convert_meuse_predict.jl` —
  R-side `.RData` → JLD2 conversion (the existing fixture pipeline
  pattern from PR-7c's `cameletti_pm10` script).
- `packages/INLASPDERasters.jl/test/oracle/fixtures/meuse_spde_predict.jld2`
  (generated artifact; check in via Git LFS per the existing
  fixture convention).
- `packages/INLASPDERasters.jl/test/oracle/test_meuse_predict.jl`
  (new):
  - Load the fixture; build the Meuse mesh with the fixture's
    `points`/`tv`; run `predict_raster(vertex_mean, mesh, template)`.
  - Assert pixel-wise agreement with R-INLA's `inla.mesh.project`
    output to **1e-10** absolute tolerance — both implementations
    are exact P1 barycentric interpolation, so this is a tight
    regression gate, not a Phase F-style 5% tolerance.
  - Same gate for `predict_raster(vertex_sd, mesh, template)`.
  - Larger gate (Phase F's 5%) for the `predict_raster(model, res,
    template; quantity = :mean)` end-to-end fit-then-project path,
    where the vertex mean carries Julia-INLA vs R-INLA fit-side
    differences.
- `packages/INLASPDERasters.jl/test/runtests.jl` — wire in the new
  oracle testset alongside the existing M1/M2/M3 testsets.
- `packages/INLASPDERasters.jl/Project.toml` — add `JLD2` to test
  extras (matching `INLASPDE.jl/Project.toml`'s pattern).

**Tests:** the fixture itself is the test surface. The 1e-10 gate
against `inla.mesh.project` is the closest thing to a "byte-for-byte
match R-INLA" assertion in the entire repo — both sides evaluate the
same closed-form barycentric formula.

### PR-5 (stretch) — Sample-based `predict_raster` for `KroneckerComponent` space-time and `SPDE2NonStationary`

If PRs 1–4 close inside week 2, ship the obvious next-step:

- `KroneckerComponent` space-time: the SPDE block in
  `model.latent_ranges[i]` has length `n_v · n_t`. The user-facing
  call needs to pick a time slice (e.g. `time_index = 7`) or
  produce a stack of time-slice rasters. Stretch goal because the
  Cameletti vignette doesn't currently render rasters either.
- `SPDE2NonStationary`: same projector geometry as `SPDE2`, but the
  user-facing component-resolution helper currently rejects it.
  Single-line dispatch addition.

If neither lands inside Phase O's 2-week budget, both defer to a
Phase O+1 follow-up under v0.4.0 (next minor) with a tracking issue.

## Test surface (Phase O close gate)

- 1 new oracle fixture (PR-4) — Meuse `inla.mesh.project` reference;
  oracle tier promotes from "Phase M's fit-side oracle only" to
  "fit-side + raster-side oracle." Counted into the JLD2 fixture
  total: 36 → 37.
- 6 new regression tests (PR-1: 3 — model overload, sd/lower/upper,
  CRS; PR-2: 3 — sample mean, exceedance, RNG reproducibility).
- 1 new doctest gate (PR-3) — the Meuse vignette `@example` blocks
  fail under `Documenter.makedocs(strict = true)` if any
  `predict_raster` signature or rendering breaks.
- Aqua + JET clean across the new code paths
  (`packages/INLASPDERasters.jl/test/quality/`).

## Out of scope for Phase O

- **Reproject-on-the-fly between mismatched CRSs.** ADR-041 lands
  the assertion-only path; reprojection requires reassembling FEM
  matrices on the reprojected mesh and is deferred to v0.4+.
- **`INLAMesh` carrying `crs` field.** Promoted only when CRS
  becomes load-bearing for upstream SPDE machinery (sphere SPDE,
  geodesic priors). Per
  [`packages/INLASPDE.jl/plans/plan.md:213-219`](../packages/INLASPDE.jl/plans/plan.md),
  3D / sphere / non-separable space-time SPDE all land at v0.3+.
- **3D raster prediction.** `Rasters.Raster` supports more than 2
  dims, but `MeshProjector` is 2D-only; 3D projection waits on the
  3D SPDE deferral.
- **`predict_raster` for non-mesh components.** IID / RW1 / Besag
  fields are areal, not coordinate-indexed; raster projection has
  no meaning for them. The `_resolve_spde_component` helper
  rejects with a user-readable pointer; "areal-on-raster" is a
  separate v0.4+ design question.
- **R-INLA `inla.posterior.sample` parity.** Sampling shapes differ
  enough between R-INLA and our `posterior_sample` that an oracle
  fixture comparing sample-based exceedance probabilities would
  need an MC tolerance (~5%, Phase L tier). Skip the comparison;
  the Gaussian-approximation oracle in PR-4 covers the
  numerically-tight gate.
- **Pixel-level exact `diag(P Σ Pᵀ)`.** The current `quantile_rasters`
  docstring already flags this as deferred until LGM exposes a
  `marginal_variances`-shaped pixel covariance helper. Continue to
  defer; the sample-based path in PR-2 is the practical workaround.
- **Raster-side parallelism / threading.** The PR-2 sparse GEMM is
  already vectorised; threaded reduction over draws is a v0.4+
  performance polish, not a Phase O scope item.

## Critical files

| Concern | Path |
|---|---|
| Vertex-vector `predict_raster` (target of PR-1, PR-2 extensions) | [`packages/INLASPDERasters.jl/src/predict.jl`](../packages/INLASPDERasters.jl/src/predict.jl) |
| `extract_at_mesh` (target of PR-1 CRS keyword + test fill-in) | [`packages/INLASPDERasters.jl/src/extract.jl`](../packages/INLASPDERasters.jl/src/extract.jl) |
| `MeshProjector` (the seam everything else composes through) | [`packages/INLASPDE.jl/src/projector.jl`](../packages/INLASPDE.jl/src/projector.jl) |
| `posterior_sample`, `posterior_predictive` (PR-2 sample-based path) | [`packages/LatentGaussianModels.jl/src/inference/diagnostics.jl`](../packages/LatentGaussianModels.jl/src/inference/diagnostics.jl) |
| `random_effects`, `model.latent_ranges` (PR-1 component slicing) | [`packages/LatentGaussianModels.jl/src/inference/accessors.jl`](../packages/LatentGaussianModels.jl/src/inference/accessors.jl) |
| `SPDE2` mesh field (PR-1 component-mesh resolution depends on Phase N PR-7a) | [`packages/INLASPDE.jl/src/components/spde2.jl`](../packages/INLASPDE.jl/src/components/spde2.jl) |
| Meuse vignette (target of PR-3 acceptance) | [`docs/src/vignettes/meuse-spde.md`](../docs/src/vignettes/meuse-spde.md) |
| Package docs page (PR-3 sync) | [`docs/src/packages/inlaspderasters.md`](../docs/src/packages/inlaspderasters.md) |
| Package README (PR-3 sync — currently stale) | [`packages/INLASPDERasters.jl/README.md`](../packages/INLASPDERasters.jl/README.md) |
| Package plan (close M1/M2/M3 + add Phase O completion) | [`packages/INLASPDERasters.jl/plans/plan.md`](../packages/INLASPDERasters.jl/plans/plan.md) |
| Fixture script (PR-4 generation) | `scripts/generate-fixtures/meuse_spde_predict.R` (new) |
| Oracle fixture (PR-4 artifact) | `packages/INLASPDERasters.jl/test/oracle/fixtures/meuse_spde_predict.jld2` (new) |
| ADR registry (Phase O entries + N housekeeping) | [`plans/decisions.md`](decisions.md) — append ADR-036–039 (PR-0 housekeeping), then ADR-040, ADR-041, optionally ADR-042 |

## Verification

Phase O closes when **all** of the following hold:

1. PRs 1–4 land green on CI; the existing 36 oracle fixtures + the
   new `meuse_spde_predict.jld2` (37 total) all pass.
2. The Meuse vignette renders the posterior-mean raster, the
   posterior credible-interval rasters, and the exceedance-probability
   raster in `Documenter.makedocs(strict = true)` mode without
   warnings.
3. `random_effects(model, res)["…"]` → `predict_raster(model, res,
   template)` agrees with the manual `predict_raster(values, mesh,
   template)` form to machine precision (PR-1 regression gate).
4. R-INLA `inla.mesh.project` agreement to 1e-10 absolute on the
   Meuse fixture (PR-4 oracle gate).
5. `INLASPDERasters.jl/README.md` no longer claims "scaffolding only."
6. ADRs 036–039 (PR-0) and 040–041 (Phase O) appended to
   `plans/decisions.md`; ADR-042 elected as inline CHANGELOG note or
   numbered ADR per PR-1 implementation outcome.
7. CHANGELOG entry for `v0.3.0` covering all five PRs (PR-0 through
   PR-4).

## Release target

`v0.3.0` — first minor bump after Phase M's `v0.2.0`. Per-package:

- `INLASPDERasters.jl` 0.2.0 → 0.3.0 (PR-1 additive contract change),
  → 0.4.0 (PR-2 new public type `Exceedance` + sample-based method).
  Fold both bumps into a single `0.3.0` release at Phase O close
  rather than tagging separately — the audit shows both shipping as
  one cohesive surface, and the user-facing CHANGELOG line is "you
  can now call `predict_raster(model, res, template)`," not two
  separate beats.
- `INLASPDE.jl` 0.3.0 → 0.3.1 (no contract changes; patch bump for
  the Phase N PR-7 ADRs landing in `decisions.md`).
- `LatentGaussianModels.jl` 0.2.0 — unchanged.
- `LGMFormula.jl` 0.4.0 — unchanged.
- `GMRFs.jl` 0.1.2 — unchanged.
- Umbrella `INLA.jl` 0.2.0 → 0.3.0 (Phase O close); CHANGELOG entry
  covers all five PRs.
- Tag at Phase O close: `v0.3.0` on `main`, release commit
  `chore(release): v0.3.0 — Phase O close (raster bridge maturity)`.
- PR-0 prerequisites: tag `v0.2.2` retroactively on `cbfae1e`
  (Phase N PR-7c close).

The replan's exit register
([`plans/replan-2026-04-28.md:604`](replan-2026-04-28.md))
predicted Phase O at `v0.4.0-alpha`. The audit-driven re-scope
(scope items 2 and 3 already shipped) means Phase O is a single
minor bump rather than a beta-tail; promote Phase O's exit tag to
`v0.3.0` and let the replan's `v0.4.0` land at the first phase that
actually requires a major-feature beta cycle (likely sphere SPDE or
3D under the deferred v0.3+ surface).

## Cadence

Solo-developer pace, ~2 weeks total.

- **Week 1**
  - Day 1: PR-0 housekeeping (ADRs 036–039 backfill + `v0.2.2` tag).
    Pure metadata; ~half a day.
  - Days 2–3: PR-1 (model-overload `predict_raster` + CRS keyword
    + `test_extract_crs.jl` fill-in). The component-resolution
    helper is the load-bearing piece; everything else is glue.
  - Days 4–5: PR-2 (sample-based path + `Exceedance` type). The
    sparse-GEMM-once design keeps the PR mechanically small;
    most of the work is regression-test scaffolding.

- **Week 2**
  - Days 1–2: PR-4 (Meuse `inla.mesh.project` oracle). R-side
    fixture script first, JLD2 conversion, oracle test. The
    1e-10 tolerance gate is tight; expect to debug a coordinate
    ordering or row-major / column-major mismatch.
  - Days 3–4: PR-3 (Meuse vignette renders rasters end-to-end +
    docs / README sync + `@lgm` form replaces the stale stub
    paragraphs). The slowest part is Makie heatmap polish under
    Documenter; budget extra for that.
  - Day 5: ADRs 040 / 041, CHANGELOG, version bumps, tag, release
    script.

- **Stretch (optional, week 3 if time)**: PR-5 — `KroneckerComponent`
  + `SPDE2NonStationary` overloads. If this slips, defer to a
  follow-up tracking issue in v0.4 — Phase O closes regardless.

The aggressive 2-week target reflects the audit's main finding:
the replan budgeted 4–6 weeks for what is now mostly publication
work. If PR-4's R-side fixture pipeline reveals a substantive
disagreement (always possible — fmesher coordinates, projector row
order, mass-matrix lumping conventions), the cadence stretches one
week.
