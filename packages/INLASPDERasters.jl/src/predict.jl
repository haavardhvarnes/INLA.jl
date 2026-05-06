"""
    predict_raster(values::AbstractVector{<:Real},
                   mesh::INLAMesh,
                   template::Raster;
                   outside::Symbol = :missing,
                   missingval::Real = NaN) -> Raster

Project a vertex-valued field onto a raster grid matching `template`.

Builds a barycentric `INLASPDE.MeshProjector` from `mesh` to the cell
centres of `template` and applies it to `values`. The returned raster
has the same dimensions, extent, resolution, and dim order as
`template`; only the underlying array is replaced.

# Arguments

- `values` — length-`num_vertices(mesh)` vector: the per-vertex field
  to project. Typically the posterior mean (or any linear functional
  of it) of an SPDE component.
- `mesh` — the SPDE mesh `values` live on.
- `template` — a 2D `Raster` with `X` and `Y` dimensions that defines
  the target grid.

# Keywords

- `outside = :missing` — policy for raster cells that fall outside the
  mesh domain: `:error` throws, `:missing` substitutes `missingval`.
- `missingval = NaN` — sentinel used when `outside = :missing`.

# Returns

A `Raster` with the same dims as `template` whose data is the
projected field. Cells outside the mesh domain carry `missingval`.
"""
function predict_raster(
        values::AbstractVector{<:Real},
        mesh::INLAMesh,
        template::Raster;
        outside::Symbol=:missing,
        missingval::Real=NaN
)
    length(values) == num_vertices(mesh) || throw(ArgumentError(
        "values has length $(length(values)) but mesh has $(num_vertices(mesh)) vertices",
    ))
    outside ∈ (:error, :missing) ||
        throw(ArgumentError("outside must be :error or :missing; got $outside"))

    xs = collect(Rasters.lookup(template, X))
    ys = collect(Rasters.lookup(template, Y))
    nx = length(xs)
    ny = length(ys)

    # Stack cell centres in a fixed (i, j) scan so we can post-index.
    n_cells = nx * ny
    locs = Matrix{Float64}(undef, n_cells, 2)
    cell_ij = Vector{Tuple{Int, Int}}(undef, n_cells)
    k = 0
    for j in 1:ny
        for i in 1:nx
            k += 1
            locs[k, 1] = xs[i]
            locs[k, 2] = ys[j]
            cell_ij[k] = (i, j)
        end
    end

    # :error wants the projector to throw on first outside cell; :missing
    # lets us mask zero rows afterwards (barycentric rows sum to 1 for
    # interior points, so an all-zero row is an unambiguous outside
    # marker under :zero).
    proj_outside = outside === :error ? :error : :zero
    P = INLASPDE.MeshProjector(mesh, locs; outside=proj_outside)
    projected = P.A * Vector{Float64}(values)

    out = similar(template, Float64)
    fill!(out, missingval)
    @inbounds for k in 1:n_cells
        i, j = cell_ij[k]
        row = @view P.A[k, :]
        if !iszero(row)
            out[X=i, Y=j] = projected[k]
        end
    end
    return out
end

"""
    quantile_rasters(mean::AbstractVector{<:Real},
                     sd::AbstractVector{<:Real},
                     mesh::INLAMesh,
                     template::Raster;
                     z::Real = 1.959963984540054,
                     outside::Symbol = :missing,
                     missingval::Real = NaN) -> NamedTuple

Project per-vertex posterior mean and standard deviation onto a raster
grid, together with symmetric Gaussian credible-interval rasters.

# Semantics

Each output raster is a linear projection of the corresponding vertex
quantity via the barycentric mesh-to-pixel projector `P`:

- `mean  = P * mean_v`
- `sd    = P * sd_v`
- `lower = P * (mean_v - z * sd_v)`
- `upper = P * (mean_v + z * sd_v)`

Projecting `sd` linearly is a convenient but approximate
"vertex-level-quantile" view: the SPDE posterior at vertices is not
diagonal, so a cell-level exact standard deviation would require
`diag(P Σ P')`. The per-vertex interval is sharp at the vertex and
interpolated linearly inside each triangle; this matches the Gaussian
summary fields callers plot on R-INLA outputs.

# Keywords

- `z = 1.96` — the Gaussian quantile half-width (default = 97.5th
  percentile, giving a 95% interval). Pass `z = 1.645` for 90%, etc.
- `outside`, `missingval` — forwarded to [`predict_raster`](@ref).

# Returns

A `NamedTuple` `(; mean, sd, lower, upper)` of four `Raster`s, all
sharing the dims of `template`.
"""
function quantile_rasters(
        mean::AbstractVector{<:Real},
        sd::AbstractVector{<:Real},
        mesh::INLAMesh,
        template::Raster;
        z::Real=1.959963984540054,
        outside::Symbol=:missing,
        missingval::Real=NaN
)
    n = num_vertices(mesh)
    length(mean) == n ||
        throw(ArgumentError("mean has length $(length(mean)) but mesh has $n vertices"))
    length(sd) == n ||
        throw(ArgumentError("sd has length $(length(sd)) but mesh has $n vertices"))
    z >= 0 ||
        throw(ArgumentError("z must be non-negative; got $z"))
    all(s -> s >= 0, sd) ||
        throw(ArgumentError("sd must be non-negative"))

    m = Vector{Float64}(mean)
    s = Vector{Float64}(sd)
    lo_v = m .- z .* s
    up_v = m .+ z .* s

    mean_r = predict_raster(m, mesh, template; outside=outside, missingval=missingval)
    sd_r = predict_raster(s, mesh, template; outside=outside, missingval=missingval)
    lower_r = predict_raster(lo_v, mesh, template; outside=outside, missingval=missingval)
    upper_r = predict_raster(up_v, mesh, template; outside=outside, missingval=missingval)

    return (mean=mean_r, sd=sd_r, lower=lower_r, upper=upper_r)
end

"""
    predict_raster(model::LatentGaussianModel, res::INLAResult,
                   template::Raster;
                   component,
                   quantity::Symbol = :mean,
                   level::Real = 0.95,
                   outside::Symbol = :missing,
                   missingval::Real = NaN,
                   mesh_crs = nothing) -> Raster

Project a posterior summary of an SPDE component onto a raster grid
matching `template`. Wraps the vertex-vector
[`predict_raster`](@ref) primitive — the user-facing entry point lifts
the per-component slice out of `LatentGaussianModels.random_effects` and forwards
the corresponding vector through the barycentric mesh→raster
projector.

# Arguments

- `model` — a fitted `LatentGaussianModel` containing at least one
  `SPDE2` component constructed via `SPDE2(mesh::INLAMesh; …)` (so the
  mesh is retained per ADR-036).
- `res` — the `INLAResult` returned by `inla(model, y; …)`.
- `template` — a 2D `Raster` defining the target grid; the returned
  raster shares its dims, extent, resolution, and dim order.

# Keywords

- `component` — selector identifying the SPDE component to project.
  Accepts `Int` (1-based index), `String` (matching the
  `random_effects` key, e.g. `"SPDE2[3]"`), or `Type{<:SPDE2}` (auto-
  locates the unique SPDE2 component, errors if there are zero or
  multiple).
- `quantity = :mean` — which posterior summary to project: one of
  `:mean`, `:sd`, `:lower`, `:upper`.
- `level = 0.95` — credible level forwarded to `random_effects` for
  `:lower` / `:upper`.
- `outside`, `missingval` — forwarded to the vertex-vector
  `predict_raster` form.
- `mesh_crs = nothing` — optional CRS assertion (ADR-041). When
  supplied, must equal `Rasters.crs(template)` or an `ArgumentError`
  is raised. Default `nothing` preserves the v0.2.x "trust the
  caller" behaviour.

# Returns

A `Raster` matching `template`'s dims with the projected per-pixel
posterior summary.
"""
function predict_raster(
        model::LatentGaussianModel,
        res::INLAResult,
        template::Raster;
        component,
        quantity::Symbol=:mean,
        level::Real=0.95,
        outside::Symbol=:missing,
        missingval::Real=NaN,
        mesh_crs=nothing
)
    quantity ∈ (:mean, :sd, :lower, :upper) || throw(ArgumentError(
        "predict_raster: quantity must be one of (:mean, :sd, :lower, :upper); got $(quantity)"
    ))
    _check_crs(Rasters.crs(template), mesh_crs)

    i, mesh = _resolve_spde_component(model, component)
    re = random_effects(model, res; level=level)
    name = LatentGaussianModels._component_name(model.components[i], i)
    haskey(re, name) || throw(ArgumentError(
        "predict_raster: component $(name) is missing from random_effects; " *
        "this typically means the component has length 1 (e.g. Intercept). " *
        "SPDE components are vector-valued and should always appear."
    ))
    nt = re[name]
    values = getfield(nt, quantity)
    return predict_raster(values, mesh, template;
        outside=outside, missingval=missingval)
end

"""
    predict_raster(rng::AbstractRNG, model::LatentGaussianModel,
                   res::INLAResult, template::Raster;
                   component, quantity,
                   n_samples::Integer = 1000,
                   outside::Symbol = :missing,
                   missingval::Real = NaN,
                   mesh_crs = nothing) -> Raster

Sample-based projection of an SPDE component onto a raster grid. Draws
`n_samples` joint posterior samples via `LatentGaussianModels.posterior_sample`,
slices the SPDE block out of each draw, projects all draws through the
barycentric mesh→raster projector in a single sparse-dense GEMM, then
reduces the resulting `n_cells × n_samples` matrix column-wise per
`quantity`.

# Quantity dispatch

- `quantity = :mean` — per-cell posterior mean. Asymptotically agrees
  with the Gaussian-approximation `predict_raster(model, res, template;
  quantity = :mean)` overload up to MC error.
- `quantity isa Real` with `0 ≤ quantity ≤ 1` — per-cell empirical
  quantile (`Statistics.quantile`). E.g. `quantity = 0.5` is the
  posterior median; `quantity = 0.025` / `0.975` are the 95% credible
  bounds.
- `quantity isa Exceedance(c)` — per-cell exceedance probability
  `P(u(s) > c | y) ≈ mean(η_samples .> c, dims = 2)`. The honest path
  for tail probabilities; cannot be derived from Gaussian-approximation
  summaries.

# Performance

`P.A * x_samples[spde_range, :]` is one sparse-dense GEMM yielding the
full `n_cells × n_samples` projected block. All reductions are
column-wise scans of this block. Memory cost is `n_cells × n_samples`
Float64s; pick `n_samples` accordingly for very large rasters.

# Reproducibility

Pass a seeded `rng` (e.g. `Xoshiro(1234)`) to obtain bit-wise
identical rasters across calls.

See also: the `predict_raster(model, res, template; …)` overload for the
Gaussian-approximation overload (cheaper, no sampling, but cannot
produce exceedance rasters).
"""
function predict_raster(
        rng::Random.AbstractRNG,
        model::LatentGaussianModel,
        res::INLAResult,
        template::Raster;
        component,
        quantity,
        n_samples::Integer=1000,
        outside::Symbol=:missing,
        missingval::Real=NaN,
        mesh_crs=nothing
)
    n_samples >= 1 ||
        throw(ArgumentError("predict_raster: n_samples must be ≥ 1; got $(n_samples)"))
    outside ∈ (:error, :missing) ||
        throw(ArgumentError("predict_raster: outside must be :error or :missing; got $(outside)"))
    _validate_sample_quantity(quantity)
    _check_crs(Rasters.crs(template), mesh_crs)

    i, mesh = _resolve_spde_component(model, component)

    xs = collect(Rasters.lookup(template, X))
    ys = collect(Rasters.lookup(template, Y))
    nx = length(xs)
    ny = length(ys)
    n_cells = nx * ny

    locs = Matrix{Float64}(undef, n_cells, 2)
    cell_ij = Vector{Tuple{Int, Int}}(undef, n_cells)
    k = 0
    for j in 1:ny
        for ii in 1:nx
            k += 1
            locs[k, 1] = xs[ii]
            locs[k, 2] = ys[j]
            cell_ij[k] = (ii, j)
        end
    end

    proj_outside = outside === :error ? :error : :zero
    P = INLASPDE.MeshProjector(mesh, locs; outside=proj_outside)

    rng_range = LatentGaussianModels.component_range(model, i)
    draws = LatentGaussianModels.posterior_sample(rng, res, model;
        n_samples=n_samples)
    x_block = draws.x[rng_range, :]
    η_samples = P.A * x_block
    cell_values = _sample_reduce(η_samples, quantity)

    out = similar(template, Float64)
    fill!(out, missingval)
    @inbounds for k in 1:n_cells
        ii, jj = cell_ij[k]
        row = @view P.A[k, :]
        if !iszero(row)
            out[X=ii, Y=jj] = cell_values[k]
        end
    end
    return out
end

_validate_sample_quantity(q::Symbol) = q === :mean || throw(ArgumentError(
    "predict_raster: quantity Symbol must be :mean (sample-based path); " *
    "got $(q). Use the Gaussian-approximation overload for :sd / :lower / :upper."
))
_validate_sample_quantity(q::Real) = (zero(q) <= q <= oneunit(q)) ||
    throw(ArgumentError(
        "predict_raster: quantile quantity must be in [0, 1]; got $(q)"
    ))
_validate_sample_quantity(::Exceedance) = nothing
_validate_sample_quantity(q) = throw(ArgumentError(
    "predict_raster: quantity must be :mean, a Real in [0, 1] (quantile), " *
    "or an Exceedance; got $(typeof(q))"
))

_sample_reduce(η::AbstractMatrix, ::Symbol) = vec(Statistics.mean(η; dims=2))
function _sample_reduce(η::AbstractMatrix, q::Real)
    return [Statistics.quantile(view(η, k, :), q) for k in axes(η, 1)]
end
_sample_reduce(η::AbstractMatrix, q::Exceedance) =
    vec(Statistics.mean(η .> q.c; dims=2))
