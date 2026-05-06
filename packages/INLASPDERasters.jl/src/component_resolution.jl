"""
    _resolve_spde_component(model::LatentGaussianModel, id; time_index = nothing)
        -> (i::Int, mesh::INLAMesh, slice::AbstractRange{Int})

Locate a raster-projectable spatial component inside `model.components`
and return its 1-based index, the retained `INLAMesh`, and the
component-local index range that picks out the spatial slice the
projector consumes. Used by the `predict_raster(model, res, template)`
overloads to bridge a `random_effects` block (SPDE-flavored vertex
vector or Kronecker-structured space-time block) to a barycentric
mesh→raster projector.

# Accepted spatial components

- `SPDE2` and `SPDE2NonStationary` constructed via the mesh-aware
  constructor (so the `mesh` field is populated per ADR-036). Slice is
  `1:n_v` (identity).
- `KroneckerComponent(spatial, temporal)` where `spatial` is one of the
  two SPDE flavours. Requires `time_index::Integer ∈ 1:length(temporal)`;
  slice is `time_index:n_t:n_v · n_t`, picking the spatial vector at
  the requested time slot under `KroneckerMapping`'s
  time-inner-index flattening (`x[(s - 1) · n_t + t]`).

# Accepted `id` shapes

- `Int` — 1-based component index.
- `String` — matches `LatentGaussianModels._component_name(c, i)`
  (e.g. `"SPDE2[3]"`, `"SPDE2NonStationary[2]"`,
  `"KroneckerComponent[1]"`).
- `Type{<:SPDE2}`, `Type{<:SPDE2NonStationary}`, or
  `Type{<:KroneckerComponent}` — auto-locate the unique component of
  that type; throws if zero or multiple match.

# Errors

Raises `ArgumentError` with a user-readable enumeration of the model's
components when the id cannot be resolved or the resolved component is
not raster-projectable. KroneckerComponent paths additionally require
`time_index`; non-Kronecker components reject `time_index ≠ nothing`.
"""
function _resolve_spde_component(
        model::LatentGaussianModel, id::Integer; time_index=nothing)
    n = length(model.components)
    (1 ≤ id ≤ n) || throw(ArgumentError(
        "predict_raster: component index $(id) is out of range; model has " *
        "$(n) components ($(_components_summary(model)))"
    ))
    return _check_spde_component(model, Int(id), time_index)
end

function _resolve_spde_component(
        model::LatentGaussianModel, id::AbstractString; time_index=nothing)
    for (i, c) in enumerate(model.components)
        LatentGaussianModels._component_name(c, i) == id || continue
        return _check_spde_component(model, i, time_index)
    end
    throw(ArgumentError(
        "predict_raster: no component named $(repr(id)); known names: " *
        "[$(join(string.([LatentGaussianModels._component_name(c, i)
                         for (i, c) in enumerate(model.components)]), ", "))]"
    ))
end

function _resolve_spde_component(
        model::LatentGaussianModel, ::Type{T}; time_index=nothing) where {T}
    matches = Int[]
    for (i, c) in enumerate(model.components)
        c isa T && push!(matches, i)
    end
    isempty(matches) && throw(ArgumentError(
        "predict_raster: no $(T) component found; model components are " *
        "($(_components_summary(model)))"
    ))
    length(matches) == 1 || throw(ArgumentError(
        "predict_raster: $(length(matches)) $(T) components found at indices " *
        "$(matches); pass an explicit `Int` index or component name to disambiguate"
    ))
    return _check_spde_component(model, matches[1], time_index)
end

function _check_spde_component(model::LatentGaussianModel, i::Int, time_index)
    return _resolve_component_geometry(model.components[i], i, time_index)
end

function _resolve_component_geometry(c::SPDE2, i::Int, time_index)
    time_index === nothing || throw(ArgumentError(
        "predict_raster: time_index is only meaningful for KroneckerComponent " *
        "with an SPDE-flavored spatial child; got time_index=$(time_index) on " *
        "an SPDE2 component at index $(i)"
    ))
    c.mesh isa INLAMesh || throw(ArgumentError(
        "predict_raster: SPDE2 component at index $(i) was constructed via " *
        "`SPDE2(points, triangles; …)` and does not retain its mesh. Use " *
        "`SPDE2(mesh; …)` (where `mesh = inla_mesh_2d(...)`) so that " *
        "predict_raster can reuse the mesh for barycentric projection."
    ))
    n_v = length(c)
    return i, c.mesh::INLAMesh, 1:1:n_v
end

function _resolve_component_geometry(c::SPDE2NonStationary, i::Int, time_index)
    time_index === nothing || throw(ArgumentError(
        "predict_raster: time_index is only meaningful for KroneckerComponent " *
        "with an SPDE-flavored spatial child; got time_index=$(time_index) on " *
        "an SPDE2NonStationary component at index $(i)"
    ))
    c.mesh isa INLAMesh || throw(ArgumentError(
        "predict_raster: SPDE2NonStationary component at index $(i) was " *
        "constructed via `SPDE2NonStationary(points, triangles; …)` and does " *
        "not retain its mesh. Use `SPDE2NonStationary(mesh; …)` (where " *
        "`mesh = inla_mesh_2d(...)`) so that predict_raster can reuse the " *
        "mesh for barycentric projection."
    ))
    n_v = length(c)
    return i, c.mesh::INLAMesh, 1:1:n_v
end

function _resolve_component_geometry(
        c::LatentGaussianModels.KroneckerComponent, i::Int, time_index)
    spatial = c.space
    spatial isa Union{SPDE2, SPDE2NonStationary} || throw(ArgumentError(
        "predict_raster: KroneckerComponent at index $(i) has spatial child of " *
        "type $(typeof(spatial).name.name); raster prediction requires SPDE2 " *
        "or SPDE2NonStationary as the spatial child."
    ))
    time_index === nothing && throw(ArgumentError(
        "predict_raster: KroneckerComponent at index $(i) requires the " *
        "`time_index` keyword to specify which time slot to project. Pass an " *
        "integer in 1:$(length(c.time))."
    ))
    n_t = length(c.time)
    n_v = length(spatial)
    (time_index isa Integer) || throw(ArgumentError(
        "predict_raster: time_index must be an Integer; got $(typeof(time_index))"
    ))
    (1 <= time_index <= n_t) || throw(ArgumentError(
        "predict_raster: time_index ($(time_index)) must be in 1:$(n_t)"
    ))
    spatial.mesh isa INLAMesh || throw(ArgumentError(
        "predict_raster: KroneckerComponent at index $(i) has a spatial child " *
        "that does not retain its mesh. Construct the spatial child with " *
        "`SPDE2(mesh; …)` or `SPDE2NonStationary(mesh; …)`."
    ))
    return i, spatial.mesh::INLAMesh,
    StepRange{Int, Int}(time_index, n_t, time_index + (n_v - 1) * n_t)
end

function _resolve_component_geometry(c, i::Int, time_index)
    throw(ArgumentError(
        "predict_raster: component $(i) is a $(typeof(c).name.name); raster " *
        "prediction requires SPDE2, SPDE2NonStationary, or " *
        "KroneckerComponent(SPDE2-flavored, time)."
    ))
end

function _components_summary(model::LatentGaussianModel)
    parts = String[]
    for (i, c) in enumerate(model.components)
        push!(parts, "[$(i)] $(LatentGaussianModels._component_name(c, i))")
    end
    return join(parts, ", ")
end

"""
    _check_crs(rs_crs, mesh_crs) -> Nothing

Assert `mesh_crs == rs_crs` when `mesh_crs !== nothing`. The keyword
defaults to `nothing` everywhere it appears — back-compat with v0.2.x
"trust the caller" CRS handling. When supplied, CRS mismatches at the
API boundary raise an informative error rather than producing silent
geographic foot-guns (ADR-041).
"""
function _check_crs(rs_crs, mesh_crs)
    mesh_crs === nothing && return nothing
    rs_crs == mesh_crs && return nothing
    throw(ArgumentError(
        "CRS mismatch: raster CRS is $(repr(rs_crs)) but mesh_crs was " *
        "supplied as $(repr(mesh_crs)). Pre-project one side or omit " *
        "`mesh_crs` to skip the check."
    ))
end
