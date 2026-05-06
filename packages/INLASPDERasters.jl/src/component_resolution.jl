"""
    _resolve_spde_component(model::LatentGaussianModel, id) -> (i::Int, mesh::INLAMesh)

Locate an SPDE2 component inside `model.components` and return its
1-based index plus the retained `INLAMesh`. Used by the
`predict_raster(model, res, template)` overload to bridge a
`random_effects` slice to a barycentric mesh→raster projector.

# Accepted `id` shapes

- `Int` — 1-based component index. Bounds-checked against
  `length(model.components)`.
- `String` — matches the component name produced by
  `LatentGaussianModels._component_name(c, i)` (the same key that
  `random_effects(model, res)` uses, e.g. `"SPDE2[3]"`).
- `Type{<:SPDE2}` — auto-locate the unique SPDE2 component; throws if
  zero or two-or-more components match.

# Errors

Raises `ArgumentError` with a user-readable enumeration of the model's
components when the id cannot be resolved or the resolved component is
not an `SPDE2` carrying a retained mesh (i.e. constructed via
`SPDE2(mesh::INLAMesh; …)` rather than the back-compat
`SPDE2(points, triangles; …)` path; ADR-036).
"""
function _resolve_spde_component(model::LatentGaussianModel, id::Integer)
    n = length(model.components)
    (1 ≤ id ≤ n) || throw(ArgumentError(
        "predict_raster: component index $(id) is out of range; model has " *
        "$(n) components ($(_components_summary(model)))"
    ))
    return _check_spde_component(model, Int(id))
end

function _resolve_spde_component(model::LatentGaussianModel, id::AbstractString)
    for (i, c) in enumerate(model.components)
        LatentGaussianModels._component_name(c, i) == id || continue
        return _check_spde_component(model, i)
    end
    throw(ArgumentError(
        "predict_raster: no component named $(repr(id)); known names: " *
        "[$(join(string.([LatentGaussianModels._component_name(c, i)
                         for (i, c) in enumerate(model.components)]), ", "))]"
    ))
end

function _resolve_spde_component(model::LatentGaussianModel, ::Type{T}) where {T <: SPDE2}
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
    return _check_spde_component(model, matches[1])
end

function _check_spde_component(model::LatentGaussianModel, i::Int)
    c = model.components[i]
    c isa SPDE2 || throw(ArgumentError(
        "predict_raster: component $(i) is a $(typeof(c).name.name); " *
        "raster prediction requires an SPDE2 component"
    ))
    c.mesh isa INLAMesh || throw(ArgumentError(
        "predict_raster: SPDE2 component at index $(i) was constructed via " *
        "`SPDE2(points, triangles; …)` and does not retain its mesh. Use " *
        "`SPDE2(mesh; …)` (where `mesh = inla_mesh_2d(...)`) so that " *
        "predict_raster can reuse the mesh for barycentric projection."
    ))
    return i, c.mesh::INLAMesh
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
