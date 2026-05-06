"""
    LGMFormulaINLASPDEExt

Weakdep extension that teaches `@lgm` how to expand tuple-coordinate
`f((east, north), spde::SPDE2)` and
`f((east, north, time), KroneckerComponent(spde, time_comp))` terms
into a barycentric design-matrix block via [`MeshProjector`](@ref).
Loaded automatically when both `LGMFormula` and `INLASPDE` are imported
(Phase N PR-7b / PR-7c, ADR-039).
"""
module LGMFormulaINLASPDEExt

using LGMFormula: LGMFormula
using INLASPDE: INLASPDE, SPDE2, MeshProjector, INLAMesh
using LatentGaussianModels: KroneckerComponent
using SparseArrays: SparseMatrixCSC, sparse, findnz
using Tables: Tables

function LGMFormula._build_spatial_block(
        spde::SPDE2,
        data_cols,
        coord_cols::Tuple{Symbol, Symbol},
        n_obs::Integer
)
    spde.mesh isa INLAMesh || throw(ArgumentError(
        "@lgm: f((east, north), spde): the `SPDE2` component was constructed via " *
        "`SPDE2(points, triangles; …)` and does not retain its mesh. Use " *
        "`SPDE2(mesh; …)` (where `mesh = inla_mesh_2d(...)`) so the macro " *
        "can build a barycentric projector at runtime."
    ))
    locs = _gather_spatial_locs(data_cols, coord_cols, n_obs)
    P = MeshProjector(spde.mesh, locs)
    return SparseMatrixCSC{Float64, Int}(P.A)
end

# `f((east, north, time), KroneckerComponent(spde, time_comp))` — ADR-038.
# Builds a sparse Khatri-Rao (row-product) design matrix with column layout
# `(s - 1) · n_t + t`, matching `KroneckerComponent`'s
# `vec(X)` (`X` of shape `(n_t × n_s)`) and `KroneckerMapping`'s flattening
# convention. Per-obs row `i` carries the spatial barycentric weights of
# observation `i` shifted to the `time[i]` slot of the temporal axis.
#
# Note (deviation from plans/phase-n-pr7.md §"PR-7c"): the plan conjectured
# pure-Kron emit `KroneckerMapping(A_space, A_time)` and a mixed RHS
# `StackedMapping(...)`. That is mis-typed for arbitrary per-obs data:
# `KroneckerMapping` is matrix-Kron with `nrows = nrows(A_s) · nrows(A_t)`,
# so feeding both factors at `n_obs` rows would yield `n_obs²` observation
# rows. `StackedMapping` is row-stacked, not column-stacked. The Cameletti
# oracle ([test_cameletti_pm10.jl:90](packages/INLASPDE.jl/test/oracle/test_cameletti_pm10.jl))
# already builds `kron(A_space_j, I(n_t))` as a `SparseMatrixCSC`, which is
# exactly the gridded-data special case of this Khatri-Rao construction.
function LGMFormula._build_spatial_block(
        c::KroneckerComponent,
        data_cols,
        coord_cols::Tuple{Symbol, Symbol, Symbol},
        n_obs::Integer
)
    c.space isa SPDE2 || throw(ArgumentError(
        "@lgm: f((east, north, time), KroneckerComponent(space, time)): the " *
        "spatial child must be an `SPDE2` component, got " *
        "$(typeof(c.space)). Phase N PR-7c only supports `SPDE2` as the " *
        "spatial factor; non-SPDE Kronecker compositions can be built " *
        "manually via the explicit `LatentGaussianModel(...)` constructor."
    ))
    spde = c.space::SPDE2
    spde.mesh isa INLAMesh || throw(ArgumentError(
        "@lgm: f((east, north, time), KroneckerComponent(spde, …)): the `SPDE2` " *
        "spatial child was constructed via `SPDE2(points, triangles; …)` and " *
        "does not retain its mesh. Use `SPDE2(mesh; …)` (where " *
        "`mesh = inla_mesh_2d(...)`) so the macro can build a barycentric " *
        "projector at runtime."
    ))
    east_sym, north_sym, time_sym = coord_cols
    locs = _gather_spatial_locs(data_cols, (east_sym, north_sym), n_obs)
    n_t = length(c.time)
    time_idx = _gather_time_index(data_cols, time_sym, n_t, n_obs)

    P = MeshProjector(spde.mesh, locs)
    A_space = SparseMatrixCSC{Float64, Int}(P.A)
    n_v = size(A_space, 2)

    rows, cols, vals = findnz(A_space)
    new_cols = similar(cols)
    @inbounds for k in eachindex(rows)
        s = cols[k]
        t = time_idx[rows[k]]
        new_cols[k] = (s - 1) * n_t + t
    end
    return sparse(rows, new_cols, vals, n_obs, n_v * n_t)
end

# Mis-pair error: `f((east, north, time), spde::SPDE2)` — point at the
# `KroneckerComponent` wrapper.
function LGMFormula._build_spatial_block(
        spde::SPDE2,
        data_cols,
        coord_cols::Tuple{Symbol, Symbol, Symbol},
        n_obs::Integer
)
    throw(ArgumentError(
        "@lgm: f((east, north, time), SPDE2(...)): a 3-tuple coordinate " *
        "(space + time) requires a `KroneckerComponent` wrapper. Use " *
        "`f((east, north, time), KroneckerComponent(spde, time_comp))` — " *
        "for example `KroneckerComponent(SPDE2(mesh), AR1(n_t))` for the " *
        "Cameletti-style separable space-time SPDE."
    ))
end

# Mis-pair error: `f((east, north), KroneckerComponent(...))` — point at the
# 3-tuple form.
function LGMFormula._build_spatial_block(
        c::KroneckerComponent,
        data_cols,
        coord_cols::Tuple{Symbol, Symbol},
        n_obs::Integer
)
    throw(ArgumentError(
        "@lgm: f((east, north), KroneckerComponent(...)): a `KroneckerComponent` " *
        "needs a 3-tuple coordinate `(east, north, time)` so the macro can " *
        "build the space-time design block. Use " *
        "`f((east, north, time), KroneckerComponent(spde, time_comp))`, or " *
        "drop the wrapper and pass just the spatial child for a pure-spatial " *
        "term."
    ))
end

function _gather_spatial_locs(data_cols, coord_cols::Tuple{Symbol, Symbol},
        n_obs::Integer)
    east_sym, north_sym = coord_cols
    e_raw = Tables.getcolumn(data_cols, east_sym)
    n_raw = Tables.getcolumn(data_cols, north_sym)
    length(e_raw) == n_obs || throw(DimensionMismatch(
        "@lgm: spatial-coord column `$(east_sym)` has length $(length(e_raw)); expected $(n_obs)"
    ))
    length(n_raw) == n_obs || throw(DimensionMismatch(
        "@lgm: spatial-coord column `$(north_sym)` has length $(length(n_raw)); expected $(n_obs)"
    ))
    locs = Matrix{Float64}(undef, n_obs, 2)
    @inbounds for i in 1:n_obs
        locs[i, 1] = Float64(e_raw[i])
        locs[i, 2] = Float64(n_raw[i])
    end
    return locs
end

function _gather_time_index(data_cols, time_sym::Symbol, n_t::Int,
        n_obs::Integer)
    t_raw = Tables.getcolumn(data_cols, time_sym)
    length(t_raw) == n_obs || throw(DimensionMismatch(
        "@lgm: time-coord column `$(time_sym)` has length $(length(t_raw)); expected $(n_obs)"
    ))
    out = Vector{Int}(undef, n_obs)
    @inbounds for i in 1:n_obs
        v = t_raw[i]
        v isa Integer || throw(ArgumentError(
            "@lgm: time-coord column `$(time_sym)` must contain integers in 1:$(n_t); " *
            "got entry of type $(typeof(v)) at row $(i)"
        ))
        (1 ≤ v ≤ n_t) || throw(ArgumentError(
            "@lgm: time-coord column `$(time_sym)` value $(v) at row $(i) is " *
            "outside 1:$(n_t) (temporal child has $(n_t) levels)"
        ))
        out[i] = Int(v)
    end
    return out
end

end # module
