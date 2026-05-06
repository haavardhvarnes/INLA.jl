"""
    LGMFormulaINLASPDEExt

Weakdep extension that teaches `@lgm` how to expand tuple-coordinate
`f((east, north), spde::SPDE2)` terms into a barycentric design-matrix
block via [`MeshProjector`](@ref). Loaded automatically when both
`LGMFormula` and `INLASPDE` are imported (Phase N PR-7b / ADR-039).
"""
module LGMFormulaINLASPDEExt

using LGMFormula: LGMFormula
using INLASPDE: INLASPDE, SPDE2, MeshProjector, INLAMesh
using SparseArrays: SparseMatrixCSC
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
    P = MeshProjector(spde.mesh, locs)
    return SparseMatrixCSC{Float64, Int}(P.A)
end

end # module
