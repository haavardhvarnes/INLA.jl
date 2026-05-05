"""
    MeshProjector1D{T, S}

Sparse linear-interpolation projector from a 1D segment mesh to
observation locations on the line. For a field `u` defined on mesh
vertices, the interpolated value at observation `i` is `(A * u)[i]`,
where each row of `A` has at most two nonzeros — the 1D barycentric
weights of observation `i` in its enclosing segment, summing to 1.

# Fields
- `mesh::INLAMesh1D` — the source mesh.
- `locations::Vector{T}` — observation coordinates.
- `A::SparseMatrixCSC{T, Int}` — `n_obs × n_vertices` sparse projector.

Construct via [`MeshProjector1D(mesh, locations)`](@ref).
"""
struct MeshProjector1D{T <: Real, M <: INLAMesh1D, S <: AbstractSparseMatrix{T, Int}}
    mesh::M
    locations::Vector{T}
    A::S
end

"""
    MeshProjector1D(mesh, locations; outside = :error, atol = 0.0) -> MeshProjector1D

Build a sparse linear-interpolation projector from `mesh` to the
points in `locations`.

# Keyword arguments

- `outside::Symbol = :error` — policy for locations falling outside
  `[first(mesh.points), last(mesh.points)]`:
    - `:error` — throw `ArgumentError` on the first outside point.
    - `:zero` — emit an empty row.
- `atol::Real = 0.0` — extra tolerance when classifying as outside;
  a location within `atol` of either endpoint is accepted and
  clamped.
"""
function MeshProjector1D(
        mesh::INLAMesh1D,
        locations::AbstractVector{<:Real};
        outside::Symbol=:error,
        atol::Real=0.0
)
    outside in (:error, :zero) ||
        throw(ArgumentError("outside must be :error or :zero; got $outside"))
    atol >= 0 ||
        throw(ArgumentError("atol must be non-negative; got $atol"))

    T = promote_type(eltype(mesh.points), float(eltype(locations)))
    locs = Vector{T}(locations)
    n_obs = length(locs)
    n_v = length(mesh.points)

    Is = Int[]
    Js = Int[]
    Vs = T[]
    sizehint!(Is, 2 * n_obs)
    sizehint!(Js, 2 * n_obs)
    sizehint!(Vs, 2 * n_obs)

    pts = mesh.points
    lo = first(pts)
    hi = last(pts)

    for i in 1:n_obs
        x = locs[i]
        if x < lo - atol || x > hi + atol
            outside === :error &&
                throw(ArgumentError("location $i at $x is outside the mesh domain [$lo, $hi]"))
            continue
        end
        x_clamped = clamp(x, lo, hi)
        k = searchsortedlast(pts, x_clamped)
        # k is the largest index with pts[k] ≤ x_clamped. Map endpoint
        # to its left segment so we always have segment [k, k+1].
        if k == n_v
            k = n_v - 1
        end
        x_left = pts[k]
        x_right = pts[k + 1]
        h = x_right - x_left
        λ_right = (x_clamped - x_left) / h
        λ_left = one(T) - λ_right

        push!(Is, i); push!(Js, k);     push!(Vs, λ_left)
        push!(Is, i); push!(Js, k + 1); push!(Vs, λ_right)
    end

    A = sparse(Is, Js, Vs, n_obs, n_v)
    return MeshProjector1D{T, typeof(mesh), typeof(A)}(mesh, locs, A)
end

# Matrix-like interface ----------------------------------------

Base.size(P::MeshProjector1D) = size(P.A)
Base.size(P::MeshProjector1D, d::Integer) = size(P.A, d)
Base.eltype(::Type{MeshProjector1D{T, M, S}}) where {T, M, S} = T
SparseArrays.sparse(P::MeshProjector1D) = P.A

Base.:*(P::MeshProjector1D, x::AbstractVector) = P.A * x
Base.:*(P::MeshProjector1D, X::AbstractMatrix) = P.A * X

function LinearAlgebra.mul!(y::AbstractVecOrMat, P::MeshProjector1D, x::AbstractVecOrMat)
    mul!(y, P.A, x)
end
function LinearAlgebra.mul!(
        y::AbstractVecOrMat, P::MeshProjector1D, x::AbstractVecOrMat,
        α::Number, β::Number
)
    mul!(y, P.A, x, α, β)
end

function Base.show(io::IO, P::MeshProjector1D)
    n_obs, n_v = size(P.A)
    return print(io, "MeshProjector1D(", n_obs, " locations → ", n_v, " vertices)")
end
