"""
    INLAMesh1D{T}

A 1D segment mesh in INLA-SPDE form: a sorted vector of vertex
coordinates plus the implicit segment topology (each row of
`segments` connects consecutive vertices). Suitable for the 1D
[`FEMMatrices`](@ref) overload and [`SPDE1D`](@ref).

# Fields
- `points::Vector{T}` — vertex coordinates in ascending order.
- `segments::Matrix{Int}` — `(n - 1) × 2` one-based vertex indices,
  `segments[k, :] == [k, k + 1]`.
- `boundary::NTuple{2, Int}` — endpoint vertex indices `(1, n)`.

Construct via [`inla_mesh_1d`](@ref).
"""
struct INLAMesh1D{T <: Real}
    points::Vector{T}
    segments::Matrix{Int}
    boundary::NTuple{2, Int}
end

function Base.show(io::IO, mesh::INLAMesh1D)
    n = length(mesh.points)
    print(
        io,
        "INLAMesh1D(",
        n,
        " vertices, ",
        n - 1,
        " segments, [",
        first(mesh.points),
        ", ",
        last(mesh.points),
        "])"
    )
end

num_vertices(mesh::INLAMesh1D) = length(mesh.points)
num_segments(mesh::INLAMesh1D) = size(mesh.segments, 1)

"""
    inla_mesh_1d(loc; max_edge, cutoff = 0.0, boundary = nothing) -> INLAMesh1D

Build a 1D segment mesh covering the interval implied by `loc` and
`boundary`, refined so every segment is no longer than `max_edge`.
The 1D analogue of [`inla_mesh_2d`](@ref).

# Arguments
- `loc::AbstractVector{<:Real}` — observation locations on the line.
  Need not be sorted; duplicates and points within `cutoff` of each
  other are collapsed (first-wins).
- `max_edge::Real` — upper bound on segment length. Required.
- `cutoff::Real = 0` — minimum allowed spacing; closer points are
  merged.
- `boundary::Union{Nothing, NTuple{2, <:Real}} = nothing` — explicit
  domain endpoints `(a, b)` with `a < b`. Must contain all `loc`
  values. If omitted, the domain defaults to `(minimum(loc),
  maximum(loc))`.

# Returns
An [`INLAMesh1D`](@ref) with sorted vertex coordinates and
consecutive-segment topology.

# Example

```julia
loc  = randn(50)
mesh = inla_mesh_1d(loc; max_edge = 0.1, cutoff = 0.01,
                    boundary = (-3.0, 3.0))
fem  = FEMMatrices(mesh)
```
"""
function inla_mesh_1d(
        loc::AbstractVector{<:Real};
        max_edge::Real,
        cutoff::Real=0.0,
        boundary::Union{Nothing, NTuple{2, <:Real}}=nothing
)
    max_edge > 0 ||
        throw(ArgumentError("inla_mesh_1d: max_edge must be positive; got $max_edge"))
    cutoff >= 0 ||
        throw(ArgumentError("inla_mesh_1d: cutoff must be non-negative; got $cutoff"))
    isempty(loc) &&
        throw(ArgumentError("inla_mesh_1d: loc must be non-empty"))

    T = float(eltype(loc))
    sorted = sort(collect(T, loc))

    if boundary !== nothing
        a, b = T(boundary[1]), T(boundary[2])
        a < b ||
            throw(ArgumentError("inla_mesh_1d: boundary endpoints must satisfy a < b; got ($a, $b)"))
        first(sorted) >= a ||
            throw(ArgumentError("inla_mesh_1d: loc contains a point $(first(sorted)) below boundary[1]=$a"))
        last(sorted) <= b ||
            throw(ArgumentError("inla_mesh_1d: loc contains a point $(last(sorted)) above boundary[2]=$b"))
        sorted = T[a; sorted; b]
        sort!(sorted)
    end

    deduped = _dedup_sorted(sorted, T(cutoff))

    refined = _refine_to_max_edge(deduped, T(max_edge))
    n = length(refined)
    segments = Matrix{Int}(undef, n - 1, 2)
    @inbounds for k in 1:(n - 1)
        segments[k, 1] = k
        segments[k, 2] = k + 1
    end
    return INLAMesh1D{T}(refined, segments, (1, n))
end

# Collapse points within `cutoff` of the previous kept point. Input
# must already be sorted ascending.
function _dedup_sorted(sorted::AbstractVector{T}, cutoff::T) where {T <: Real}
    isempty(sorted) && return T[]
    cutoff == 0 && return _dedup_exact(sorted)
    out = T[first(sorted)]
    for x in @view sorted[2:end]
        if x - last(out) > cutoff
            push!(out, x)
        end
    end
    return out
end

function _dedup_exact(sorted::AbstractVector{T}) where {T <: Real}
    out = T[first(sorted)]
    for x in @view sorted[2:end]
        if x != last(out)
            push!(out, x)
        end
    end
    return out
end

# Subdivide each gap > `max_edge` into equal-length sub-segments
# inserting `ceil(gap / max_edge) - 1` interior points.
function _refine_to_max_edge(sorted::AbstractVector{T}, max_edge::T) where {T <: Real}
    length(sorted) == 1 && return collect(T, sorted)
    out = T[first(sorted)]
    for k in 1:(length(sorted) - 1)
        a = sorted[k]
        b = sorted[k + 1]
        gap = b - a
        n_sub = max(1, Int(ceil(gap / max_edge)))
        if n_sub == 1
            push!(out, b)
        else
            step = gap / n_sub
            for j in 1:(n_sub - 1)
                push!(out, a + j * step)
            end
            push!(out, b)
        end
    end
    return out
end

"""
    FEMMatrices(mesh::INLAMesh1D)

Convenience overload: assemble FEM matrices from a 1D mesh.
"""
FEMMatrices(mesh::INLAMesh1D) = FEMMatrices(mesh.points, mesh.segments)
