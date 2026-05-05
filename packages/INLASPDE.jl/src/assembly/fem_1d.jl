"""
    assemble_fem_matrices_1d(points, segments) -> (C, G1)

Assemble the P1 finite-element mass matrix `C` and stiffness matrix
`G₁` on a 1D segment mesh.

The mass matrix has entries `C[i,j] = ∫_Ω φ_i φ_j dx` and the
stiffness matrix `G₁[i,j] = ∫_Ω (dφ_i/dx) · (dφ_j/dx) dx`, where
`φ_i` is the piecewise-linear hat basis function at vertex `i`.

# Arguments
- `points::AbstractVector{<:Real}` — vertex coordinates (length `n`).
- `segments::AbstractMatrix{<:Integer}` — `m × 2` one-based vertex
  indices, one row per segment. Orientation is irrelevant; segment
  length `|x_{i_2} - x_{i_1}|` is used.

# Returns
- `C::SparseMatrixCSC` — `n × n` symmetric positive-definite mass
  matrix.
- `G1::SparseMatrixCSC` — `n × n` symmetric positive-semidefinite
  stiffness matrix. `G₁ · 1 = 0` (constant-preserving).

# Element-level formulas

For a segment of length `h` with endpoint indices `(i, j)`:

- mass:        `C_loc = (h/6) · [2 1; 1 2]`
- stiffness:   `G_loc = (1/h) · [1 -1; -1 1]`
"""
function assemble_fem_matrices_1d(
        points::AbstractVector{<:Real},
        segments::AbstractMatrix{<:Integer}
)
    size(segments, 2) == 2 ||
        throw(ArgumentError("segments must be m × 2 for P1 elements; got size $(size(segments))"))

    n_vertices = length(points)
    n_segments = size(segments, 1)
    T = float(eltype(points))

    nnz_upper_bound = 4 * n_segments
    Is = Vector{Int}(undef, nnz_upper_bound)
    Js = Vector{Int}(undef, nnz_upper_bound)
    Vc = Vector{T}(undef, nnz_upper_bound)
    Vg = Vector{T}(undef, nnz_upper_bound)

    k = 0
    for s in axes(segments, 1)
        i1 = Int(segments[s, 1])
        i2 = Int(segments[s, 2])
        (1 <= i1 <= n_vertices && 1 <= i2 <= n_vertices) ||
            throw(ArgumentError("segment $s references vertex out of range 1:$n_vertices"))
        i1 == i2 &&
            throw(ArgumentError("segment $s is degenerate (i1 == i2)"))

        h = abs(T(points[i2]) - T(points[i1]))
        h > 0 ||
            throw(ArgumentError("segment $s has zero length"))

        idx = (i1, i2)
        # Local mass: (h/6) [2 1; 1 2]
        mass_diag = h / 3
        mass_off = h / 6
        # Local stiffness: (1/h) [1 -1; -1 1]
        stiff_diag = inv(h)
        stiff_off = -inv(h)

        for a in 1:2, b in 1:2
            k += 1
            Is[k] = idx[a]
            Js[k] = idx[b]
            Vc[k] = a == b ? mass_diag : mass_off
            Vg[k] = a == b ? stiff_diag : stiff_off
        end
    end

    C = sparse(Is, Js, Vc, n_vertices, n_vertices)
    G1 = sparse(Is, Js, Vg, n_vertices, n_vertices)
    return C, G1
end

"""
    FEMMatrices(points::AbstractVector, segments::AbstractMatrix)

1D variant: assemble `C`, `G₁`, `C̃`, and `G₂` from a 1D segment mesh
given as raw arrays. See [`assemble_fem_matrices_1d`](@ref) for the
argument conventions.

For α = 2 in 1D the smoothness is `ν = 3/2` (well-defined under PC-Matérn);
the same `G₂ = G₁ · C̃⁻¹ · G₁` is used as in 2D.
"""
function FEMMatrices(
        points::AbstractVector{<:Real},
        segments::AbstractMatrix{<:Integer}
)
    C, G1 = assemble_fem_matrices_1d(points, segments)
    C_lumped = lumped_mass(C)
    G2 = stiffness_squared(G1, C_lumped)
    return FEMMatrices(C, G1, C_lumped, G2)
end
