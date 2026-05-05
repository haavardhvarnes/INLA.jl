"""
    FEMMatrices{T, SC, SG}

Container for the precomputed FEM matrices `(C, G₁, C̃, G₂)` on a fixed
triangular mesh. These matrices depend only on mesh geometry, not on the
SPDE hyperparameters `(τ, κ)`, so they are assembled once and reused for
every precision evaluation.

# Fields
- `C::SC`        — full P1 mass matrix.
- `G1::SG`       — P1 stiffness matrix (`∫ ∇φ_i · ∇φ_j`).
- `C_lumped::SC` — diagonal mass matrix from `lumped_mass(C)`.
- `G2::SG`       — `G₁ · C̃⁻¹ · G₁`, used for α = 2.

Construct via `FEMMatrices(points, triangles)`.
"""
struct FEMMatrices{T, SC <: AbstractSparseMatrix{T}, SG <: AbstractSparseMatrix{T}}
    C::SC
    G1::SG
    C_lumped::SC
    G2::SG
end

"""
    FEMMatrices(points, triangles)

Assemble `C`, `G₁`, `C̃`, and `G₂` from a 2D triangular mesh given as raw
arrays. See [`assemble_fem_matrices`](@ref) for the argument conventions.
"""
function FEMMatrices(
        points::AbstractMatrix{<:Real},
        triangles::AbstractMatrix{<:Integer}
)
    C, G1 = assemble_fem_matrices(points, triangles)
    C_lumped = lumped_mass(C)
    G2 = stiffness_squared(G1, C_lumped)
    return FEMMatrices(C, G1, C_lumped, G2)
end

"""
    stiffness_squared(G1, C_lumped) -> G2

Construct `G₂ = G₁ · C̃⁻¹ · G₁`, where `C̃` is a lumped (diagonal) mass
matrix. This is the sparse approximation of `G₁ · C⁻¹ · G₁` used for α = 2
SPDE precision.

Throws `ArgumentError` if any diagonal entry of `C_lumped` is zero — this
indicates a vertex with zero associated area, i.e. a degenerate mesh.
"""
function stiffness_squared(G1::AbstractSparseMatrix, C_lumped::AbstractSparseMatrix)
    d = diag(C_lumped)
    any(iszero, d) &&
        throw(ArgumentError("C_lumped has a zero diagonal entry; mesh is degenerate"))
    D_inv = Diagonal(inv.(d))
    return G1 * D_inv * G1
end

"""
    spde_precision(fem::FEMMatrices, α, τ, κ) -> Q

Assemble the SPDE-Matérn precision matrix on user-scale parameters
`(τ, κ)`. Supported orders are `α ∈ {1, 2}`.

- α = 1: `Q = τ² · (κ² C + G₁)` (Matérn smoothness ν = 0).
- α = 2: `Q = τ² · (κ⁴ C̃ + 2κ² G₁ + G₂)` (ν = 1), using the lumped mass
  matrix `C̃` — this matches R-INLA's implementation per
  Lindgren-Rue-Lindström (2011, Appendix C).

Fractional α is deferred to v0.3 (Bolin–Kirchner 2020 rational
approximation).
"""
function spde_precision(fem::FEMMatrices, α::Integer, τ::Real, κ::Real)
    # Parametric (τ, κ) checks raise DomainError so the LGM safety net
    # in `_neg_log_posterior_θ` can recover when LBFGS overshoots into
    # `exp(very_negative_log_κ) → 0.0`. The structural `α ∈ {1, 2}` check
    # stays as ArgumentError because no θ-step can produce it (ADR-031).
    τ > 0 ||
        throw(DomainError(τ, "τ must be positive; got τ=$τ"))
    κ > 0 ||
        throw(DomainError(κ, "κ must be positive; got κ=$κ"))
    if α == 1
        return τ^2 * (κ^2 * fem.C + fem.G1)
    elseif α == 2
        return τ^2 * (κ^4 * fem.C_lumped + 2 * κ^2 * fem.G1 + fem.G2)
    end
    throw(ArgumentError("α must be 1 or 2; got α=$α. Fractional α deferred to v0.3."))
end

"""
    spde_precision_nonstationary(fem::FEMMatrices, α, τ_v, κ_v) -> Q

Sparse SPDE-Matérn precision with per-vertex `τ_v`, `κ_v` (vectors of
length `n_vertices`). Implements R-INLA's `inla.spde2.matern` formula:

- α = 2: `Q = D_τ · (D_κ²·C̃·D_κ² + D_κ²·G₁ + G₁·D_κ² + G₂) · D_τ`,
  where `D_τ = diag(τ_v)`, `D_κ² = diag(κ²_v)`, and `C̃` is the lumped
  (diagonal) mass matrix. `D_κ²·C̃·D_κ²` is just
  `Diagonal(κ⁴_v ⊙ c̃_v)`; the cross-terms `D_κ²·G₁ + G₁·D_κ²`
  symmetrise to `G₁_{ij}·(κ²_i + κ²_j)`.

In the stationary limit (`τ_v ≡ τ`, `κ_v ≡ κ`) this reduces to
`τ² · (κ⁴ C̃ + 2κ² G₁ + G₂)`, matching the `SPDE2` stationary
precision. Only `α = 2` is supported for the non-stationary form;
`α = 1` and fractional `α` are deferred (no R-INLA-parity oracle for
either yet).
"""
function spde_precision_nonstationary(
        fem::FEMMatrices, α::Integer,
        τ_v::AbstractVector{<:Real}, κ_v::AbstractVector{<:Real}
)
    n = size(fem.C, 1)
    length(τ_v) == n || throw(ArgumentError(
        "spde_precision_nonstationary: length(τ_v) ($(length(τ_v))) " *
        "must equal n_vertices ($n)"
    ))
    length(κ_v) == n || throw(ArgumentError(
        "spde_precision_nonstationary: length(κ_v) ($(length(κ_v))) " *
        "must equal n_vertices ($n)"
    ))
    all(>(0), τ_v) || throw(DomainError(
        τ_v, "spde_precision_nonstationary: all τ_v must be positive"
    ))
    all(>(0), κ_v) || throw(DomainError(
        κ_v, "spde_precision_nonstationary: all κ_v must be positive"
    ))
    α == 2 || throw(ArgumentError(
        "spde_precision_nonstationary: only α = 2 is supported in v0.2; " *
        "got α = $α. α = 1 and fractional α deferred."
    ))
    κ²_v = κ_v .^ 2
    c̃ = diag(fem.C_lumped)
    inner = spdiagm(0 => (κ²_v .^ 2) .* c̃) +
            spdiagm(0 => κ²_v) * fem.G1 +
            fem.G1 * spdiagm(0 => κ²_v) +
            fem.G2
    D_τ = spdiagm(0 => τ_v)
    return D_τ * inner * D_τ
end

"""
    spde_precision(α, τ, κ, C, G1[, C_lumped, G2]) -> Q

Stateless form: assemble `Q(α, τ, κ)` directly from the raw FEM matrices.
Missing `C_lumped` and `G2` are derived on the fly. Prefer
[`spde_precision(::FEMMatrices, ...)`](@ref) in hot loops — the
`FEMMatrices` constructor precomputes `C̃` and `G₂` once.
"""
function spde_precision(
        α::Integer, τ::Real, κ::Real,
        C::AbstractSparseMatrix, G1::AbstractSparseMatrix,
        C_lumped::Union{Nothing, AbstractSparseMatrix}=nothing,
        G2::Union{Nothing, AbstractSparseMatrix}=nothing
)
    if α == 1
        τ > 0 && κ > 0 ||
            throw(DomainError((τ, κ),
                "τ and κ must be positive; got τ=$τ, κ=$κ"))
        return τ^2 * (κ^2 * C + G1)
    end
    Cl = C_lumped === nothing ? lumped_mass(C) : C_lumped
    G2_ = G2 === nothing ? stiffness_squared(G1, Cl) : G2
    fem = FEMMatrices(C, G1, Cl, G2_)
    return spde_precision(fem, α, τ, κ)
end
