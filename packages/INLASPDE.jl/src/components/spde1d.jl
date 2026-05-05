"""
    SPDE1D{α}(points, segments; pc = PCMatern{1}())
    SPDE1D{α}(mesh::INLAMesh1D; pc = PCMatern{1}())

Integer-α SPDE–Matérn component on a 1D segment mesh, implementing
the `AbstractLatentComponent` contract from `LatentGaussianModels`.

The latent field has length `num_vertices(mesh)` and follows a GMRF
with sparse precision `Q(τ, κ)` assembled from the 1D FEM matrices
`(C, G₁, C̃, G₂)`. Internal hyperparameters are `θ = [log τ, log κ]`
per ADR-013; the user-facing scale is the Matérn `(ρ, σ)` pair with
mapping

    ρ = √(8ν) / κ,
    σ² = Γ(ν) / (Γ(ν + 0.5) · (4π)^(1/2) · κ^(2ν) · τ²)         (d = 1)

where `ν = α - d/2`. Both `α = 1` (⇒ `ν = 0.5`, exponential covariance)
and `α = 2` (⇒ `ν = 1.5`) are supported; the 1D PC-Matern density is
non-degenerate at both.

# Fields

- `fem::FEMMatrices` — the 1D FEM matrices, assembled once at
  construction.
- `graph::GMRFGraph` — mesh adjacency, derived from the off-diagonal
  pattern of `C`. Used for diagnostics and component-machinery
  compatibility.
- `pc::PCMatern{1}` — joint PC prior on `(ρ, σ)` with `D = 1`.

# Example

```julia
mesh = inla_mesh_1d([0.0, 1.0, 2.0]; max_edge = 0.1)
spde = SPDE1D(mesh; α = 2,
              pc = PCMatern{1}(range_U=0.5, range_α=0.05,
                               sigma_U=1.0, sigma_α=0.01))
Q = precision_matrix(spde, [0.0, 0.0])   # (log τ, log κ) = (0, 0)
```
"""
struct SPDE1D{
        α, T, FE <: FEMMatrices{T}, G <: GMRFs.AbstractGMRFGraph, PR <: PCMatern{1}
} <: AbstractLatentComponent
    fem::FE
    graph::G
    pc::PR
end

function SPDE1D(
        points::AbstractVector{<:Real},
        segments::AbstractMatrix{<:Integer};
        α::Integer=2,
        pc::PCMatern{1}=PCMatern{1}()
)
    α in (1, 2) ||
        throw(ArgumentError("SPDE1D: only α ∈ {1, 2} is supported; got α=$α. " *
                            "Fractional α is deferred to v0.3 (Bolin–Kirchner)."))
    fem = FEMMatrices(points, segments)
    graph = _mesh_graph_from_C(fem.C)
    T = eltype(fem.C)
    return SPDE1D{Int(α), T, typeof(fem), typeof(graph), typeof(pc)}(fem, graph, pc)
end

SPDE1D(mesh::INLAMesh1D; kwargs...) = SPDE1D(mesh.points, mesh.segments; kwargs...)

# --- AbstractLatentComponent contract ---------------------------------

Base.length(c::SPDE1D) = size(c.fem.C, 1)
LatentGaussianModels.nhyperparameters(::SPDE1D) = 2
LatentGaussianModels.initial_hyperparameters(::SPDE1D) = [0.0, 0.0]   # (log τ, log κ) = (0, 0)

"""
    precision_matrix(c::SPDE1D{α}, θ) -> SparseMatrixCSC

Sparse SPDE precision at `θ = [log τ, log κ]`, delegating to
`spde_precision(fem, α, τ, κ)`.
"""
function LatentGaussianModels.precision_matrix(c::SPDE1D{α}, θ) where {α}
    τ = exp(θ[1])
    κ = exp(θ[2])
    return spde_precision(c.fem, α, τ, κ)
end

"""
    log_hyperprior(c::SPDE1D{α}, θ) -> Real

PC-Matern log-prior density evaluated at `θ = [log τ, log κ]`. Maps
`(log τ, log κ) → (log ρ, log σ)` per the d=1 formulas and calls
`pc_matern_log_density`. The (log τ, log κ) ↔ (log ρ, log σ) Jacobian
has absolute determinant 1 (ADR-013).
"""
function LatentGaussianModels.log_hyperprior(c::SPDE1D{1}, θ)
    log_τ, log_κ = θ[1], θ[2]
    # α = 1, d = 1 ⇒ ν = 0.5
    #   ρ = 2/κ                ⇒ log ρ = log 2 - log κ
    #   σ² = 1 / (2 κ τ²)      ⇒ log σ = -0.5 log(2κ) - log τ
    log_ρ = log(2.0) - log_κ
    log_σ = -0.5 * (log(2.0) + log_κ) - log_τ
    return pc_matern_log_density(c.pc, log_ρ, log_σ)
end

function LatentGaussianModels.log_hyperprior(c::SPDE1D{2}, θ)
    log_τ, log_κ = θ[1], θ[2]
    # α = 2, d = 1 ⇒ ν = 1.5
    #   ρ = √12 / κ            ⇒ log ρ = 0.5 log 12 - log κ
    #   σ² = 1 / (4 κ³ τ²)     ⇒ log σ = -log 2 - 1.5 log κ - log τ
    log_ρ = 0.5 * log(12.0) - log_κ
    log_σ = -log(2.0) - 1.5 * log_κ - log_τ
    return pc_matern_log_density(c.pc, log_ρ, log_σ)
end

"""
    GMRFs.constraints(::SPDE1D) -> NoConstraint

The SPDE precision is strictly positive-definite for all `(τ, κ) > 0`
on a non-degenerate 1D mesh, so no hard linear constraint is required.
"""
GMRFs.constraints(::SPDE1D) = GMRFs.NoConstraint()

"""
    log_normalizing_constant(c::SPDE1D, θ) -> Real

R-INLA-style log normalizing constant for the proper Gaussian SPDE
prior, evaluated at internal `θ = [log τ, log κ]`:

    log NC = -½ d log(2π) + ½ log|Q(θ)|

where `d = length(c)` is the number of mesh vertices and `Q(θ)` is the
SPDE precision.
"""
function LatentGaussianModels.log_normalizing_constant(c::SPDE1D, θ)
    Q = LatentGaussianModels.precision_matrix(c, θ)
    F = cholesky(Symmetric(Q))
    return -0.5 * length(c) * log(2π) + 0.5 * logdet(F)
end

"""
    spde_user_scale(c::SPDE1D{α}, θ) -> (ρ, σ)

Convert internal `θ = [log τ, log κ]` to the user-facing Matérn pair
`(ρ, σ)` for `d = 1`, `α ∈ {1, 2}`.
"""
function spde_user_scale(::SPDE1D{1}, θ)
    log_τ, log_κ = θ[1], θ[2]
    log_ρ = log(2.0) - log_κ
    log_σ = -0.5 * (log(2.0) + log_κ) - log_τ
    return exp(log_ρ), exp(log_σ)
end

function spde_user_scale(::SPDE1D{2}, θ)
    log_τ, log_κ = θ[1], θ[2]
    log_ρ = 0.5 * log(12.0) - log_κ
    log_σ = -log(2.0) - 1.5 * log_κ - log_τ
    return exp(log_ρ), exp(log_σ)
end

"""
    spde_internal_scale(c::SPDE1D{α}, ρ, σ) -> (log τ, log κ)

Inverse of [`spde_user_scale`](@ref) for `d = 1`, `α ∈ {1, 2}`.
"""
function spde_internal_scale(::SPDE1D{1}, ρ::Real, σ::Real)
    ρ > 0 || throw(ArgumentError("ρ must be positive; got $ρ"))
    σ > 0 || throw(ArgumentError("σ must be positive; got $σ"))
    log_κ = log(2.0) - log(ρ)
    log_τ = -0.5 * (log(2.0) + log_κ) - log(σ)
    return log_τ, log_κ
end

function spde_internal_scale(::SPDE1D{2}, ρ::Real, σ::Real)
    ρ > 0 || throw(ArgumentError("ρ must be positive; got $ρ"))
    σ > 0 || throw(ArgumentError("σ must be positive; got $σ"))
    log_κ = 0.5 * log(12.0) - log(ρ)
    log_τ = -log(2.0) - 1.5 * log_κ - log(σ)
    return log_τ, log_κ
end
