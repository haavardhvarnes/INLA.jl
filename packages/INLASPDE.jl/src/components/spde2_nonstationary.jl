"""
    SPDE2NonStationary(points, triangles; α=2, B_τ, B_κ, prior=GaussianBasisPrior(p_τ + p_κ))
    SPDE2NonStationary(mesh::INLAMesh; α=2, B_τ, B_κ, prior=...)

Non-stationary integer-α SPDE-Matérn component on a 2D triangular
mesh. Per-vertex `log τ_v = (B_τ θ_τ)_v` and `log κ_v = (B_κ θ_κ)_v`
let the Matérn precision and range vary across the domain — the
parameterisation behind R-INLA's `inla.spde2.matern(..., B.tau, B.kappa,
theta.prior.mean, theta.prior.prec)`.

Currently restricted to `α = 2` (Matérn smoothness `ν = 1` in 2D);
`α = 1` is deferred to a follow-up.

# Hyperparameter layout

The internal hyperparameter vector concatenates τ-coefficients and
κ-coefficients in that order:

    θ = [θ_τ_1, …, θ_τ_{p_τ}, θ_κ_1, …, θ_κ_{p_κ}]

with `length(θ) = p_τ + p_κ = size(B_τ, 2) + size(B_κ, 2)`. The
[`GaussianBasisPrior`](@ref) on `θ` (ADR-028) factorises
independently across coefficients with per-coefficient `mean` and
`prec`, matching R-INLA's `theta.prior.mean` / `theta.prior.prec`.

# Fields

- `fem::FEMMatrices` — precomputed M1 FEM matrices.
- `graph::GMRFGraph` — mesh adjacency.
- `B_τ::Matrix{T}` — `(n_v, p_τ)` τ-basis matrix.
- `B_κ::Matrix{T}` — `(n_v, p_κ)` κ-basis matrix.
- `prior::GaussianBasisPrior` — independent Gaussian on basis
  coefficients.

# Stationary limit

If `B_τ = ones(n_v, 1)` and `B_κ = ones(n_v, 1)` and `θ = [log τ, log
κ]`, the precision reduces to `SPDE2`'s stationary form
`τ²·(κ⁴ C̃ + 2κ² G₁ + G₂)`. This makes stationary-recovery a clean
sub-test of the LRL §3.2 oracle fixture.

# Example

```julia
mesh = inla_mesh_2d(loc; max_edge = 0.1)
n_v = num_vertices(mesh)

# Two-region piecewise-constant κ, constant τ:
region = [x < 0.5 ? 1 : 2 for (x, _) in eachrow(mesh.points)]
B_κ = hcat([region .== 1, region .== 2]...) .|> Float64
B_τ = ones(n_v, 1)

prior = GaussianBasisPrior(
    mean = [0.0, 0.0, 0.0],
    prec = [1.0e-3, 1.0, 1.0]            # wide on intercept, tight on regions
)
spde = SPDE2NonStationary(mesh.points, mesh.triangles;
                          B_τ = B_τ, B_κ = B_κ, prior = prior)
```
"""
struct SPDE2NonStationary{α, T,
        FE <: FEMMatrices{T},
        G <: GMRFs.AbstractGMRFGraph,
        PR <: GaussianBasisPrior} <: AbstractLatentComponent
    fem::FE
    graph::G
    B_τ::Matrix{T}
    B_κ::Matrix{T}
    prior::PR
end

function SPDE2NonStationary(
        points::AbstractMatrix{<:Real},
        triangles::AbstractMatrix{<:Integer};
        α::Integer=2,
        B_τ::AbstractMatrix{<:Real},
        B_κ::AbstractMatrix{<:Real},
        prior::Union{GaussianBasisPrior, Nothing}=nothing
)
    α == 2 || throw(ArgumentError(
        "SPDE2NonStationary: only α = 2 is supported in v0.2; got α = $α. " *
        "α = 1 and fractional α deferred."
    ))
    fem = FEMMatrices(points, triangles)
    n_v = size(fem.C, 1)
    size(B_τ, 1) == n_v || throw(ArgumentError(
        "SPDE2NonStationary: size(B_τ, 1) ($(size(B_τ, 1))) must equal " *
        "n_vertices ($n_v)"
    ))
    size(B_κ, 1) == n_v || throw(ArgumentError(
        "SPDE2NonStationary: size(B_κ, 1) ($(size(B_κ, 1))) must equal " *
        "n_vertices ($n_v)"
    ))
    p_τ = size(B_τ, 2)
    p_κ = size(B_κ, 2)
    p_τ >= 1 || throw(ArgumentError("SPDE2NonStationary: B_τ needs ≥ 1 column"))
    p_κ >= 1 || throw(ArgumentError("SPDE2NonStationary: B_κ needs ≥ 1 column"))

    pr = prior === nothing ? GaussianBasisPrior(p_τ + p_κ) : prior
    length(pr) == p_τ + p_κ || throw(ArgumentError(
        "SPDE2NonStationary: length(prior) ($(length(pr))) must equal " *
        "p_τ + p_κ ($(p_τ + p_κ))"
    ))

    T = eltype(fem.C)
    graph = _mesh_graph_from_C(fem.C)
    return SPDE2NonStationary{Int(α), T, typeof(fem), typeof(graph), typeof(pr)}(
        fem, graph, Matrix{T}(B_τ), Matrix{T}(B_κ), pr
    )
end

function SPDE2NonStationary(
        mesh::INLAMesh;
        α::Integer=2,
        B_τ::AbstractMatrix{<:Real},
        B_κ::AbstractMatrix{<:Real},
        prior::Union{GaussianBasisPrior, Nothing}=nothing
)
    return SPDE2NonStationary(mesh.points, mesh.triangles;
        α = α, B_τ = B_τ, B_κ = B_κ, prior = prior)
end

# --- AbstractLatentComponent contract --------------------------------

Base.length(c::SPDE2NonStationary) = size(c.fem.C, 1)

function LatentGaussianModels.nhyperparameters(c::SPDE2NonStationary)
    return size(c.B_τ, 2) + size(c.B_κ, 2)
end

"""
    initial_hyperparameters(c::SPDE2NonStationary) -> Vector

Default initial point for the BFGS θ-mode finder: the prior mean.
This is the natural starting point for a Gaussian prior and matches
R-INLA's `theta.initial = theta.prior.mean` default.
"""
function LatentGaussianModels.initial_hyperparameters(c::SPDE2NonStationary)
    return copy(c.prior.mean)
end

"""
    precision_matrix(c::SPDE2NonStationary, θ) -> SparseMatrixCSC

Sparse non-stationary precision at `θ = [θ_τ; θ_κ]` of length
`p_τ + p_κ`. Computes per-vertex `τ_v = exp(B_τ θ_τ)` and
`κ_v = exp(B_κ θ_κ)`, then dispatches to
[`spde_precision_nonstationary`](@ref).
"""
function LatentGaussianModels.precision_matrix(c::SPDE2NonStationary{α}, θ) where {α}
    p_τ = size(c.B_τ, 2)
    p_κ = size(c.B_κ, 2)
    length(θ) == p_τ + p_κ || throw(ArgumentError(
        "SPDE2NonStationary: length(θ) ($(length(θ))) must equal " *
        "p_τ + p_κ ($(p_τ + p_κ))"
    ))
    θ_τ = @view θ[1:p_τ]
    θ_κ = @view θ[(p_τ + 1):(p_τ + p_κ)]
    τ_v = exp.(c.B_τ * θ_τ)
    κ_v = exp.(c.B_κ * θ_κ)
    return spde_precision_nonstationary(c.fem, α, τ_v, κ_v)
end

"""
    log_hyperprior(c::SPDE2NonStationary, θ) -> Real

Independent-Gaussian log density on `θ = [θ_τ; θ_κ]` via
[`log_basis_prior_density`](@ref).
"""
function LatentGaussianModels.log_hyperprior(c::SPDE2NonStationary, θ)
    return log_basis_prior_density(c.prior, θ)
end

"""
    GMRFs.constraints(::SPDE2NonStationary) -> NoConstraint

The non-stationary SPDE precision is SPD as long as `τ_v, κ_v > 0`,
which the construction enforces.
"""
GMRFs.constraints(::SPDE2NonStationary) = GMRFs.NoConstraint()

"""
    log_normalizing_constant(c::SPDE2NonStationary, θ) -> Real

Proper-Gaussian normalizing constant `-½ d log(2π) + ½ log|Q(θ)|`,
matching `SPDE2`. Required for marginal-likelihood comparison with
R-INLA.
"""
function LatentGaussianModels.log_normalizing_constant(
        c::SPDE2NonStationary, θ
)
    Q = LatentGaussianModels.precision_matrix(c, θ)
    F = cholesky(Symmetric(Q))
    return -0.5 * length(c) * log(2π) + 0.5 * logdet(F)
end
