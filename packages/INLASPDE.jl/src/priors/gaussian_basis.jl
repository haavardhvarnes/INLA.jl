"""
    GaussianBasisPrior(; mean, prec)
    GaussianBasisPrior(p::Integer; mean = 0.0, prec = 1.0)

Independent Gaussian prior on the basis-coefficient vector
`θ = [θ_τ; θ_κ]` of [`SPDE2NonStationary`](@ref). Each coefficient has
its own mean and precision:

    log π(θ) = ∑_k -½ λ_k (θ_k - μ_k)^2 + ½ log(λ_k / 2π)

The normalising constant is included so `log_basis_prior_density`
returns the proper log density — required for marginal-likelihood
comparison against R-INLA's `inla.spde2.matern`, which uses the same
parameterisation under `theta.prior.mean` / `theta.prior.prec`.

# Fields

- `mean::Vector{T}` — prior means `μ_k`, one per basis coefficient.
- `prec::Vector{T}` — prior precisions `λ_k > 0`, one per basis
  coefficient. `prec[k] = 0` is rejected; users wanting an improper
  flat prior on a coefficient should choose a small but positive value
  (e.g. `1.0e-3`), mirroring R-INLA practice.

# Length contract

`length(prior.mean) == length(prior.prec)` is enforced at construction
and equals the total number of basis coefficients
(`p_τ + p_κ` for [`SPDE2NonStationary`](@ref)).

# Defaults

`mean = 0`, `prec = 1` for every coefficient — R-INLA's documented
`inla.spde2.matern` defaults. Override per-coefficient for canonical
non-stationary fits (a wide intercept prior + tighter spline-basis
priors is the standard pattern).

# Example

```julia
# Two τ-coefficients (intercept + spline) and three κ-coefficients
# (intercept + two region-indicators):
prior = GaussianBasisPrior(
    mean = [0.0, 0.0, 0.0, 0.0, 0.0],
    prec = [1.0e-3, 1.0, 1.0e-3, 1.0, 1.0]   # wide on intercepts
)
```

See ADR-028 for the design rationale (R-INLA-parity choice over the
predecessor's unit-Gaussian default and over PC-on-basis-norm).
"""
struct GaussianBasisPrior{T <: Real}
    mean::Vector{T}
    prec::Vector{T}
end

function GaussianBasisPrior(;
        mean::AbstractVector{<:Real},
        prec::AbstractVector{<:Real}
)
    length(mean) == length(prec) || throw(ArgumentError(
        "GaussianBasisPrior: length(mean) ($(length(mean))) must equal " *
        "length(prec) ($(length(prec)))"
    ))
    all(>(0), prec) || throw(ArgumentError(
        "GaussianBasisPrior: all entries of prec must be > 0; got $prec. " *
        "Use a small positive value (e.g. 1.0e-3) for a wide prior."
    ))
    T = promote_type(eltype(float.(mean)), eltype(float.(prec)))
    return GaussianBasisPrior{T}(Vector{T}(mean), Vector{T}(prec))
end

function GaussianBasisPrior(p::Integer; mean::Real=0.0, prec::Real=1.0)
    p >= 1 || throw(ArgumentError(
        "GaussianBasisPrior: number of coefficients must be ≥ 1; got $p"
    ))
    return GaussianBasisPrior(; mean=fill(float(mean), p), prec=fill(float(prec), p))
end

Base.length(prior::GaussianBasisPrior) = length(prior.mean)

"""
    log_basis_prior_density(prior::GaussianBasisPrior, θ) -> Real

Log density of the independent-Gaussian basis prior at coefficient
vector `θ`, including the normalising constant
`½ ∑_k log(λ_k / 2π)`.

Throws `ArgumentError` on length mismatch.
"""
function log_basis_prior_density(prior::GaussianBasisPrior, θ::AbstractVector{<:Real})
    length(θ) == length(prior) || throw(ArgumentError(
        "GaussianBasisPrior: length(θ) ($(length(θ))) must equal " *
        "length(prior) ($(length(prior)))"
    ))
    s = zero(promote_type(eltype(prior.mean), eltype(θ)))
    for k in eachindex(θ)
        d = θ[k] - prior.mean[k]
        λ = prior.prec[k]
        s += -0.5 * λ * d * d + 0.5 * (log(λ) - log(2π))
    end
    return s
end
