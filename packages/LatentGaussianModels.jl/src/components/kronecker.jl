"""
    KroneckerComponent(space::AbstractLatentComponent,
                       time::AbstractLatentComponent)

Separable space-time latent component composing two child
[`AbstractLatentComponent`](@ref)s into a Kronecker-structured GMRF
with precision

```
Q = Q_space ⊗ Q_time.
```

R-INLA's `f(field, model = spatial, group = time, control.group =
list(model = "ar1"))` builds the same separable form for SPDE-Matérn
× AR(1); the Cameletti et al. (2013) PM₁₀ case study is the canonical
oracle. Pairs with the [`KroneckerMapping`](@ref) projector so the
flattening conventions agree end-to-end.

# Hyperparameter layout

`θ = [θ_space; θ_time]` concatenated. The composer adds *no*
hyperparameters of its own:

```julia
nhyperparameters(c) == nhyperparameters(c.space) + nhyperparameters(c.time)
```

`log_hyperprior(c, θ)` sums the two children's log-priors.

# Latent flattening — matches `KroneckerMapping`

The latent vector `x` is `vec(X)` where `X` has shape
`(length(time) × length(space))`. Equivalently:
`x[(s - 1) · n_t + t]` is the value at space slot `s`, time slot `t`.
This matches the `KroneckerMapping(A_space, A_time)` storage convention
(`(A_space ⊗ A_time) vec(X) = vec(A_time · X · A_space')`), so a
`KroneckerComponent(space, time)` paired with
`KroneckerMapping(A_space, A_time)` projects correctly without any
transposition.

# Constraints

- Both children proper (`NoConstraint`): composer is proper.
- One child constrained, the other proper: the composer's constraint
  is the child's constraint lifted via Kronecker with the identity on
  the unconstrained axis.
- Both children constrained: rejected at construction. The
  joint-constraint case (rank, conditioning-by-kriging interaction
  with the Kronecker factor) is deferred — see ADR-029.

# Example

```julia
using LatentGaussianModels: KroneckerComponent, AR1, IID

# AR1 in space × AR1 in time, sharing nothing.
space = AR1(50)
time  = AR1(20)
spt   = KroneckerComponent(space, time)
length(spt) == 50 * 20  # 1000

# Compose with a KroneckerMapping projector for observation-side
# Kronecker structure (e.g. SPDE2 ⊗ AR1 for Cameletti).
```
"""
struct KroneckerComponent{S <: AbstractLatentComponent,
    T <: AbstractLatentComponent} <: AbstractLatentComponent
    space::S
    time::T

    function KroneckerComponent{S, T}(space::S, time::T) where {
            S <: AbstractLatentComponent, T <: AbstractLatentComponent}
        cs_s = GMRFs.constraints(space)
        cs_t = GMRFs.constraints(time)
        if !(cs_s isa NoConstraint) && !(cs_t isa NoConstraint)
            throw(ArgumentError(
                "KroneckerComponent: both children carry linear constraints " *
                "(spatial and temporal). The doubly-constrained Kronecker " *
                "case is deferred — see ADR-029. Wrap at most one " *
                "intrinsic child for now."
            ))
        end
        return new{S, T}(space, time)
    end
end

function KroneckerComponent(
        space::AbstractLatentComponent, time::AbstractLatentComponent)
    return KroneckerComponent{typeof(space), typeof(time)}(space, time)
end

Base.length(c::KroneckerComponent) = length(c.space) * length(c.time)

function nhyperparameters(c::KroneckerComponent)
    return nhyperparameters(c.space) + nhyperparameters(c.time)
end

function initial_hyperparameters(c::KroneckerComponent)
    return vcat(initial_hyperparameters(c.space), initial_hyperparameters(c.time))
end

# Slice the concatenated θ into the spatial and temporal child views.
function _kron_split_theta(c::KroneckerComponent, θ)
    p_s = nhyperparameters(c.space)
    p_t = nhyperparameters(c.time)
    length(θ) == p_s + p_t || throw(ArgumentError(
        "KroneckerComponent: θ has length $(length(θ)) but expected " *
        "$(p_s + p_t) = $(p_s) (space) + $(p_t) (time)"
    ))
    @views θ_s = θ[1:p_s]
    @views θ_t = θ[(p_s + 1):(p_s + p_t)]
    return θ_s, θ_t
end

function precision_matrix(c::KroneckerComponent, θ)
    θ_s, θ_t = _kron_split_theta(c, θ)
    Q_s = SparseMatrixCSC{Float64, Int}(precision_matrix(c.space, θ_s))
    Q_t = SparseMatrixCSC{Float64, Int}(precision_matrix(c.time, θ_t))
    return kron(Q_s, Q_t)
end

function log_hyperprior(c::KroneckerComponent, θ)
    θ_s, θ_t = _kron_split_theta(c, θ)
    return log_hyperprior(c.space, θ_s) + log_hyperprior(c.time, θ_t)
end

# Kronecker logdet identity: log|Q_s ⊗ Q_t| = n_t · log|Q_s| + n_s · log|Q_t|.
# Each child's `log_normalizing_constant` is the proper Gaussian log-NC
# `-½ d log(2π) + ½ log|Q_d|`, so we can recover ½ log|Q_d| as
# `log_normc_d + ½ d log(2π)` without needing to factor the dense
# Kronecker product.
function log_normalizing_constant(c::KroneckerComponent, θ)
    θ_s, θ_t = _kron_split_theta(c, θ)
    n_s = length(c.space)
    n_t = length(c.time)
    log_normc_s = log_normalizing_constant(c.space, θ_s)
    log_normc_t = log_normalizing_constant(c.time, θ_t)
    return 0.5 * n_s * n_t * log(2π) +
           n_t * log_normc_s + n_s * log_normc_t
end

# Outer-product mean: μ[(s - 1) n_t + t] = μ_s[s] · μ_t[t]. When either
# child has the default zero mean, kron returns a zero vector; both
# zero is the common case.
function prior_mean(c::KroneckerComponent, θ)
    θ_s, θ_t = _kron_split_theta(c, θ)
    μ_s = prior_mean(c.space, θ_s)
    μ_t = prior_mean(c.time, θ_t)
    return kron(μ_s, μ_t)
end

# At most one child carries a constraint (the doubly-constrained case
# is rejected at construction). The constrained child's constraint
# `(A_c, e_c)` lifts via Kronecker with the identity on the
# unconstrained axis, matching the latent flattening convention.
function GMRFs.constraints(c::KroneckerComponent)
    cs_s = GMRFs.constraints(c.space)
    cs_t = GMRFs.constraints(c.time)
    if cs_s isa NoConstraint && cs_t isa NoConstraint
        return NoConstraint()
    end
    n_s = length(c.space)
    n_t = length(c.time)
    if cs_s isa NoConstraint
        # Time-only constraint: A = I_{n_s} ⊗ A_t lifted to (k_t · n_s) rows.
        # `(I_s ⊗ A_t) vec(X) = vec(A_t X)` gives the temporal constraint
        # applied to each spatial slot independently.
        A_t = GMRFs.constraint_matrix(cs_t)
        e_t = GMRFs.constraint_rhs(cs_t)
        I_s = sparse(1.0I, n_s, n_s)
        A = Matrix(kron(I_s, sparse(Float64.(A_t))))
        e = repeat(Float64.(e_t), n_s)
        return LinearConstraint(A, e)
    else
        # Space-only constraint: A = A_s ⊗ I_{n_t} lifted to (k_s · n_t) rows.
        # `(A_s ⊗ I_t) vec(X) = vec(X A_s')` gives the spatial constraint
        # applied to each time slot independently.
        A_s = GMRFs.constraint_matrix(cs_s)
        e_s = GMRFs.constraint_rhs(cs_s)
        I_t = sparse(1.0I, n_t, n_t)
        A = Matrix(kron(sparse(Float64.(A_s)), I_t))
        e = repeat(Float64.(e_s), inner = n_t)
        return LinearConstraint(A, e)
    end
end

# Rank-deficiency under Kronecker: the eigenvalues of `Q_s ⊗ Q_t` are
# `λ_i μ_j` over all pairs, so the null-space dimension is
# `rd_s · n_t + n_s · rd_t − rd_s · rd_t` (inclusion-exclusion). Both
# children proper ⇒ rd = 0.
function gmrf(c::KroneckerComponent, θ)
    θ_s, θ_t = _kron_split_theta(c, θ)
    Q = precision_matrix(c, θ)
    rd_s = GMRFs.rankdef(gmrf(c.space, θ_s))
    rd_t = GMRFs.rankdef(gmrf(c.time, θ_t))
    n_s = length(c.space)
    n_t = length(c.time)
    rd = rd_s * n_t + n_s * rd_t - rd_s * rd_t
    return GMRFs.Generic0GMRF(Q; τ=1.0, rankdef=rd)
end
