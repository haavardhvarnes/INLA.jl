"""
    AbstractIIDND <: AbstractLatentComponent

Umbrella type for the multivariate-IID family. Two concrete subtypes
(see ADR-022):

- `IIDND_Sep{N}` — separable parameterisation: `N` marginal precisions
  + `N(N-1)/2` correlations, each carrying its own scalar
  `AbstractHyperPrior`. This matches R-INLA's
  `f(., model = "2diid"/"iid3d", ...)` defaults.
- `IIDND_Joint{N}` — joint Λ parameterisation with a single
  `AbstractJointHyperPrior` (Wishart / InvWishart). Lands in PR-1c.

The latent vector has length `n · N` and stores the `N` blocks of length
`n` consecutively: `x = (x^{(1)}, x^{(2)}, …, x^{(N)})` where each
`x^{(j)} ∈ ℝ^n`. Joint precision is `Q = Λ ⊗ I_n` for both subtypes.
"""
abstract type AbstractIIDND <: AbstractLatentComponent end

"""
    IIDND_Sep{N, PP, PC} <: AbstractIIDND

Separable parameterisation of an `N`-variate IID component. `PP` is the
type of `precpriors` (an `NTuple{N, AbstractHyperPrior}`); `PC` is the
type of `corrpriors` (an `NTuple{N(N-1)/2, AbstractHyperPrior}`).

Hyperparameter vector layout (length `N + N(N-1)/2`):

    θ = (log τ_1, …, log τ_N,  atanh z_{2,1}, atanh z_{3,1}, atanh z_{3,2}, …,  atanh z_{N,N-1})

with the strictly-lower-triangular entries enumerated in row-major order
(`(i, j)` with `j < i`). The `z_{i,j}` are the LKJ canonical partial
correlations of the correlation matrix `R` (Lewandowski-Kurowicka-Joe
2009): for the first column they coincide with the direct correlations
(`z_{i,1} = ρ_{i,1}`); for entries with `i, j > 1` they are partial
correlations of variables `i` and `j` controlling for variables
`1, …, j − 1`. The correlation matrix is reconstructed from the `z`'s
via stick-breaking, which guarantees `R` is positive-definite for any
`(z_{i,j}) ∈ (-1, 1)^{N(N-1)/2}`.

PR-1a ships the `N = 2` case (where `z_{2,1} = ρ_{1,2}` directly);
PR-1b lifts the `N = 3` case via the LKJ stick-breaking step.
"""
struct IIDND_Sep{N, PP <: Tuple, PC <: Tuple} <: AbstractIIDND
    n::Int
    precpriors::PP
    corrpriors::PC
end

"""
    IIDND(n, N = 2; hyperprior_precs, hyperprior_corr, hyperprior_corrs)

Multivariate IID latent component of dimension `n × N`. See ADR-022 for
the parameterisation choice; PR-1a/PR-1b ship `N ∈ {2, 3}`. ADR-022
caps the separable form at `N ≤ 3` (the `atanh-of-each-pairwise-corr`
parameterisation is not injective onto positive-definite correlation
matrices for `N ≥ 4`).

For `N = 2`, pass `hyperprior_corr` (a single `AbstractHyperPrior` on
`atanh ρ`). For `N ≥ 3`, pass `hyperprior_corrs` as an `NTuple` of
length `N(N-1)/2`.

Defaults: each marginal precision uses `PCPrecision()`; each correlation
uses `PCCor0()` (R-INLA's `pc.cor0(0.5, 0.5)` default).
"""
function IIDND(n::Integer, N::Integer=2;
        hyperprior_precs::Union{Nothing, NTuple}=nothing,
        hyperprior_corr::Union{Nothing, AbstractHyperPrior}=nothing,
        hyperprior_corrs::Union{Nothing, NTuple}=nothing)
    n > 0 || throw(ArgumentError("IIDND: n must be positive"))
    N ≥ 2 || throw(ArgumentError("IIDND: N must be ≥ 2"))
    N ≤ 3 ||
        throw(ArgumentError("IIDND: N ≥ 4 deferred to a successor ADR per ADR-022 — separable form caps at N = 3"))

    if hyperprior_corr !== nothing && hyperprior_corrs !== nothing
        throw(ArgumentError("IIDND: pass exactly one of `hyperprior_corr` (N=2) or `hyperprior_corrs` (N≥3)"))
    end

    pps = hyperprior_precs === nothing ?
          ntuple(_ -> PCPrecision(), N) :
          hyperprior_precs
    length(pps) == N ||
        throw(ArgumentError("IIDND: hyperprior_precs must have length $N, got $(length(pps))"))
    all(p -> p isa AbstractHyperPrior, pps) ||
        throw(ArgumentError("IIDND: every entry of hyperprior_precs must be an AbstractHyperPrior"))

    K = N * (N - 1) ÷ 2
    cps = if hyperprior_corr !== nothing
        N == 2 ||
            throw(ArgumentError("IIDND: `hyperprior_corr` is for N = 2; use `hyperprior_corrs` for N ≥ 3"))
        (hyperprior_corr,)
    elseif hyperprior_corrs !== nothing
        length(hyperprior_corrs) == K ||
            throw(ArgumentError("IIDND: hyperprior_corrs must have length $K, got $(length(hyperprior_corrs))"))
        all(p -> p isa AbstractHyperPrior, hyperprior_corrs) ||
            throw(ArgumentError("IIDND: every entry of hyperprior_corrs must be an AbstractHyperPrior"))
        hyperprior_corrs
    else
        ntuple(_ -> PCCor0(), K)
    end

    return IIDND_Sep{N, typeof(pps), typeof(cps)}(Int(n), pps, cps)
end

"""
    IID2D(n; hyperprior_precs = (PCPrecision(), PCPrecision()),
              hyperprior_corr  = PCCor0())

Bivariate IID component — ergonomic alias for `IIDND(n, 2; …)`. Mirrors
R-INLA's `f(idx, model = "2diid", …)` defaults.
"""
function IID2D(n::Integer;
        hyperprior_precs::NTuple{2, AbstractHyperPrior}=(PCPrecision(), PCPrecision()),
        hyperprior_corr::AbstractHyperPrior=PCCor0())
    return IIDND(n, 2;
        hyperprior_precs=hyperprior_precs,
        hyperprior_corr=hyperprior_corr)
end

"""
    IID3D(n; hyperprior_precs = (PCPrecision(), PCPrecision(), PCPrecision()),
              hyperprior_corrs = (PCCor0(), PCCor0(), PCCor0()))

Trivariate IID component — ergonomic alias for `IIDND(n, 3; …)`. The
three correlation slots are LKJ canonical partial correlations on
`atanh` scale (see `IIDND_Sep`'s docstring): θ entries `4` and `5` are
direct correlations `ρ_{1,2}` and `ρ_{1,3}` on `atanh` scale; entry `6`
is the partial correlation of variables 2 and 3 controlling for
variable 1, also on `atanh` scale.

The closest R-INLA analogue is `f(idx, model = "iid3d", …)`, which
defaults to a Wishart prior on the joint precision; the separable
default here uses `PCCor0` per pairwise CPC and so will not match
R-INLA's defaults exactly. Per ADR-022, the Wishart alternative lands
in PR-1c.
"""
function IID3D(n::Integer;
        hyperprior_precs::NTuple{3, AbstractHyperPrior}=(
            PCPrecision(), PCPrecision(), PCPrecision()),
        hyperprior_corrs::NTuple{3, AbstractHyperPrior}=(PCCor0(), PCCor0(), PCCor0()))
    return IIDND(n, 3;
        hyperprior_precs=hyperprior_precs,
        hyperprior_corrs=hyperprior_corrs)
end

Base.length(c::IIDND_Sep{N}) where {N} = c.n * N

nhyperparameters(::IIDND_Sep{N}) where {N} = N + N * (N - 1) ÷ 2

initial_hyperparameters(c::IIDND_Sep{N}) where {N} = zeros(Float64, nhyperparameters(c))

GMRFs.constraints(::IIDND_Sep) = NoConstraint()

# ---------------------------------------------------------------------
# N = 2 specialisations (PR-1a). N ≥ 3 lifted in PR-1b.
# ---------------------------------------------------------------------

# Internal: bivariate precision-matrix entries `(Λ_11, Λ_12, Λ_22)` from
# `(τ_1, τ_2, t)` with `t = atanh(ρ)`. The closed form
#
#     Λ = (1/(1 - ρ²)) · [τ_1            -ρ √(τ_1 τ_2);
#                          -ρ √(τ_1 τ_2)   τ_2         ]
#
# is rewritten via `1 - ρ² = sech²(t)` and `ρ/(1 - ρ²) = sinh(t) cosh(t)`
# to stay finite where `1 - tanh(t)²` underflows to `0` in float64
# (|t| ≳ 19) — `precision_matrix(IIDND_Sep{2}, θ)` is called inside the
# outer-θ LBFGS line search, which can probe |θ[3]| ≫ 19 before the
# objective is evaluated to discover the region is bad.
#
# Returned as a 3-tuple so `precision_matrix` can build the sparse
# `Λ ⊗ I_n` directly without allocating a 2×2 matrix.
function _iid2d_lambda(τ1, τ2, t)
    cosh_t = cosh(t)
    sinh_t = sinh(t)
    cosh²_t = cosh_t * cosh_t
    Λ11 = τ1 * cosh²_t
    Λ22 = τ2 * cosh²_t
    Λ12 = -sqrt(τ1 * τ2) * sinh_t * cosh_t
    return Λ11, Λ12, Λ22
end

function precision_matrix(c::IIDND_Sep{2}, θ)
    n = c.n
    τ1 = exp(θ[1])
    τ2 = exp(θ[2])
    Λ11, Λ12, Λ22 = _iid2d_lambda(τ1, τ2, θ[3])
    diag_main = vcat(fill(Λ11, n), fill(Λ22, n))
    off_diag = fill(Λ12, n)
    return spdiagm(0 => diag_main, n => off_diag, -n => off_diag)
end

function log_hyperprior(c::IIDND_Sep{2}, θ)
    return log_prior_density(c.precpriors[1], θ[1]) +
           log_prior_density(c.precpriors[2], θ[2]) +
           log_prior_density(c.corrpriors[1], θ[3])
end

# Proper, full-rank Q = Λ ⊗ I_n with `det(Q) = det(Λ)^n` and
# `det(Λ) = τ_1 τ_2 / (1 - ρ²)`. Total dimension is `2n`. The
# `-log(1 - ρ²)` term is computed as `2·logcosh(θ[3])` to remain finite
# at saturation; see `_iid2d_lambda` and `_logcosh` for context.
#
#     log NC = -½ · 2n · log(2π) + ½ · n · (θ_1 + θ_2 + 2·logcosh(θ_3)).
function log_normalizing_constant(c::IIDND_Sep{2}, θ)
    n = c.n
    return -n * log(2π) + 0.5 * n * (θ[1] + θ[2] + 2 * _logcosh(θ[3]))
end

function gmrf(c::IIDND_Sep{2}, θ)
    Q = precision_matrix(c, θ)
    return GMRFs.Generic0GMRF(Q; τ=1.0, rankdef=0, scale_model=false)
end

# ---------------------------------------------------------------------
# N = 3 specialisations (PR-1b). LKJ stick-breaking with canonical
# partial correlations `z_{2,1}, z_{3,1}, z_{3,2} = tanh(θ[4..6])`.
#
# The 3×3 correlation matrix `R` has Cholesky factor `L` with
#   L_{1,1} = 1
#   L_{2,1} = z_{2,1},                    L_{2,2} = √(1 − z_{2,1}²)
#   L_{3,1} = z_{3,1}, L_{3,2} = z_{3,2}·√(1 − z_{3,1}²),
#                                          L_{3,3} = √((1 − z_{3,1}²)(1 − z_{3,2}²))
# so `R = L L'` is positive-definite by construction. The joint
# precision is `Λ = D_τ^{1/2} R^{-1} D_τ^{1/2} = G' G` with
# `G = L^{-1} D_τ^{1/2}`.
#
# Closed-form `L^{-1}` entries simplify cleanly under
# `1/√(1 − tanh²(t)) = cosh(t)`:
#   L^{-1}_{2,1} = -sinh(a),
#   L^{-1}_{2,2} =  cosh(a),
#   L^{-1}_{3,1} =  sinh(c) sinh(a) − sinh(b) cosh(c),
#   L^{-1}_{3,2} = -sinh(c) cosh(a),
#   L^{-1}_{3,3} =  cosh(b) cosh(c),
# where (a, b, c) = (θ[4], θ[5], θ[6]). The cosh/sinh form stays finite
# for any θ — the same saturation guarantee `_iid2d_lambda` provides for
# N = 2.
function _iid3d_lambda(τ1, τ2, τ3, a, b, c)
    s1 = sqrt(τ1)
    s2 = sqrt(τ2)
    s3 = sqrt(τ3)
    ca = cosh(a)
    sha = sinh(a)
    cb = cosh(b)
    shb = sinh(b)
    cc = cosh(c)
    shc = sinh(c)
    # G = L^{-1} D_τ^{1/2}; Λ = G' G via columnwise inner products.
    g11 = s1
    g21 = -s1 * sha
    g22 = s2 * ca
    g31 = s1 * (shc * sha - shb * cc)
    g32 = -s2 * shc * ca
    g33 = s3 * cb * cc
    Λ11 = g11 * g11 + g21 * g21 + g31 * g31
    Λ22 = g22 * g22 + g32 * g32
    Λ33 = g33 * g33
    Λ12 = g21 * g22 + g31 * g32
    Λ13 = g31 * g33
    Λ23 = g32 * g33
    return Λ11, Λ22, Λ33, Λ12, Λ13, Λ23
end

function precision_matrix(c::IIDND_Sep{3}, θ)
    n = c.n
    τ1 = exp(θ[1])
    τ2 = exp(θ[2])
    τ3 = exp(θ[3])
    Λ11, Λ22, Λ33, Λ12, Λ13, Λ23 = _iid3d_lambda(τ1, τ2, τ3, θ[4], θ[5], θ[6])
    diag_main = vcat(fill(Λ11, n), fill(Λ22, n), fill(Λ33, n))
    off_n = vcat(fill(Λ12, n), fill(Λ23, n))   # length 2n
    off_2n = fill(Λ13, n)                       # length n
    return spdiagm(0 => diag_main,
        n => off_n, -n => off_n,
        2n => off_2n, -2n => off_2n)
end

function log_hyperprior(c::IIDND_Sep{3}, θ)
    return log_prior_density(c.precpriors[1], θ[1]) +
           log_prior_density(c.precpriors[2], θ[2]) +
           log_prior_density(c.precpriors[3], θ[3]) +
           log_prior_density(c.corrpriors[1], θ[4]) +
           log_prior_density(c.corrpriors[2], θ[5]) +
           log_prior_density(c.corrpriors[3], θ[6])
end

# Proper, full-rank `Q = Λ ⊗ I_n` with `det(Q) = det(Λ)^n` and
#   det(Λ) = τ_1 τ_2 τ_3 · cosh²(a) · cosh²(b) · cosh²(c).
# Total dimension is `3n`. The cosh terms come from `det R⁻¹ = det L⁻²`;
# `2·logcosh(θ_k)` is the saturation-stable form.
#
#   log NC = -½ · 3n · log(2π)
#          + ½ · n · (θ_1 + θ_2 + θ_3 + 2·logcosh(θ_4) + 2·logcosh(θ_5) + 2·logcosh(θ_6)).
function log_normalizing_constant(c::IIDND_Sep{3}, θ)
    n = c.n
    return -1.5 * n * log(2π) +
           0.5 * n *
           (θ[1] + θ[2] + θ[3] +
            2 * (_logcosh(θ[4]) + _logcosh(θ[5]) + _logcosh(θ[6])))
end

function gmrf(c::IIDND_Sep{3}, θ)
    Q = precision_matrix(c, θ)
    return GMRFs.Generic0GMRF(Q; τ=1.0, rankdef=0, scale_model=false)
end
