# Posterior marginal densities for latent coordinates and hyperparameters.
#
# For `x_i | y` the strategies live under `AbstractMarginalStrategy` (see
# ADR-026):
#
# - `Gaussian()` — mixture of per-θ Gaussians with
#   `(mean = mode_i(θ_k), var = diag(H(θ_k)⁻¹)_i)`, weighted by `w_k`.
# - `SimplifiedLaplace()` — same mixture, but each Gaussian component is
#   multiplied by an Edgeworth first-order skewness correction
#   `(1 + γ/6 · H₃(s))` where `γ = κ_3/σ³` is the posterior skewness of
#   `x_i` under the Laplace at `θ_k` and `H₃(s) = s³ - 3s` is the third
#   Hermite polynomial. This corresponds to Rue-Martino-Chopin (2009) §4.2
#   with the Gaussian marginal augmented by the leading skew term.
# - `FullLaplace()` — per-`x_i` refitted Laplace via constraint injection.
#   At each output grid point `a`, runs `laplace_mode_fixed_xi` with the
#   augmented constraint `e_i^T x = a` and reads the constrained
#   log-marginal; mixes across θ. R-INLA's `strategy = "laplace"`. See
#   `inference/full_laplace.jl`.
#
# For `θ_j | y` two methods are available (ADR-046):
#
# - `:gaussian` — Gaussian at (θ̂_j, Σθ[j,j]), the pre-ADR-046 behaviour.
# - `:integrated` — numerically integrated density. For dim(θ) = 1 the
#   density is read directly off the retained design points
#   (`INLAResult.log_π`) and interpolated; for dim(θ) ≥ 2 it is a
#   conditional-mode profile slice of `log π̂(θ | y)` along the
#   Gaussian-conditional path through θ̂ (fresh `laplace_mode` calls,
#   requires `model` and `y`).

"""
    posterior_marginal_x(res::INLAResult, i::Int;
                         strategy = Gaussian(),
                         model = nothing, y = nothing,
                         grid_size = 75, span = 5.0,
                         grid = nothing) -> @NamedTuple{x::Vector, pdf::Vector}

Posterior density of the `i`-th latent coordinate, evaluated on a grid.

Returns a named tuple `(x, pdf)`. The density is the θ-mixture

    p(x_i | y) ≈ Σ_k w_k · π_k(x_i)

where `π_k` depends on `strategy` (an `AbstractMarginalStrategy`
instance, or one of the legacy symbols `:gaussian` /
`:simplified_laplace` — see ADR-026):

- `Gaussian()` / `:gaussian` (default) — `π_k = φ(x_i; mode_k[i], var_k[i])`.
- `SimplifiedLaplace()` / `:simplified_laplace` — `π_k = φ · (1 + γ_k / 6 · H₃(s))`
  with `s = (x_i - mode_k[i]) / σ_k`, `H₃(s) = s³ - 3s`, and posterior
  skewness `γ_k = κ_3(x_i|θ_k) / σ_k³`. The third cumulant is
  assembled from `∇³_η_log_density` and the Laplace precision at
  `θ_k`; for a Gaussian likelihood this collapses to the Gaussian
  strategy. Requires `model` and `y` to be supplied.
- `FullLaplace()` / `:full_laplace` — at each grid point `a` in `xs`,
  refit the joint Newton mode under the augmented constraint
  `e_i^T x = a` (via [`laplace_mode_fixed_xi`](@ref)) and read the
  constrained Laplace log-marginal. Per-θ density is the renormalised
  ratio `exp(log p̂(y | θ_k, x_i = a) - log p̂(y | θ_k))`; the mixture
  is again renormalised on the grid. Requires `model` and `y`.

If `grid` is supplied it is used verbatim; otherwise a grid of `grid_size`
equally spaced points spanning `±span · √posterior_var` about the posterior
mean is generated.
"""
function posterior_marginal_x(res::INLAResult, i::Integer;
        strategy::Union{Symbol, AbstractMarginalStrategy}=Gaussian(),
        model::Union{Nothing, LatentGaussianModel}=nothing,
        y=nothing,
        grid_size::Integer=75,
        span::Real=5.0,
        grid::Union{Nothing, AbstractVector{<:Real}}=nothing)
    1 ≤ i ≤ length(res.x_mean) ||
        throw(ArgumentError("posterior_marginal_x: index $i out of bounds (1:$(length(res.x_mean)))"))
    s = _resolve_marginal_strategy(strategy)
    if (s isa SimplifiedLaplace || s isa FullLaplace) &&
       (model === nothing || y === nothing)
        throw(ArgumentError("strategy = $(typeof(s).name.name)() requires " *
                            "keyword arguments `model` and `y`"))
    end

    μ = res.x_mean[i]
    σ = sqrt(max(res.x_var[i], 0.0))
    xs = grid === nothing ? _default_grid(μ, σ, grid_size, span) : collect(Float64, grid)

    if s isa FullLaplace
        return (x=xs, pdf=_full_laplace_pdf(res, i, model::LatentGaussianModel, y, xs))
    end

    # Per-θ conditional mean and variance (constraint-corrected).
    m_k = [lp.mode[i] for lp in res.laplaces]
    v_k = [_constrained_marginal_variances(lp.factor, lp.constraint)[i]
           for lp in res.laplaces]

    # Precompute per-θ skewness if requested. `_density_skewness` returns
    # zero for any strategy that does not need skew correction; the inner
    # loop branches on `γ == 0` so the Gaussian path is bit-for-bit
    # equivalent regardless of which strategy is requested.
    γ_k = [_density_skewness(s, res.laplaces[k], model, y, i, v_k[k])
           for k in eachindex(res.laplaces)]

    pdf = zeros(Float64, length(xs))
    @inbounds for k in eachindex(res.laplaces)
        w = res.θ_weights[k]
        w == 0 && continue
        σk = sqrt(max(v_k[k], 0.0))
        σk == 0 && continue
        γ = γ_k[k]
        if γ == 0
            for (j, x) in pairs(xs)
                pdf[j] += w * _normal_pdf(x, m_k[k], σk)
            end
        else
            for (j, x) in pairs(xs)
                t = (x - m_k[k]) / σk
                # Edgeworth first-order skewness correction.
                # H_3(t) = t^3 - 3t.  ∫ φ(t) H_3(t) dt = 0 ⇒ density
                # integrates to 1 without renormalisation.
                c = 1 + γ / 6 * (t^3 - 3 * t)
                # Floor to zero — Edgeworth densities can go slightly
                # negative in the deep tails when |γ| is large. Clamping
                # preserves non-negativity without destroying mass.
                c = max(c, 0.0)
                pdf[j] += w * c * _normal_pdf(x, m_k[k], σk)
            end
        end
    end
    return (x=xs, pdf=pdf)
end

# Per-coordinate density-skew hook. `Gaussian` returns 0 (collapses to
# the unweighted Gaussian mixture). `SimplifiedLaplace` returns the
# Edgeworth coefficient assembled from `∇³_η_log_density`.
_density_skewness(::Gaussian, lp::LaplaceResult, model, y, i, var_i) = 0.0

function _density_skewness(::SimplifiedLaplace,
        lp::LaplaceResult,
        model::LatentGaussianModel,
        y,
        i::Integer,
        var_i::Real)
    return _latent_skewness(lp, model, y, i, var_i)
end

"""
    posterior_marginal_θ(res::INLAResult, j::Int;
                         method = :auto,
                         model = nothing, y = nothing,
                         laplace = Laplace(),
                         grid_size = 75, span = 5.0,
                         grid = nothing) -> @NamedTuple{θ::Vector, pdf::Vector}

Posterior density of the `j`-th hyperparameter on the internal scale,
evaluated on a grid.

### `method`

- `:gaussian` — Gaussian centred at `res.θ̂[j]` with standard deviation
  `√res.Σθ[j,j]`. Bit-for-bit the pre-ADR-046 behaviour; keeps working
  without `model`/`y`.
- `:integrated` — density numerically integrated over the INLA design
  (the analogue of R-INLA's integrated `marginals.hyperpar`):
  - `dim(θ) == 1` — the density is read directly off the retained design
    points: `exp.(res.log_π)`, trapezoid-normalised and interpolated onto
    the output grid (shape-preserving cubic in log-density, linear
    log-density extrapolation beyond the design span). No new Laplace
    fits; requires ≥ 3 retained design points — use
    [`refine_hyperposterior`](@ref) for a denser design.
  - `dim(θ) ≥ 2` — conditional-mode profile slice: `log π̂(θ | y)` is
    evaluated along the Gaussian-conditional path
    `θ(t) = θ̂ + Σθ[:, j] / Σθ[j, j] · (t − θ̂[j])` on an internal profile
    grid and renormalised. Under the local-Gaussian approximation the
    conditional normaliser is `t`-independent, so the renormalised slice
    is the Laplace-approximate marginal. Requires `model` and `y`; each
    profile point costs one warm-started [`laplace_mode`](@ref) fit
    (≤ 21 by default, fewer with tail truncation).
- `:auto` (default) — `:integrated` whenever it is free or enabled
  (`dim(θ) == 1` with ≥ 3 design points; `dim(θ) ≥ 2` with `model` and
  `y` supplied), `:gaussian` otherwise.

The integrated pdf is renormalised on the output grid (trapezoid), so it
integrates to 1 over the returned `θ` range; pass a wide `span` or an
explicit `grid` if tail mass matters. On the internal scale skewed
hyperparameter posteriors (e.g. log-precisions) can differ visibly from
the `:gaussian` method — that difference, plus the Jacobian to the user
scale, is exactly what R-INLA's `summary.hyperpar` reflects.

!!! note "Known limitation: non-identified ridges"
    The `dim(θ) ≥ 2` profile slice captures the skew of `log π̂` along
    the straight conditional-mode path but cannot see posterior mass
    far off that path. For weakly identified models whose
    hyperparameter posterior is a curved ridge (e.g. classical BYM's
    `(τ_v, τ_b)`, Eberly & Carlin 2000), the resulting marginal is
    faithful to the local Laplace picture but can differ materially
    from a full re-integration — R-INLA shows the same sensitivity on
    such models. Prefer BYM2-style identifiable parameterisations when
    the hyperparameter marginal itself is the quantity of interest.

If `grid` is supplied it is used verbatim; otherwise a grid of `grid_size`
equally spaced points spanning `±span · √Σθ[j,j]` about `θ̂[j]` is
generated.
"""
function posterior_marginal_θ(res::INLAResult, j::Integer;
        method::Symbol=:auto,
        model::Union{Nothing, LatentGaussianModel}=nothing,
        y=nothing,
        laplace::Laplace=Laplace(),
        grid_size::Integer=75,
        span::Real=5.0,
        grid::Union{Nothing, AbstractVector{<:Real}}=nothing)
    mθ = length(res.θ̂)
    1 ≤ j ≤ mθ ||
        throw(ArgumentError("posterior_marginal_θ: index $j out of bounds (1:$mθ)"))
    method in (:auto, :integrated, :gaussian) ||
        throw(ArgumentError("posterior_marginal_θ: unknown method :$method; " *
                            "use :auto, :integrated, or :gaussian"))

    μ = res.θ̂[j]
    σ = sqrt(max(res.Σθ[j, j], 0.0))
    θs = grid === nothing ? _default_grid(μ, σ, grid_size, span) : collect(Float64, grid)

    # The 1-D reuse path needs ≥ 3 retained design points for a usable
    # interpolant (tail points can be dropped by `_inla_integrate`); the
    # profile-slice path needs the model and data for fresh Laplace fits.
    can_reuse_1d = mθ == 1 && length(res.θ_points) ≥ 3
    can_profile = mθ ≥ 2 && model !== nothing && y !== nothing
    resolved = method
    if method === :auto
        resolved = (can_reuse_1d || can_profile) ? :integrated : :gaussian
    elseif method === :integrated
        if mθ == 1 && !can_reuse_1d
            throw(ArgumentError("posterior_marginal_θ: method = :integrated " *
                                "needs ≥ 3 retained design points for dim(θ) == 1 " *
                                "(got $(length(res.θ_points))); refit with a denser " *
                                "scheme or use refine_hyperposterior"))
        elseif mθ ≥ 2 && !can_profile
            throw(ArgumentError("posterior_marginal_θ: method = :integrated " *
                                "requires keyword arguments `model` and `y` " *
                                "for dim(θ) ≥ 2"))
        end
    end

    if resolved === :gaussian
        pdf = [_normal_pdf(θ, μ, σ) for θ in θs]
        return (θ=θs, pdf=pdf)
    end

    if mθ == 1
        ld = _integrated_theta_log_density_1d(res, θs)
    else
        # Unreachable when `resolved === :integrated` (checked above);
        # repeated here as a `=== nothing` branch so inference narrows
        # `model` to `LatentGaussianModel` for the call below.
        model === nothing &&
            throw(ArgumentError("posterior_marginal_θ: method = :integrated " *
                                "requires keyword arguments `model` and `y` " *
                                "for dim(θ) ≥ 2"))
        ld = _profile_slice_log_density(res, j, model, y, laplace, θs)
    end
    return (θ=θs, pdf=_normalise_log_density(θs, ld))
end

# 1-D integrated marginal: the retained design points carry the
# unnormalised log posterior directly (`INLAResult.log_π`); sort along
# the single axis and interpolate. Skew-corrected Grid designs are
# asymmetric but still collinear, so sorting is all that's needed.
function _integrated_theta_log_density_1d(res::INLAResult,
        θs::AbstractVector{<:Real})
    ts = [p[1] for p in res.θ_points]
    perm = sortperm(ts)
    return _interp_log_density(ts[perm], res.log_π[perm], θs)
end

# Conditional-mode profile slice for dim(θ) ≥ 2 (ADR-046). Evaluates
# `log π̂(θ | y) = log p̂(y | θ) + log π(θ)` along the Gaussian-conditional
# path through θ̂, sweeping outward from the centre with warm-started
# Newton fits and the same 25-nat tail truncation as the FullLaplace
# per-`x_i` sweep (`full_laplace.jl`). Failure paths mirror
# `_inla_integrate`: bad-θ exceptions and non-finite log-marginals drop
# the point (it stays at `-Inf` and the interpolant bridges past it).
function _profile_slice_log_density(res::INLAResult, j::Integer,
        model::LatentGaussianModel, y, laplace::Laplace,
        θs::AbstractVector{<:Real};
        n_profile::Integer=21)
    lo, hi = extrema(θs)
    ts = collect(range(lo, hi; length=n_profile))
    ld = fill(-Inf, n_profile)

    θ̂ = res.θ̂
    Σjj = res.Σθ[j, j]
    # Direction of the Gaussian-conditional mode path E[θ | θ_j = t]. Σθ
    # is PD by construction (`_safe_inverse_hessian` floors eigenvalues),
    # but guard the division anyway: fall back to varying θ_j alone.
    dir = if Σjj > 0
        res.Σθ[:, j] ./ Σjj
    else
        e = zeros(Float64, length(θ̂))
        e[j] = 1.0
        e
    end

    k_center = argmin(abs.(ts .- θ̂[j]))
    mode0 = _try_profile_fit!(ld, k_center, model, y,
        θ̂ .+ dir .* (ts[k_center] - θ̂[j]), laplace, nothing)
    x_right = mode0
    x_left = mode0 === nothing ? nothing : copy(mode0)

    running_max = ld[k_center]
    for k in (k_center + 1):n_profile
        new_mode = _try_profile_fit!(ld, k, model, y,
            θ̂ .+ dir .* (ts[k] - θ̂[j]), laplace, x_right)
        new_mode === nothing || (x_right = new_mode)
        isfinite(ld[k]) || continue
        running_max = max(running_max, ld[k])
        ld[k] - running_max < -_TRUNCATE_THRESHOLD && break
    end
    running_max = ld[k_center]
    for k in (k_center - 1):-1:1
        new_mode = _try_profile_fit!(ld, k, model, y,
            θ̂ .+ dir .* (ts[k] - θ̂[j]), laplace, x_left)
        new_mode === nothing || (x_left = new_mode)
        isfinite(ld[k]) || continue
        running_max = max(running_max, ld[k])
        ld[k] - running_max < -_TRUNCATE_THRESHOLD && break
    end

    return _interp_log_density(ts, ld, θs)
end

# Single profile-point fit + log-density write. Returns the converged
# mode on success (next warm-start), `nothing` on a dropped point.
function _try_profile_fit!(ld::Vector{Float64}, k::Integer,
        model::LatentGaussianModel, y, θ_t::AbstractVector{Float64},
        laplace::Laplace,
        x_warm::Union{Nothing, AbstractVector{<:Real}})
    local lp
    try
        lp = laplace_mode(model, y, θ_t; strategy=laplace, x0=x_warm)
    catch err
        _is_bad_theta_failure(err) || rethrow(err)
        return nothing
    end
    isfinite(lp.log_marginal) || return nothing
    ld[k] = lp.log_marginal + log_hyperprior(model, θ_t)
    return lp.mode
end

# Normalise a log-density evaluated on `θs` to a trapezoid-unit pdf.
function _normalise_log_density(θs::AbstractVector{<:Real},
        ld::AbstractVector{Float64})
    mx = maximum(ld)
    isfinite(mx) || return zeros(Float64, length(θs))
    pdf = exp.(ld .- mx)
    Z = _trapz(θs, pdf)
    Z > 0 || return zeros(Float64, length(θs))
    return pdf ./ Z
end

# Interpolate a log-density known at scattered knots onto `θs`:
# shape-preserving cubic Hermite (Fritsch-Carlson slopes) inside the
# knot range, linear log-density extrapolation outside — the boundary
# slope is negative for any decaying density, so the extrapolated tails
# decay exponentially instead of ringing. Non-finite knots (dropped or
# truncated evaluations) carry no information and are excluded; the
# interpolant bridges or extrapolates past them.
#
# Knot arguments are concrete `Vector{Float64}` (both callers construct
# them) so `searchsortedlast` dispatch resolves statically — with an
# abstract container JET union-splits it into unrelated container
# methods.
function _interp_log_density(ts::Vector{Float64},
        ld::Vector{Float64},
        θs::AbstractVector{<:Real})
    keep = findall(isfinite, ld)
    length(keep) ≥ 2 || return fill(-Inf, length(θs))
    tsf = ts[keep]
    ldf = ld[keep]
    d = _pchip_slopes(tsf, ldf)
    nf = length(tsf)
    out = zeros(Float64, length(θs))
    for (i, x) in pairs(θs)
        if x ≤ tsf[1]
            out[i] = ldf[1] + d[1] * (x - tsf[1])
        elseif x ≥ tsf[nf]
            out[i] = ldf[nf] + d[nf] * (x - tsf[nf])
        else
            k = searchsortedlast(tsf, x)
            hk = tsf[k + 1] - tsf[k]
            t = (x - tsf[k]) / hk
            h00 = (1 + 2t) * (1 - t)^2
            h10 = t * (1 - t)^2
            h01 = t^2 * (3 - 2t)
            h11 = t^2 * (t - 1)
            out[i] = h00 * ldf[k] + h10 * hk * d[k] +
                     h01 * ldf[k + 1] + h11 * hk * d[k + 1]
        end
    end
    return out
end

# Fritsch-Carlson (PCHIP) slopes on strictly increasing knots: local
# monotonicity is preserved between knots, so the interpolated
# log-density cannot overshoot — no negative-density ringing after exp.
function _pchip_slopes(xs::Vector{Float64}, ys::Vector{Float64})
    n = length(xs)
    d = zeros(Float64, n)
    n == 1 && return d
    h = diff(xs)
    δ = diff(ys) ./ h
    if n == 2
        d[1] = δ[1]
        d[2] = δ[1]
        return d
    end
    for i in 2:(n - 1)
        if δ[i - 1] * δ[i] > 0
            w1 = 2 * h[i] + h[i - 1]
            w2 = h[i] + 2 * h[i - 1]
            d[i] = (w1 + w2) / (w1 / δ[i - 1] + w2 / δ[i])
        end
    end
    d[1] = _pchip_endpoint(h[1], h[2], δ[1], δ[2])
    d[n] = _pchip_endpoint(h[n - 1], h[n - 2], δ[n - 1], δ[n - 2])
    return d
end

# One-sided three-point endpoint slope with the standard monotonicity
# clamps (Fritsch-Carlson).
function _pchip_endpoint(h1::Float64, h2::Float64, δ1::Float64, δ2::Float64)
    d = ((2 * h1 + h2) * δ1 - h1 * δ2) / (h1 + h2)
    if sign(d) != sign(δ1)
        return 0.0
    elseif sign(δ1) != sign(δ2) && abs(d) > 3 * abs(δ1)
        return 3 * δ1
    end
    return d
end

function _default_grid(μ::Real, σ::Real, n::Integer, span::Real)
    σe = σ > 0 ? σ : 1.0
    lo = μ - span * σe
    hi = μ + span * σe
    return collect(range(lo, hi; length=n))
end

_normal_pdf(x::Real, μ::Real, σ::Real) = exp(-0.5 * ((x - μ) / σ)^2) / (σ * sqrt(2π))

# Posterior skewness of x_i under the Laplace at θ (including constraint
# correction). Returns `γ = κ_3(x_i) / σ_i³`.
#
# Derivation: log posterior = log π(x|θ) + log p(y|A x, θ_ℓ). The prior is
# Gaussian so its third derivative is zero. The likelihood factorises over
# observations with η_j = (A x)_j, so
#
#   ∂³ log p / ∂x_a ∂x_b ∂x_c = Σ_j c³_j · A_{ja} A_{jb} A_{jc}
#
# where c³_j = ∇³_η log p evaluated at η̂_j. Contracting with H^{-1} e_i in
# each slot (H^{-1} is the posterior precision inverse, constraint-corrected):
#
#   κ_3(x_i) = Σ_j c³_j · (A u_i)_j³   with   u_i = H^{-1} e_i (constrained)
#
# and σ_i² = (u_i)_i = constraint-corrected marginal variance.
function _latent_skewness(lp::LaplaceResult,
        model::LatentGaussianModel,
        y,
        i::Integer,
        var_i::Real)
    σ_i = sqrt(max(var_i, 0.0))
    σ_i == 0 && return 0.0

    A = as_matrix(model.mapping)
    # `J` is the effective Jacobian `dη/dx` — equals `A` for non-Copy
    # models; otherwise includes per-block β-rows. The skewness contracts
    # `c³` with `(J u)_j³`, so the Copy contribution must enter J here.
    J = joint_effective_jacobian(model, lp.θ)
    η̂ = A * lp.mode
    joint_apply_copy_contributions!(η̂, model, lp.mode, lp.θ)
    c³ = joint_∇³_η_log_density(model, y, η̂, lp.θ)
    all(iszero, c³) && return 0.0

    # u = H⁻¹ e_i (unconstrained). Sparse triangular solve against a unit
    # RHS — one forward + one back substitution.
    e_i = zeros(Float64, length(lp.mode))
    e_i[i] = 1.0
    u = lp.factor \ e_i

    if lp.constraint !== nothing
        U = lp.constraint.U
        W_fact = lp.constraint.W_fact
        C = lp.constraint.C
        u .-= U * (W_fact \ (C * u))
    end

    Ju = J * u
    κ3 = zero(Float64)
    @inbounds for j in eachindex(c³)
        κ3 += c³[j] * Ju[j]^3
    end
    return κ3 / σ_i^3
end
