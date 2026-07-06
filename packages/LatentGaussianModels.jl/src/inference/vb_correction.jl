# Low-rank variational-Bayes mean correction (ADR-048) — van Niekerk &
# Rue (2024, JMLR 25(62)) and van Niekerk, Krainski, Rustand & Rue
# (2023, CSDA 181) §2.3. R-INLA's `control.vb` mean strategy.
#
# At a fixed θ (one integration design point) the Laplace approximation
# `N(μ, Q_X⁻¹)` carries an accurate conditional precision but a
# mode-based mean that is biased for skewed likelihoods. The VBC
# improves the mean as `μ* = μ + Mλ`, where `M` holds `p` selected
# columns of the (constraint-corrected) posterior covariance —
# correcting `p` influential coordinates and letting the GMRF graph
# propagate the change to the full field — and `λ` minimises the
# variational objective
#
#     F(λ) = E_{x ~ N(μ + Mλ, Q_X⁻¹)}[−log p(y | x)]
#            + ½ (μ + Mλ − m₀)ᵀ Q(θ) (μ + Mλ − m₀)
#
# (CSDA eq. 16; `m₀` the prior mean, `Q(θ)` the prior precision). The
# expected log-likelihood factorises over observations and is evaluated
# with probabilists' Gauss-Hermite quadrature against
# `η_i ~ N(η_i^cur, S_i²)`; the Gaussian integration-by-parts (Stein)
# identities give its gradient and curvature in a mean shift `δ_i`
# without any likelihood-derivative evaluations at the nodes:
#
#     ∂/∂δ_i  E[ℓ_i] = Σ_r ω_r (z_r / S_i)        ℓ_i(η_i + S_i z_r)
#     ∂²/∂δ_i² E[ℓ_i] = Σ_r ω_r ((z_r² − 1)/S_i²) ℓ_i(η_i + S_i z_r)
#
# (the `B_i` / `C_i` of JMLR eq. 18). With `δ_η = J M λ` the objective
# is locally quadratic in `λ`; each Newton step is one `p × p` solve,
# re-expanded at the updated mean until the shift converges.
#
# For a Gaussian likelihood `ℓ'` is linear, so `E[ℓ'] = ℓ'(η̂)` exactly
# and the Newton-mode stationarity gives a zero initial gradient:
# `λ = 0` and the correction vanishes identically — the Laplace mean is
# already the exact posterior mean.

"""
    VBMeanCorrection(; n_gh = 15, max_block = 30, indices = nothing)

Configuration for the low-rank variational-Bayes mean correction
(ADR-048; R-INLA's `control.vb` mean strategy), applied per
integration design point when passed to
`INLA(vb_correction = VBMeanCorrection(...))` (or the shorthand
`vb_correction = :mean` for the defaults).

- `n_gh` — Gauss-Hermite nodes per observation for the expected
  log-likelihood (JMLR eq. 17).
- `max_block` — components with `length(c) ≤ max_block` are included
  in the correction set, alongside all `Intercept`/`FixedEffects`
  coordinates ("fixed effects and short random effects", the R-INLA
  default heuristic).
- `indices` — explicit latent coordinates to correct; overrides the
  component heuristic entirely when supplied.

See [`INLA`](@ref) and ADR-048.
"""
Base.@kwdef struct VBMeanCorrection
    n_gh::Int = 15
    max_block::Int = 30
    indices::Union{Nothing, Vector{Int}} = nothing

    function VBMeanCorrection(n_gh::Int, max_block::Int,
            indices::Union{Nothing, Vector{Int}})
        n_gh ≥ 3 ||
            throw(ArgumentError("VBMeanCorrection: n_gh must be ≥ 3, got $n_gh"))
        max_block ≥ 0 ||
            throw(ArgumentError("VBMeanCorrection: max_block must be ≥ 0, got $max_block"))
        if indices !== nothing
            isempty(indices) &&
                throw(ArgumentError("VBMeanCorrection: explicit `indices` must be non-empty"))
            allunique(indices) ||
                throw(ArgumentError("VBMeanCorrection: `indices` must be unique"))
        end
        return new(n_gh, max_block, indices)
    end
end

"""
    _resolve_vb_correction(v) -> Union{Nothing, VBMeanCorrection}

Symbol-or-type shim for the `vb_correction` keyword (mirrors
`_resolve_marginal_strategy`): `:none`/`nothing` disable the
correction, `:mean` selects the defaults, a `VBMeanCorrection` is
passed through.
"""
_resolve_vb_correction(::Nothing) = nothing
_resolve_vb_correction(v::VBMeanCorrection) = v
function _resolve_vb_correction(v::Symbol)
    v === :none && return nothing
    v === :mean && return VBMeanCorrection()
    throw(ArgumentError("unknown vb_correction :$v; " *
                        "use :none, :mean, or a VBMeanCorrection instance"))
end

# Default correction set: every coordinate of fixed-effect-shaped
# components plus any component short enough to be "connected to many
# datapoints" in the R-INLA sense.
function _default_vb_indices(m::LatentGaussianModel, max_block::Int)
    idx = Int[]
    for (i, c) in enumerate(m.components)
        if c isa Intercept || c isa FixedEffects || length(c) ≤ max_block
            append!(idx, m.latent_ranges[i])
        end
    end
    return idx
end

function _vb_indices(vb::VBMeanCorrection, m::LatentGaussianModel)
    vb.indices === nothing && return _default_vb_indices(m, vb.max_block)
    all(i -> 1 ≤ i ≤ m.n_x, vb.indices) ||
        throw(ArgumentError("VBMeanCorrection: indices out of bounds (1:$(m.n_x))"))
    return vb.indices
end

# Gauss-Hermite gradient/curvature of the expected log-likelihood in a
# per-observation mean shift (JMLR eq. 18's B and C). Observations with
# S_i ≈ 0 (fully pinned linear predictors) fall back to the plain
# pointwise derivatives — the S → 0 limit of the Stein identities.
# Observations producing non-finite node values (GH tail nodes pushing
# the likelihood out of domain) are dropped from the correction
# information; the remaining observations carry the solve.
function _vb_gh_moments(m::LatentGaussianModel, y,
        η::Vector{Float64}, θ::AbstractVector{<:Real},
        S::Vector{Float64},
        z_nodes::Vector{Float64}, ω_nodes::Vector{Float64})
    n = length(η)
    B = zeros(Float64, n)
    C = zeros(Float64, n)
    η_node = similar(η)
    small = 1.0e-8
    ℓ_node = Vector{Vector{Float64}}(undef, length(z_nodes))
    for (r, zr) in enumerate(z_nodes)
        @. η_node = η + S * zr
        ℓ_node[r] = Vector{Float64}(joint_pointwise_log_density(m, y, η_node, θ))
    end
    @inbounds for i in 1:n
        if S[i] ≤ small
            continue    # handled by the derivative fallback below
        end
        ok = true
        b = 0.0
        c = 0.0
        for (r, zr) in enumerate(z_nodes)
            ℓ = ℓ_node[r][i]
            if !isfinite(ℓ)
                ok = false
                break
            end
            b += ω_nodes[r] * (zr / S[i]) * ℓ
            c += ω_nodes[r] * ((zr^2 - 1.0) / S[i]^2) * ℓ
        end
        if ok
            B[i] = b
            C[i] = c
        end
    end
    if any(≤(small), S)
        g1 = joint_∇_η_log_density(m, y, η, θ)
        g2 = joint_∇²_η_log_density(m, y, η, θ)
        @inbounds for i in 1:n
            if S[i] ≤ small
                B[i] = g1[i]
                C[i] = g2[i]
            end
        end
    end
    return B, C
end

"""
    _vb_mean_shift(vb::VBMeanCorrection, lp::LaplaceResult,
                   model::LatentGaussianModel, y) -> Vector{Float64}

Low-rank VB mean shift `Mλ` at the Laplace fit `lp` (one integration
design point). `M` is the kriging-projected posterior covariance
restricted to the correction columns, so hard linear constraints are
preserved exactly: `C (x̂ + Mλ) = C x̂ = e`. Returns a zero vector when
the correction set is empty.
"""
function _vb_mean_shift(vb::VBMeanCorrection, lp::LaplaceResult,
        m::LatentGaussianModel, y)
    I_set = _vb_indices(vb, m)
    n_x = length(lp.mode)
    isempty(I_set) && return zeros(Float64, n_x)

    A = as_matrix(m.mapping)
    J = joint_effective_jacobian(m, lp.θ)
    η̂ = A * lp.mode
    joint_apply_copy_contributions!(η̂, m, lp.mode, lp.θ)

    # σ²_η under the constraint-corrected Laplace — the same multi-RHS
    # solve as `_sla_mean_shift`.
    Z = lp.factor \ Matrix(transpose(J))
    if lp.constraint !== nothing
        U = lp.constraint.U
        W_fact = lp.constraint.W_fact
        C_c = lp.constraint.C
        Z .-= U * (W_fact \ (C_c * Z))
    end
    σ²_η = vec(sum(J .* transpose(Z), dims=2))
    S = sqrt.(max.(σ²_η, 0.0))

    # Propagation matrix M: constrained covariance columns of the
    # correction set.
    E_I = zeros(Float64, n_x, length(I_set))
    for (k, i) in enumerate(I_set)
        E_I[i, k] = 1.0
    end
    M = lp.factor \ E_I
    if lp.constraint !== nothing
        U = lp.constraint.U
        W_fact = lp.constraint.W_fact
        C_c = lp.constraint.C
        M .-= U * (W_fact \ (C_c * M))
    end
    JM = J * M

    # Prior pieces. Q(θ) may be intrinsic/singular — the likelihood
    # curvature keeps the Newton system usable, and a tiny scaled ridge
    # guards the corner cases.
    Q = joint_precision(m, lp.θ)
    m₀ = joint_prior_mean(m, lp.θ)
    G_prior = Matrix(transpose(M) * (Q * M))

    gh = FastGaussQuadrature.gausshermite(vb.n_gh)
    z_nodes = gh[1] .* sqrt(2)
    ω_nodes = gh[2] ./ sqrt(π)

    p = length(I_set)
    λ_tot = zeros(Float64, p)
    shift = zeros(Float64, n_x)
    η_cur = copy(η̂)
    scale = 1.0 + maximum(abs, lp.mode)
    tol = 1.0e-9 * scale
    max_iter = 20
    for _ in 1:max_iter
        B, C = _vb_gh_moments(m, y, η_cur, lp.θ, S, z_nodes, ω_nodes)
        r = lp.mode .+ shift .- m₀
        g = transpose(JM) * (.-B) .+ transpose(M) * (Q * r)
        H = transpose(JM) * (Diagonal(max.(.-C, 0.0)) * JM) .+ G_prior
        ridge = 1.0e-10 * max(1.0, maximum(abs, view(H, diagind(H))))
        H[diagind(H)] .+= ridge
        local step
        try
            step = -(Symmetric(H) \ g)
        catch err
            _is_bad_theta_failure(err) || rethrow(err)
            break    # keep the shift accumulated so far
        end
        all(isfinite, step) || break
        δ = M * step
        λ_tot .+= step
        shift .+= δ
        η_cur .+= JM * step
        maximum(abs, δ) < tol && break
    end

    return shift
end
