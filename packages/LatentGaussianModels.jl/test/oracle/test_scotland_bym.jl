# Oracle test: classical (non-reparametrised) BYM on Scotland lip
# cancer vs R-INLA. Companion to test_scotland_bym2.jl that exercises
# the `model = "bym"` pathway with separate τ_v (iid) and τ_b (besag).
#
# Checked quantities:
#   - posterior means of fixed effects, 7% relative tolerance
#   - τ_b (Besag precision): Julia's point estimate (`exp(θ̂[2])`) must
#     fall inside R-INLA's 95% credible interval, and the integrated
#     θ-marginal (ADR-046) must be two-sided-consistent with R-INLA's
#     summary (Julia median inside R's CI; R's mean inside Julia's 95%
#     band). Tight point parity is impossible on this fixture — the
#     classical-BYM (τ_v, τ_b) posterior is a non-identified ridge and
#     the implementations distribute ridge mass differently; see the
#     comment above the @testset.
#   - τ_v (IID precision) only sanity-checked: known weakly identified
#     in classical BYM, posterior is heavy-tailed.
#   - mlik: 2% relative tolerance (matches R-INLA's `extra()` for F_BYM).
#
# Fixture: scripts/generate-fixtures/lgm/scotland_bym.R.

include("load_fixture.jl")

using Test
using SparseArrays
using LinearAlgebra: I
using LatentGaussianModels: PoissonLikelihood, Intercept, FixedEffects,
                            BYM, LatentGaussianModel, inla, PCPrecision,
                            fixed_effects, hyperparameters, log_marginal_likelihood,
                            posterior_marginal_θ
using GMRFs: GMRFGraph

const BYM_FIXTURE = "scotland_bym"

const BYM_FIXED_EFFECT_TOL = 0.07
# Classical BYM is non-identified (Eberly & Carlin 2000): only τ_b/τ_v
# is constrained by data, posteriors on each are heavy-tailed. On
# Scotland (K=4 connected components: 53-node main + 3 island singletons)
# the τ_b summaries Julia's *mode-based* point estimate and R-INLA
# expose are not the same statistic: `exp(θ̂[2])` is the value at the
# mode of the *internal-scale* (log τ_b) marginal, while R-INLA's
# `summary.hyperpar` columns (`mode`, `0.5quant`, `mean`) are all
# summaries of the *user-scale* density `p(τ_b)`, related by the
# Jacobian `p(τ_b) = (1/τ_b) p_int(log τ_b)`. The Jacobian alone shifts
# a user-scale mode below `exp(mode_int)` by a factor near `exp(-σ²)`,
# so even a perfectly-correct Julia fit shows a 50–80 % gap on a
# direct mode-to-mean comparison.
#
# The integrated θ-marginal accessor (ADR-046) removes the
# statistic-mismatch problem but exposed the deeper one: the (τ_v, τ_b)
# posterior is a non-identified ridge, and the two implementations
# legitimately distribute ridge mass differently. The fixture's own
# τ_v row (mean 1005, sd 7756, mode 15) shows R-INLA integrating far
# into the τ_v → ∞ / τ_b-small end of the ridge; Julia's design stays
# nearer the mode, leaving its τ_b marginal centred higher (~2×). Fixed
# effects, mlik, and the BYM2 fixtures (identified τ) are unaffected —
# tight τ parity gates live there. Here the assertions are the
# statistic-independent, tail-robust pair:
#   1. `exp(θ̂[2])` inside R-INLA's 95 % credible interval, and
#   2. two-sided consistency of the *integrated user-scale marginal*:
#      Julia's median inside R-INLA's 95 % CI, and R-INLA's mean inside
#      Julia's central 95 % band.
# Per-CC Sørbye-Rue scaling (Freni-Sterrantino et al. 2018) is
# implemented and is c-invariant under PCPrecision priors, so it
# doesn't shift τ_b on its own.
const BYM_MLIK_REL_TOL = 0.02

_rel_bym(a, b) = abs(a - b) / max(abs(b), 1.0)

function _bym_row(frame, name::AbstractString, col::AbstractString)
    rn_raw = frame["rownames"]
    rn = rn_raw isa AbstractString ? [String(rn_raw)] : String.(rn_raw)
    idx = findfirst(==(name), rn)
    idx === nothing && error("row '$name' not found (have: $(join(rn, "; ")))")
    col_raw = frame[col]
    return col_raw isa Real ? Float64(col_raw) : Float64(col_raw[idx])
end

@testset "scotland_bym vs R-INLA" begin
    if !has_oracle_fixture(BYM_FIXTURE)
        @test_skip "oracle fixture $BYM_FIXTURE not generated (see scripts/generate-fixtures/)"
    else
        fx = load_oracle_fixture(BYM_FIXTURE)
        @test fx["name"] == BYM_FIXTURE

        if !haskey(fx, "input")
            @test_skip "fixture has no `input` field — regenerate with the current R script"
        else
            inp = fx["input"]
            y = Int.(inp["cases"])
            E = Float64.(inp["expected"])
            x = Float64.(inp["x"])
            W = inp["W"]
            n = length(y)

            # Latent layout: [α; β; v; b]. Predictor η_i = α + β x_i + v_i + b_i.
            ℓ = PoissonLikelihood(; E=E)
            c_int = Intercept()
            c_beta = FixedEffects(1)
            c_bym = BYM(GMRFGraph(W);
                hyperprior_iid=PCPrecision(1.0, 0.01),
                hyperprior_besag=PCPrecision(1.0, 0.01))
            A = sparse(hcat(
                ones(n),                        # intercept → α
                reshape(x, n, 1),               # AFF slope → β
                Matrix{Float64}(I, n, n),       # v_i contribution
                Matrix{Float64}(I, n, n)       # b_i contribution
            ))
            model = LatentGaussianModel(ℓ, (c_int, c_beta, c_bym), A)

            res = inla(model, y; int_strategy=:grid)

            # --- Fixed effects --------------------------------------------
            sf = fx["summary_fixed"]
            α_R = _bym_row(sf, "(Intercept)", "mean")
            β_R = _bym_row(sf, "x", "mean")
            fe = fixed_effects(model, res)
            @test length(fe) == 2
            @test _rel_bym(fe[1].mean, α_R) < BYM_FIXED_EFFECT_TOL
            @test _rel_bym(fe[2].mean, β_R) < BYM_FIXED_EFFECT_TOL

            # --- Hyperparameters ------------------------------------------
            # Internal θ = [log τ_v, log τ_b]. R-INLA reports user-scale
            # precisions. τ_v is weakly identified in classical BYM —
            # only sanity-check finiteness. For τ_b, the ridge-robust
            # checks described in the BYM_TAU_B notes above: CI
            # containment of the mode-based point estimate, plus the
            # two-sided consistency of the integrated θ marginal
            # (ADR-046).
            sh = fx["summary_hyperpar"]
            τ_b_R_lo = _bym_row(
                sh, "Precision for region (spatial component)", "0.025quant")
            τ_b_R_hi = _bym_row(
                sh, "Precision for region (spatial component)", "0.975quant")
            τ_b_J = exp(res.θ̂[2])
            @test τ_b_R_lo ≤ τ_b_J ≤ τ_b_R_hi

            τ_b_mean_R = _bym_row(
                sh, "Precision for region (spatial component)", "mean")
            mθ = posterior_marginal_θ(res, 2; model=model, y=y)
            τs = exp.(mθ.θ)
            pτ = mθ.pdf ./ τs
            τ_b_med_J = marginal_quantile(τs, pτ, 0.5)
            τ_b_lo_J = marginal_quantile(τs, pτ, 0.025)
            τ_b_hi_J = marginal_quantile(τs, pτ, 0.975)
            @test τ_b_R_lo ≤ τ_b_med_J ≤ τ_b_R_hi
            @test τ_b_lo_J ≤ τ_b_mean_R ≤ τ_b_hi_J

            hp = hyperparameters(model, res)
            @test length(hp) == 2
            @test all(isfinite(r.mean) && r.sd > 0 for r in hp)

            # --- Marginal log-likelihood ----------------------------------
            mlik_R = Float64(fx["mlik"][1])
            mlik_J = log_marginal_likelihood(res)
            @test _rel_bym(mlik_J, mlik_R) < BYM_MLIK_REL_TOL
        end
    end
end
