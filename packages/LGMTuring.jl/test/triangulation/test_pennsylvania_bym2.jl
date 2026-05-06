# Triangulation tier — Pennsylvania BYM2 (Phase P PR-1 v1.0).
#
# Sibling of `test_scotland_bym2.jl`. Same Poisson-BYM2 architecture,
# larger n (67 counties), proper covariate (standardised smoking rate),
# and indirect-standardised expected counts. The role of this test is
# to certify that NUTS recovers the same θ posterior INLA's grid sees
# on a *second* canonical Poisson-BYM2 dataset — Scotland alone could
# fail to surface a regression that only manifests on a larger graph
# or a non-trivial fixed effect.
#
# Tolerances (Phase P PR-1 / ADR-044 — uniform v1.0 levels):
#   tol_mean = 1.5 SDs
#   tol_sd   = 0.60 (relative; see Scotland test for the structural
#                    rationale — NUTS sees 30–45% wider posterior on
#                    the BYM2 logit-φ mixing weight than INLA's grid).
# Chain length: 1000 post-warmup after 200 warmup. Wall-clock ~40 s on
# `Apple M2 Pro / 16 GB` (≈ 1.3× Scotland; 67 counties vs 56 + denser W).
#
# Fixture is the LGM-side JLD2 (`packages/LatentGaussianModels.jl/test/
# oracle/fixtures/pennsylvania_bym2.jld2`). If absent, skip transparently
# rather than block the suite for users without the R toolchain.

using Test
using LinearAlgebra: I
using SparseArrays: sparse
using Random
using LatentGaussianModels: PoissonLikelihood, Intercept, FixedEffects,
                            BYM2, LatentGaussianModel, inla, PCPrecision
using GMRFs: GMRFGraph
using LGMTuring: nuts_sample, compare_posteriors

const PA_FIXTURE_PATH = joinpath(@__DIR__, "..", "..", "..",
    "LatentGaussianModels.jl", "test", "oracle", "fixtures",
    "pennsylvania_bym2.jld2")

@testset "triangulation — Pennsylvania BYM2 (INLA vs NUTS)" begin
    if !isfile(PA_FIXTURE_PATH)
        @test_skip "Pennsylvania BYM2 fixture missing at $PA_FIXTURE_PATH"
    else
        using JLD2
        fx = jldopen(PA_FIXTURE_PATH, "r") do f
            f["fixture"]
        end
        if !haskey(fx, "input")
            @test_skip "fixture has no `input` field — regenerate"
        else
            inp = fx["input"]
            y = Int.(inp["cases"])
            E = Float64.(inp["expected"])
            x = Float64.(inp["x"])
            W = inp["W"]
            n = length(y)

            ℓ = PoissonLikelihood(; E=E)
            c_int = Intercept()
            c_beta = FixedEffects(1)
            c_bym2 = BYM2(GMRFGraph(W); hyperprior_prec=PCPrecision(1.0, 0.01))
            A = sparse(hcat(
                ones(n),
                reshape(x, n, 1),
                Matrix{Float64}(I, n, n),
                zeros(n, n)
            ))
            model = LatentGaussianModel(ℓ, (c_int, c_beta, c_bym2), A)

            inla_fit = inla(model, y; int_strategy=:grid)

            chain = nuts_sample(model, y, 1000;
                n_adapts=200,
                init_from_inla=inla_fit,
                rng=Random.Xoshiro(20260506),
                progress=false)

            rows = compare_posteriors(inla_fit, chain;
                model=model,
                tol_mean=1.5,
                tol_sd=0.60)
            @test length(rows) == 2  # log τ, logit φ
            for r in rows
                @test isfinite(r.inla_mean) && isfinite(r.nuts_mean)
                @test isfinite(r.inla_sd) && isfinite(r.nuts_sd)
                if r.flagged
                    @info "Pennsylvania NUTS-vs-INLA disagreement" name=r.name inla_mean=r.inla_mean inla_sd=r.inla_sd nuts_mean=r.nuts_mean nuts_sd=r.nuts_sd mean_abs_diff=r.mean_abs_diff sd_rel_diff=r.sd_rel_diff
                end
                @test !r.flagged
            end
        end
    end
end
