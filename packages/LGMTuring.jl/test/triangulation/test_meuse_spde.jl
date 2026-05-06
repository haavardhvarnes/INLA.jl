# Triangulation tier — Meuse SPDE (Phase P PR-1 v1.0).
#
# The geostatistics flagship: Gaussian likelihood + Intercept +
# FixedEffects + SPDE2 (α = 2). 3-D θ posterior on
# (log τ_noise, log τ_spde, log κ_spde). Latent dim ≈ 355
# (n_obs = 155 + n_v ≈ 200). The role of this test is to certify that
# NUTS recovers INLA's grid posterior on the SPDE precision and the
# noise precision *jointly* — the LGM oracle (`test_meuse_spde.jl`)
# already confirms agreement with R-INLA, but only on the marginal
# means; the cross-engine HMC test verifies the joint posterior shape.
#
# Tolerances (Phase P PR-1 / ADR-044 — uniform v1.0 levels):
#   tol_mean = 1.5 SDs
#   tol_sd   = 0.60 (relative; uniform with Scotland / PA — see
#                    `test_scotland_bym2.jl` for the structural
#                    rationale on weakly-identified hyperparameters)
# Chain length: 1000 post-warmup after 200 warmup. Wall-clock ~10 min
# on `Apple M2 Pro / 16 GB` — each leapfrog step costs a full Laplace
# fit of a 355-dim latent. **This test is gated behind
# `--triangulation`** in `runtests.jl` for that reason; CI runs it
# nightly, not on every PR.
#
# Fixture is the INLASPDE-side JLD2
# (`packages/INLASPDE.jl/test/oracle/fixtures/meuse_spde.jld2`).
# If absent, skip transparently.

using Test
using LinearAlgebra: I
using SparseArrays: SparseMatrixCSC
using Random
using LatentGaussianModels: GaussianLikelihood, Intercept, FixedEffects,
                            LatentGaussianModel, inla, PCPrecision
using INLASPDE: SPDE2, PCMatern
using LGMTuring: nuts_sample, compare_posteriors

const MEUSE_FIXTURE_PATH = joinpath(@__DIR__, "..", "..", "..",
    "INLASPDE.jl", "test", "oracle", "fixtures", "meuse_spde.jld2")

@testset "triangulation — Meuse SPDE (INLA vs NUTS)" begin
    if !isfile(MEUSE_FIXTURE_PATH)
        @test_skip "Meuse SPDE fixture missing at $MEUSE_FIXTURE_PATH"
    else
        using JLD2
        fxt = jldopen(MEUSE_FIXTURE_PATH, "r") do f
            f["fixture"]
        end

        y = Float64.(fxt["input"]["y"])
        dist_cov = Float64.(fxt["input"]["dist"])
        points = fxt["mesh"]["loc"]::Matrix{Float64}
        tv = fxt["mesh"]["tv"]::Matrix{Int}
        A_field = SparseMatrixCSC{Float64, Int}(fxt["A_field"])

        n_obs = length(y)
        n_v = size(points, 1)

        spde = SPDE2(points, tv; α=2,
            pc=PCMatern(
                range_U=0.5, range_α=0.5,
                sigma_U=1.0, sigma_α=0.5
            ))
        # Match the LGM oracle: proper N(0, 1000) intercept (R-INLA's
        # `prec.intercept = 1e-3`), proper N(0, 1000) covariate slope.
        intercept = Intercept(prec=1.0e-3, improper=false)
        beta_dist = FixedEffects(1; prec=1.0e-3)

        A_intercept = ones(n_obs, 1)
        A_dist = reshape(dist_cov, n_obs, 1)
        A = hcat(A_intercept, A_dist, A_field)

        like = GaussianLikelihood(hyperprior=PCPrecision(1.0, 0.01))
        model = LatentGaussianModel(like, (intercept, beta_dist, spde), A)

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
        @test length(rows) == 3  # log τ_noise, log τ_spde, log κ_spde
        for r in rows
            @test isfinite(r.inla_mean) && isfinite(r.nuts_mean)
            @test isfinite(r.inla_sd) && isfinite(r.nuts_sd)
            if r.flagged
                @info "Meuse NUTS-vs-INLA disagreement" name=r.name inla_mean=r.inla_mean inla_sd=r.inla_sd nuts_mean=r.nuts_mean nuts_sd=r.nuts_sd mean_abs_diff=r.mean_abs_diff sd_rel_diff=r.sd_rel_diff
            end
            @test !r.flagged
        end
    end
end
