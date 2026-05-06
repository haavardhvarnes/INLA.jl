# Triangulation tier — Scotland BYM2 (Phase P PR-1 v1.0 tightening).
#
# Same model R-INLA fits in `scripts/generate-fixtures/lgm/scotland_bym2.R`,
# fit two ways:
#   1. Julia INLA via `inla(model, y; int_strategy = :grid)`.
#   2. NUTS over the same INLA hyperparameter posterior via
#      `nuts_sample(model, y, n_samples; init_from_inla = inla_fit)`.
#
# We then assert that the two posterior summaries on θ agree —
# tier-3 cross-validation against the upstream R-INLA fit is covered by
# the existing oracle in LGM (`test/oracle/test_scotland_bym2.jl`); this
# test only certifies that NUTS on the *same* `INLALogDensity` recovers
# the same θ posterior INLA's grid integration sees. That's enough to
# catch a regression in either side: `INLALogDensity` evaluation, the
# Laplace mode finder, the inner Newton, the BYM2 precision build, or
# the NUTS bridge itself.
#
# Tolerances (Phase P PR-1 / ADR-044 — uniform v1.0 levels):
#   tol_mean = 1.5 SDs (envelope INLA / NUTS within 1.5σ)
#   tol_sd   = 0.60 (relative)
#
# `tol_sd` is set wide enough to absorb a structural feature, not MC
# noise: NUTS finds 30–45% wider posteriors than INLA on the BYM2
# mixing weight `logit φ` (`BYM2[3][2]`) on both Scotland and PA. The
# means agree within 0.06–0.23 SDs — well below `tol_mean`. The SD gap
# is inherent to grid integration vs. full HMC: INLA's grid is bounded
# by the Hessian-explored region around the Laplace mode, while NUTS
# samples the full posterior including the heavy tails of weakly
# identified mixing parameters. Tightening `tol_sd` would force per-
# parameter tolerances or n=10k+ chains — neither is the right tier-3
# contract. The 0.60 envelope still flags genuine regressions
# (gradient bug → SDs blow up; precision build bug → both shift).
#
# Chain length: 1000 post-warmup after 200 warmup. Wall-clock ~30 s on
# `Apple M2 Pro / 16 GB`.
#
# Fixture is the LGM-side JLD2 (`packages/LatentGaussianModels.jl/test/
# oracle/fixtures/scotland_bym2.jld2`). If absent, skip transparently
# rather than block the suite for users without the R toolchain.

using Test
using LinearAlgebra: I
using SparseArrays: sparse
using Random
using LatentGaussianModels: PoissonLikelihood, Intercept, FixedEffects,
                            BYM2, LatentGaussianModel, inla, PCPrecision
using GMRFs: GMRFGraph
using LGMTuring: nuts_sample, compare_posteriors

const SCOTLAND_FIXTURE_PATH = joinpath(@__DIR__, "..", "..", "..",
    "LatentGaussianModels.jl", "test", "oracle", "fixtures",
    "scotland_bym2.jld2")

@testset "triangulation — Scotland BYM2 (INLA vs NUTS)" begin
    if !isfile(SCOTLAND_FIXTURE_PATH)
        @test_skip "Scotland BYM2 fixture missing at $SCOTLAND_FIXTURE_PATH"
    else
        using JLD2
        fx = jldopen(SCOTLAND_FIXTURE_PATH, "r") do f
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
                    @info "Scotland NUTS-vs-INLA disagreement" name=r.name inla_mean=r.inla_mean inla_sd=r.inla_sd nuts_mean=r.nuts_mean nuts_sd=r.nuts_sd mean_abs_diff=r.mean_abs_diff sd_rel_diff=r.sd_rel_diff
                end
                @test !r.flagged
            end
        end
    end
end
