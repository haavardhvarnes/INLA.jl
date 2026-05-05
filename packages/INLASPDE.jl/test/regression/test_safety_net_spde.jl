# Phase M PR-4 — SPDE2 triangulation against the LGM safety net.
#
# The safety-net code lives in `LatentGaussianModels.jl` and has its own
# closed-form regression suite. This file is the SPDE2-side
# triangulation: it confirms that the safety net activates correctly
# when the bad-θ trigger is the actual numerical pathology that
# motivated the patch — Cholesky failing on a numerically-singular
# `Q_SPDE(τ, κ)` at extreme `(log τ, log κ)`.
#
# Pre-PR-4 this scenario crashed the outer LBFGS line search with a
# `PosDefException` from inside the inner Newton; the closure now
# returns the smooth quadratic penalty `1e10 + 1e3·‖θ‖²` instead.

using LatentGaussianModels
using LatentGaussianModels: GaussianLikelihood, Intercept,
                            LatentGaussianModel, Laplace, INLA, fit
using SparseArrays
using Random

const LGM = LatentGaussianModels

@testset "SPDE2 — extreme (log τ, log κ) activates LGM safety net" begin
    # Tiny SPDE2 fixture from the FEM unit-square mesh + a synthetic
    # Gaussian outcome on the four vertices.
    pc = PCMatern(range_U=0.5, range_α=0.05,
        sigma_U=1.0, sigma_α=0.01)
    spde = SPDE2(POINTS_SQ, TRIS_SQ; pc=pc)

    rng = Random.Xoshiro(20260505)
    n_obs = 8
    A_obs = sparse(rand(rng, n_obs, 4))
    intercept = Intercept()
    A = hcat(ones(n_obs, 1), A_obs)
    y = randn(rng, n_obs)

    like = GaussianLikelihood()
    model = LatentGaussianModel(like, (intercept, spde), A)

    f = LGM._neg_log_posterior_θ(model, y, Laplace())

    # Sane θ — finite objective inside the normal range.
    val_normal = f([0.0, 0.0, 0.0], nothing)
    @test isfinite(val_normal)
    @test val_normal < 1.0e9

    # Extreme `log κ = 200` overflows `κ⁴ ≈ exp(800)` inside the SPDE
    # precision; cholesky on the resulting non-finite Q raises
    # `PosDefException`. The closure must return the smooth penalty.
    θ_extreme = [0.0, 0.0, 200.0]   # [log τ_obs, log τ_field, log κ_field]
    val_extreme = f(θ_extreme, nothing)
    @test isfinite(val_extreme)
    @test val_extreme > 1.0e9
    @test val_extreme≈1.0e10 + 1.0e3 * sum(abs2, θ_extreme) atol=1.0e-6

    # End-to-end: starting LBFGS from the bad-θ region, the safety net
    # must let the optimizer escape and the integration stage must
    # produce a finite `log_marginal`.
    res = fit(model, y, INLA(; θ0=[0.0, 0.0, 200.0],
        int_strategy=:grid))
    @test isfinite(res.log_marginal)
    @test all(isfinite, res.x_mean)
    @test all(isfinite, res.x_var)
    @test all(≥(0), res.x_var)
    @test abs(res.θ̂[3]) < 50.0   # `log κ` recovered well inside feasible region
end
