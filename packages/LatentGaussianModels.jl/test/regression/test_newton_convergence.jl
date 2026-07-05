using LatentGaussianModels: PoissonLikelihood, GaussianLikelihood, Intercept, IID,
                            LatentGaussianModel, inla, laplace_mode, Laplace,
                            PCPrecision
using Distributions: Poisson
using SparseArrays
using LinearAlgebra: I
using Random

# Regression guard for the inner Newton solver (see ADR-010, laplace.jl).
# The point is not statistical accuracy — other tests cover that — but that
# every Laplace fit the INLA pipeline performs actually reaches its
# convergence flag, and that the solver stays robust at the extreme θ the
# integration grid can visit. A regression that broke Newton (bad step,
# wrong stopping rule, constraint mishandling) would surface here.

@testset "Newton convergence — every design-point Laplace converges" begin
    rng = Random.Xoshiro(20260705)

    maxiter = Laplace().maxiter

    # 1. Gaussian + Intercept + IID (2 hyperparameters).
    n = 150
    y_g = 1.0 .+ 0.5 .* randn(rng, n)
    A_g = sparse([ones(n) Matrix{Float64}(I, n,n)])
    model_g = LatentGaussianModel(GaussianLikelihood(),
        (Intercept(), IID(n; hyperprior=PCPrecision(1.0, 0.01))), A_g)

    # 2. Poisson + Intercept + IID.
    τ_true = 4.0
    u = randn(rng, n) ./ sqrt(τ_true)
    y_p = [rand(rng, Poisson(exp(0.5 + u[i]))) for i in 1:n]
    model_p = LatentGaussianModel(PoissonLikelihood(),
        (Intercept(), IID(n; hyperprior=PCPrecision(1.0, 0.01))), A_g)

    for (name, model, y) in (("Gaussian", model_g, y_g), ("Poisson", model_p, y_p))
        res = inla(model, y; int_strategy=:grid)
        @test !isempty(res.laplaces)
        # Every stored per-θ-point Laplace must have converged within budget.
        @test all(lp.converged for lp in res.laplaces)
        @test all(lp.iterations <= maxiter for lp in res.laplaces)
        @test all(lp.iterations >= 1 for lp in res.laplaces)
    end
end

@testset "Newton convergence — robust at extreme fixed θ" begin
    # Newton must converge across the whole θ range the integration scheme
    # can probe, including near-degenerate precisions. We fit at fixed θ
    # far above and below the mode and require convergence each time.
    rng = Random.Xoshiro(11)
    n = 120
    u = randn(rng, n) ./ 2
    y = [rand(rng, Poisson(exp(0.3 + u[i]))) for i in 1:n]
    A = sparse([ones(n) Matrix{Float64}(I, n,n)])
    model = LatentGaussianModel(PoissonLikelihood(),
        (Intercept(), IID(n; hyperprior=PCPrecision(1.0, 0.01))), A)

    # Poisson + Intercept + IID has a single hyperparameter (the IID
    # log-precision); the likelihood and intercept carry none. Sweep a
    # wide, demanding range for that one coordinate.
    for logτ in (-6.0, -3.0, 0.0, 3.0, 6.0, 9.0)
        lp = laplace_mode(model, y, [logτ])
        @test lp.converged
        @test isfinite(lp.log_marginal)
    end
end
