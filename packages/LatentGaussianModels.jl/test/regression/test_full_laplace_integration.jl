# FullLaplace integration-stage summaries (ADR-026 "PR-4" follow-up,
# Tier-4 item 16.2): `INLA(latent_strategy = FullLaplace())` replaces
# `INLAResult.x_mean` / `x_var` with per-coordinate trapezoid moments
# of the refitted-Laplace mixture.
#
# Tier-1 coverage:
#  - struct parameters + validation,
#  - Gaussian-likelihood collapse onto the Gaussian strategy (the joint
#    posterior is exactly Gaussian, so the refit moments must reproduce
#    the pilot summaries up to grid/trapezoid error),
#  - θ-stage invariance (the post-pass must not touch θ̂, weights, mlik),
#  - skewed Poisson posterior: FL shifts the summary, and the summary
#    equals the moments of the `posterior_marginal_x` FullLaplace
#    density on the same grid,
#  - 0-hyperparameter fast path.
#
# Accuracy against R-INLA `strategy = "laplace"` is tier-2 territory —
# see the integration-stage block in `test/oracle/test_synthetic_brunei.jl`.

using LatentGaussianModels: GaussianLikelihood, PoissonLikelihood, Intercept,
                            IID, LatentGaussianModel, inla,
                            posterior_marginal_x, Gaussian, FullLaplace
using Distributions: Poisson
using SparseArrays
using LinearAlgebra: I
using Random

function _fl_trapz(xs::AbstractVector, ys::AbstractVector)
    s = 0.0
    @inbounds for i in 1:(length(xs) - 1)
        s += 0.5 * (xs[i + 1] - xs[i]) * (ys[i] + ys[i + 1])
    end
    return s
end

@testset "FullLaplace struct parameters" begin
    fl = FullLaplace()
    @test fl.n_grid == 51
    @test fl.span == 5.0
    fl2 = FullLaplace(n_grid=31, span=4.0)
    @test fl2.n_grid == 31
    @test fl2.span == 4.0
    @test_throws ArgumentError FullLaplace(n_grid=2)
    @test_throws ArgumentError FullLaplace(span=0.0)
end

@testset "Gaussian likelihood — FL summaries collapse to Gaussian" begin
    rng = Random.Xoshiro(20260706)
    n = 16
    y = 0.4 .+ randn(rng, n)
    A = sparse(I, n, n)
    model = LatentGaussianModel(GaussianLikelihood(), (IID(n),), A)

    res_g = inla(model, y; int_strategy=:grid)
    res_f = inla(model, y; int_strategy=:grid,
        latent_strategy=FullLaplace())

    # θ-stage outputs must be untouched by the post-pass (the two fits
    # are deterministic and identical up to the summary replacement).
    @test res_f.θ̂ == res_g.θ̂
    @test res_f.θ_weights == res_g.θ_weights
    @test res_f.log_marginal == res_g.log_marginal
    @test res_f.log_π == res_g.log_π

    # Joint posterior exactly Gaussian → refit moments reproduce the
    # pilot summaries up to grid-truncation/trapezoid error (±5σ grid
    # clips ~6e-7 of the mass; moments move by O(1e-5 σ), well under
    # the 2%/3% bounds asserted here).
    σ = sqrt.(max.(res_g.x_var, 0.0))
    @test maximum(abs.(res_f.x_mean .- res_g.x_mean) ./ σ) < 0.02
    @test maximum(abs.(sqrt.(res_f.x_var) .- σ) ./ σ) < 0.03
end

@testset "Poisson low counts — FL shifts and matches the accessor" begin
    rng = Random.Xoshiro(20260707)
    n = 40
    k = 8
    grp = repeat(1:k, inner=n ÷ k)
    u_true = 0.9 .* randn(rng, k)
    y = [rand(rng, Poisson(exp(-0.5 + u_true[g]))) for g in grp]

    Acols = [Float64(j == g) for g in grp, j in 1:k]
    A = sparse([ones(n) Acols])
    model = LatentGaussianModel(PoissonLikelihood(), (Intercept(), IID(k)), A)

    res_g = inla(model, y; int_strategy=:grid)
    fl = FullLaplace()
    res_f = inla(model, y; int_strategy=:grid, latent_strategy=fl)

    # Low-count Poisson latents are measurably non-Gaussian: the FL
    # summary must move away from the Newton-mode accumulation.
    @test maximum(abs, res_f.x_mean .- res_g.x_mean) > 1.0e-3
    @test all(res_f.x_var .> 0)

    # Self-consistency: the replaced summary is exactly the trapezoid
    # moment of the FullLaplace density `posterior_marginal_x` exposes,
    # evaluated on the same grid the post-pass used (anchored at the
    # pilot summaries of the Gaussian-strategy fit, which are identical
    # to the pilot pass inside `res_f`'s fit).
    i = argmax(abs.(res_f.x_mean .- res_g.x_mean))
    σ0 = sqrt(max(res_g.x_var[i], 0.0))
    σ0 = σ0 > 0 ? σ0 : 1.0
    xs = collect(range(res_g.x_mean[i] - fl.span * σ0,
        res_g.x_mean[i] + fl.span * σ0; length=fl.n_grid))
    mfl = posterior_marginal_x(res_f, i; strategy=FullLaplace(),
        model=model, y=y, grid=xs)
    Z = _fl_trapz(xs, mfl.pdf)
    μ = _fl_trapz(xs, xs .* mfl.pdf) / Z
    σ² = _fl_trapz(xs, (xs .- μ) .^ 2 .* mfl.pdf) / Z
    @test isapprox(μ, res_f.x_mean[i]; atol=1.0e-8)
    @test isapprox(σ², res_f.x_var[i]; atol=1.0e-8)
end

@testset "0-hyperparameter fast path" begin
    rng = Random.Xoshiro(7)
    n = 30
    y = [rand(rng, Poisson(exp(0.3))) for _ in 1:n]
    model = LatentGaussianModel(PoissonLikelihood(), (Intercept(),),
        sparse(ones(n, 1)))

    res_g = inla(model, y)
    res_f = inla(model, y; latent_strategy=FullLaplace())
    @test isfinite(res_f.x_mean[1])
    @test res_f.x_var[1] > 0
    # Same posterior, different summary construction: agreement within
    # a fraction of the posterior sd (n = 30 Poisson intercept is only
    # mildly skewed).
    @test abs(res_f.x_mean[1] - res_g.x_mean[1]) <
          0.5 * sqrt(res_g.x_var[1])
end
