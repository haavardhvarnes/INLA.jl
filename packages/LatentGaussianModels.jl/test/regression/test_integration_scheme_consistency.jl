using LatentGaussianModels: GaussianLikelihood, Intercept, IID,
                            LatentGaussianModel, inla, PCPrecision
using SparseArrays
using LinearAlgebra: I
using Random

# On a well-behaved (near-Gaussian hyperposterior) model, the three
# integration schemes — tensor Grid, CCD, and Gauss–Hermite — must agree
# on the posterior hyperparameter mean, the latent posterior mean, and the
# log marginal likelihood to within their quadrature differences. A bug in
# design-point placement or quadrature weights (a real risk in CCD assembly
# and GH weight scaling) would break the agreement even though each scheme
# still "runs". This is the cross-check that single-scheme tests can't give.

@testset "Integration schemes agree — Gaussian + Intercept + IID (dim θ = 2)" begin
    rng = Random.Xoshiro(20260705)
    n = 200
    α_true = 1.0
    σ = 0.4
    τ_u = 6.0
    u = randn(rng, n) ./ sqrt(τ_u)
    y = α_true .+ u .+ σ .* randn(rng, n)

    A = sparse([ones(n) Matrix{Float64}(I, n, n)])
    model = LatentGaussianModel(GaussianLikelihood(),
        (Intercept(), IID(n; hyperprior=PCPrecision(1.0, 0.01))), A)

    res_grid = inla(model, y; int_strategy=:grid)
    res_ccd = inla(model, y; int_strategy=:ccd)
    res_gh = inla(model, y; int_strategy=:gauss_hermite)

    # Posterior hyperparameter mean (internal log scale, dim 2): schemes
    # differ by quadrature order, so compare with a modest absolute band.
    @test isapprox(res_grid.θ_mean, res_ccd.θ_mean; atol=0.05)
    @test isapprox(res_grid.θ_mean, res_gh.θ_mean; atol=0.05)

    # Latent posterior mean: the dominant summary; should agree tightly.
    @test isapprox(res_grid.x_mean, res_ccd.x_mean; rtol=1.0e-3, atol=1.0e-3)
    @test isapprox(res_grid.x_mean, res_gh.x_mean; rtol=1.0e-3, atol=1.0e-3)

    # Log marginal likelihood: schemes differ by their quadrature order
    # (CCD is a low-order design), so allow up to ~half a nat — still tight
    # enough to catch a gross weighting/normalisation bug worth several nats.
    @test isapprox(res_grid.log_marginal, res_ccd.log_marginal; atol=0.5)
    @test isapprox(res_grid.log_marginal, res_gh.log_marginal; atol=0.5)
end
