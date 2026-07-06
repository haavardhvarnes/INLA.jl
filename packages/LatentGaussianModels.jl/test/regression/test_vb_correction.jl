# Low-rank variational-Bayes mean correction (ADR-048).
#
# Tier-1 coverage for `INLA(vb_correction = :mean)`:
#
#  - resolution shim + struct validation + `:none` bit-identity,
#  - exact λ = 0 collapse on Gaussian likelihoods (the Laplace mean is
#    the exact posterior mean; the Newton-mode stationarity zeroes the
#    initial variational gradient, and GH quadrature is exact for the
#    linear ℓ′),
#  - brute-force posterior-mean agreement on a tiny low-count Poisson
#    model (dense 1-D quadrature of the exact conditional posterior),
#  - internal triangulation on a skewed Poisson-IID model: the VB
#    corrected mean must land closer to the FullLaplace
#    integration-stage mean (item 16.2, reference quality) than the
#    uncorrected Gaussian mean does,
#  - hard-constraint preservation (sum-to-zero) of the corrected mean.
#
# Parity against modern (VB-corrected) R-INLA is tier-2 territory —
# see the vb block in `test/oracle/test_synthetic_brunei.jl`.

using LatentGaussianModels: GaussianLikelihood, PoissonLikelihood, Intercept,
                            IID, Generic0, GammaPrecision, LatentGaussianModel,
                            inla, INLA, VBMeanCorrection, FullLaplace
using GMRFs: LinearConstraint
using Distributions: Poisson
using SparseArrays
using LinearAlgebra: I, SymTridiagonal
using Random

@testset "resolution + validation" begin
    @test VBMeanCorrection().n_gh == 15
    @test VBMeanCorrection().max_block == 30
    @test VBMeanCorrection().indices === nothing
    @test VBMeanCorrection(n_gh=9, max_block=5, indices=[1, 3]).indices == [1, 3]
    @test_throws ArgumentError VBMeanCorrection(n_gh=2)
    @test_throws ArgumentError VBMeanCorrection(max_block=-1)
    @test_throws ArgumentError VBMeanCorrection(indices=Int[])
    @test_throws ArgumentError VBMeanCorrection(indices=[1, 1])
    @test_throws ArgumentError INLA(vb_correction=:bogus)
    @test INLA(vb_correction=:none).vb === nothing
    @test INLA(vb_correction=:mean).vb == VBMeanCorrection()
end

@testset "Gaussian likelihood — correction is exactly zero" begin
    rng = Random.Xoshiro(20260706)
    n = 40
    k = 8
    grp = repeat(1:k, inner=n ÷ k)
    u_true = 0.7 .* randn(rng, k)
    y = [0.3 + u_true[g] + 0.5 * randn(rng) for g in grp]

    Acols = [Float64(j == g) for g in grp, j in 1:k]
    A = sparse([ones(n) Acols])
    model = LatentGaussianModel(GaussianLikelihood(), (Intercept(), IID(k)), A)

    res_none = inla(model, y; int_strategy=:grid)
    res_vb = inla(model, y; int_strategy=:grid, vb_correction=:mean)

    # θ-stage identical (the correction is a summary-stage shift).
    @test res_vb.θ̂ == res_none.θ̂
    @test res_vb.log_marginal == res_none.log_marginal
    # λ = 0 exactly (measured at machine zero; the tolerance leaves
    # headroom for platform variation in the Newton residual).
    @test maximum(abs, res_vb.x_mean .- res_none.x_mean) < 1.0e-7
    # x_var accumulates mode_k² and inherits last-bit noise from the
    # (machine-zero) shift — approximate equality, not bit identity.
    @test isapprox(res_vb.x_var, res_none.x_var; rtol=1.0e-10, atol=1.0e-14)
end

@testset "tiny Poisson — brute-force posterior mean" begin
    # One latent coordinate (intercept only), a handful of low counts:
    # p(α | y) ∝ exp(Σᵢ (yᵢ α − exp(α))) · N(α; 0, τ₀⁻¹) is 1-D, so the
    # exact posterior mean at fixed θ comes from dense quadrature. The
    # model has no hyperparameters (Poisson likelihood, improper-free
    # intercept prior handled by the fit at θ = []), so INLAResult is a
    # single Laplace fit — the corrected x_mean must beat the
    # uncorrected one against the dense truth.
    y = [0, 1, 0, 0, 2, 0, 1, 0]
    n = length(y)
    model = LatentGaussianModel(PoissonLikelihood(), (Intercept(),),
        sparse(ones(n, 1)))

    res_none = inla(model, y)
    res_vb = inla(model, y; vb_correction=:mean)

    # Dense truth: the intercept prior is improper flat (R-INLA
    # convention, prec.intercept = 0), so the posterior is just the
    # Poisson likelihood in α.
    αs = range(-6.0, 3.0; length=20001)
    logp = [sum(y) * α - n * exp(α) for α in αs]
    w = exp.(logp .- maximum(logp))
    μ_true = sum(w .* αs) / sum(w)

    err_none = abs(res_none.x_mean[1] - μ_true)
    err_vb = abs(res_vb.x_mean[1] - μ_true)
    @test err_vb < err_none
    # Measured: the correction closes ~96% of the 0.13 Laplace-mean gap
    # (residual 5.4e-3 — the Gaussian-family variational optimum is the
    # E_q[∇ log p] = 0 point, not the exact mean, and this posterior has
    # sd ≈ 0.5). Assert both the absolute residual and the relative
    # gap reduction with slack.
    @test err_vb < 1.0e-2
    @test err_vb < 0.2 * err_none
end

@testset "skewed Poisson-IID — triangulation against FullLaplace" begin
    rng = Random.Xoshiro(99)
    n = 48
    k = 8
    grp = repeat(1:k, inner=n ÷ k)
    u_true = 0.9 .* randn(rng, k)
    y = [rand(rng, Poisson(exp(-0.4 + u_true[g]))) for g in grp]

    Acols = [Float64(j == g) for g in grp, j in 1:k]
    A = sparse([ones(n) Acols])
    model = LatentGaussianModel(PoissonLikelihood(), (Intercept(), IID(k)), A)

    res_g = inla(model, y; int_strategy=:grid)
    res_vb = inla(model, y; int_strategy=:grid, vb_correction=:mean)
    res_fl = inla(model, y; int_strategy=:grid,
        latent_strategy=FullLaplace())

    # The correction must do real work on a skewed posterior…
    @test maximum(abs, res_vb.x_mean .- res_g.x_mean) > 1.0e-4
    # …and move the summary toward the reference-quality FullLaplace
    # mean, coordinate-wise in aggregate.
    gap_g = sum(abs, res_g.x_mean .- res_fl.x_mean)
    gap_vb = sum(abs, res_vb.x_mean .- res_fl.x_mean)
    @test gap_vb < gap_g
end

@testset "sum-to-zero constraint preserved" begin
    rng = Random.Xoshiro(20260504)
    n = 12
    R = sparse(SymTridiagonal(fill(2.0, n), fill(-1.0, n - 1)))
    R[1, 1] = 1.0
    R[n, n] = 1.0
    constraint = LinearConstraint(ones(1, n), zeros(1))
    g0 = Generic0(R; rankdef=1, constraint=constraint,
        hyperprior=GammaPrecision(10.0, 10.0))
    u_true = randn(rng, n)
    u_true .-= sum(u_true) / n
    y = [rand(rng, Poisson(exp(0.2 + 0.8 * u_true[i]))) for i in 1:n]
    A = sparse([ones(n) Matrix{Float64}(I, n,n)])
    model = LatentGaussianModel(PoissonLikelihood(), (Intercept(), g0), A)

    res_vb = inla(model, y; int_strategy=:grid, vb_correction=:mean)
    # The Generic0 block occupies coordinates 2:(n+1); the kriging
    # projection inside the shift must keep the sum-to-zero constraint
    # on the corrected mean.
    @test abs(sum(res_vb.x_mean[2:(n + 1)])) < 1.0e-8
end
