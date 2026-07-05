# Integrated hyperparameter marginals (ADR-046).
#
# Tier-1 coverage for `posterior_marginal_θ(method = :integrated)`:
#
#  - dim(θ) = 1 design-point reuse: normalisation, self-consistency of
#    the density mean against `res.θ_mean` (same posterior, different
#    quadrature), near-Gaussian collapse onto `:gaussian`, and
#    composition with `refine_hyperposterior`.
#  - dim(θ) = 2 profile slice: normalisation, near-Gaussian moment
#    agreement with `:gaussian`, and the `model`/`y` requirement.
#  - `:auto` resolution rules and error paths.
#
# Accuracy against R-INLA's `marginals.hyperpar` is tier-2 territory —
# see the integrated-marginal blocks in `test/oracle/`.

using LatentGaussianModels: GaussianLikelihood, PoissonLikelihood, Intercept,
                            IID, LatentGaussianModel, inla,
                            posterior_marginal_θ, refine_hyperposterior
using SparseArrays
using Random
using Distributions: Poisson

# Trapezoid moments of a gridded density, renormalised on the grid so
# that narrow grids still yield the conditional mean/sd of the visible
# mass (mirrors the oracle-test idiom).
function _theta_pdf_moments(xs::AbstractVector, pdf::AbstractVector)
    Z = 0.0
    μ = 0.0
    @inbounds for i in 1:(length(xs) - 1)
        h = xs[i + 1] - xs[i]
        Z += 0.5 * h * (pdf[i] + pdf[i + 1])
        μ += 0.5 * h * (xs[i] * pdf[i] + xs[i + 1] * pdf[i + 1])
    end
    μ /= Z
    σ² = 0.0
    @inbounds for i in 1:(length(xs) - 1)
        h = xs[i + 1] - xs[i]
        σ² += 0.5 * h * ((xs[i] - μ)^2 * pdf[i] + (xs[i + 1] - μ)^2 * pdf[i + 1])
    end
    σ² /= Z
    return (mean=μ, sd=sqrt(max(σ², 0.0)), mass=Z)
end

@testset "dim(θ) = 1 — design-point reuse" begin
    rng = Random.Xoshiro(20260705)
    n = 60
    k = 12
    grp = repeat(1:k, inner=n ÷ k)
    u_true = 0.8 .* randn(rng, k)
    y = [rand(rng, Poisson(exp(0.5 + u_true[g]))) for g in grp]

    Acols = [Float64(j == g) for g in grp, j in 1:k]
    A = sparse([ones(n) Acols])
    model = LatentGaussianModel(PoissonLikelihood(), (Intercept(), IID(k)), A)

    res = inla(model, y; int_strategy=:grid)
    @test length(res.log_π) == length(res.θ_points)
    @test all(isfinite, res.log_π)

    m_auto = posterior_marginal_θ(res, 1)
    m_int = posterior_marginal_θ(res, 1; method=:integrated)
    @test m_auto.pdf == m_int.pdf          # :auto resolves to :integrated

    @test all(≥(0), m_int.pdf)
    mom = _theta_pdf_moments(m_int.θ, m_int.pdf)
    # Renormalised on the returned grid by construction.
    @test isapprox(mom.mass, 1.0; atol=1.0e-8)

    # Self-consistency: the density mean and `θ_mean` estimate the same
    # posterior mean by different quadratures (interpolated trapezoid vs
    # IS-weighted design points). On a skewed log-precision posterior
    # both sit off θ̂; they must agree with each other.
    σ_g = sqrt(res.Σθ[1, 1])
    @test abs(mom.mean - res.θ_mean[1]) < 0.1 * σ_g

    # Custom grid is used verbatim.
    gr = collect(range(res.θ̂[1] - 3σ_g, res.θ̂[1] + 3σ_g; length=41))
    m_gr = posterior_marginal_θ(res, 1; grid=gr)
    @test m_gr.θ == gr

    # :gaussian stays bit-for-bit the closed form.
    m_g = posterior_marginal_θ(res, 1; method=:gaussian)
    @test m_g.pdf[38] ≈
          exp(-0.5 * ((m_g.θ[38] - res.θ̂[1]) / σ_g)^2) / (σ_g * sqrt(2π))
end

@testset "dim(θ) = 1 — near-Gaussian collapse + refine composition" begin
    rng = Random.Xoshiro(42)
    n = 150
    y = 1.2 .+ 0.6 .* randn(rng, n)
    model = LatentGaussianModel(GaussianLikelihood(), (Intercept(),),
        sparse(ones(n, 1)))

    res = inla(model, y; int_strategy=:grid)
    m_int = posterior_marginal_θ(res, 1; method=:integrated)
    m_g = posterior_marginal_θ(res, 1; method=:gaussian)

    # log τ posterior at n = 150 is close to Gaussian: the integrated
    # density must collapse onto the Gaussian one. 10% of the peak
    # covers the residual skewness plus FD-Hessian σ error.
    peak = maximum(m_g.pdf)
    @test maximum(abs, m_int.pdf .- m_g.pdf) < 0.1 * peak

    # A denser design (refine_hyperposterior) must reproduce the same
    # density, not shift it: sup-norm agreement at the few-% level.
    res_ref = refine_hyperposterior(res, model, y; n_grid=15,
        skewness_correction=false)
    m_ref = posterior_marginal_θ(res_ref, 1; method=:integrated,
        grid=m_int.θ)
    @test maximum(abs, m_ref.pdf .- m_int.pdf) < 0.05 * peak
end

@testset "dim(θ) = 2 — profile slice" begin
    rng = Random.Xoshiro(99)
    n = 120
    k = 8
    grp = repeat(1:k, inner=n ÷ k)
    u_true = 0.7 .* randn(rng, k)
    y = [0.3 + u_true[g] + 0.5 * randn(rng) for g in grp]

    Acols = [Float64(j == g) for g in grp, j in 1:k]
    A = sparse([ones(n) Acols])
    model = LatentGaussianModel(GaussianLikelihood(), (Intercept(), IID(k)), A)

    res = inla(model, y; int_strategy=:grid)
    @test length(res.θ̂) == 2

    # :auto without model/y falls back to :gaussian exactly.
    m_auto = posterior_marginal_θ(res, 1)
    m_g = posterior_marginal_θ(res, 1; method=:gaussian)
    @test m_auto.pdf == m_g.pdf

    # Explicit :integrated without model/y is a user error.
    @test_throws ArgumentError posterior_marginal_θ(res, 1;
        method=:integrated)
    @test_throws ArgumentError posterior_marginal_θ(res, 1;
        method=:integrated, model=model)   # y still missing

    for j in 1:2
        m_int = posterior_marginal_θ(res, j; model=model, y=y)
        @test all(≥(0), m_int.pdf)
        mom = _theta_pdf_moments(m_int.θ, m_int.pdf)
        @test isapprox(mom.mass, 1.0; atol=1.0e-8)
        # Near-Gaussian posterior (Gaussian likelihood, n = 120): the
        # profile-slice moments must land near the Gaussian ones.
        mom_g = _theta_pdf_moments(posterior_marginal_θ(res, j;
                method=:gaussian).θ,
            posterior_marginal_θ(res, j;
                method=:gaussian).pdf)
        σ_g = sqrt(res.Σθ[j, j])
        @test abs(mom.mean - mom_g.mean) < 0.5 * σ_g
        @test 0.7 < mom.sd / mom_g.sd < 1.4
    end
end

@testset "error paths" begin
    rng = Random.Xoshiro(7)
    n = 40
    y = 0.3 .+ 0.7 .* randn(rng, n)
    model = LatentGaussianModel(GaussianLikelihood(), (Intercept(),),
        sparse(ones(n, 1)))
    res = inla(model, y)

    @test_throws ArgumentError posterior_marginal_θ(res, 0)
    @test_throws ArgumentError posterior_marginal_θ(res, 2)
    @test_throws ArgumentError posterior_marginal_θ(res, 1; method=:bogus)
end
