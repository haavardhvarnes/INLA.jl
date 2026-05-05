# GaussianBasisPrior — independent-Gaussian prior on
# `SPDE2NonStationary` basis coefficients. Anchors:
# (1) the log density matches the closed-form ∑_k −½ λ_k (θ−μ)² + NC;
# (2) defaults reproduce the unit-Gaussian baseline (μ=0, λ=1);
# (3) length contract is enforced; (4) prior totals integrate to 1
# (sanity check via two-point density ratio).

using INLASPDE
using INLASPDE: log_basis_prior_density

@testset "GaussianBasisPrior — construction + length contract" begin
    prior = GaussianBasisPrior(mean = [0.0, 1.0, -1.0], prec = [1.0, 0.5, 2.0])
    @test length(prior) == 3
    @test prior.mean == [0.0, 1.0, -1.0]
    @test prior.prec == [1.0, 0.5, 2.0]

    # `p`-arg form with scalar defaults
    p = GaussianBasisPrior(4)
    @test length(p) == 4
    @test all(p.mean .== 0.0)
    @test all(p.prec .== 1.0)

    p2 = GaussianBasisPrior(3; mean = -0.5, prec = 1.0e-3)
    @test all(p2.mean .== -0.5)
    @test all(p2.prec .== 1.0e-3)
end

@testset "GaussianBasisPrior — invalid input rejected" begin
    @test_throws ArgumentError GaussianBasisPrior(
        mean = [0.0, 0.0], prec = [1.0]
    )
    @test_throws ArgumentError GaussianBasisPrior(
        mean = [0.0], prec = [-1.0]
    )
    @test_throws ArgumentError GaussianBasisPrior(
        mean = [0.0], prec = [0.0]
    )
    @test_throws ArgumentError GaussianBasisPrior(0)
end

@testset "GaussianBasisPrior — log density matches closed form" begin
    prior = GaussianBasisPrior(
        mean = [0.0, 1.0, -1.0],
        prec = [1.0, 0.5, 2.0]
    )
    for θ in ([0.0, 0.0, 0.0], [1.5, 0.5, -0.7], [-2.3, 1.2, 0.4])
        expected = sum(
            -0.5 * prior.prec[k] * (θ[k] - prior.mean[k])^2 +
            0.5 * (log(prior.prec[k]) - log(2π)) for k in eachindex(θ)
        )
        @test log_basis_prior_density(prior, θ) ≈ expected rtol = 1.0e-12
    end
end

@testset "GaussianBasisPrior — unit-Gaussian default ≈ -d/2 log(2π) at θ=0" begin
    p = GaussianBasisPrior(5)             # μ=0, λ=1
    θ = zeros(5)
    @test log_basis_prior_density(p, θ) ≈ -2.5 * log(2π) rtol = 1.0e-12
end

@testset "GaussianBasisPrior — length-mismatch raises" begin
    prior = GaussianBasisPrior(3)
    @test_throws ArgumentError log_basis_prior_density(prior, [0.0, 0.0])
end

@testset "GaussianBasisPrior — density factorises across coefficients" begin
    # Pull two independent coordinates: log p(θ) = log p(θ_1) + log p(θ_-1)
    # where each marginal is N(μ_k, 1/λ_k). Equivalent test: density
    # difference between two θ-vectors equals the sum of per-coord
    # quadratic differences (no cross terms).
    prior = GaussianBasisPrior(mean = [0.0, 0.5, -0.3], prec = [1.5, 0.8, 2.2])
    θa = [0.1, 0.2, 0.3]
    θb = [-0.4, 0.6, 0.0]
    direct = log_basis_prior_density(prior, θa) - log_basis_prior_density(prior, θb)
    quad = sum(
        -0.5 * prior.prec[k] * ((θa[k] - prior.mean[k])^2 -
                                (θb[k] - prior.mean[k])^2)
        for k in eachindex(θa)
    )
    @test direct ≈ quad rtol = 1.0e-12
end
