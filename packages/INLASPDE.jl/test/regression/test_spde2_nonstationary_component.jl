# SPDE2NonStationary — the AbstractLatentComponent contract,
# stationary recovery against SPDE2 (constant B_τ, B_κ → stationary
# precision), and the per-vertex non-stationary precision build via
# spde_precision_nonstationary. Anchors:
#   1. AbstractLatentComponent surface
#   2. α ∈ {2} accepted, α ∉ {2} rejected
#   3. Stationary recovery: B_τ = B_κ = ones(n_v, 1) reproduces SPDE2
#      precision exactly
#   4. spde_precision_nonstationary stationary limit equals SPDE2
#   5. Length contract on θ
#   6. log_hyperprior matches GaussianBasisPrior log density
#   7. log_normalizing_constant equals -½ d log(2π) + ½ log|Q|
#   8. Two-region piecewise-constant κ → SPD precision

using GMRFs
using LatentGaussianModels
using LatentGaussianModels: initial_hyperparameters, nhyperparameters,
                            precision_matrix, log_hyperprior,
                            log_normalizing_constant
using LinearAlgebra
using SparseArrays
using INLASPDE: log_basis_prior_density

@testset "SPDE2NonStationary — construction + AbstractLatentComponent contract" begin
    n_v = 4
    B_τ = ones(n_v, 1)
    B_κ = ones(n_v, 1)
    spde = SPDE2NonStationary(POINTS_SQ, TRIS_SQ; B_τ=B_τ, B_κ=B_κ)

    @test spde isa AbstractLatentComponent
    @test length(spde) == 4
    @test nhyperparameters(spde) == 2     # 1 + 1 = 2
    @test initial_hyperparameters(spde) == [0.0, 0.0]   # default prior μ = 0

    @test GMRFs.constraints(spde) isa GMRFs.NoConstraint
end

@testset "SPDE2NonStationary — α=2 accepted; α∉{2} rejected" begin
    n_v = 4
    B = ones(n_v, 1)
    spde = SPDE2NonStationary(POINTS_SQ, TRIS_SQ; α=2, B_τ=B, B_κ=B)
    @test spde isa SPDE2NonStationary{2}

    @test_throws ArgumentError SPDE2NonStationary(POINTS_SQ, TRIS_SQ;
        α=1, B_τ=B, B_κ=B)
    @test_throws ArgumentError SPDE2NonStationary(POINTS_SQ, TRIS_SQ;
        α=3, B_τ=B, B_κ=B)
end

@testset "SPDE2NonStationary — stationary recovery matches SPDE2" begin
    # B_τ = B_κ = ones(n_v, 1) and θ = [log τ, log κ] should reproduce
    # the SPDE2 stationary precision τ²·(κ⁴ C̃ + 2κ² G₁ + G₂).
    n_v = 4
    B_τ = ones(n_v, 1)
    B_κ = ones(n_v, 1)
    spde_ns = SPDE2NonStationary(POINTS_SQ, TRIS_SQ; B_τ=B_τ, B_κ=B_κ)
    spde_st = SPDE2(POINTS_SQ, TRIS_SQ)

    for (log_τ, log_κ) in ((0.0, 0.0), (-0.7, 0.4), (1.3, -0.6))
        Q_ns = precision_matrix(spde_ns, [log_τ, log_κ])
        Q_st = precision_matrix(spde_st, [log_τ, log_κ])
        @test Matrix(Q_ns)≈Matrix(Q_st) rtol=1.0e-12
    end
end

@testset "SPDE2NonStationary — multi-coefficient stationary recovery" begin
    # Replicate the intercept p_τ + 1 = 2 times; the same θ_τ shared
    # across the two columns gives 2·θ_τ as the effective intercept.
    # Equivalent constant intercept = 2·θ_τ_each.
    n_v = 4
    B_τ = hcat(ones(n_v), ones(n_v))   # n_v × 2
    B_κ = ones(n_v, 1)
    spde_ns = SPDE2NonStationary(POINTS_SQ, TRIS_SQ; B_τ=B_τ, B_κ=B_κ)
    spde_st = SPDE2(POINTS_SQ, TRIS_SQ)

    # log τ_v = θ_τ_1 + θ_τ_2 = 2·a   ⇒ effective log τ = 2a
    a = 0.3
    log_κ = 0.5
    Q_ns = precision_matrix(spde_ns, [a, a, log_κ])
    Q_st = precision_matrix(spde_st, [2 * a, log_κ])
    @test Matrix(Q_ns)≈Matrix(Q_st) rtol=1.0e-12
end

@testset "SPDE2NonStationary — spde_precision_nonstationary stationary limit" begin
    # Pin a constant τ_v, κ_v and check the build matches SPDE2's
    # stationary precision via spde_precision(fem, 2, τ, κ).
    fem = INLASPDE.FEMMatrices(POINTS_SQ, TRIS_SQ)
    n_v = 4
    τ = 1.7
    κ = 0.6
    Q_const = spde_precision_nonstationary(fem, 2,
        fill(τ, n_v), fill(κ, n_v))
    Q_stat = spde_precision(fem, 2, τ, κ)
    @test Matrix(Q_const)≈Matrix(Q_stat) rtol=1.0e-12
end

@testset "SPDE2NonStationary — length-mismatch raises" begin
    n_v = 4
    spde = SPDE2NonStationary(POINTS_SQ, TRIS_SQ;
        B_τ=ones(n_v, 1), B_κ=ones(n_v, 1))
    # nhyperparameters = 2; passing 3 should throw.
    @test_throws ArgumentError precision_matrix(spde, [0.0, 0.0, 0.0])
    @test_throws ArgumentError precision_matrix(spde, [0.0])
end

@testset "SPDE2NonStationary — log_hyperprior matches GaussianBasisPrior" begin
    n_v = 4
    p_τ, p_κ = 1, 2
    prior = GaussianBasisPrior(
        mean=[0.5, -0.3, 0.7],
        prec=[1.0, 2.0, 0.5]
    )
    spde = SPDE2NonStationary(POINTS_SQ, TRIS_SQ;
        B_τ=ones(n_v, p_τ), B_κ=ones(n_v, p_κ), prior=prior)
    for θ in ([0.0, 0.0, 0.0], [0.5, -0.3, 0.7], [1.2, -0.4, 0.1])
        @test log_hyperprior(spde, θ)≈log_basis_prior_density(prior, θ) rtol=1.0e-12
    end

    # Default prior (μ=0, λ=1) for comparison
    spde_default = SPDE2NonStationary(POINTS_SQ, TRIS_SQ;
        B_τ=ones(n_v, 1), B_κ=ones(n_v, 1))
    @test log_hyperprior(spde_default, [0.0, 0.0])≈-log(2π) rtol=1.0e-12
end

@testset "SPDE2NonStationary — initial_hyperparameters = prior mean" begin
    n_v = 4
    prior = GaussianBasisPrior(mean=[0.5, -0.3, 0.7], prec=[1.0, 1.0, 1.0])
    spde = SPDE2NonStationary(POINTS_SQ, TRIS_SQ;
        B_τ=ones(n_v, 1), B_κ=ones(n_v, 2), prior=prior)
    @test initial_hyperparameters(spde) == [0.5, -0.3, 0.7]
    # Ensure it's a copy, not a view, to prevent caller mutation.
    init = initial_hyperparameters(spde)
    init[1] = 999.0
    @test prior.mean[1] == 0.5
end

@testset "SPDE2NonStationary — log_normalizing_constant finite + ½ log|Q|" begin
    n_v = 4
    spde = SPDE2NonStationary(POINTS_SQ, TRIS_SQ;
        B_τ=ones(n_v, 1), B_κ=ones(n_v, 1))
    θ = [0.4, -0.2]
    Q = precision_matrix(spde, θ)
    expected = -0.5 * length(spde) * log(2π) + 0.5 * logdet(cholesky(Symmetric(Q)))
    @test log_normalizing_constant(spde, θ)≈expected rtol=1.0e-12
end

@testset "SPDE2NonStationary — invalid input rejected" begin
    n_v = 4
    # B_τ, B_κ size mismatch
    @test_throws ArgumentError SPDE2NonStationary(POINTS_SQ, TRIS_SQ;
        B_τ=ones(3, 1), B_κ=ones(n_v, 1))
    @test_throws ArgumentError SPDE2NonStationary(POINTS_SQ, TRIS_SQ;
        B_τ=ones(n_v, 1), B_κ=ones(3, 1))
    # prior length mismatch
    bad_prior = GaussianBasisPrior(5)
    @test_throws ArgumentError SPDE2NonStationary(POINTS_SQ, TRIS_SQ;
        B_τ=ones(n_v, 1), B_κ=ones(n_v, 1), prior=bad_prior)
    # zero-column basis
    @test_throws ArgumentError SPDE2NonStationary(POINTS_SQ, TRIS_SQ;
        B_τ=zeros(n_v, 0), B_κ=ones(n_v, 1))
    # negative or zero precision in prior
    @test_throws ArgumentError GaussianBasisPrior(mean=[0.0], prec=[0.0])
end

@testset "SPDE2NonStationary — spde_precision_nonstationary input validation" begin
    fem = INLASPDE.FEMMatrices(POINTS_SQ, TRIS_SQ)
    # Length mismatches and the structural α-validity check stay as
    # ArgumentError (programming-bug tier); parametric per-vertex (τ_v,
    # κ_v) positivity checks raise DomainError so the LGM safety net
    # catches them under LBFGS overshoot (ADR-031).
    @test_throws ArgumentError spde_precision_nonstationary(fem, 2, ones(3), ones(4))
    @test_throws ArgumentError spde_precision_nonstationary(fem, 2, ones(4), ones(3))
    @test_throws DomainError spde_precision_nonstationary(fem, 2, fill(-1.0, 4), ones(4))
    @test_throws DomainError spde_precision_nonstationary(fem, 2, ones(4), fill(0.0, 4))
    @test_throws ArgumentError spde_precision_nonstationary(fem, 1, ones(4), ones(4))
end

@testset "SPDE2NonStationary — non-stationary κ produces SPD precision" begin
    # Use the larger fmesher unit-square mesh from M3 to get a real
    # non-stationary precision rather than the trivial 4-vertex case.
    pts = [Float64(x) for x in 0.0:0.2:1.0]
    grid_pts = collect(Iterators.product(pts, pts))
    coords = reduce(vcat, [collect(p)' for p in grid_pts])
    mesh = inla_mesh_2d(coords; max_edge=0.25)
    n_v = num_vertices(mesh)
    points = mesh.points
    triangles = mesh.triangles

    # Two-region piecewise-constant κ_v: x < 0.5 → region 1, else region 2.
    region1 = points[:, 1] .< 0.5
    region2 = .!region1
    B_κ = hcat(Float64.(region1), Float64.(region2))
    B_τ = ones(n_v, 1)
    spde = SPDE2NonStationary(points, triangles;
        B_τ=B_τ, B_κ=B_κ)

    # θ_τ = 0 (so τ_v = 1 everywhere); θ_κ = [log 1, log 4]
    # Region 1 has κ = 1, region 2 has κ = 4 — large contrast.
    θ = [0.0, 0.0, log(4.0)]
    Q = precision_matrix(spde, θ)
    @test size(Q) == (n_v, n_v)
    @test maximum(abs, Q - Q') < 1.0e-9
    @test isposdef(Symmetric(Matrix(Q)))
end
