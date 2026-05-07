# SPDE1D component — the AbstractLatentComponent contract, the
# (log τ, log κ) ↔ (ρ, σ) mapping for d = 1, and interop with the 1D
# FEM assembly. Anchors: graph reflects 1D adjacency, ρ and σ user
# scales match the closed-form formulas, both α ∈ {1, 2} are valid.

using GMRFs
using LatentGaussianModels
using LatentGaussianModels: initial_hyperparameters, nhyperparameters,
                            precision_matrix, log_hyperprior,
                            log_normalizing_constant
using LinearAlgebra
using SparseArrays

@testset "SPDE1D — construction + AbstractLatentComponent contract" begin
    mesh = inla_mesh_1d([0.0, 1.0, 2.0]; max_edge=0.5)
    pc = PCMatern{1}(range_U=0.5, range_α=0.05,
        sigma_U=1.0, sigma_α=0.01)
    spde = SPDE1D(mesh; α=2, pc=pc)

    @test spde isa AbstractLatentComponent
    @test length(spde) == num_vertices(mesh)
    @test nhyperparameters(spde) == 2
    @test initial_hyperparameters(spde) == [0.0, 0.0]

    @test GMRFs.constraints(spde) isa GMRFs.NoConstraint
end

@testset "SPDE1D — both α ∈ {1, 2} accepted" begin
    mesh = inla_mesh_1d([0.0, 1.0]; max_edge=0.25)
    spde1 = SPDE1D(mesh; α=1)
    spde2 = SPDE1D(mesh; α=2)
    @test spde1 isa SPDE1D{1}
    @test spde2 isa SPDE1D{2}

    @test_throws ArgumentError SPDE1D(mesh; α=0)
    @test_throws ArgumentError SPDE1D(mesh; α=3)
end

@testset "SPDE1D — graph matches 1D vertex adjacency" begin
    # 5-vertex chain: edges (1,2), (2,3), (3,4), (4,5) → 4 edges.
    mesh = inla_mesh_1d([0.0, 1.0, 2.0, 3.0, 4.0]; max_edge=1.0)
    spde = SPDE1D(mesh)
    g = spde.graph
    @test g isa GMRFs.AbstractGMRFGraph
    @test GMRFs.num_nodes(g) == 5

    A = GMRFs.adjacency_matrix(g)
    @test A == A'
    @test sum(A) == 2 * 4
    for k in 1:4
        @test A[k, k + 1] == true
    end
    @test A[1, 3] == false   # not adjacent on a 1D chain
end

@testset "SPDE1D — precision_matrix agrees with spde_precision(fem, …)" begin
    mesh = inla_mesh_1d(collect(0.0:0.5:3.0); max_edge=0.5)
    spde1 = SPDE1D(mesh; α=1)
    spde2 = SPDE1D(mesh; α=2)
    for (log_τ, log_κ) in ((0.0, 0.0), (-1.0, 0.7), (2.0, -0.5))
        Q1c = precision_matrix(spde1, [log_τ, log_κ])
        Q1f = spde_precision(spde1.fem, 1, exp(log_τ), exp(log_κ))
        @test Matrix(Q1c)≈Matrix(Q1f) rtol=1.0e-12

        Q2c = precision_matrix(spde2, [log_τ, log_κ])
        Q2f = spde_precision(spde2.fem, 2, exp(log_τ), exp(log_κ))
        @test Matrix(Q2c)≈Matrix(Q2f) rtol=1.0e-12
    end
end

@testset "SPDE1D — (log τ, log κ) ↔ (ρ, σ) mapping for α = 1 (ν = 0.5)" begin
    mesh = inla_mesh_1d([0.0, 1.0]; max_edge=1.0)
    spde = SPDE1D(mesh; α=1)

    # At (τ, κ) = (1, 1): ρ = 2, σ² = 1/2 ⇒ σ = 1/√2.
    ρ, σ = spde_user_scale(spde, [0.0, 0.0])
    @test ρ≈2.0 rtol=1.0e-12
    @test σ≈inv(sqrt(2.0)) rtol=1.0e-12

    # σ² = 1 / (2 κ τ²) — direct check at several points.
    for (log_τ, log_κ) in ((0.5, -0.3), (-1.0, 1.2), (2.0, 0.0))
        τ, κ = exp(log_τ), exp(log_κ)
        ρ, σ = spde_user_scale(spde, [log_τ, log_κ])
        @test ρ≈2.0 / κ rtol=1.0e-12
        @test σ^2≈1.0 / (2.0 * κ * τ^2) rtol=1.0e-12
    end

    # Round-trip user → internal → user.
    for (ρ_target, σ_target) in ((0.5, 0.3), (2.7, 1.9), (1.0, 1.0))
        log_τ, log_κ = spde_internal_scale(spde, ρ_target, σ_target)
        ρ_back, σ_back = spde_user_scale(spde, [log_τ, log_κ])
        @test ρ_back≈ρ_target rtol=1.0e-12
        @test σ_back≈σ_target rtol=1.0e-12
    end
end

@testset "SPDE1D — (log τ, log κ) ↔ (ρ, σ) mapping for α = 2 (ν = 1.5)" begin
    mesh = inla_mesh_1d([0.0, 1.0]; max_edge=1.0)
    spde = SPDE1D(mesh; α=2)

    # At (τ, κ) = (1, 1): ρ = √12, σ² = 1/4 ⇒ σ = 0.5.
    ρ, σ = spde_user_scale(spde, [0.0, 0.0])
    @test ρ≈sqrt(12.0) rtol=1.0e-12
    @test σ≈0.5 rtol=1.0e-12

    # σ² = 1 / (4 κ³ τ²)  —  direct check.
    for (log_τ, log_κ) in ((0.5, -0.3), (-1.0, 1.2), (2.0, 0.0))
        τ, κ = exp(log_τ), exp(log_κ)
        ρ, σ = spde_user_scale(spde, [log_τ, log_κ])
        @test ρ≈sqrt(12.0) / κ rtol=1.0e-12
        @test σ^2≈1.0 / (4.0 * κ^3 * τ^2) rtol=1.0e-12
    end

    # Round-trip.
    for (ρ_target, σ_target) in ((0.5, 0.3), (2.7, 1.9), (1.0, 1.0))
        log_τ, log_κ = spde_internal_scale(spde, ρ_target, σ_target)
        ρ_back, σ_back = spde_user_scale(spde, [log_τ, log_κ])
        @test ρ_back≈ρ_target rtol=1.0e-12
        @test σ_back≈σ_target rtol=1.0e-12
    end
end

@testset "SPDE1D — log_hyperprior matches PC-Matern{1} evaluator" begin
    mesh = inla_mesh_1d([0.0, 1.0]; max_edge=1.0)
    pc = PCMatern{1}(range_U=0.5, range_α=0.05,
        sigma_U=1.0, sigma_α=0.01)

    spde1 = SPDE1D(mesh; α=1, pc=pc)
    for θ in ([0.0, 0.0], [-1.3, 0.7], [2.1, -0.4])
        log_τ, log_κ = θ
        log_ρ = log(2.0) - log_κ
        log_σ = -0.5 * (log(2.0) + log_κ) - log_τ
        @test log_hyperprior(spde1, θ)≈pc_matern_log_density(pc, log_ρ, log_σ) rtol=1.0e-12
    end

    spde2 = SPDE1D(mesh; α=2, pc=pc)
    for θ in ([0.0, 0.0], [-1.3, 0.7], [2.1, -0.4])
        log_τ, log_κ = θ
        log_ρ = 0.5 * log(12.0) - log_κ
        log_σ = -log(2.0) - 1.5 * log_κ - log_τ
        @test log_hyperprior(spde2, θ)≈pc_matern_log_density(pc, log_ρ, log_σ) rtol=1.0e-12
    end
end

@testset "SPDE1D — invalid input rejected" begin
    mesh = inla_mesh_1d([0.0, 1.0]; max_edge=1.0)
    @test_throws ArgumentError spde_internal_scale(SPDE1D(mesh; α=1), -1.0, 1.0)
    @test_throws ArgumentError spde_internal_scale(SPDE1D(mesh; α=1), 1.0, -1.0)
    @test_throws ArgumentError spde_internal_scale(SPDE1D(mesh; α=2), 0.0, 1.0)
    @test_throws ArgumentError spde_internal_scale(SPDE1D(mesh; α=2), 1.0, 0.0)
end

@testset "SPDE1D — precision SPD on a fine mesh" begin
    mesh = inla_mesh_1d([0.0, 5.0]; max_edge=0.05)
    spde1 = SPDE1D(mesh; α=1)
    spde2 = SPDE1D(mesh; α=2)
    n = num_vertices(mesh)

    Q1 = precision_matrix(spde1, [0.0, 0.0])
    Q2 = precision_matrix(spde2, [0.0, 0.0])
    # α = 1: sum of symmetric matrices, exactly symmetric.
    @test issymmetric(Q1)
    # α = 2: G₂ = G₁·C̃⁻¹·G₁ via sparse matmul can have tiny FP
    # asymmetries on fine meshes; require approximate symmetry only.
    @test maximum(abs, Q2 - Q2') < 1.0e-10
    @test isposdef(Symmetric(Matrix(Q1)))
    @test isposdef(Symmetric(Matrix(Q2)))
    @test size(Q1) == (n, n)
    @test size(Q2) == (n, n)
    @test nnz(Q1) < n^2
    @test nnz(Q2) < n^2
end

@testset "SPDE1D — log_normalizing_constant finite + ½ log|Q|" begin
    mesh = inla_mesh_1d([0.0, 2.0]; max_edge=0.2)
    spde = SPDE1D(mesh; α=2)
    θ = [0.4, -0.2]
    Q = precision_matrix(spde, θ)
    expected = -0.5 * length(spde) * log(2π) + 0.5 * logdet(cholesky(Symmetric(Q)))
    @test log_normalizing_constant(spde, θ)≈expected rtol=1.0e-12
end
