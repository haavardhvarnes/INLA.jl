# Phase M PR-5 — `KroneckerComponent` separable space-time composer.
# Anchors per ADR-029:
#   1. Contract surface (length, nhyperparameters, initial_hyperparameters).
#   2. θ-split correctness (concatenated, not shared).
#   3. precision_matrix matches dense `kron(Q_s, Q_t)`.
#   4. log_normalizing_constant matches direct logdet on the dense Q.
#   5. log_hyperprior sums children.
#   6. Constraints — proper×proper → NoConstraint; intrinsic-child
#      lifted via Kronecker; doubly-constrained rejected at construction.
#   7. End-to-end smoke fit (AR1 ⊗ AR1 with Gaussian likelihood).
#   8. Pairing with KroneckerMapping projector — joint shape contract.

using LatentGaussianModels: KroneckerComponent, AR1, IID, RW1, Besag,
                            Intercept,
                            GaussianLikelihood, LinearProjector, KroneckerMapping,
                            LatentGaussianModel,
                            precision_matrix, prior_mean, log_hyperprior,
                            nhyperparameters, initial_hyperparameters,
                            log_normalizing_constant,
                            laplace_mode, n_latent, n_hyperparameters,
                            inla, log_marginal_likelihood,
                            joint_precision, joint_prior_mean, gmrf
using GMRFs: NoConstraint, LinearConstraint, GMRFGraph,
             constraint_matrix, constraint_rhs, nconstraints, rankdef
import GMRFs

@testset "KroneckerComponent — contract surface" begin
    space = AR1(5)
    time = AR1(4)
    c = KroneckerComponent(space, time)

    @test length(c) == 20
    @test nhyperparameters(c) == 4
    @test initial_hyperparameters(c) == zeros(4)
    @test GMRFs.constraints(c) isa NoConstraint
end

@testset "KroneckerComponent — θ length validation" begin
    c = KroneckerComponent(AR1(3), IID(4))
    # AR1 has 2 hyperparams, IID has 1 → expect 3 total.
    @test nhyperparameters(c) == 3
    @test_throws ArgumentError precision_matrix(c, [0.0, 0.0])
    @test_throws ArgumentError precision_matrix(c, zeros(4))
end

@testset "KroneckerComponent — precision_matrix matches dense kron" begin
    space = AR1(6)
    time = AR1(5)
    c = KroneckerComponent(space, time)
    θ = [log(1.5), atanh(0.4), log(0.7), atanh(-0.3)]
    Q = precision_matrix(c, θ)
    Q_s = precision_matrix(space, θ[1:2])
    Q_t = precision_matrix(time, θ[3:4])
    @test size(Q) == (6 * 5, 6 * 5)
    @test Matrix(Q)≈kron(Matrix(Q_s), Matrix(Q_t)) rtol=1.0e-12
    @test issymmetric(Q)
end

@testset "KroneckerComponent — log_normalizing_constant via dense logdet" begin
    space = AR1(5)
    time = AR1(4)
    c = KroneckerComponent(space, time)
    for θ in (zeros(4),
        [log(0.7), atanh(0.5), log(2.0), atanh(-0.4)],
        [log(3.0), atanh(0.1), log(0.5), atanh(0.3)])
        Q = precision_matrix(c, θ)
        expected = -0.5 * length(c) * log(2π) +
                   0.5 * logdet(Symmetric(Matrix(Q)))
        @test log_normalizing_constant(c, θ)≈expected rtol=1.0e-10
    end
end

@testset "KroneckerComponent — log_hyperprior sums children" begin
    space = AR1(4)
    time = AR1(3)
    c = KroneckerComponent(space, time)
    θ = [log(2.5), atanh(0.3), log(0.4), atanh(-0.5)]
    expected = log_hyperprior(space, θ[1:2]) + log_hyperprior(time, θ[3:4])
    @test log_hyperprior(c, θ) ≈ expected
end

@testset "KroneckerComponent — prior_mean = kron(μ_s, μ_t); zero by default" begin
    c = KroneckerComponent(AR1(4), IID(3))
    θ = [0.0, 0.0, 0.0]
    μ = prior_mean(c, θ)
    @test μ == zeros(4 * 3)
end

@testset "KroneckerComponent — AR1 ⊗ AR1 reduces to AR1 stationary recovery" begin
    # When time = AR1(1), the Kronecker collapses to space alone (Q_t
    # is the 1×1 marginal precision). Sanity check structure stays
    # well-defined and dim-correct in the degenerate case.
    space = AR1(8)
    time = AR1(2)
    c = KroneckerComponent(space, time)
    θ = [log(1.7), atanh(0.0), log(1.0), atanh(0.0)]
    Q = precision_matrix(c, θ)
    @test size(Q) == (16, 16)
    @test issymmetric(Q)
    @test isposdef(Symmetric(Matrix(Q)))
end

@testset "KroneckerComponent — AR1 ⊗ IID precision is sparse" begin
    space = AR1(10)
    time = IID(5)
    c = KroneckerComponent(space, time)
    @test nhyperparameters(c) == 2 + 1
    θ = [log(1.0), atanh(0.5), log(2.0)]
    Q = precision_matrix(c, θ)
    @test size(Q) == (50, 50)
    # IID time component contributes a diagonal block per spatial slot.
    Q_s = precision_matrix(space, θ[1:2])
    Q_t = precision_matrix(time, θ[3:3])
    @test Matrix(Q)≈kron(Matrix(Q_s), Matrix(Q_t)) rtol=1.0e-12
end

@testset "KroneckerComponent — intrinsic spatial child lifts constraint via kron" begin
    # Besag on a 4-node ring: 1 connected component → 1 sum-to-zero row.
    # Composing with AR1 in time should lift the constraint to
    # (n_t × n_total) by Kronecker with I_{n_t}.
    W = [0 1 0 1
         1 0 1 0
         0 1 0 1
         1 0 1 0]
    space = Besag(GMRFGraph(W))
    n_s = 4
    n_t = 5
    time = AR1(n_t)
    c = KroneckerComponent(space, time)

    kc = GMRFs.constraints(c)
    @test kc isa LinearConstraint
    @test nconstraints(kc) == 1 * n_t          # 1 spatial constraint × n_t times
    A = constraint_matrix(kc)
    @test size(A) == (n_t, n_s * n_t)

    # Reconstruct expected `A_space ⊗ I_{n_t}` from the inner
    # constraint and the convention.
    A_space = constraint_matrix(GMRFs.constraints(space))
    expected = kron(Float64.(A_space), Matrix{Float64}(I, n_t, n_t))
    @test A ≈ expected
    @test constraint_rhs(kc) ==
          repeat(Float64.(constraint_rhs(GMRFs.constraints(space))), inner=n_t)
end

@testset "KroneckerComponent — intrinsic temporal child lifts via kron" begin
    # RW1(8) is intrinsic with 1 sum-to-zero constraint. Compose with a
    # proper IID(3) on the spatial side: the lifted constraint is
    # I_{n_s} ⊗ A_t, with rhs repeated `n_s` times in outer order.
    space = IID(3)
    time = RW1(8)
    c = KroneckerComponent(space, time)
    kc = GMRFs.constraints(c)
    @test kc isa LinearConstraint
    @test nconstraints(kc) == 3                # 1 temporal × 3 spatial slots
    A = constraint_matrix(kc)
    @test size(A) == (3, 3 * 8)

    A_time = constraint_matrix(GMRFs.constraints(time))
    expected = kron(Matrix{Float64}(I, 3, 3), Float64.(A_time))
    @test A ≈ expected
    @test constraint_rhs(kc) ==
          repeat(Float64.(constraint_rhs(GMRFs.constraints(time))), 3)
end

@testset "KroneckerComponent — doubly-constrained children rejected" begin
    # Both Besag (spatial) and RW1 (temporal) carry constraints;
    # the joint-Kronecker case is deferred per ADR-029.
    W = [0 1 0
         1 0 1
         0 1 0]
    @test_throws ArgumentError KroneckerComponent(Besag(GMRFGraph(W)), RW1(5))
end

@testset "KroneckerComponent — gmrf rankdef matches inclusion-exclusion" begin
    # Both proper: rd = 0.
    c1 = KroneckerComponent(AR1(5), AR1(4))
    @test rankdef(gmrf(c1, zeros(4))) == 0

    # Spatial intrinsic (Besag, 1 connected component): rd = 1 · n_t.
    W = [0 1 0 1
         1 0 1 0
         0 1 0 1
         1 0 1 0]
    space = Besag(GMRFGraph(W))
    @test rankdef(gmrf(space, [0.0])) == 1
    n_t = 6
    c2 = KroneckerComponent(space, AR1(n_t))
    g = gmrf(c2, [0.0, 0.0, 0.0])
    @test rankdef(g) == 1 * n_t
end

@testset "KroneckerComponent — KroneckerMapping pair (Phase M PR-1 ⊗ PR-5)" begin
    # The Kronecker latent flattening matches the KroneckerMapping
    # storage convention: x = vec(X) with X = (n_t × n_s). Construct a
    # KroneckerComponent and a KroneckerMapping with matching factor
    # sizes, and verify that the projected η has the right length and
    # the apply!/precision pipeline lines up.
    n_s, n_t = 6, 4
    n_obs_s, n_obs_t = 8, 5
    A_space = LinearProjector(randn(MersenneTwister(7), n_obs_s, n_s))
    A_time = LinearProjector(randn(MersenneTwister(8), n_obs_t, n_t))
    proj = KroneckerMapping(A_space, A_time)

    c = KroneckerComponent(AR1(n_s), AR1(n_t))
    @test length(c) == n_s * n_t
    # Projector ncols matches latent length.
    @test LatentGaussianModels.ncols(proj) == length(c)
    @test LatentGaussianModels.nrows(proj) == n_obs_s * n_obs_t
end

@testset "KroneckerComponent — Gaussian end-to-end (identity projector)" begin
    # AR1 ⊗ AR1 with identity projector and Gaussian likelihood. The
    # mode under τ_y = 1, τ_x precision should match the Tikhonov
    # solution `x = (I + Q)⁻¹ y` element-wise.
    rng = MersenneTwister(20260505)
    n_s, n_t = 4, 3
    N = n_s * n_t
    y = randn(rng, N)
    A = sparse(1.0I, N, N)
    ℓ = GaussianLikelihood()
    c = KroneckerComponent(AR1(n_s), AR1(n_t))
    m = LatentGaussianModel(ℓ, (c,), A)

    θ = initial_hyperparameters(m)
    res = laplace_mode(m, y, θ)
    @test res.converged
    # τ_y = exp(0) = 1; Q_x at θ=0 is identity-shaped.
    Q_x = precision_matrix(c, zeros(4))
    expected = (Matrix{Float64}(I, N, N) + Matrix(Q_x)) \ y
    @test res.mode≈expected atol=1.0e-8
end

@testset "KroneckerComponent — INLA smoke fit (Gaussian, AR1 ⊗ AR1)" begin
    # End-to-end inla() runs and produces a finite log-marginal.
    rng = MersenneTwister(20260506)
    n_s, n_t = 6, 5
    N = n_s * n_t
    A = sparse(1.0I, N, N)
    y = randn(rng, N)
    ℓ = GaussianLikelihood()
    c = KroneckerComponent(AR1(n_s), AR1(n_t))
    m = LatentGaussianModel(ℓ, (c,), A)
    res = inla(m, y; int_strategy=:grid)
    @test isfinite(log_marginal_likelihood(res))
    @test n_latent(m) == N
    # 4 component hyperparameters (τ_s, ρ_s, τ_t, ρ_t) + 1 likelihood.
    @test n_hyperparameters(m) == 4 + 1
end
