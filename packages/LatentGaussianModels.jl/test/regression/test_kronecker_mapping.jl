using LatentGaussianModels
using LatentGaussianModels: AbstractObservationMapping, IdentityMapping,
                            LinearProjector, StackedMapping, KroneckerMapping,
                            apply!, apply_adjoint!, nrows, ncols, as_matrix
using LinearAlgebra
using SparseArrays
using Random

# Phase M PR-1: `KroneckerMapping` operator. The struct ships in
# Phase G with `apply!` / `apply_adjoint!` stubbed; this test
# exercises the real implementations introduced here. Correctness is
# anchored to `as_matrix(m) * x` / `as_matrix(m)' * y` (which uses
# `kron` and is closed-form), and we verify that the operator path
# avoids materializing the dense Kronecker by using small block
# sizes that would still trip dimension errors if the algebra were
# wrong.

@testset "KroneckerMapping — dimensions" begin
    A_space = LinearProjector(randn(MersenneTwister(1), 5, 3))
    A_time = LinearProjector(randn(MersenneTwister(2), 7, 4))
    m = KroneckerMapping(A_space, A_time)

    @test nrows(m) == 5 * 7
    @test ncols(m) == 3 * 4
    @test size(m) == (35, 12)
end

@testset "KroneckerMapping — apply! matches as_matrix" begin
    rng = MersenneTwister(11)
    A_space = LinearProjector(randn(rng, 5, 3))
    A_time = LinearProjector(randn(rng, 7, 4))
    m = KroneckerMapping(A_space, A_time)

    x = randn(rng, ncols(m))
    η = similar(x, nrows(m))
    apply!(η, m, x)

    expected = as_matrix(m) * x
    @test η ≈ expected atol = 1.0e-12
end

@testset "KroneckerMapping — apply_adjoint! matches as_matrix'" begin
    rng = MersenneTwister(12)
    A_space = LinearProjector(randn(rng, 5, 3))
    A_time = LinearProjector(randn(rng, 7, 4))
    m = KroneckerMapping(A_space, A_time)

    r = randn(rng, nrows(m))
    g = similar(r, ncols(m))
    apply_adjoint!(g, m, r)

    expected = as_matrix(m)' * r
    @test g ≈ expected atol = 1.0e-12
end

@testset "KroneckerMapping — adjoint dot-product identity" begin
    # ⟨M x, r⟩ = ⟨x, Mᵀ r⟩ for any x, r — a stronger correctness
    # check than just matching `as_matrix`, because both branches
    # share that reference.
    rng = MersenneTwister(13)
    A_space = LinearProjector(randn(rng, 5, 3))
    A_time = LinearProjector(randn(rng, 7, 4))
    m = KroneckerMapping(A_space, A_time)

    x = randn(rng, ncols(m))
    r = randn(rng, nrows(m))

    η = similar(x, nrows(m))
    g = similar(r, ncols(m))
    apply!(η, m, x)
    apply_adjoint!(g, m, r)

    @test dot(η, r) ≈ dot(x, g) atol = 1.0e-12
end

@testset "KroneckerMapping — identity factor short-circuit" begin
    # `IdentityMapping ⊗ A_time` should equal a block-diagonal
    # repeated `A_time`; verify the operator path agrees with
    # `as_matrix`.
    rng = MersenneTwister(14)
    A_time = LinearProjector(randn(rng, 4, 2))
    m = KroneckerMapping(IdentityMapping(3), A_time)

    @test nrows(m) == 3 * 4
    @test ncols(m) == 3 * 2

    x = randn(rng, ncols(m))
    η = similar(x, nrows(m))
    apply!(η, m, x)
    @test η ≈ as_matrix(m) * x atol = 1.0e-12

    r = randn(rng, nrows(m))
    g = similar(r, ncols(m))
    apply_adjoint!(g, m, r)
    @test g ≈ as_matrix(m)' * r atol = 1.0e-12
end

@testset "KroneckerMapping — sparse factors" begin
    # Real SPDE space-time cases use a sparse mesh projector
    # (LinearProjector wrapping a SparseMatrixCSC) crossed with a
    # small dense AR1 selector. Ensure the mixed-density case works.
    rng = MersenneTwister(15)
    A_sparse = sprandn(rng, 6, 4, 0.3)
    A_dense  = randn(rng, 5, 3)
    m = KroneckerMapping(LinearProjector(A_sparse), LinearProjector(A_dense))

    x = randn(rng, ncols(m))
    η = similar(x, nrows(m))
    apply!(η, m, x)
    @test η ≈ as_matrix(m) * x atol = 1.0e-12

    r = randn(rng, nrows(m))
    g = similar(r, ncols(m))
    apply_adjoint!(g, m, r)
    @test g ≈ as_matrix(m)' * r atol = 1.0e-12
end

@testset "KroneckerMapping — dimension errors" begin
    A_space = LinearProjector(randn(MersenneTwister(20), 5, 3))
    A_time = LinearProjector(randn(MersenneTwister(21), 7, 4))
    m = KroneckerMapping(A_space, A_time)

    @test_throws DimensionMismatch apply!(zeros(nrows(m)), m, zeros(ncols(m) + 1))
    @test_throws DimensionMismatch apply!(zeros(nrows(m) - 1), m, zeros(ncols(m)))
    @test_throws DimensionMismatch apply_adjoint!(zeros(ncols(m)), m, zeros(nrows(m) + 1))
    @test_throws DimensionMismatch apply_adjoint!(zeros(ncols(m) - 1), m, zeros(nrows(m)))
end

@testset "KroneckerMapping — composability with StackedMapping" begin
    # Two response blocks sharing the same latent: each block is a
    # KroneckerMapping. The Stacked container should propagate
    # `apply!` / `apply_adjoint!` through to the Kronecker
    # implementations transparently.
    rng = MersenneTwister(22)
    n_x = 6                                     # shared latent length
    A1_s = LinearProjector(randn(rng, 2, 3))
    A1_t = LinearProjector(randn(rng, 3, 2))
    A2_s = LinearProjector(randn(rng, 4, 3))
    A2_t = LinearProjector(randn(rng, 5, 2))
    K1 = KroneckerMapping(A1_s, A1_t)
    K2 = KroneckerMapping(A2_s, A2_t)
    @assert ncols(K1) == n_x
    @assert ncols(K2) == n_x

    rows = [1:nrows(K1), nrows(K1) .+ (1:nrows(K2))]
    stacked = StackedMapping((K1, K2), rows)

    @test ncols(stacked) == n_x
    @test nrows(stacked) == nrows(K1) + nrows(K2)

    x = randn(rng, n_x)
    η = similar(x, nrows(stacked))
    apply!(η, stacked, x)

    expected_top = as_matrix(K1) * x
    expected_bot = as_matrix(K2) * x
    @test η[rows[1]] ≈ expected_top atol = 1.0e-12
    @test η[rows[2]] ≈ expected_bot atol = 1.0e-12

    r = randn(rng, nrows(stacked))
    g = similar(r, n_x)
    apply_adjoint!(g, stacked, r)
    expected_g = as_matrix(K1)' * r[rows[1]] + as_matrix(K2)' * r[rows[2]]
    @test g ≈ expected_g atol = 1.0e-12
end

@testset "KroneckerMapping — recursive Kronecker" begin
    # KroneckerMapping(KroneckerMapping(A, B), C) — exercise the
    # generic apply! recursion on a non-trivial nesting. Concrete
    # use case: separable space-time with a depth ⊗ time-of-day
    # structure.
    rng = MersenneTwister(23)
    A = LinearProjector(randn(rng, 2, 2))
    B = LinearProjector(randn(rng, 3, 2))
    C = LinearProjector(randn(rng, 2, 3))
    inner = KroneckerMapping(A, B)
    outer = KroneckerMapping(inner, C)

    @test nrows(outer) == 2 * 3 * 2
    @test ncols(outer) == 2 * 2 * 3

    x = randn(rng, ncols(outer))
    η = similar(x, nrows(outer))
    apply!(η, outer, x)
    @test η ≈ as_matrix(outer) * x atol = 1.0e-12

    r = randn(rng, nrows(outer))
    g = similar(r, ncols(outer))
    apply_adjoint!(g, outer, r)
    @test g ≈ as_matrix(outer)' * r atol = 1.0e-12
end
