using LatentGaussianModels: _symmetrize!, _is_structurally_symmetric
using SparseArrays
using LinearAlgebra: I, Diagonal
using Random

# `_symmetrize!` cleans up assembly-order float asymmetry in the posterior
# precision `H = Qᵀ + JᵀDJ`. The in-place path must be bit-identical to the
# reference `(H + H') ./ 2` and must not allocate; a structurally
# asymmetric input must fall back to the allocating form.

@testset "_symmetrize! — in-place, bit-identical to (H + H')/2" begin
    rng = Random.Xoshiro(20260705)
    for _ in 1:100
        n = rand(rng, 3:60)
        J = sprandn(rng, rand(rng, n:(2n)), n, 0.2)
        D = Diagonal(abs.(randn(rng, size(J, 1))) .+ 0.1)
        A = sprandn(rng, n, n, 0.1)
        H = SparseMatrixCSC(A'A + I + J' * D * J)      # structurally symmetric
        @test _is_structurally_symmetric(H)
        ref = (H + H') ./ 2
        got = _symmetrize!(copy(H))
        @test got == ref                               # exact, not approximate
    end
end

@testset "_symmetrize! — zero allocation on the in-place path" begin
    rng = Random.Xoshiro(1)
    A = sprandn(rng, 300, 300, 0.03)
    H = SparseMatrixCSC(A'A + I)
    _symmetrize!(copy(H))                              # warm up
    Hc = copy(H)
    @test @allocated(_symmetrize!(Hc)) == 0
end

@testset "_symmetrize! — structurally asymmetric input falls back" begin
    # (1,2) stored with no (2,1): not structurally symmetric.
    Ha = sparse([1, 1, 2], [1, 2, 2], [1.0, 2.0, 3.0], 2, 2)
    @test !_is_structurally_symmetric(Ha)
    @test _symmetrize!(copy(Ha)) == (Ha + Ha') ./ 2
end
