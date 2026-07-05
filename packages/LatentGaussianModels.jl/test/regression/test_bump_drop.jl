using LatentGaussianModels: GaussianLikelihood, Besag, LatentGaussianModel, laplace_mode
using GMRFs: GMRFGraph, num_nodes, nconnected_components
using SparseArrays
using LinearAlgebra: I
using Random

# ADR-045: for a constrained intrinsic component the null-space bump
# V Vᵀ = Cᵀ(CCᵀ)⁻¹C is dense (O(n²)) but rank-k, and is unnecessary whenever
# H_s = Q + JᵀDJ is already positive definite — every constraint-corrected
# quantity is invariant to precision changes confined to range(Cᵀ). A
# Gaussian likelihood with an identity projector identifies every latent
# coordinate, so H_s is PD and the bump must be dropped: the stored
# posterior precision stays sparse instead of densifying to n² nonzeros.
# (End-to-end result equivalence is covered by the constrained oracle and
# regression suites; here we guard the sparsity that the win depends on.)

# Rook-adjacency grid. Grids are used (rather than tiny hand graphs) so the
# sparse Besag precision — O(n) nonzeros — is unambiguously far below the n²
# the dense null-space bump would produce, giving the sparsity assertion a
# comfortable margin.
function _grid_adj(nr, nc)
    n = nr * nc
    W = spzeros(Int, n, n)
    idx = (i, j) -> (i - 1) * nc + j
    for i in 1:nr, j in 1:nc
        k = idx(i, j)
        i < nr && (W[k, idx(i + 1, j)] = 1; W[idx(i + 1, j), k] = 1)
        j < nc && (W[k, idx(i, j + 1)] = 1; W[idx(i, j + 1), k] = 1)
    end
    return W
end

@testset "ADR-045 — constrained precision stays sparse (bump dropped)" begin
    rng = Random.Xoshiro(20260705)
    W_disc = let g6 = _grid_adj(6, 6), z = spzeros(Int, 36, 36)
        [g6 z; z g6]                                       # two disconnected 6×6 grids
    end
    for (name, W) in (("connected 6×6 grid", _grid_adj(6, 6)),
        ("disconnected two-block", W_disc))
        g = GMRFGraph(W)
        n = num_nodes(g)
        model = LatentGaussianModel(GaussianLikelihood(), (Besag(g),),
            sparse(1.0I, n, n))
        y = randn(rng, n)
        lp = laplace_mode(model, y, [log(3.0), log(0.9)])

        @testset "$name" begin
            @test lp.constraint !== nothing
            @test size(lp.constraint.C, 1) == nconnected_components(g)
            # Bump dropped: a grid Besag precision has O(n) nonzeros, far
            # below the (block-)dense n² the null-space bump would produce.
            @test nnz(lp.precision) < n^2 ÷ 4
            # Still a valid constrained fit: C x̂ = 0 per component.
            @test maximum(abs, lp.constraint.C * lp.mode) < 1.0e-8
            @test lp.converged
        end
    end
end
