using LatentGaussianModels: Besag, BYM, precision_matrix
using GMRFs: GMRFs, GMRFGraph, BesagGMRF, num_nodes
using SparseArrays
using LinearAlgebra: I

# The Sørbye-Rue per-node scale constants depend only on the graph (they are
# θ-independent) but computing them builds a dense inv(Qperp) per connected
# component. Besag/BYM cache them at construction and reuse them in every
# precision_matrix call. These tests guard two things: the cached result is
# bit-identical to recomputing from scratch, and the per-call cost no longer
# includes the dense scale-factor computation.

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

@testset "scale_factor caching — bit-identical, allocation-light" begin
    W_conn = _grid_adj(6, 6)
    W_disc = let g = _grid_adj(5, 5), z = spzeros(Int, 25, 25)
        [g z; z g]
    end

    @testset "Besag precision == uncached recompute" begin
        for (name, W) in (("connected", W_conn), ("disconnected", W_disc)),
            sm in (true, false)

            g = GMRFGraph(W)
            θ = [log(2.3)]
            cached = precision_matrix(Besag(g; scale_model=sm), θ)
            # Uncached BesagGMRF (no `c`): recomputes scale factors from scratch.
            recompute = GMRFs.precision_matrix(
                BesagGMRF(g; τ=exp(θ[1]), scale_model=sm))
            @test cached == recompute
        end
    end

    @testset "BYM precision finite + symmetric" begin
        g = GMRFGraph(W_conn)
        Q = precision_matrix(BYM(g), [log(1.4), log(0.8)])
        @test size(Q, 1) == 2 * num_nodes(g)
        @test Q == permutedims(Q)
        @test all(isfinite, nonzeros(Q))
    end

    @testset "cached precision_matrix is allocation-light (625 nodes)" begin
        # The dense scale-factor recompute allocated ~18 MiB per call; the
        # cached path is well under 1 MiB. A generous ceiling keeps the guard
        # robust across platforms while still catching a regression to the
        # per-call recompute.
        g = GMRFGraph(_grid_adj(25, 25))
        c = Besag(g)
        θ = [log(2.0)]
        precision_matrix(c, θ)                    # warm
        @test @allocated(precision_matrix(c, θ)) < 4 * 2^20
    end
end
