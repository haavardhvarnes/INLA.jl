using LatentGaussianModels: GaussianLikelihood, Besag, LatentGaussianModel,
                            laplace_mode, _constrained_marginal_variances
using GMRFs: GMRFGraph, num_nodes, nconnected_components
using SparseArrays
using LinearAlgebra: I, diag, Symmetric
using Random

# `_constrained_marginal_variances` implements the sparse-selinv + kriging
# correction Var(x_i | Cx = e) = (H⁻¹)_ii − (U W⁻¹ Uᵀ)_ii. The existing
# constrained-Laplace test checks a single sum-to-zero constraint on a
# diagonal precision against a closed form. Here we validate the general
# case — a non-diagonal (Besag/Laplacian) precision with multiple
# constraint rows — against a full dense oracle, and cover the disconnected
# graph where the constraint is one sum-to-zero row per connected component.

# Dense reference: constrained covariance of x ~ N(·, H⁻¹) under Cx = e is
# Σ − ΣCᵀ(CΣCᵀ)⁻¹CΣ with Σ = H⁻¹. Its diagonal is the marginal variance.
function _dense_constrained_variances(H::AbstractMatrix, C::AbstractMatrix)
    Σ = inv(Symmetric(Matrix(H)))
    ΣCt = Σ * C'
    Σc = Σ - ΣCt * ((C * ΣCt) \ (C * Σ))
    return diag(Σc)
end

# Two-component graph: {1-2-3} ⊔ {4-5-6}.
const _W_DISCONNECTED = [0 1 0 0 0 0;
                         1 0 1 0 0 0;
                         0 1 0 0 0 0;
                         0 0 0 0 1 0;
                         0 0 0 1 0 1;
                         0 0 0 0 1 0]

# Connected path graph: 1-2-3-4.
const _W_PATH = [0 1 0 0;
                 1 0 1 0;
                 0 1 0 1;
                 0 0 1 0]

@testset "Constrained marginal variances vs dense oracle" begin
    rng = Random.Xoshiro(20260705)

    @testset "connected Besag (1 constraint, non-diagonal H)" begin
        g = GMRFGraph(_W_PATH)
        n = num_nodes(g)
        model = LatentGaussianModel(GaussianLikelihood(), (Besag(g),),
            sparse(1.0I, n, n))
        y = randn(rng, n)
        lp = laplace_mode(model, y, [log(4.0), log(1.0)])

        @test lp.constraint !== nothing
        @test size(lp.constraint.C, 1) == nconnected_components(g) == 1

        got = _constrained_marginal_variances(lp.precision, lp.constraint)
        ref = _dense_constrained_variances(lp.precision, lp.constraint.C)
        @test maximum(abs, got .- ref) < 1.0e-8
        @test all(got .> 0)
    end

    @testset "disconnected Besag (one sum-to-zero per component)" begin
        g = GMRFGraph(_W_DISCONNECTED)
        n = num_nodes(g)
        @test nconnected_components(g) == 2
        model = LatentGaussianModel(GaussianLikelihood(), (Besag(g),),
            sparse(1.0I, n, n))
        y = randn(rng, n)
        lp = laplace_mode(model, y, [log(3.0), log(0.7)])

        # One constraint row per connected component.
        @test lp.constraint !== nothing
        @test size(lp.constraint.C, 1) == 2

        got = _constrained_marginal_variances(lp.precision, lp.constraint)
        ref = _dense_constrained_variances(lp.precision, lp.constraint.C)
        @test maximum(abs, got .- ref) < 1.0e-8
        @test all(got .> 0)
    end

    @testset "cached-factor overload is bit-faithful to re-factorisation" begin
        # The hot path reuses LaplaceResult.factor for selected inversion
        # instead of re-factorising lp.precision. The two must agree
        # exactly (same matrix, same Cholesky), so equality — not just
        # approximate — is the right assertion here.
        for W in (_W_PATH, _W_DISCONNECTED)
            g = GMRFGraph(W)
            n = num_nodes(g)
            model = LatentGaussianModel(GaussianLikelihood(), (Besag(g),),
                sparse(1.0I, n, n))
            y = randn(rng, n)
            lp = laplace_mode(model, y, [log(2.5), log(0.9)])
            v_refactor = _constrained_marginal_variances(lp.precision, lp.constraint)
            v_cached = _constrained_marginal_variances(lp.factor, lp.constraint)
            @test v_cached == v_refactor
        end
    end
end
