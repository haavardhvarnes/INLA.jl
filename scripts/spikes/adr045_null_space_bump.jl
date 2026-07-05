# ADR-045 validation spike — de-densify the constrained-Laplace null-space bump.
#
# Validates option B against the current dense-bump path: when the sparse
# H_s = Q + JᵀDJ is positive definite, is the dense bump V Vᵀ unnecessary?
# Tests three quantities for bump-invariance on real Besag models:
#   H1: constrained marginal variances       cvar(H_s) == cvar(H)
#   H2: constrained log-det (_log_det_HC)     logdetHC(H_s) == logdetHC(H)
#   H3: the constrained Newton mode is a fixed point of the bump-free step
# plus a Woodbury solve cross-check. Reference `H = H_s + bump` is what
# `laplace_mode` currently factorizes.
#
# Run from the repo root:  julia --startup-file=no scripts/spikes/adr045_null_space_bump.jl
using Pkg
const REPO = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(mktempdir())
Pkg.develop([PackageSpec(path=joinpath(REPO, "packages", "GMRFs.jl")),
    PackageSpec(path=joinpath(REPO, "packages", "LatentGaussianModels.jl"))])
Pkg.add("Random")

using LatentGaussianModels, GMRFs, Random, SparseArrays, LinearAlgebra
const LGM = LatentGaussianModels

# dense constrained marginal variances of N(·, M⁻¹) under C x = e
function cvar_dense(M, C)
    Σ = inv(Symmetric(Matrix(M)))
    ΣCt = Σ * C'
    return diag(Σ - ΣCt * ((C * ΣCt) \ (C * Σ)))
end

# constrained log-det via the code's formula (_log_det_HC), from an explicit M
function logdetHC(M, C)
    F = cholesky(Symmetric(SparseMatrixCSC(M)))
    HinvCt = F \ Matrix(C')
    return logdet(F) + logdet(Symmetric(C * HinvCt)) - logdet(Symmetric(Matrix(C * C')))
end

function grid_adj(nr, nc)
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

function run_case(name, W)
    rng = Xoshiro(4)
    n = size(W, 1)
    model = LatentGaussianModel(
        GaussianLikelihood(), (Besag(GMRFGraph(W)),), sparse(1.0I, n, n))
    y = randn(rng, n)
    θ = [log(3.0), log(0.8)]
    lp = laplace_mode(model, y, θ)                       # reference (dense bump)

    Q = LGM.joint_precision(model, θ)
    J = LGM.joint_effective_jacobian(model, θ)
    A = LGM.as_matrix(model.mapping)
    η = A * lp.mode
    D = Diagonal(-LGM.joint_∇²_η_log_density(model, y, η, θ))
    Hs = SparseMatrixCSC(Q + J' * D * J)                 # sparse, no bump
    C = GMRFs.constraint_matrix(LGM.model_constraints(model))
    H = SparseMatrixCSC(Hs + LGM._null_bump(C))          # = what the code factorizes

    println("\n=== ", name, "  n=", n, " k=", size(C, 1),
        " nnz(Hs)=", nnz(Hs), " nnz(H)=", nnz(H), " ===")
    println("  Hs positive definite:            ", isposdef(Symmetric(Matrix(Hs))))
    println("  reconstruction == lp.precision:  ",
        maximum(abs, Matrix(H) - Matrix(lp.precision)) < 1e-9)

    v_lp = LGM._constrained_marginal_variances(lp.precision, lp.constraint)
    println("  H1 cvar(Hs) vs cvar(H)   max|Δ|: ",
        maximum(abs, cvar_dense(Hs, C) - cvar_dense(H, C)))
    println("     cvar(Hs) vs lp result max|Δ|: ", maximum(abs, cvar_dense(Hs, C) - v_lp))
    println("  H2 logdetHC(Hs) vs (H)      |Δ|: ", abs(logdetHC(Hs, C) - logdetHC(H, C)))

    ∇η = LGM.joint_∇_η_log_density(model, y, η, θ)
    g = J' * ∇η - Q * lp.mode
    Fs = cholesky(Symmetric(Hs))
    Δ = Fs \ g
    U = Fs \ Matrix(C')
    Δ .-= U * ((C * U) \ (C * Δ))                        # project onto null(C)
    println("  H3 ||projected Hs step||_inf:    ", norm(Δ, Inf))

    b = randn(rng, n)
    L = cholesky(Symmetric(Matrix(C * C'))).L
    Vlr = Matrix(C') / L'                                # Vlr Vlrᵀ = bump
    HsInvV = Fs \ Vlr
    Mw = I + Vlr' * HsInvV
    HsInvb = Fs \ b
    x_wood = HsInvb - HsInvV * (Mw \ (Vlr' * HsInvb))
    println("  Woodbury vs dense H\\b    max|Δ|: ",
        maximum(abs, x_wood - Symmetric(Matrix(H)) \ b))
end

run_case("connected 12x12 Besag", grid_adj(12, 12))
run_case("disconnected two-block",
    [grid_adj(5, 5) spzeros(Int, 25, 25); spzeros(Int, 25, 25) grid_adj(5, 5)])
