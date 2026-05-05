# SPDE–Matérn covariance reproduction in 1D on a fine regular mesh.
#
# Theory: for the SPDE in 1D with smoothness `ν = α - d/2 = α - 0.5`,
#   C(r; κ, σ²) = σ² · 2^(1-ν) / Γ(ν) · (κr)^ν · K_ν(κr).
#
# Closed forms used here:
#   α = 1, ν = 0.5:  C(r) = σ² · exp(-κr),                σ² = 1/(2κτ²)
#   α = 2, ν = 1.5:  C(r) = σ² · (1 + κr) · exp(-κr),     σ² = 1/(4κ³τ²)
#
# This test builds a fine regular 1D mesh, computes a column of Q⁻¹
# via sparse Cholesky, and checks against the closed-form covariance
# at known lags. Tolerances are generous — finite-element error and
# finite-domain (boundary) effects both contribute bias.

using LinearAlgebra
using SparseArrays

"1D Matérn ν = 0.5 (exponential): C(r) = σ² exp(-κr)."
_matern_1d_nu_half(r, κ, σ²) = σ² * exp(-κ * r)

"1D Matérn ν = 1.5: C(r) = σ² (1 + κr) exp(-κr)."
_matern_1d_nu_three_halves(r, κ, σ²) = σ² * (1 + κ * r) * exp(-κ * r)

@testset "Matérn 1D reproduction — α = 1, ν = 0.5" begin
    L = 30.0                      # long domain to push boundaries far
    n = 600                       # 601 vertices, h = 0.05
    κ = 1.0
    σ² = 1.0
    τ = inv(sqrt(2 * κ * σ²))     # σ² = 1 / (2κτ²)

    h = L / n
    points = collect(-L/2:h:L/2)
    segments = hcat(1:n, 2:(n + 1))
    fem = FEMMatrices(points, segments)
    Q = spde_precision(fem, 1, τ, κ)

    F = cholesky(Symmetric(Q))
    c = cld(length(points), 2)         # central vertex
    e_c = zeros(length(points))
    e_c[c] = 1.0
    q_inv_col = F \ e_c

    # Marginal variance at the centre.
    @test q_inv_col[c] ≈ σ² rtol = 0.05

    for di in (1, 5, 10, 20)
        r = di * h
        theoretical = _matern_1d_nu_half(r, κ, σ²)
        empirical = q_inv_col[c + di]
        @test isapprox(empirical, theoretical; rtol = 0.10, atol = 0.02 * σ²)
    end

    # Symmetry around the centre — interior pair, well away from
    # the boundary.
    @test q_inv_col[c + 5] ≈ q_inv_col[c - 5] rtol = 1.0e-2

    # Monotone decay.
    radial = [q_inv_col[c + di] for di in 0:30]
    @test all(diff(radial) .< 0)
end

@testset "Matérn 1D reproduction — α = 2, ν = 1.5" begin
    L = 30.0
    n = 600
    κ = 1.0
    σ² = 1.0
    τ = inv(sqrt(4 * κ^3 * σ²))    # σ² = 1 / (4κ³τ²)

    h = L / n
    points = collect(-L/2:h:L/2)
    segments = hcat(1:n, 2:(n + 1))
    fem = FEMMatrices(points, segments)
    Q = spde_precision(fem, 2, τ, κ)

    F = cholesky(Symmetric(Q))
    c = cld(length(points), 2)
    e_c = zeros(length(points))
    e_c[c] = 1.0
    q_inv_col = F \ e_c

    @test q_inv_col[c] ≈ σ² rtol = 0.05

    for di in (1, 5, 10, 20)
        r = di * h
        theoretical = _matern_1d_nu_three_halves(r, κ, σ²)
        empirical = q_inv_col[c + di]
        @test isapprox(empirical, theoretical; rtol = 0.10, atol = 0.02 * σ²)
    end

    @test q_inv_col[c + 5] ≈ q_inv_col[c - 5] rtol = 1.0e-2
    radial = [q_inv_col[c + di] for di in 0:30]
    @test all(diff(radial) .< 0)
end
