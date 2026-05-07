# 1D FEM mass + stiffness assembly correctness on hand-meshed
# intervals. The element formulas for a segment of length `h` are:
#
#   C_loc = (h/6) [2 1; 1 2]      G_loc = (1/h) [1 -1; -1 1]
#
# These are the only error-prone pieces; everything downstream
# (`spde_precision`, lumped mass, G₂) is dimension-agnostic and
# already tested.

using LinearAlgebra
using SparseArrays

@testset "assemble_fem_matrices_1d — 3-vertex hand reference" begin
    # Vertices at x = 0, 1, 3 → segments [1, 2] (h = 1) and [2, 3] (h = 2).
    # By-hand assembly of the global matrices:
    #
    #   C = [1/3   1/6   0;
    #        1/6   1     1/3;
    #        0     1/3   2/3]
    #
    #   G₁ = [1     -1     0;
    #         -1    3/2    -1/2;
    #         0     -1/2   1/2]
    points = [0.0, 1.0, 3.0]
    segments = [1 2; 2 3]
    C, G1 = assemble_fem_matrices_1d(points, segments)

    expected_C = [1/3 1/6 0.0;
                  1/6 1.0 1/3;
                  0.0 1/3 2/3]
    expected_G1 = [1.0 -1.0 0.0;
                   -1.0 3/2 -1/2;
                   0.0 -1/2 1/2]
    @test Array(C)≈expected_C atol=1.0e-12
    @test Array(G1)≈expected_G1 atol=1.0e-12

    # Symmetry.
    @test Array(C) ≈ Array(C')
    @test Array(G1) ≈ Array(G1')

    # Mass row sums recover ∫ φ_i dx.
    @test sum(C; dims=2) ≈ [0.5, 1.5, 1.0]

    # Stiffness is constant-preserving: G₁ · 1 = 0.
    @test G1 * ones(3)≈zeros(3) atol=1.0e-12
end

@testset "assemble_fem_matrices_1d — uniform grid scaling" begin
    # On a uniform mesh of step h with n+1 vertices, the interior
    # rows are: G₁[i,i] = 2/h, G₁[i,i±1] = -1/h.
    n = 10
    h = 0.5
    points = collect(0.0:h:(n * h))
    @test length(points) == n + 1
    segments = hcat(1:n, 2:(n + 1))
    C, G1 = assemble_fem_matrices_1d(points, segments)

    @test G1[5, 5]≈2 / h atol=1.0e-12
    @test G1[5, 4]≈-1 / h atol=1.0e-12
    @test G1[5, 6]≈-1 / h atol=1.0e-12
    @test G1[5, 7] == 0          # local support

    # Boundary rows: half-stencil.
    @test G1[1, 1]≈1 / h atol=1.0e-12
    @test G1[n + 1, n + 1]≈1 / h atol=1.0e-12

    # Sum of mass entries = total domain length.
    @test sum(C)≈n * h atol=1.0e-12
end

@testset "FEMMatrices — 1D constructor + α∈{1,2} precision" begin
    points = collect(0.0:0.5:5.0)
    segments = hcat(1:(length(points) - 1), 2:length(points))
    fem = FEMMatrices(points, segments)

    @test size(fem.C) == (length(points), length(points))
    @test size(fem.G1) == (length(points), length(points))
    @test size(fem.C_lumped) == (length(points), length(points))
    @test size(fem.G2) == (length(points), length(points))

    # `spde_precision` already supports α ∈ {1, 2} dimension-agnostically.
    Q1 = spde_precision(fem, 1, 1.0, 1.0)
    Q2 = spde_precision(fem, 2, 1.0, 1.0)
    @test issymmetric(Symmetric(Q1))
    @test issymmetric(Symmetric(Q2))
    @test isposdef(Symmetric(Q1))
    @test isposdef(Symmetric(Q2))

    # Closed forms:
    #   α = 1: Q = τ²(κ²C + G₁)
    #   α = 2: Q = τ²(κ⁴C̃ + 2κ²G₁ + G₂)
    @test Q1≈fem.C + fem.G1 atol=1.0e-12
    @test Q2≈fem.C_lumped + 2 * fem.G1 + fem.G2 atol=1.0e-12
end

@testset "assemble_fem_matrices_1d — argument validation" begin
    @test_throws ArgumentError assemble_fem_matrices_1d([0.0, 1.0], [1 2 3])
    @test_throws ArgumentError assemble_fem_matrices_1d([0.0, 1.0], [1 3])     # out of range
    @test_throws ArgumentError assemble_fem_matrices_1d([0.0, 1.0], [2 2])     # degenerate
    @test_throws ArgumentError assemble_fem_matrices_1d([0.5, 0.5], [1 2])     # zero length
end
