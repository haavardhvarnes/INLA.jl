# 1D mesh projector: linear interpolation onto a sorted segment mesh.
# Anchors:
#   - Each row sums to 1 (partition of unity).
#   - Linear functions u(x) = a + b·x interpolate exactly.
#   - Outside-domain locations honour the `outside` policy.

using LinearAlgebra
using SparseArrays

@testset "MeshProjector1D — partition of unity" begin
    mesh = inla_mesh_1d([0.0, 1.0, 2.0]; max_edge=0.5)
    locs = [0.0, 0.25, 0.5, 0.7, 1.0, 1.5, 2.0]
    P = MeshProjector1D(mesh, locs)

    @test size(P) == (length(locs), num_vertices(mesh))
    # Each row sums to 1 — barycentric weights for an interior or
    # endpoint location are a partition of unity.
    @test all(isapprox.(sum(P.A; dims=2), 1.0; atol=1.0e-12))
    # At most two nonzeros per row.
    @test maximum(count(!iszero, view(P.A, i, :)) for i in 1:size(P.A, 1)) == 2
end

@testset "MeshProjector1D — interpolates linears exactly" begin
    mesh = inla_mesh_1d(collect(0.0:0.5:5.0); max_edge=0.5)
    locs = [0.123, 0.5, 1.4567, 2.0, 3.0001, 4.999]
    P = MeshProjector1D(mesh, locs)

    # Linear function on mesh vertices.
    a, b = 0.7, -1.3
    u = @. a + b * mesh.points
    interpolated = P * u
    expected = @. a + b * locs
    @test interpolated≈expected atol=1.0e-12
end

@testset "MeshProjector1D — endpoint clamping" begin
    mesh = inla_mesh_1d([0.0, 1.0]; max_edge=1.0)
    P = MeshProjector1D(mesh, [0.0, 1.0])
    @test P.A[1, 1] ≈ 1.0
    @test P.A[1, 2] ≈ 0.0
    @test P.A[2, end - 1] ≈ 0.0
    @test P.A[2, end] ≈ 1.0
end

@testset "MeshProjector1D — outside policy" begin
    mesh = inla_mesh_1d([0.0, 1.0, 2.0]; max_edge=0.5)

    @test_throws ArgumentError MeshProjector1D(mesh, [-0.1])
    @test_throws ArgumentError MeshProjector1D(mesh, [2.5])

    P = MeshProjector1D(mesh, [-0.1, 1.0, 2.5]; outside=:zero)
    @test sum(P.A[1, :]) == 0
    @test sum(P.A[2, :]) ≈ 1.0
    @test sum(P.A[3, :]) == 0

    # atol allows small overshoots; the location is clamped.
    P_atol = MeshProjector1D(mesh, [-1.0e-9]; atol=1.0e-6)
    @test P_atol.A[1, 1]≈1.0 atol=1.0e-6
end

@testset "MeshProjector1D — argument validation" begin
    mesh = inla_mesh_1d([0.0, 1.0]; max_edge=1.0)
    @test_throws ArgumentError MeshProjector1D(mesh, [0.5]; outside=:foo)
    @test_throws ArgumentError MeshProjector1D(mesh, [0.5]; atol=-1.0)
end
