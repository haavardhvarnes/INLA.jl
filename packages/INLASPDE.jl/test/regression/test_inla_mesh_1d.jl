# 1D mesh generation: sorting, deduplication, max_edge subdivision,
# explicit boundary expansion. Mirrors `test_inla_mesh_2d.jl` for the
# bits that have a 1D analogue.

using LinearAlgebra

@testset "inla_mesh_1d — basic refinement" begin
    loc = [0.0, 1.0, 2.0]
    mesh = inla_mesh_1d(loc; max_edge=0.4)

    @test mesh isa INLAMesh1D
    @test issorted(mesh.points)
    @test first(mesh.points) == 0.0
    @test last(mesh.points) == 2.0

    # Each gap should be ≤ max_edge.
    diffs = diff(mesh.points)
    @test all(diffs .<= 0.4 + 1.0e-12)

    # 2.0 / 0.4 = 5; ceil(1/0.4) = 3 ⇒ each unit segment becomes 3
    # subsegments. Total: 6 + 1 = 7 vertices.
    @test length(mesh.points) == 7
    @test num_segments(mesh) == 6
    @test mesh.boundary == (1, 7)

    # `segments` is consecutive-pair topology.
    @test mesh.segments == hcat(1:6, 2:7)
end

@testset "inla_mesh_1d — unsorted + duplicate input" begin
    loc = [2.0, 0.0, 1.0, 1.0, 0.0]
    mesh = inla_mesh_1d(loc; max_edge=1.0, cutoff=0.0)

    @test issorted(mesh.points)
    @test first(mesh.points) == 0.0
    @test last(mesh.points) == 2.0
    # cutoff = 0 strips exact duplicates only, leaving (0, 1, 2) →
    # max_edge 1 → no further refinement.
    @test mesh.points == [0.0, 1.0, 2.0]
end

@testset "inla_mesh_1d — cutoff merges close points" begin
    loc = [0.0, 0.005, 0.01, 1.0]
    mesh = inla_mesh_1d(loc; max_edge=1.0, cutoff=0.05)

    # 0.005 and 0.01 are within cutoff of 0.0 → merged.
    @test mesh.points[1] == 0.0
    @test 1.0 in mesh.points
    @test all(diff(mesh.points) .> 0.05)
end

@testset "inla_mesh_1d — explicit boundary widens the domain" begin
    loc = [0.5, 1.5]
    mesh = inla_mesh_1d(loc; max_edge=0.3, boundary=(0.0, 2.0))
    @test first(mesh.points) == 0.0
    @test last(mesh.points) == 2.0
    @test 0.5 in mesh.points
    @test 1.5 in mesh.points
    @test all(diff(mesh.points) .<= 0.3 + 1.0e-12)
end

@testset "inla_mesh_1d — boundary inconsistency rejected" begin
    @test_throws ArgumentError inla_mesh_1d([0.0, 5.0]; max_edge=1.0,
        boundary=(1.0, 4.0))
    @test_throws ArgumentError inla_mesh_1d([0.0, 1.0]; max_edge=1.0,
        boundary=(1.0, 0.0))
end

@testset "inla_mesh_1d — argument validation" begin
    @test_throws ArgumentError inla_mesh_1d([0.0, 1.0]; max_edge=0.0)
    @test_throws ArgumentError inla_mesh_1d([0.0, 1.0]; max_edge=-1.0)
    @test_throws ArgumentError inla_mesh_1d([0.0, 1.0]; max_edge=1.0,
        cutoff=-0.1)
    @test_throws ArgumentError inla_mesh_1d(Float64[]; max_edge=1.0)
end

@testset "FEMMatrices(mesh::INLAMesh1D) convenience" begin
    mesh = inla_mesh_1d([0.0, 1.0, 2.0]; max_edge=0.5)
    fem = FEMMatrices(mesh)
    @test size(fem.C, 1) == num_vertices(mesh)
    Q = spde_precision(fem, 2, 1.0, 1.0)
    @test isposdef(Symmetric(Q))
end
