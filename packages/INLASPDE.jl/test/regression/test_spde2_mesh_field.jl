# Phase N PR-7a (ADR-036) — `SPDE2` retains `INLAMesh`.
#
# Constructing via `SPDE2(mesh; …)` stores the mesh in `spde.mesh`;
# constructing via `SPDE2(points, triangles; …)` leaves it `nothing`.
# `MeshProjector(spde.mesh, locs)` on the mesh-aware constructor must
# agree with `MeshProjector(mesh, locs)` built from the source mesh.

using Random: MersenneTwister

@testset "SPDE2 — mesh field round-trips through SPDE2(::INLAMesh; …)" begin
    sq = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
    mesh = inla_mesh_2d(; boundary=sq, max_edge=0.25, min_angle=25.0)

    pc = PCMatern(range_U=0.5, range_α=0.05,
        sigma_U=1.0, sigma_α=0.01)
    spde = SPDE2(mesh; pc=pc)

    @test spde.mesh === mesh
    @test spde.mesh isa INLAMesh
    @test num_vertices(spde.mesh) == length(spde)

    # Precision built via the mesh constructor matches the one built
    # via the raw (points, triangles) path on the same mesh — same FEM,
    # same graph, same PC prior.
    spde_raw = SPDE2(mesh.points, mesh.triangles; pc=pc)
    @test spde_raw.mesh === nothing
    Q_mesh = LatentGaussianModels.precision_matrix(spde, [0.0, 0.0])
    Q_raw = LatentGaussianModels.precision_matrix(spde_raw, [0.0, 0.0])
    @test Matrix(Q_mesh)≈Matrix(Q_raw) rtol=1.0e-12
end

@testset "SPDE2 — MeshProjector(spde.mesh, locs) ≡ MeshProjector(mesh, locs)" begin
    sq = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
    mesh = inla_mesh_2d(; boundary=sq, max_edge=0.2, min_angle=25.0)
    spde = SPDE2(mesh)

    rng = MersenneTwister(1)
    locs = 0.1 .+ 0.8 .* rand(rng, 30, 2)

    P_via_mesh = MeshProjector(mesh, locs)
    P_via_spde = MeshProjector(spde.mesh, locs)
    # Same mesh, same locs — projectors agree to fp roundoff. Exact
    # equality is not guaranteed because `find_triangle` uses
    # randomised seeding internally; the resulting barycentric weights
    # differ by O(eps) between independent calls.
    @test size(P_via_spde.A) == size(P_via_mesh.A)
    @test P_via_spde.locations == P_via_mesh.locations
    rng_field = MersenneTwister(2)
    x = randn(rng_field, num_vertices(mesh))
    @test (P_via_spde * x)≈(P_via_mesh * x) rtol=1.0e-12
end

@testset "SPDE2 — α=2 invariant on the mesh constructor" begin
    sq = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
    mesh = inla_mesh_2d(; boundary=sq, max_edge=0.5, min_angle=25.0)
    @test_throws ArgumentError SPDE2(mesh; α=1)
    @test_throws ArgumentError SPDE2(mesh; α=3)
end
