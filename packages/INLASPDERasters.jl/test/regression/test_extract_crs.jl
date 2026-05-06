# Phase O PR-1 — `extract_at_mesh` `mesh_crs` keyword (ADR-041).
# Closes the CRS-test gap promised in `packages/INLASPDERasters.jl/CLAUDE.md`.

using Test
using INLASPDE: inla_mesh_2d
using INLASPDERasters: extract_at_mesh
using Rasters: Rasters, Raster, X, Y, EPSG

@testset "extract_at_mesh mesh_crs keyword" begin
    sq = [0.1 0.1; 0.9 0.1; 0.9 0.9; 0.1 0.9]
    mesh = inla_mesh_2d(; boundary=sq, max_edge=0.25, min_angle=25.0)
    xs = 0.0:0.1:1.0
    ys = 0.0:0.1:1.0

    @testset "default mesh_crs=nothing is back-compat (no check)" begin
        plain = Raster(rand(length(xs), length(ys)),
            (X(collect(xs)), Y(collect(ys))))
        v = extract_at_mesh(plain, mesh; outside=:missing)
        @test length(v) == size(mesh.points, 1)
    end

    @testset "matching CRS passes" begin
        r = Raster(rand(length(xs), length(ys)),
            (X(collect(xs)), Y(collect(ys))); crs=EPSG(31370))
        v = extract_at_mesh(r, mesh; mesh_crs=EPSG(31370), outside=:missing)
        @test length(v) == size(mesh.points, 1)
    end

    @testset "mismatched CRS errors at API boundary" begin
        r = Raster(rand(length(xs), length(ys)),
            (X(collect(xs)), Y(collect(ys))); crs=EPSG(4326))
        err = try
            extract_at_mesh(r, mesh; mesh_crs=EPSG(31370))
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("CRS mismatch", msg)
        @test occursin("31370", msg)
        @test occursin("4326", msg)
    end

    @testset "mesh_crs supplied but raster has no CRS errors" begin
        plain = Raster(rand(length(xs), length(ys)),
            (X(collect(xs)), Y(collect(ys))))
        # Rasters.crs(plain) === nothing, so EPSG(31370) ≠ nothing →
        # mismatch error.
        @test_throws ArgumentError extract_at_mesh(plain, mesh;
            mesh_crs=EPSG(31370))
    end
end
