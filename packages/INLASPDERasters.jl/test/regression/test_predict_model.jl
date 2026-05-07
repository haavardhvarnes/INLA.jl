# Phase O PR-1 — `predict_raster(model, res, template)` Gaussian-
# approximation overload. The wrapper lifts the per-component slice
# out of `random_effects(model, res)` and forwards through the
# vertex-vector primitive; tests assert exact agreement with the
# manual form, plus error coverage for component resolution + CRS.

using Test
using SparseArrays: SparseMatrixCSC, sparse
using LatentGaussianModels: LatentGaussianModel, GaussianLikelihood,
                            Intercept, INLAResult, LaplaceResult,
                            random_effects
using INLASPDE: SPDE2, PCMatern, MeshProjector, INLAMesh, inla_mesh_2d,
                num_vertices
using INLASPDERasters: predict_raster
using Rasters: Rasters, Raster, X, Y, EPSG

# Build an `(Intercept, SPDE2[, SPDE2])` model on a small square mesh
# plus a deterministic mock `INLAResult` whose `x_mean` / `x_var` slot
# values are known, so we can compare against the vertex-vector
# `predict_raster` form bit-for-bit.
function _build_test_model(; with_two_spde::Bool=false, no_mesh::Bool=false)
    sq = [0.1 0.1; 0.9 0.1; 0.9 0.9; 0.1 0.9]
    mesh = inla_mesh_2d(; boundary=sq, max_edge=0.25, min_angle=25.0)
    n_v = num_vertices(mesh)
    pc = PCMatern(range_U=0.3, range_α=0.5, sigma_U=1.0, sigma_α=0.5)

    spde_mesh = no_mesh ? SPDE2(mesh.points, mesh.triangles; pc=pc) :
                SPDE2(mesh; pc=pc)

    # Five interior observation points spread across the mesh.
    locs = [0.3 0.3; 0.5 0.3; 0.7 0.3; 0.5 0.5; 0.5 0.7]
    n_obs = size(locs, 1)
    P = MeshProjector(mesh, locs)
    A_spde = SparseMatrixCSC{Float64, Int}(P.A)
    A_int = sparse(1:n_obs, ones(Int, n_obs), 1.0, n_obs, 1)

    if with_two_spde
        spde2 = SPDE2(mesh; pc=pc)
        P2 = MeshProjector(mesh, locs)
        A_spde2 = SparseMatrixCSC{Float64, Int}(P2.A)
        A = hcat(A_int, A_spde, A_spde2)
        components = (Intercept(), spde_mesh, spde2)
    else
        A = hcat(A_int, A_spde)
        components = (Intercept(), spde_mesh)
    end

    model = LatentGaussianModel(GaussianLikelihood(), components, A)
    return model, mesh, n_v
end

# Construct a mock INLAResult with deterministic `x_mean` / `x_var` —
# `random_effects` only reads those two fields, so the rest can be
# trivial. Latent layout: [intercept; spde_block; (optional) spde2_block].
function _mock_inla_result(model::LatentGaussianModel; seed::Int=1)
    nx = model.n_x
    # Reproducible deterministic values for x_mean / x_var.
    x_mean = [Float64(seed) + Float64(i) for i in 1:nx]
    x_var = [0.5 + 0.01 * i for i in 1:nx]
    return INLAResult(
        Float64[], zeros(0, 0), Vector{Float64}[], Float64[],
        LaplaceResult[], x_mean, x_var, Float64[], NaN, nothing
    )
end

@testset "PR-1 Gaussian-approximation predict_raster(model, res, template)" begin
    model, mesh, n_v = _build_test_model()
    res = _mock_inla_result(model)
    re = random_effects(model, res)
    @test haskey(re, "SPDE2[2]")

    xs = 0.2:0.05:0.8
    ys = 0.2:0.05:0.8
    template = Raster(zeros(length(xs), length(ys)),
        (X(collect(xs)), Y(collect(ys))))

    @testset ":mean / :sd / :lower / :upper agree with vertex-vector form" begin
        for q in (:mean, :sd, :lower, :upper)
            r_model = predict_raster(model, res, template;
                component=2, quantity=q)
            r_manual = predict_raster(getfield(re["SPDE2[2]"], q), mesh,
                template)
            @test parent(r_model)≈parent(r_manual) rtol=1e-12
        end
    end

    @testset "component resolution: Int / String / Type{<:SPDE2}" begin
        r_int = predict_raster(model, res, template; component=2)
        r_str = predict_raster(model, res, template; component="SPDE2[2]")
        r_typ = predict_raster(model, res, template; component=SPDE2)
        @test parent(r_int)≈parent(r_str) rtol=1e-12
        @test parent(r_str)≈parent(r_typ) rtol=1e-12
    end

    @testset "level keyword propagates to random_effects" begin
        r_default = predict_raster(model, res, template;
            component=2, quantity=:lower)
        r_tight = predict_raster(model, res, template;
            component=2, quantity=:lower, level=0.50)
        # Tighter level → lower bound is closer to the mean →
        # values are larger (less negative) on average.
        @test all(parent(r_tight) .>= parent(r_default))
    end
end

@testset "PR-1 component-resolution error paths" begin
    model, _, _ = _build_test_model()
    res = _mock_inla_result(model)
    template = Raster(zeros(4, 4), (X(0.3:0.1:0.6), Y(0.3:0.1:0.6)))

    @testset "out-of-range Int index" begin
        @test_throws ArgumentError predict_raster(model, res, template;
            component=99)
        @test_throws ArgumentError predict_raster(model, res, template;
            component=0)
    end

    @testset "unknown component name" begin
        @test_throws ArgumentError predict_raster(model, res, template;
            component="Foo[7]")
    end

    @testset "wrong component type (Intercept at index 1)" begin
        @test_throws ArgumentError predict_raster(model, res, template;
            component=1)
    end

    @testset "ambiguous Type{<:SPDE2} when two SPDE2 components present" begin
        m2, _, _ = _build_test_model(with_two_spde=true)
        r2 = _mock_inla_result(m2)
        @test_throws ArgumentError predict_raster(m2, r2, template;
            component=SPDE2)
        # But by-Int / by-name still resolve unambiguously.
        @test predict_raster(m2, r2, template; component=2) isa Raster
        @test predict_raster(m2, r2, template; component="SPDE2[3]") isa
              Raster
    end

    @testset "SPDE2 without retained mesh" begin
        m_nomesh, _, _ = _build_test_model(no_mesh=true)
        r_nomesh = _mock_inla_result(m_nomesh)
        @test_throws ArgumentError predict_raster(m_nomesh, r_nomesh,
            template; component=2)
    end

    @testset "invalid quantity" begin
        @test_throws ArgumentError predict_raster(model, res, template;
            component=2, quantity=:foo)
    end
end

@testset "PR-1 mesh_crs keyword (ADR-041)" begin
    model, _, _ = _build_test_model()
    res = _mock_inla_result(model)
    xs = 0.3:0.1:0.6
    ys = 0.3:0.1:0.6

    @testset "no CRS metadata + mesh_crs=nothing (default) is back-compat" begin
        plain = Raster(zeros(length(xs), length(ys)),
            (X(collect(xs)), Y(collect(ys))))
        r = predict_raster(model, res, plain; component=2)
        @test r isa Raster
    end

    @testset "matching CRS passes" begin
        ras = Raster(zeros(length(xs), length(ys)),
            (X(collect(xs)), Y(collect(ys))); crs=EPSG(31370))
        r = predict_raster(model, res, ras;
            component=2, mesh_crs=EPSG(31370))
        @test r isa Raster
    end

    @testset "mismatched CRS errors with both sides named" begin
        ras = Raster(zeros(length(xs), length(ys)),
            (X(collect(xs)), Y(collect(ys))); crs=EPSG(4326))
        err = try
            predict_raster(model, res, ras;
                component=2, mesh_crs=EPSG(31370))
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
end
