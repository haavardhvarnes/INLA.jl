# Phase O PR-5 — predict_raster reach extended past stationary SPDE2:
#  • SPDE2NonStationary now retains its source mesh under the
#    `…(mesh::INLAMesh; …)` constructor (mirroring SPDE2 / ADR-036), so
#    `predict_raster(model, res, template; component=…)` works without
#    threading a mesh through the call.
#  • KroneckerComponent(spatial::SPDE2-flavored, temporal) gets a
#    raster-projector dispatch with a 1-based `time_index` keyword that
#    slices a single time slot out of the per-component vector before
#    barycentric projection, matching the time-inner-index flattening
#    `x[(s - 1) · n_t + t]` used by KroneckerMapping.
#
# The tests pin both paths against the manual vertex-vector primitive
# `predict_raster(values, mesh, template)`, which v0.2.0 already
# regression-tests bit-for-bit. Equality-against-manual is the sharp
# correctness gate.

using Test
using SparseArrays: SparseMatrixCSC, sparse
using LatentGaussianModels: LatentGaussianModel, GaussianLikelihood,
                            Intercept, INLAResult, LaplaceResult,
                            random_effects, KroneckerComponent, AR1
using INLASPDE: SPDE2, SPDE2NonStationary, PCMatern, GaussianBasisPrior,
                MeshProjector, INLAMesh, inla_mesh_2d, num_vertices
using INLASPDERasters: predict_raster
using Rasters: Rasters, Raster, X, Y

# Reusable square mesh + dense observation block. Mirrors the layout in
# test_predict_model.jl so the two test files stay in lockstep.
function _build_pr5_mesh()
    sq = [0.1 0.1; 0.9 0.1; 0.9 0.9; 0.1 0.9]
    mesh = inla_mesh_2d(; boundary=sq, max_edge=0.25, min_angle=25.0)
    return mesh
end

# Mock INLAResult exactly like test_predict_model.jl — `random_effects`
# touches only `x_mean` / `x_var`, so the rest is filler.
function _mock_inla_result_pr5(model::LatentGaussianModel; seed::Int=1)
    nx = model.n_x
    x_mean = [Float64(seed) + Float64(i) for i in 1:nx]
    x_var = [0.5 + 0.01 * i for i in 1:nx]
    return INLAResult(
        Float64[], zeros(0, 0), Vector{Float64}[], Float64[], Float64[],
        LaplaceResult[], x_mean, x_var, Float64[], NaN, nothing
    )
end

@testset "PR-5 — SPDE2NonStationary retains mesh + predict_raster wrapper" begin
    mesh = _build_pr5_mesh()
    n_v = num_vertices(mesh)
    locs = [0.3 0.3; 0.5 0.3; 0.7 0.3; 0.5 0.5; 0.5 0.7]
    n_obs = size(locs, 1)
    P = MeshProjector(mesh, locs)
    A_int = sparse(1:n_obs, ones(Int, n_obs), 1.0, n_obs, 1)
    A_spde = SparseMatrixCSC{Float64, Int}(P.A)

    # Stationary basis (B_τ = B_κ = ones) so the prior reduces to two
    # scalar coefficients and the precision matches a SPDE2 fit with
    # the same θ. The mesh-aware constructor is what we're testing.
    B_τ = ones(n_v, 1)
    B_κ = ones(n_v, 1)
    prior = GaussianBasisPrior(mean=[0.0, 0.0], prec=[1.0, 1.0])
    spde_ns = SPDE2NonStationary(mesh; B_τ=B_τ, B_κ=B_κ, prior=prior)

    @testset "mesh field populated by mesh-aware constructor" begin
        @test spde_ns.mesh isa INLAMesh
        @test spde_ns.mesh === mesh   # stored by reference, not rebuilt
    end

    @testset "raw points/triangles constructor leaves mesh = nothing" begin
        spde_raw = SPDE2NonStationary(mesh.points, mesh.triangles;
            B_τ=B_τ, B_κ=B_κ, prior=prior)
        @test spde_raw.mesh === nothing
    end

    model = LatentGaussianModel(GaussianLikelihood(),
        (Intercept(), spde_ns), hcat(A_int, A_spde))
    res = _mock_inla_result_pr5(model)

    xs = 0.2:0.05:0.8
    ys = 0.2:0.05:0.8
    template = Raster(zeros(length(xs), length(ys)),
        (X(collect(xs)), Y(collect(ys))))

    @testset ":mean / :sd / :lower / :upper agree with vertex-vector form" begin
        re = random_effects(model, res)
        @test haskey(re, "SPDE2NonStationary[2]")
        for q in (:mean, :sd, :lower, :upper)
            r_model = predict_raster(model, res, template;
                component=2, quantity=q)
            r_manual = predict_raster(
                getfield(re["SPDE2NonStationary[2]"], q), mesh, template)
            @test parent(r_model)≈parent(r_manual) rtol=1e-12
        end
    end

    @testset "component resolution: Int / String / Type{<:SPDE2NonStationary}" begin
        r_int = predict_raster(model, res, template; component=2)
        r_str = predict_raster(model, res, template;
            component="SPDE2NonStationary[2]")
        r_typ = predict_raster(model, res, template;
            component=SPDE2NonStationary)
        @test parent(r_int)≈parent(r_str) rtol=1e-12
        @test parent(r_str)≈parent(r_typ) rtol=1e-12
    end

    @testset "raw-constructor SPDE2NonStationary errors at predict_raster" begin
        spde_raw = SPDE2NonStationary(mesh.points, mesh.triangles;
            B_τ=B_τ, B_κ=B_κ, prior=prior)
        m_raw = LatentGaussianModel(GaussianLikelihood(),
            (Intercept(), spde_raw), hcat(A_int, A_spde))
        r_raw = _mock_inla_result_pr5(m_raw)
        err = try
            predict_raster(m_raw, r_raw, template; component=2)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("SPDE2NonStationary", msg)
        @test occursin("does not retain its mesh", msg)
    end
end

@testset "PR-5 — KroneckerComponent space-time predict_raster + time_index" begin
    mesh = _build_pr5_mesh()
    n_v = num_vertices(mesh)
    n_t = 4
    pc = PCMatern(range_U=0.3, range_α=0.5, sigma_U=1.0, sigma_α=0.5)
    spde = SPDE2(mesh; pc=pc)
    ar1 = AR1(n_t)
    spt = KroneckerComponent(spde, ar1)
    @test length(spt) == n_v * n_t

    # Single observation block — five obs at t = 1 only is enough to
    # build a valid LGM. The actual A-matrix doesn't enter the
    # predict_raster path; only the latent layout / random_effects
    # vector matter.
    locs = [0.3 0.3; 0.5 0.3; 0.7 0.3; 0.5 0.5; 0.5 0.7]
    n_obs = size(locs, 1)
    P = MeshProjector(mesh, locs)
    A_int = sparse(1:n_obs, ones(Int, n_obs), 1.0, n_obs, 1)
    # Crude Kronecker-mapping placeholder: pad each obs row with zeros
    # for the (n_t - 1) remaining time slots. The numerical content
    # doesn't matter for predict_raster's projection-side test.
    A_spt_dense = zeros(n_obs, n_v * n_t)
    A_dense_spde = Matrix(P.A)
    for k in 1:n_obs, s in 1:n_v
        A_spt_dense[k, (s - 1) * n_t + 1] = A_dense_spde[k, s]
    end
    A_spt = sparse(A_spt_dense)
    A = hcat(A_int, A_spt)

    model = LatentGaussianModel(GaussianLikelihood(),
        (Intercept(), spt), A)
    res = _mock_inla_result_pr5(model)

    xs = 0.2:0.05:0.8
    ys = 0.2:0.05:0.8
    template = Raster(zeros(length(xs), length(ys)),
        (X(collect(xs)), Y(collect(ys))))

    re = random_effects(model, res)
    @test haskey(re, "KroneckerComponent[2]")
    full_vec_mean = re["KroneckerComponent[2]"].mean
    @test length(full_vec_mean) == n_v * n_t

    @testset "time_index slice picks rows k:n_t:end (time inner)" begin
        for q in (:mean, :sd, :lower, :upper)
            full_vec = getfield(re["KroneckerComponent[2]"], q)
            for k in 1:n_t
                r_model = predict_raster(model, res, template;
                    component=2, quantity=q, time_index=k)
                r_manual = predict_raster(full_vec[k:n_t:end], mesh, template)
                @test parent(r_model)≈parent(r_manual) rtol=1e-12
            end
        end
    end

    @testset "component resolution: Int / String / Type{<:KroneckerComponent}" begin
        r_int = predict_raster(model, res, template;
            component=2, time_index=1)
        r_str = predict_raster(model, res, template;
            component="KroneckerComponent[2]", time_index=1)
        r_typ = predict_raster(model, res, template;
            component=KroneckerComponent, time_index=1)
        @test parent(r_int)≈parent(r_str) rtol=1e-12
        @test parent(r_str)≈parent(r_typ) rtol=1e-12
    end

    @testset "missing time_index errors with helpful message" begin
        err = try
            predict_raster(model, res, template; component=2)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("time_index", msg)
        @test occursin("KroneckerComponent", msg)
    end

    @testset "out-of-range time_index errors" begin
        @test_throws ArgumentError predict_raster(model, res, template;
            component=2, time_index=0)
        @test_throws ArgumentError predict_raster(model, res, template;
            component=2, time_index=n_t + 1)
        @test_throws ArgumentError predict_raster(model, res, template;
            component=2, time_index=-1)
    end

    @testset "non-integer time_index errors" begin
        @test_throws ArgumentError predict_raster(model, res, template;
            component=2, time_index=1.5)
    end

    @testset "time_index on stationary SPDE2 errors" begin
        # Build a plain SPDE2 model and pass time_index=2; the
        # SPDE2 dispatch should reject the kwarg outright.
        n_v_local = num_vertices(mesh)
        spde_only = SPDE2(mesh; pc=pc)
        P2 = MeshProjector(mesh, locs)
        A_only = hcat(A_int, SparseMatrixCSC{Float64, Int}(P2.A))
        m_only = LatentGaussianModel(GaussianLikelihood(),
            (Intercept(), spde_only), A_only)
        r_only = _mock_inla_result_pr5(m_only)
        err = try
            predict_raster(m_only, r_only, template;
                component=2, time_index=1)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("time_index", msg)
    end

    @testset "KroneckerComponent without retained spatial mesh errors" begin
        spde_raw = SPDE2(mesh.points, mesh.triangles; pc=pc)
        spt_raw = KroneckerComponent(spde_raw, ar1)
        A_raw = hcat(A_int, A_spt)
        m_raw = LatentGaussianModel(GaussianLikelihood(),
            (Intercept(), spt_raw), A_raw)
        r_raw = _mock_inla_result_pr5(m_raw)
        err = try
            predict_raster(m_raw, r_raw, template;
                component=2, time_index=1)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("does not retain its mesh", msg)
    end

    @testset "KroneckerComponent with non-SPDE spatial child errors" begin
        # AR1 ⊗ AR1 — both children are temporal-flavored, no mesh.
        spt_no_spde = KroneckerComponent(AR1(n_v), ar1)
        A_no = hcat(A_int, A_spt)   # right shape doesn't matter for predict
        m_no = LatentGaussianModel(GaussianLikelihood(),
            (Intercept(), spt_no_spde), A_no)
        r_no = _mock_inla_result_pr5(m_no)
        err = try
            predict_raster(m_no, r_no, template;
                component=2, time_index=1)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("KroneckerComponent", msg)
        @test occursin("AR1", msg)
    end
end
