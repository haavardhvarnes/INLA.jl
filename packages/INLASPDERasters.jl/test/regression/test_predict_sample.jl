# Phase O PR-2 — sample-based `predict_raster` (rng, model, res, template)
# overload, plus the `Exceedance{T}` quantity wrapper.

using Test
using Random: Xoshiro
using Statistics: mean
using SparseArrays: SparseMatrixCSC, sparse
using Distributions: Normal, cdf
using LatentGaussianModels: LatentGaussianModel, GaussianLikelihood,
                            Intercept, inla, random_effects
using INLASPDE: SPDE2, PCMatern, MeshProjector, inla_mesh_2d, num_vertices
using INLASPDERasters: predict_raster, Exceedance
using Rasters: Raster, X, Y, EPSG

# Build an end-to-end (Intercept + SPDE2) Gaussian model fit on
# synthetic data. Small mesh + few observations keeps the test fast.
function _build_and_fit_synthetic_spde()
    rng = Xoshiro(20260506)
    sq = [0.1 0.1; 0.9 0.1; 0.9 0.9; 0.1 0.9]
    mesh = inla_mesh_2d(; boundary=sq, max_edge=0.25, min_angle=25.0)
    pc = PCMatern(range_U=0.3, range_α=0.5, sigma_U=1.0, sigma_α=0.5)
    spde = SPDE2(mesh; pc=pc)

    locs = [0.3 0.3; 0.5 0.3; 0.7 0.3; 0.3 0.5; 0.5 0.5;
            0.7 0.5; 0.3 0.7; 0.5 0.7; 0.7 0.7]
    n_obs = size(locs, 1)
    P = MeshProjector(mesh, locs)
    A_spde = SparseMatrixCSC{Float64, Int}(P.A)
    A_int = sparse(1:n_obs, ones(Int, n_obs), 1.0, n_obs, 1)
    A = hcat(A_int, A_spde)

    model = LatentGaussianModel(GaussianLikelihood(), (Intercept(), spde), A)
    y = 0.5 .+ 0.3 .* randn(rng, n_obs)
    res = inla(model, y; int_strategy=:grid)
    return model, res, mesh
end

@testset "PR-2 sample-based predict_raster — :mean agrees with Gaussian approx" begin
    model, res, mesh = _build_and_fit_synthetic_spde()
    xs = 0.2:0.05:0.8
    ys = 0.2:0.05:0.8
    template = Raster(zeros(length(xs), length(ys)),
        (X(collect(xs)), Y(collect(ys))))

    rng = Xoshiro(11)
    r_sample = predict_raster(rng, model, res, template;
        component=2, quantity=:mean, n_samples=5000)
    r_gauss = predict_raster(model, res, template;
        component=2, quantity=:mean)

    # Compare only on cells inside the mesh (NaN in either side is
    # outside the domain).
    a = parent(r_sample)
    b = parent(r_gauss)
    interior = .!isnan.(a) .& .!isnan.(b)
    @test any(interior)
    # MC tolerance: per-cell sd of η is bounded by the SPDE prior sd
    # (~1.0). Five-sigma band on n=5000 draws is ≈ 5/√5000 ≈ 0.07.
    @test all(abs.(a[interior] .- b[interior]) .< 0.1)
end

@testset "PR-2 quantile dispatch — :mean ≈ median (0.5) for symmetric Gaussian posterior" begin
    model, res, mesh = _build_and_fit_synthetic_spde()
    xs = 0.3:0.1:0.7
    ys = 0.3:0.1:0.7
    template = Raster(zeros(length(xs), length(ys)),
        (X(collect(xs)), Y(collect(ys))))

    rng_a = Xoshiro(22)
    rng_b = Xoshiro(22)
    r_mean = predict_raster(rng_a, model, res, template;
        component=2, quantity=:mean, n_samples=5000)
    r_med = predict_raster(rng_b, model, res, template;
        component=2, quantity=0.5, n_samples=5000)

    a = parent(r_mean)
    b = parent(r_med)
    interior = .!isnan.(a) .& .!isnan.(b)
    @test any(interior)
    # Same draws, so MC mean / median agree more tightly than across
    # independent calls.
    @test all(abs.(a[interior] .- b[interior]) .< 0.1)
end

@testset "PR-2 Exceedance — extreme thresholds" begin
    model, res, mesh = _build_and_fit_synthetic_spde()
    xs = 0.3:0.1:0.7
    ys = 0.3:0.1:0.7
    template = Raster(zeros(length(xs), length(ys)),
        (X(collect(xs)), Y(collect(ys))))

    rng = Xoshiro(33)
    r_lo = predict_raster(rng, model, res, template;
        component=2, quantity=Exceedance(-1.0e6), n_samples=200)
    rng = Xoshiro(33)
    r_hi = predict_raster(rng, model, res, template;
        component=2, quantity=Exceedance(1.0e6), n_samples=200)

    interior_lo = .!isnan.(parent(r_lo))
    interior_hi = .!isnan.(parent(r_hi))
    @test any(interior_lo)
    @test all(parent(r_lo)[interior_lo] .== 1.0)
    @test all(parent(r_hi)[interior_hi] .== 0.0)
end

@testset "PR-2 Exceedance — at c = cell-mean → P ≈ 0.5 (Gaussian symmetry)" begin
    # The joint posterior of (η_cell, ...) is Gaussian, so for any
    # cell the marginal posterior of η_cell is symmetric Gaussian
    # around its mean. P(η_cell > μ_cell) = 0.5 exactly. Sampling
    # estimates this to MC tolerance √(0.25 / n) ≈ 0.007 at n=5000.
    model, res, mesh = _build_and_fit_synthetic_spde()
    xs = 0.3:0.1:0.7
    ys = 0.3:0.1:0.7
    template = Raster(zeros(length(xs), length(ys)),
        (X(collect(xs)), Y(collect(ys))))

    r_mean = predict_raster(model, res, template;
        component=2, quantity=:mean)
    interior = .!isnan.(parent(r_mean))
    interior_means = parent(r_mean)[interior]
    @test any(interior)

    # For each interior cell, compute the empirical exceedance at
    # c = that cell's mean. Average over cells should be ≈ 0.5.
    rng = Xoshiro(44)
    n_samples = 5000
    p_at_means = Float64[]
    for c in interior_means
        rng_local = Xoshiro(44)
        r_ex = predict_raster(rng_local, model, res, template;
            component=2, quantity=Exceedance(c), n_samples=n_samples)
        # Pick the cell whose mean is `c` — bit-wise equality on
        # Float64 keys works because we just read the same array.
        idx = findfirst(parent(r_mean) .== c)
        idx === nothing && continue
        push!(p_at_means, parent(r_ex)[idx])
    end
    @test length(p_at_means) >= 1
    # Five-sigma band on each cell: 5 * 0.5/√n_samples ≈ 0.035.
    @test all(abs.(p_at_means .- 0.5) .< 0.05)
end

@testset "PR-2 RNG reproducibility" begin
    model, res, mesh = _build_and_fit_synthetic_spde()
    xs = 0.3:0.1:0.7
    ys = 0.3:0.1:0.7
    template = Raster(zeros(length(xs), length(ys)),
        (X(collect(xs)), Y(collect(ys))))

    rng_a = Xoshiro(99)
    rng_b = Xoshiro(99)
    r_a = predict_raster(rng_a, model, res, template;
        component=2, quantity=:mean, n_samples=200)
    r_b = predict_raster(rng_b, model, res, template;
        component=2, quantity=:mean, n_samples=200)

    a = parent(r_a)
    b = parent(r_b)
    # Identical to working precision. Tightly bounded — BLAS thread
    # ordering can perturb the trailing bits of sparse-dense GEMM.
    interior = .!isnan.(a) .& .!isnan.(b)
    @test a[interior]≈b[interior] rtol=1.0e-12
    @test all(isnan.(a) .== isnan.(b))
end

@testset "PR-2 quantity validation" begin
    model, res, mesh = _build_and_fit_synthetic_spde()
    template = Raster(zeros(4, 4), (X(0.3:0.1:0.6), Y(0.3:0.1:0.6)))
    rng = Xoshiro(55)

    # :sd / :lower / :upper are Gaussian-approx-only; sample-based
    # only accepts :mean.
    @test_throws ArgumentError predict_raster(rng, model, res, template;
        component=2, quantity=:sd, n_samples=10)

    # Real quantile out of [0, 1].
    @test_throws ArgumentError predict_raster(rng, model, res, template;
        component=2, quantity=1.5, n_samples=10)
    @test_throws ArgumentError predict_raster(rng, model, res, template;
        component=2, quantity=-0.1, n_samples=10)

    # n_samples must be positive.
    @test_throws ArgumentError predict_raster(rng, model, res, template;
        component=2, quantity=:mean, n_samples=0)

    # Wholly unknown quantity type.
    @test_throws ArgumentError predict_raster(rng, model, res, template;
        component=2, quantity="median", n_samples=10)
end

@testset "PR-2 mesh_crs keyword (ADR-041)" begin
    model, res, mesh = _build_and_fit_synthetic_spde()
    xs = 0.3:0.1:0.6
    ys = 0.3:0.1:0.6
    rng = Xoshiro(66)

    @testset "matching CRS passes" begin
        ras = Raster(zeros(length(xs), length(ys)),
            (X(collect(xs)), Y(collect(ys))); crs=EPSG(31370))
        r = predict_raster(rng, model, res, ras;
            component=2, quantity=:mean, n_samples=10, mesh_crs=EPSG(31370))
        @test r isa Raster
    end

    @testset "mismatched CRS errors" begin
        ras = Raster(zeros(length(xs), length(ys)),
            (X(collect(xs)), Y(collect(ys))); crs=EPSG(4326))
        @test_throws ArgumentError predict_raster(rng, model, res, ras;
            component=2, quantity=:mean, n_samples=10, mesh_crs=EPSG(31370))
    end
end
