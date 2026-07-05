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

@testset "PR-2 Exceedance — at c = cell-median → P = 0.5 (posterior near-symmetric)" begin
    # Two facts are checked separately here — a split that removes a
    # historical intermittent flake.
    #
    #  (1) Exceedance-reduction correctness. The `= 0.5` statement is
    #      anchored at each cell's *sample median*, not its mean: for
    #      ANY distribution P(η > median) = 0.5 by construction, and for
    #      an even sample size the interpolated median sits strictly
    #      between order stats n/2 and n/2+1, so exactly half the draws
    #      exceed it. Because the median and the exceedance read the
    #      *same* seeded draws, this holds to the ±1/n discreteness of
    #      the empirical CDF — deterministic, and immune to the platform
    #      BLAS/CHOLMOD nondeterminism in the upstream fit.
    #
    #  (2) Near-symmetry. The grid-integrated posterior of η_cell is a
    #      *mixture* of the per-configuration Gaussians, so it is only
    #      approximately symmetric: the per-cell |mean − median| skew
    #      runs to ≈0.16·sd here. Anchoring the `= 0.5` assertion at the
    #      mean (as an earlier version did) folded that skew into the
    #      exceedance, leaving a ~0.04 systematic offset that a fixed
    #      0.05 band could not reliably contain across platforms. We
    #      instead bound the skew directly, with headroom.
    model, res, mesh = _build_and_fit_synthetic_spde()
    xs = 0.3:0.1:0.7
    ys = 0.3:0.1:0.7
    template = Raster(zeros(length(xs), length(ys)),
        (X(collect(xs)), Y(collect(ys))))

    # Same seed + n_samples ⇒ the median and the exceedance rasters are
    # reductions of the identical draw matrix.
    n_samples = 5000
    r_mean = predict_raster(Xoshiro(44), model, res, template;
        component=2, quantity=:mean, n_samples=n_samples)
    r_med = predict_raster(Xoshiro(44), model, res, template;
        component=2, quantity=0.5, n_samples=n_samples)
    r_sd = predict_raster(model, res, template; component=2, quantity=:sd)

    interior = .!isnan.(parent(r_mean))
    @test any(interior)

    # (1) Exceedance at the per-cell median is 0.5 by construction.
    p_at_med = Float64[]
    for k in findall(interior)
        c = parent(r_med)[k]
        r_ex = predict_raster(Xoshiro(44), model, res, template;
            component=2, quantity=Exceedance(c), n_samples=n_samples)
        push!(p_at_med, parent(r_ex)[k])
    end
    @test length(p_at_med) >= 1
    @test all(abs.(p_at_med .- 0.5) .<= 1 / n_samples)

    # (2) Near-symmetry: per-cell mean and median agree to a modest
    # fraction of the local posterior sd. Observed skew ≈0.16·sd; the
    # 0.3·sd bound leaves ~2× headroom for fit nondeterminism.
    means = parent(r_mean)[interior]
    meds = parent(r_med)[interior]
    sds = parent(r_sd)[interior]
    @test all(abs.(means .- meds) .< 0.3 .* sds)
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
