"""
Julia-side benchmark harness: builds the Scotland BYM2, Pennsylvania BYM2, and
Meuse SPDE models from the same oracle fixtures used by `test/oracle/` and
times `inla(model, y)` end-to-end with `BenchmarkTools`.

Returns a `Vector{NamedTuple}` (one per dataset) with the median wall-clock,
allocations, and posterior log-marginal so the comparison stage can verify
both sides are fitting the same problem.
"""

using BenchmarkTools
using JLD2
using LinearAlgebra
using SparseArrays
using Statistics
using GMRFs
using LatentGaussianModels
using INLASPDE
using LatentGaussianModels: PoissonLikelihood, GaussianLikelihood, Intercept,
                            FixedEffects, BYM2, PCPrecision, LatentGaussianModel,
                            inla, log_marginal_likelihood
using INLASPDE: SPDE2, PCMatern

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

const FIXTURE_PATHS = (
    scotland_bym2=joinpath(REPO_ROOT, "packages", "LatentGaussianModels.jl",
        "test", "oracle", "fixtures", "scotland_bym2.jld2"),
    pennsylvania_bym2=joinpath(REPO_ROOT, "packages", "LatentGaussianModels.jl",
        "test", "oracle", "fixtures", "pennsylvania_bym2.jld2"),
    meuse_spde=joinpath(REPO_ROOT, "packages", "INLASPDE.jl",
        "test", "oracle", "fixtures", "meuse_spde.jld2")
)

function _bym2_model(fx)
    inp = fx["input"]
    y = Int.(inp["cases"])
    E = Float64.(inp["expected"])
    x = Float64.(inp["x"])
    W = inp["W"]
    n = length(y)

    ℓ = PoissonLikelihood(; E=E)
    c_int = Intercept()
    c_beta = FixedEffects(1)
    c_bym2 = BYM2(GMRFGraph(W); hyperprior_prec=PCPrecision(1.0, 0.01))
    A = sparse(hcat(ones(n), reshape(x, n, 1),
        Matrix{Float64}(I, n, n), zeros(n, n)))
    return LatentGaussianModel(ℓ, (c_int, c_beta, c_bym2), A), y
end

function _meuse_model(fx)
    y = Float64.(fx["input"]["y"])
    dist = Float64.(fx["input"]["dist"])
    points = fx["mesh"]["loc"]::Matrix{Float64}
    tv = fx["mesh"]["tv"]::Matrix{Int}
    A_field = SparseMatrixCSC{Float64, Int}(fx["A_field"])
    n_obs = length(y)

    spde = SPDE2(points, tv; α=2,
        pc=PCMatern(; range_U=0.5, range_α=0.5,
            sigma_U=1.0, sigma_α=0.5))

    c_int = Intercept(prec=1.0e-3)
    c_dist = FixedEffects(1; prec=1.0e-3)
    A = hcat(ones(n_obs, 1), reshape(dist, n_obs, 1), A_field)

    ℓ = GaussianLikelihood(hyperprior=PCPrecision(1.0, 0.01))
    return LatentGaussianModel(ℓ, (c_int, c_dist, spde), A), y
end

"""
    build_problem(name::Symbol)

Load the named oracle fixture and return `(model, y, fixture)`.
"""
function build_problem(name::Symbol)
    path = FIXTURE_PATHS[name]
    fx = jldopen(path, "r") do f
        f["fixture"]
    end
    if name === :scotland_bym2 || name === :pennsylvania_bym2
        model, y = _bym2_model(fx)
    elseif name === :meuse_spde
        model, y = _meuse_model(fx)
    else
        error("Unknown benchmark dataset: $name")
    end
    return model, y, fx
end

"""
    run_julia_benchmark(name::Symbol; samples = 5, seconds = 60.0)

Run `inla(model, y)` once for warm-up, then time it under `BenchmarkTools`
with `BenchmarkTools.Trial` median + IQR semantics. Returns a NamedTuple
with `(name, median_seconds, iqr_seconds, allocs_bytes, mlik_julia, mlik_R)`.
"""
function run_julia_benchmark(name::Symbol; samples::Int=5, seconds::Float64=60.0)
    model, y, fx = build_problem(name)

    # Warm-up — discard. Trigger compilation of the inla path on this problem.
    res_warm = inla(model, y)
    mlik_julia = log_marginal_likelihood(res_warm)
    mlik_R = Float64(fx["mlik"][1])

    # Timed runs. `evals = 1` because each fit is many seconds — no inner loop.
    bench = @benchmark inla($model, $y) samples=samples seconds=seconds evals=1
    times = bench.times                    # nanoseconds
    med_s = median(times) / 1e9
    q25, q75 = quantile(times, 0.25) / 1e9, quantile(times, 0.75) / 1e9
    return (
        name=name,
        median_seconds=med_s,
        iqr_seconds=q75 - q25,
        min_seconds=minimum(times) / 1e9,
        max_seconds=maximum(times) / 1e9,
        allocs_bytes=bench.memory,
        n_samples=length(times),
        mlik_julia=mlik_julia,
        mlik_R=mlik_R
    )
end

"""
    run_all_julia(; samples = 5, seconds = 60.0)

Bench all three datasets in fixture-load order. Single-threaded BLAS is the
caller's responsibility — `run.jl` sets it via `BLAS.set_num_threads(1)`
before calling in.
"""
function run_all_julia(; samples::Int=5, seconds::Float64=60.0)
    return [
        run_julia_benchmark(:scotland_bym2; samples=samples, seconds=seconds),
        run_julia_benchmark(:pennsylvania_bym2; samples=samples, seconds=seconds),
        run_julia_benchmark(:meuse_spde; samples=samples, seconds=seconds)
    ]
end
