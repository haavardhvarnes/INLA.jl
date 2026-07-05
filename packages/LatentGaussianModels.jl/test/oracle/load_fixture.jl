# Shared loader for R-INLA oracle fixtures in LatentGaussianModels.jl.
#
# Fixtures are JLD2 files produced by `scripts/generate-fixtures/`.
# Each fixture is stored under the key "fixture" as a `Dict{String, Any}`.
#
# This loader is a no-op if the fixture file is absent — that lets the
# oracle testset degrade to a `@test_skip` so contributors without R can
# still run the suite.

using JLD2

const ORACLE_FIXTURE_DIR = joinpath(@__DIR__, "fixtures")

"""
    oracle_fixture_path(name) -> String

Absolute path to `fixtures/<name>.jld2`.
"""
oracle_fixture_path(name::AbstractString) = joinpath(
    ORACLE_FIXTURE_DIR, string(name, ".jld2"))

"""
    has_oracle_fixture(name) -> Bool
"""
has_oracle_fixture(name::AbstractString) = isfile(oracle_fixture_path(name))

"""
    load_oracle_fixture(name) -> Dict{String, Any}

Load a fixture dict written by `scripts/generate-fixtures/convert_to_jld2.jl`.
Raises if the file is missing — callers should check `has_oracle_fixture`
first and skip.
"""
function load_oracle_fixture(name::AbstractString)
    path = oracle_fixture_path(name)
    isfile(path) ||
        error("oracle fixture not found: $path (run scripts/generate-fixtures/)")
    return jldopen(path, "r") do f
        return f["fixture"]
    end
end

"""
    fixed_summary_mean(fx, rowname) -> Float64

Lookup helper: pull the posterior mean of a fixed effect by row name.
"""
function fixed_summary_mean(fx::AbstractDict, rowname::AbstractString)
    sf = fx["summary_fixed"]
    rn = String.(sf["rownames"])
    idx = findfirst(==(rowname), rn)
    idx === nothing && error("row '$rowname' not found in summary_fixed (have: $rn)")
    return Float64(sf["mean"][idx])
end

"""
    marginal_grid(m) -> (xs::Vector{Float64}, ys::Vector{Float64})

Extract the `(x, y)` density grid from a fixture `marginals_*` entry.
JLD2 conversion stores these as NamedTuples; older fixtures may carry
Dicts — accept both.
"""
function marginal_grid(m)
    if m isa NamedTuple
        return Float64.(m.x), Float64.(m.y)
    end
    return Float64.(m["x"]), Float64.(m["y"])
end

"""
    oracle_pdf_moments(xs, ys) -> @NamedTuple{mean, sd, mass}

Trapezoid mean/sd of a gridded density, renormalised on the grid. Used
by the integrated θ-marginal comparisons (ADR-046): the R-INLA
`marginals_hyperpar` grids and the Julia `posterior_marginal_θ` output
are both reduced to moments through this single quadrature so the
comparison is not sensitive to grid layout.
"""
function oracle_pdf_moments(xs::AbstractVector, ys::AbstractVector)
    Z = 0.0
    μ = 0.0
    for i in 1:(length(xs) - 1)
        h = Float64(xs[i + 1] - xs[i])
        Z += 0.5 * h * (ys[i] + ys[i + 1])
        μ += 0.5 * h * (xs[i] * ys[i] + xs[i + 1] * ys[i + 1])
    end
    Z > 0 || return (mean=NaN, sd=NaN, mass=Z)
    μ /= Z
    σ² = 0.0
    for i in 1:(length(xs) - 1)
        h = Float64(xs[i + 1] - xs[i])
        σ² += 0.5 * h * ((xs[i] - μ)^2 * ys[i] + (xs[i + 1] - μ)^2 * ys[i + 1])
    end
    return (mean=μ, sd=sqrt(max(σ² / Z, 0.0)), mass=Z)
end

"""
    precision_marginal_moments(m) -> @NamedTuple{mean, sd, mass}

User-scale (τ = exp θ) trapezoid moments of an internal-scale
θ-marginal `(θ, pdf)` returned by `posterior_marginal_θ` for a
log-precision hyperparameter: `p_τ(τ) = p_θ(log τ) / τ`.
"""
function precision_marginal_moments(m::NamedTuple)
    τ = exp.(m.θ)
    return oracle_pdf_moments(τ, m.pdf ./ τ)
end

"""
    marginal_quantile(xs, ys, p) -> Float64

`p`-quantile of a gridded density `(xs, ys)` (ascending `xs`) via the
cumulative trapezoid, linearly interpolated. Tail-robust companion to
[`oracle_pdf_moments`](@ref) for heavy-tailed hyperparameter marginals.
"""
function marginal_quantile(xs::AbstractVector, ys::AbstractVector, p::Real)
    n = length(xs)
    c = zeros(Float64, n)
    for i in 2:n
        c[i] = c[i - 1] + 0.5 * Float64(xs[i] - xs[i - 1]) * (ys[i] + ys[i - 1])
    end
    c[end] > 0 || return NaN
    c ./= c[end]
    k = findfirst(≥(p), c)
    k === nothing && return Float64(xs[end])
    k == 1 && return Float64(xs[1])
    t = (p - c[k - 1]) / (c[k] - c[k - 1])
    return Float64(xs[k - 1] + t * (xs[k] - xs[k - 1]))
end
