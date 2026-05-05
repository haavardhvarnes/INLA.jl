# Runtime helpers used by the `@lgm` expansion: bind column symbols to
# `Tables.jl`-compatible data sources, validate column existence, and
# build the per-term blocks of the linear projector matrix.

"""
    _check_columns(data, lhs, covariates, randoms)

Validate that `data` is a `Tables.jl` source and that every referenced
column is present. Errors refer to user-visible names, not table
internals.
"""
function _check_columns(data, lhs::Symbol, covariates::AbstractVector,
        randoms::AbstractVector)
    Tables.istable(data) ||
        throw(ArgumentError("@lgm: `data` is not a Tables.jl-compatible source (got $(typeof(data)))"))
    names = Tuple(Tables.columnnames(Tables.columns(data)))
    lhs in names ||
        throw(ArgumentError("@lgm: outcome column `$(lhs)` not found in `data`. Available columns: $(names)"))
    for c in covariates
        c isa Symbol ||
            throw(ArgumentError("@lgm: covariate must be a column name (Symbol), got $(typeof(c))"))
        c in names ||
            throw(ArgumentError("@lgm: covariate column `$(c)` not found in `data`. Available columns: $(names)"))
    end
    for term in randoms
        col_name = first(term)
        col_name isa Symbol ||
            throw(ArgumentError("@lgm: f-term column must be a Symbol, got $(typeof(col_name))"))
        col_name in names ||
            throw(ArgumentError("@lgm: f-term column `$(col_name)` not found in `data`. Available columns: $(names)"))
    end
    return nothing
end

"""
    _build_design_matrix(data, lhs, has_intercept, covariates, randoms)
        -> SparseMatrixCSC{Float64, Int}

Assemble the linear projector matrix `A` for the formula. Columns are
ordered as `[intercept | covariates | random-effect indicators...]`.

Each random-effect term contributes `length(comp)` columns; the
indicator at row `i` is 1 in column `idx[i]` (within that block). The
input column must contain integers in `1:length(comp)`.

PR-1 restriction: a column for an `f(col, Component)` term must be an
integer index column. Categorical / string-keyed columns ship in PR-3.
"""
function _build_design_matrix(data, lhs::Symbol, has_intercept::Bool,
        covariates::AbstractVector, randoms::AbstractVector)
    _check_columns(data, lhs, covariates, randoms)
    cols = Tables.columns(data)
    y = Tables.getcolumn(cols, lhs)
    n_obs = length(y)
    blocks = SparseMatrixCSC{Float64, Int}[]

    if has_intercept
        push!(blocks, sparse(ones(Float64, n_obs, 1)))
    end

    if !isempty(covariates)
        X = Matrix{Float64}(undef, n_obs, length(covariates))
        for (j, name) in enumerate(covariates)
            col = Tables.getcolumn(cols, name)
            length(col) == n_obs ||
                throw(DimensionMismatch("@lgm: covariate column `$(name)` has length $(length(col)); outcome `$(lhs)` has length $(n_obs)"))
            @inbounds for i in 1:n_obs
                X[i, j] = Float64(col[i])
            end
        end
        push!(blocks, sparse(X))
    end

    for term in randoms
        col_name, comp = first(term), last(term)
        comp isa LatentGaussianModels.AbstractLatentComponent ||
            throw(ArgumentError("@lgm: f-term `$(col_name)`: second argument must be an `AbstractLatentComponent`, got $(typeof(comp))"))
        col = Tables.getcolumn(cols, col_name)
        length(col) == n_obs ||
            throw(DimensionMismatch("@lgm: f-term column `$(col_name)` has length $(length(col)); outcome `$(lhs)` has length $(n_obs)"))
        n_levels = length(comp)
        idx = _validate_index_column(col, col_name, n_levels)
        push!(blocks, sparse(1:n_obs, idx, 1.0, n_obs, n_levels))
    end

    isempty(blocks) &&
        throw(ArgumentError("@lgm: model has no terms — formula must include at least `1`, a covariate, or `f(...)`"))
    return length(blocks) == 1 ? blocks[1] : reduce(hcat, blocks)
end

function _validate_index_column(col, col_name::Symbol, n_levels::Int)
    idx = Vector{Int}(undef, length(col))
    @inbounds for i in eachindex(col)
        v = col[i]
        v isa Integer ||
            throw(ArgumentError("@lgm: f-term column `$(col_name)` must contain integers in 1:$(n_levels); got entry of type $(typeof(v)) at row $(i)"))
        (1 ≤ v ≤ n_levels) ||
            throw(ArgumentError("@lgm: f-term column `$(col_name)` value $(v) at row $(i) is outside 1:$(n_levels) (component has $(n_levels) levels)"))
        idx[i] = Int(v)
    end
    return idx
end

"""
    lgmformula(data; lhs, intercept = true, covariates = Symbol[],
               randoms = [], family) -> LatentGaussianModel

Function form of [`@lgm`](@ref). Accepts a structured description of
the formula and returns the same `LatentGaussianModel` the macro
would produce.

# Arguments

- `data`: a `Tables.jl`-compatible source.
- `lhs::Symbol`: outcome column name.
- `intercept::Bool = true`: whether to include `Intercept()`.
- `covariates::Vector{Symbol} = Symbol[]`: scalar fixed-effect column
  names. Becomes `FixedEffects(length(covariates))` if non-empty.
- `randoms::Vector{<:Tuple{Symbol, AbstractLatentComponent}} = []`:
  list of `(col, component)` pairs for `f(...)` terms.
- `family::AbstractLikelihood`: observation likelihood.

PR-3 will add a `lgmformula(formula::FormulaTerm, data; ...)` method
that decomposes a `StatsModels.@formula` into the structured form.
"""
function lgmformula(data;
        lhs::Symbol,
        intercept::Bool = true,
        covariates::Vector{Symbol} = Symbol[],
        randoms::AbstractVector = Tuple{Symbol, LatentGaussianModels.AbstractLatentComponent}[],
        family::LatentGaussianModels.AbstractLikelihood)
    A = _build_design_matrix(data, lhs, intercept, covariates, randoms)
    components = _build_components(intercept, length(covariates), randoms)
    return LatentGaussianModels.LatentGaussianModel(family, components, A)
end

function _build_components(has_intercept::Bool, n_covariates::Int,
        randoms::AbstractVector)
    parts = LatentGaussianModels.AbstractLatentComponent[]
    has_intercept && push!(parts, LatentGaussianModels.Intercept())
    n_covariates > 0 &&
        push!(parts, LatentGaussianModels.FixedEffects(n_covariates))
    for term in randoms
        push!(parts, last(term))
    end
    isempty(parts) &&
        throw(ArgumentError("@lgm: model has no components — formula must include at least `1`, a covariate, or `f(...)`"))
    return Tuple(parts)
end
