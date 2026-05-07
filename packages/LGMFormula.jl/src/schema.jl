# Runtime helpers used by the `@lgm` expansion: bind column symbols to
# `Tables.jl`-compatible data sources, validate column existence, and
# build the per-term blocks of the linear projector matrix.

"""
    _wrap_term(comp_or_factory, data, replicate_col, group_col) ->
        AbstractLatentComponent

PR-5 runtime helper. Wraps an `f(col, Comp; replicate = id)` or
`f(col, Factory; group = grp)` term against the actual `data` table:

- `replicate_col::Symbol`: returns `Replicate(comp, R)` where
  `R = maximum(data.\$replicate_col)`. `comp_or_factory` must be an
  `AbstractLatentComponent` instance.
- `group_col::Symbol`: returns `Group(factory, data.\$group_col)`.
  `comp_or_factory` is the factory (e.g. `IID`, `AR1`); the
  per-group inner components are constructed by the LGM core
  `Group(factory, group_id)` constructor — see
  `packages/LatentGaussianModels.jl/src/components/group.jl`.
- both `nothing`: returns `comp_or_factory` unchanged (must be a
  component instance — caller's responsibility).
"""
function _wrap_term(comp_or_factory, data, replicate_col, group_col)
    replicate_col === nothing || group_col === nothing ||
        throw(ArgumentError("@lgm: `f(...)` cannot have both `replicate` and `group` set."))
    if replicate_col === nothing && group_col === nothing
        comp_or_factory isa LatentGaussianModels.AbstractLatentComponent ||
            throw(ArgumentError("@lgm: f-term: second argument must be an `AbstractLatentComponent`, got $(typeof(comp_or_factory))"))
        return comp_or_factory
    elseif replicate_col !== nothing
        comp_or_factory isa LatentGaussianModels.AbstractLatentComponent ||
            throw(ArgumentError("@lgm: f-term with `replicate = $(replicate_col)`: second argument must be an `AbstractLatentComponent` instance (e.g. `AR1(n)`), got $(typeof(comp_or_factory))"))
        rep_col = _checked_index_column(data, replicate_col)
        R = maximum(rep_col)
        R ≥ 1 ||
            throw(ArgumentError("@lgm: replicate column `$(replicate_col)` must contain integers ≥ 1; got max $(R)"))
        return LatentGaussianModels.Replicate(comp_or_factory, R)
    else
        grp_col = _checked_index_column(data, group_col)
        return LatentGaussianModels.Group(comp_or_factory, grp_col)
    end
end

function _checked_index_column(data, col_name::Symbol)
    Tables.istable(data) ||
        throw(ArgumentError("@lgm: `data` is not a Tables.jl-compatible source (got $(typeof(data)))"))
    cols = Tables.columns(data)
    names = Tuple(Tables.columnnames(cols))
    col_name in names ||
        throw(ArgumentError("@lgm: column `$(col_name)` not found in `data`. Available columns: $(names)"))
    raw = Tables.getcolumn(cols, col_name)
    out = Vector{Int}(undef, length(raw))
    @inbounds for i in eachindex(raw)
        v = raw[i]
        v isa Integer ||
            throw(ArgumentError("@lgm: column `$(col_name)` must contain integers; got entry of type $(typeof(v)) at row $(i)"))
        v ≥ 1 ||
            throw(ArgumentError("@lgm: column `$(col_name)` must contain integers ≥ 1; got $(v) at row $(i)"))
        out[i] = Int(v)
    end
    return out
end

# Normalise heterogeneous `randoms` inputs to 4-tuples
# `(col, comp, replicate_or_nothing, group_or_nothing)`. Accepts both the
# new 4-tuple/NamedTuple form and the PR-1..PR-4 2-tuple form for
# `lgmformula` backward-compat.
function _normalize_term(t)
    if t isa NamedTuple
        col = t.col
        comp = t.comp_expr
        rep = haskey(t, :replicate) ? t.replicate : nothing
        grp = haskey(t, :group) ? t.group : nothing
        return (col, comp, rep, grp)
    elseif t isa Tuple && length(t) == 2
        return (t[1], t[2], nothing, nothing)
    elseif t isa Tuple && length(t) == 4
        return (t[1], t[2], t[3], t[4])
    else
        throw(ArgumentError("@lgm: f-term must be `(col, comp)` or `(col, comp, replicate, group)`, got $(typeof(t))"))
    end
end

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
        col_name, _, rep, grp = _normalize_term(term)
        if col_name isa Symbol
            col_name in names ||
                throw(ArgumentError("@lgm: f-term column `$(col_name)` not found in `data`. Available columns: $(names)"))
        elseif col_name isa Tuple
            for c in col_name
                c isa Symbol ||
                    throw(ArgumentError("@lgm: f-term tuple-coordinate entry must be a Symbol, got $(typeof(c))"))
                c in names ||
                    throw(ArgumentError("@lgm: f-term tuple-coordinate column `$(c)` not found in `data`. Available columns: $(names)"))
            end
        else
            throw(ArgumentError("@lgm: f-term column must be a Symbol or tuple of Symbols, got $(typeof(col_name))"))
        end
        if rep !== nothing
            rep in names ||
                throw(ArgumentError("@lgm: f-term `replicate = $(rep)` column not found in `data`. Available columns: $(names)"))
        end
        if grp !== nothing
            grp in names ||
                throw(ArgumentError("@lgm: f-term `group = $(grp)` column not found in `data`. Available columns: $(names)"))
        end
    end
    return nothing
end

"""
    _build_design_matrix(data, lhs, has_intercept, covariates, randoms)
        -> SparseMatrixCSC{Float64, Int}

Assemble the linear projector matrix `A` for the formula. Columns are
ordered as `[intercept | covariates | random-effect indicators...]`.

Each random-effect term contributes one block:

- Plain `f(col, Comp)`: `length(comp)` columns; row i is 1 in column
  `idx[i]` (within the block). Input column must contain integers in
  `1:length(comp)`.
- `f(col, Comp; replicate = id_col)`: `R · length(comp)` columns,
  laid out `[x⁽¹⁾; x⁽²⁾; …; x⁽ᴿ⁾]`. Row i is 1 in column
  `(id[i] - 1) · length(comp) + col[i]`.
- `f(col, Factory; group = grp_col)`: `Σ_g s_g` columns where `s_g`
  is the per-group size (number of rows with `grp == g`). Row i is
  1 in column `offset[grp[i]] + col[i]`. `col[i]` is required to lie
  in `1:s_{grp[i]}`.
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
        col_name, comp_or_factory, rep_col, grp_col = _normalize_term(term)
        if col_name isa Tuple
            rep_col === nothing && grp_col === nothing ||
                throw(ArgumentError("@lgm: tuple-coordinate f-term `$(col_name)` cannot also have `replicate` or `group`"))
            push!(blocks, _build_spatial_block(comp_or_factory, cols, col_name, n_obs))
        else
            col = Tables.getcolumn(cols, col_name)
            length(col) == n_obs ||
                throw(DimensionMismatch("@lgm: f-term column `$(col_name)` has length $(length(col)); outcome `$(lhs)` has length $(n_obs)"))
            push!(blocks,
                _build_term_block(comp_or_factory, col, col_name,
                    rep_col, grp_col, cols, n_obs))
        end
    end

    isempty(blocks) &&
        throw(ArgumentError("@lgm: model has no terms — formula must include at least `1`, a covariate, or `f(...)`"))
    return length(blocks) == 1 ? blocks[1] : reduce(hcat, blocks)
end

function _build_term_block(comp_or_factory, col, col_name::Symbol,
        rep_col::Union{Symbol, Nothing}, grp_col::Union{Symbol, Nothing},
        cols, n_obs::Int)
    if rep_col === nothing && grp_col === nothing
        comp_or_factory isa LatentGaussianModels.AbstractLatentComponent ||
            throw(ArgumentError("@lgm: f-term `$(col_name)`: second argument must be an `AbstractLatentComponent`, got $(typeof(comp_or_factory))"))
        n_levels = length(comp_or_factory)
        idx = _validate_index_column(col, col_name, n_levels)
        return sparse(1:n_obs, idx, 1.0, n_obs, n_levels)
    elseif rep_col !== nothing
        comp_or_factory isa LatentGaussianModels.AbstractLatentComponent ||
            throw(ArgumentError("@lgm: f-term `$(col_name)` with `replicate = $(rep_col)`: second argument must be an `AbstractLatentComponent` instance, got $(typeof(comp_or_factory))"))
        inner_n = length(comp_or_factory)
        idx = _validate_index_column(col, col_name, inner_n)
        rep_raw = Tables.getcolumn(cols, rep_col)
        length(rep_raw) == n_obs ||
            throw(DimensionMismatch("@lgm: replicate column `$(rep_col)` has length $(length(rep_raw)); expected $(n_obs)"))
        R = _max_index(rep_raw, rep_col)
        rep_idx = _validate_index_column(rep_raw, rep_col, R)
        block_cols = (rep_idx .- 1) .* inner_n .+ idx
        return sparse(1:n_obs, block_cols, 1.0, n_obs, R * inner_n)
    else
        grp_raw = Tables.getcolumn(cols, grp_col)
        length(grp_raw) == n_obs ||
            throw(DimensionMismatch("@lgm: group column `$(grp_col)` has length $(length(grp_raw)); expected $(n_obs)"))
        G = _max_index(grp_raw, grp_col)
        grp_idx = _validate_index_column(grp_raw, grp_col, G)
        sizes = zeros(Int, G)
        @inbounds for g in grp_idx
            sizes[g] += 1
        end
        all(>(0), sizes) ||
            throw(ArgumentError("@lgm: group column `$(grp_col)` must have every label in 1:$(G) present; missing groups: $(findall(==(0), sizes))"))
        offsets = [0; cumsum(sizes)[1:(end - 1)]]
        idx = _validate_within_group_column(col, col_name, grp_idx, sizes)
        block_cols = offsets[grp_idx] .+ idx
        return sparse(1:n_obs, block_cols, 1.0, n_obs, sum(sizes))
    end
end

function _max_index(raw, col_name::Symbol)
    isempty(raw) &&
        throw(ArgumentError("@lgm: column `$(col_name)` is empty"))
    m = 0
    @inbounds for i in eachindex(raw)
        v = raw[i]
        v isa Integer ||
            throw(ArgumentError("@lgm: column `$(col_name)` must contain integers; got entry of type $(typeof(v)) at row $(i)"))
        v ≥ 1 ||
            throw(ArgumentError("@lgm: column `$(col_name)` must contain integers ≥ 1; got $(v) at row $(i)"))
        v > m && (m = Int(v))
    end
    return m
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

function _validate_within_group_column(col, col_name::Symbol,
        grp_idx::AbstractVector{Int}, sizes::AbstractVector{Int})
    idx = Vector{Int}(undef, length(col))
    @inbounds for i in eachindex(col)
        v = col[i]
        v isa Integer ||
            throw(ArgumentError("@lgm: f-term column `$(col_name)` must contain integers; got entry of type $(typeof(v)) at row $(i)"))
        s = sizes[grp_idx[i]]
        (1 ≤ v ≤ s) ||
            throw(ArgumentError("@lgm: f-term column `$(col_name)` value $(v) at row $(i) is outside 1:$(s) for its group (group size = $(s))"))
        idx[i] = Int(v)
    end
    return idx
end

"""
    _build_multi_likelihood_mapping(data, lhs, has_intercept, covariates, randoms)
        -> StackedMapping

Build a row-partitioned `StackedMapping` for tuple-LHS multi-likelihood
models. The shared RHS produces a single sparse `A`; each likelihood
block wraps the same `LinearProjector(A)`. Observation rows are
partitioned contiguously: block `k` owns rows `((k-1)·n + 1):(k·n)`.

All columns in `lhs` must have equal length `n` (wide-format only —
long-format with a `type` column is left for a follow-up).
"""
function _build_multi_likelihood_mapping(data, lhs::AbstractVector{Symbol},
        has_intercept::Bool, covariates::AbstractVector,
        randoms::AbstractVector)
    Tables.istable(data) ||
        throw(ArgumentError("@lgm: `data` is not a Tables.jl-compatible source (got $(typeof(data)))"))
    cols = Tables.columns(data)
    names = Tuple(Tables.columnnames(cols))
    n = nothing
    for name in lhs
        name in names ||
            throw(ArgumentError("@lgm: outcome column `$(name)` not found in `data`. Available columns: $(names)"))
        col_len = length(Tables.getcolumn(cols, name))
        if n === nothing
            n = col_len
        elseif col_len != n
            throw(DimensionMismatch("@lgm: tuple-LHS columns must have equal length; column `$(name)` has length $(col_len), but `$(first(lhs))` has length $(n)"))
        end
    end
    A = _build_design_matrix(data, first(lhs), has_intercept, covariates, randoms)
    proj = LatentGaussianModels.LinearProjector(A)
    k = length(lhs)
    blocks = ntuple(_ -> proj, k)
    rows = [((i - 1) * n + 1):(i * n) for i in 1:k]
    return LatentGaussianModels.StackedMapping(blocks, rows)
end

"""
    lgmformula(data; lhs, intercept = true, covariates = Symbol[],
               randoms = [], family) -> LatentGaussianModel

Function form of [`@lgm`](@ref). Accepts a structured description of
the formula and returns the same `LatentGaussianModel` the macro
would produce.

# Arguments

- `data`: a `Tables.jl`-compatible source.
- `lhs::Union{Symbol, AbstractVector{Symbol}}`: outcome column name(s).
  A vector triggers tuple-LHS multi-likelihood; `family` must then be
  a tuple of likelihoods of matching length.
- `intercept::Bool = true`: whether to include `Intercept()`.
- `covariates::Vector{Symbol} = Symbol[]`: scalar fixed-effect column
  names. Becomes `FixedEffects(length(covariates))` if non-empty.
- `randoms::AbstractVector = []`: list of f-term specifications. Each
  entry may be:
  - `(col::Symbol, comp::AbstractLatentComponent)` — plain f-term.
  - `(col, comp, replicate::Symbol, group::Nothing)` — replicated
    component; runtime wraps as `Replicate(comp, R)`.
  - `(col, factory, replicate::Nothing, group::Symbol)` — grouped
    component; runtime wraps as `Group(factory, grp_col_values)`.
  - A `NamedTuple{(:col, :comp_expr, :replicate, :group)}` — internal
    form emitted by the macro.
- `family`: observation likelihood (single-LHS) or tuple of
  likelihoods (multi-LHS).
"""
function lgmformula(data;
        lhs::Union{Symbol, AbstractVector{Symbol}},
        intercept::Bool=true,
        covariates::Vector{Symbol}=Symbol[],
        randoms::AbstractVector=Tuple{
            Symbol, LatentGaussianModels.AbstractLatentComponent}[],
        family)
    components = _build_components(intercept, length(covariates), randoms, data)
    if lhs isa Symbol
        family isa LatentGaussianModels.AbstractLikelihood ||
            throw(ArgumentError("@lgm: single-LHS `lgmformula` requires `family::AbstractLikelihood`, got $(typeof(family))"))
        A = _build_design_matrix(data, lhs, intercept, covariates, randoms)
        return LatentGaussianModels.LatentGaussianModel(family, components, A)
    else
        family isa Tuple ||
            throw(ArgumentError("@lgm: tuple-LHS `lgmformula` requires `family::Tuple` of likelihoods, got $(typeof(family))"))
        length(lhs) == length(family) ||
            throw(ArgumentError("@lgm: tuple-LHS has $(length(lhs)) columns but `family` has $(length(family)) likelihoods — must match"))
        all(ℓ -> ℓ isa LatentGaussianModels.AbstractLikelihood, family) ||
            throw(ArgumentError("@lgm: every entry of `family` must be an `AbstractLikelihood`; got $(map(typeof, family))"))
        mapping = _build_multi_likelihood_mapping(data, lhs, intercept, covariates, randoms)
        return LatentGaussianModels.LatentGaussianModel(family, components, mapping)
    end
end

function _build_components(has_intercept::Bool, n_covariates::Int,
        randoms::AbstractVector, data)
    parts = LatentGaussianModels.AbstractLatentComponent[]
    has_intercept && push!(parts, LatentGaussianModels.Intercept())
    n_covariates > 0 &&
        push!(parts, LatentGaussianModels.FixedEffects(n_covariates))
    for term in randoms
        col, comp_or_factory, rep, grp = _normalize_term(term)
        if col isa Tuple
            comp_or_factory isa LatentGaussianModels.AbstractLatentComponent ||
                throw(ArgumentError("@lgm: f-term `$(col)`: tuple-coordinate path requires an `AbstractLatentComponent` instance, got $(typeof(comp_or_factory))"))
            push!(parts, comp_or_factory)
        else
            push!(parts, _wrap_term(comp_or_factory, data, rep, grp))
        end
    end
    isempty(parts) &&
        throw(ArgumentError("@lgm: model has no components — formula must include at least `1`, a covariate, or `f(...)`"))
    return Tuple(parts)
end

"""
    _build_spatial_block(component, data_cols, coord_cols::Tuple, n_obs::Int)
        -> SparseMatrixCSC{Float64, Int}

Build the design-matrix block for a tuple-coordinate `f((cols...),
component)` term. The default method throws — concrete implementations
live in package extensions. `LGMFormulaINLASPDEExt` overloads this for
`SPDE2` to build a barycentric `MeshProjector`.
"""
function _build_spatial_block(component, data_cols, coord_cols, n_obs)
    throw(ArgumentError(
        "@lgm: component $(typeof(component)) does not support " *
        "tuple-coordinate `f((cols...), Component)` syntax. Spatial " *
        "SPDE components (e.g. `SPDE2`) require `using INLASPDE` to " *
        "load the `LGMFormulaINLASPDEExt` extension."
    ))
end
