# Expression-tree parsing for `@lgm`. PR-1 scope: walk a formula `lhs ~
# rhs`, separate the RHS into intercept marker / fixed-effects column
# symbols / `f(col, Component)` random-effect calls. No StatsModels
# integration here yet — that's PR-3 (`terms.jl`).

"""
    _parse_args(args) -> (formula_expr, opts::Dict{Symbol,Any})

Split macro args into the formula expression and `key = value` options.
Accepts both bare-form (`@lgm y ~ 1 data=df`) and parenthesised-form
(`@lgm(y ~ 1, data=df)`).
"""
function _parse_args(args)
    isempty(args) &&
        error("@lgm: no arguments — expected a formula `lhs ~ rhs` and `data=...`, `family=...`")
    formula_expr = nothing
    opts = Dict{Symbol, Any}()
    for arg in args
        if arg isa Expr && arg.head === :parameters
            for kw in arg.args
                _record_kw!(opts, kw)
            end
        elseif arg isa Expr && (arg.head === :(=) || arg.head === :kw)
            _record_kw!(opts, arg)
        elseif formula_expr === nothing
            formula_expr = arg
        else
            error("@lgm: unexpected argument `$arg`. Expected `formula, data=df, family=Likelihood()`.")
        end
    end
    formula_expr === nothing &&
        error("@lgm: missing formula `lhs ~ rhs`")
    return formula_expr, opts
end

function _record_kw!(opts::Dict{Symbol, Any}, kw)
    if kw isa Expr && (kw.head === :(=) || kw.head === :kw) && kw.args[1] isa Symbol
        opts[kw.args[1]::Symbol] = kw.args[2]
    else
        error("@lgm: malformed keyword argument `$kw`")
    end
    return opts
end

"""
    _parse_formula(expr) -> (lhs::Symbol, rhs)

Split `lhs ~ rhs` into LHS symbol and RHS expression.

PR-1 restricts the LHS to a single `Symbol`; tuple-LHS for joint
likelihoods ships in PR-4 (ADR-033).
"""
function _parse_formula(expr)
    if !(expr isa Expr && expr.head === :call && length(expr.args) == 3 &&
         expr.args[1] === :~)
        error("@lgm: expected a formula `lhs ~ rhs`, got `$expr`")
    end
    lhs = expr.args[2]
    lhs isa Symbol ||
        error("@lgm: LHS must be a column name (Symbol). Joint-likelihood tuple-LHS is not yet supported (PR-4).")
    rhs = expr.args[3]
    return lhs, rhs
end

"""
    _split_rhs(rhs) -> (has_intercept::Bool,
                       covariates::Vector{Symbol},
                       randoms::Vector{Tuple{Symbol, Any}})

Walk the RHS, splitting at `+`. Each summand is one of:

- `1` — explicit intercept marker (default if no marker present).
- `0` or `-1` — explicit "no intercept".
- bare `Symbol` — fixed-effects covariate column.
- `f(col, Component(...))` — random-effect term; `col` is a column
  symbol, the second arg is the (un-evaluated) component expression.

Other forms (transformations, interactions, etc.) raise an error
referring to user concepts.
"""
function _split_rhs(rhs)
    has_intercept = true
    covariates = Symbol[]
    randoms = Tuple{Symbol, Any}[]
    for s in _flatten_plus(rhs)
        if s === 1
            has_intercept = true
        elseif s === 0 || s === -1 || _is_neg_one(s)
            has_intercept = false
        elseif s isa Symbol
            push!(covariates, s)
        elseif s isa Expr && s.head === :call && s.args[1] === :f
            length(s.args) == 3 ||
                error("@lgm: malformed `f(...)` term `$s`. Expected `f(column, Component(...))`.")
            col = s.args[2]
            col isa Symbol ||
                error("@lgm: `f(...)`: first argument must be a column name (Symbol), got `$col`.")
            comp_expr = s.args[3]
            push!(randoms, (col, comp_expr))
        else
            error("@lgm: unsupported RHS term `$s`. PR-1 supports `1`, `0`, `-1`, bare column symbols, and `f(col, Component(...))` only.")
        end
    end
    return has_intercept, covariates, randoms
end

"""
    _flatten_plus(expr) -> Vector{Any}

Flatten an `a + b + c` chain into `[a, b, c]`. Anything else returns
as a single-element list.
"""
function _flatten_plus(expr)
    if expr isa Expr && expr.head === :call && expr.args[1] === :+
        result = Any[]
        for arg in expr.args[2:end]
            append!(result, _flatten_plus(arg))
        end
        return result
    else
        return Any[expr]
    end
end

_is_neg_one(s) = s isa Expr && s.head === :call && length(s.args) == 2 &&
                 s.args[1] === :- && s.args[2] === 1
