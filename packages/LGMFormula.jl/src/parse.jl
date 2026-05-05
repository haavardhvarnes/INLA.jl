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
    _parse_formula(expr) -> (lhs::Vector{Symbol}, rhs)

Split `lhs ~ rhs` into LHS column names and RHS expression. The LHS
is always returned as a `Vector{Symbol}` (length 1 for single-likelihood,
length k > 1 for multi-likelihood tuple-LHS `(y1, y2, ...) ~ rhs`,
ADR-033).
"""
function _parse_formula(expr)
    if !(expr isa Expr && expr.head === :call && length(expr.args) == 3 &&
         expr.args[1] === :~)
        error("@lgm: expected a formula `lhs ~ rhs`, got `$expr`")
    end
    lhs_expr = expr.args[2]
    rhs = expr.args[3]
    lhs = _parse_lhs(lhs_expr)
    return lhs, rhs
end

function _parse_lhs(lhs_expr)
    if lhs_expr isa Symbol
        return Symbol[lhs_expr]
    elseif lhs_expr isa Expr && lhs_expr.head === :tuple
        names = Symbol[]
        for a in lhs_expr.args
            a isa Symbol ||
                error("@lgm: tuple-LHS must contain column names (Symbol), got `$a`")
            push!(names, a)
        end
        length(names) ≥ 1 ||
            error("@lgm: tuple-LHS must contain at least one column name")
        return names
    else
        error("@lgm: LHS must be a column name `y` or a tuple of column names `(y1, y2, ...)`, got `$lhs_expr`")
    end
end

"""
    _split_rhs(rhs) -> (has_intercept::Bool,
                       covariates::Vector{Symbol},
                       randoms::Vector{<:NamedTuple})

Walk the RHS, splitting at `+`. Each summand is one of:

- `1` — explicit intercept marker (default if no marker present).
- `0` or `-1` — explicit "no intercept".
- bare `Symbol` — fixed-effects covariate column.
- `f(col, Component(...))` — random-effect term; `col` is a column
  symbol, the second arg is the (un-evaluated) component expression.
- `f(col, Component(...); replicate = id_col)` —
  R-INLA-style replicated component (PR-5). The macro emits a runtime
  call that wraps the inner component as `Replicate(comp, R)` where
  `R = maximum(id_col)`.
- `f(col, Component; group = grp_col)` — R-INLA-style grouped
  component (PR-5, factory form). The second positional argument is
  the *factory* (a `Symbol` or callable, not an instance); the macro
  emits a runtime `Group(factory, grp_col)` wrap.

Each `f(...)` term lowers to a `NamedTuple{(:col, :comp_expr,
:replicate, :group)}` where `replicate` / `group` carry the keyword-
argument column symbols (or `nothing`).

Other forms (transformations, interactions, etc.) raise an error
referring to user concepts.
"""
function _split_rhs(rhs)
    has_intercept = true
    covariates = Symbol[]
    randoms = NamedTuple{
        (:col, :comp_expr, :replicate, :group),
        Tuple{Symbol, Any, Union{Symbol, Nothing}, Union{Symbol, Nothing}},
    }[]
    for s in _flatten_plus(rhs)
        if s === 1
            has_intercept = true
        elseif s === 0 || s === -1 || _is_neg_one(s)
            has_intercept = false
        elseif s isa Symbol
            push!(covariates, s)
        elseif s isa Expr && s.head === :call && s.args[1] === :f
            push!(randoms, _parse_f_term(s))
        else
            error("@lgm: unsupported RHS term `$s`. Supported: `1`, `0`, `-1`, bare column symbols, and `f(col, Component(...)[; replicate=…|group=…])`.")
        end
    end
    return has_intercept, covariates, randoms
end

"""
    _parse_f_term(s::Expr) -> NamedTuple

Pull the column symbol, component expression, and `replicate`/`group`
keyword arguments out of an `f(...)` call. Accepts both
`f(col, comp; replicate=id)` (semicolon-style) and
`f(col, comp, replicate=id)` (trailing-kw style); rejects unsupported
keywords with a user-visible error.
"""
function _parse_f_term(s::Expr)
    kw_pairs = Pair{Symbol, Any}[]
    positional = Any[]
    for a in s.args[2:end]
        if a isa Expr && a.head === :parameters
            for p in a.args
                _record_f_kw!(kw_pairs, p, s)
            end
        elseif a isa Expr && (a.head === :kw || a.head === :(=)) && a.args[1] isa Symbol
            _record_f_kw!(kw_pairs, a, s)
        else
            push!(positional, a)
        end
    end
    length(positional) == 2 ||
        error("@lgm: malformed `f(...)` term `$s`. Expected `f(column, Component[; replicate=…, group=…])` — got $(length(positional)) positional arg(s).")
    col, comp_expr = positional
    col isa Symbol ||
        error("@lgm: `f(...)`: first argument must be a column name (Symbol), got `$col`.")
    replicate = nothing
    group = nothing
    for (k, v) in kw_pairs
        if k === :replicate
            replicate === nothing ||
                error("@lgm: `f(...)`: `replicate` keyword given more than once.")
            v isa Symbol ||
                error("@lgm: `f(...; replicate = …)` expects a column name (Symbol), got `$v`.")
            replicate = v
        elseif k === :group
            group === nothing ||
                error("@lgm: `f(...)`: `group` keyword given more than once.")
            v isa Symbol ||
                error("@lgm: `f(...; group = …)` expects a column name (Symbol), got `$v`.")
            group = v
        else
            error("@lgm: `f(...)`: unsupported keyword `$k`. Supported keywords: `replicate`, `group`.")
        end
    end
    replicate === nothing || group === nothing ||
        error("@lgm: `f(...)`: `replicate` and `group` are mutually exclusive — got both in `$s`.")
    return (; col = col, comp_expr = comp_expr,
        replicate = replicate, group = group)
end

function _record_f_kw!(kw_pairs::Vector{Pair{Symbol, Any}}, p, s)
    p isa Expr && (p.head === :kw || p.head === :(=)) && p.args[1] isa Symbol ||
        error("@lgm: malformed `f(...)` keyword argument in `$s`.")
    push!(kw_pairs, p.args[1] => p.args[2])
    return kw_pairs
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
