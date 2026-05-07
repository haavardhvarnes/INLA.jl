# Build the AST returned by `@lgm`. The expansion is an explicit
# `LatentGaussianModel(...)` constructor call: components and
# likelihood appear literally; only the design matrix / mapping
# construction (and `Replicate`/`Group` wrap calls per PR-5) is
# deferred to runtime, since those depend on the data table.

"""
    _build_expansion(lhs::Vector{Symbol}, has_intercept, covariates,
                     randoms, data_expr, family_expr) -> Expr

Return the `Expr` that `@lgm` expands to. Module references use
absolute interpolation (`\$LatentGaussianModels.Intercept()`,
`\$LGMFormula._build_design_matrix(...)`) so the expansion resolves
regardless of the caller's `using` imports.

Single-LHS expands to a `LinearProjector`-equivalent `SparseMatrixCSC`
(LGM auto-wraps). Multi-LHS expands to a `StackedMapping` with one
`LinearProjector(A)` block per likelihood, sharing the RHS-built `A`.

PR-5: `f(col, Comp; replicate = id_col)` and
`f(col, Factory; group = grp_col)` lower to `_wrap_term(...)` runtime
calls in the components tuple — the wrapper resolves `id_col`/`grp_col`
against `data` to construct `Replicate(comp, R)` / `Group(factory,
grp)`. The macro itself does no I/O; the AST is still data-free.
"""
function _build_expansion(lhs::Vector{Symbol}, has_intercept::Bool,
        covariates::Vector{Symbol}, randoms::AbstractVector,
        data_expr, family_expr)
    components_expr = _components_tuple_expr(has_intercept, covariates,
        randoms, data_expr)
    randoms_vec_expr = _randoms_vec_expr(randoms)
    covariates_vec_expr = _covariates_vec_expr(covariates)

    if length(lhs) == 1
        mapping_expr = :(
            $LGMFormula._build_design_matrix(
            $data_expr,
            $(QuoteNode(first(lhs))),
            $has_intercept,
            $covariates_vec_expr,
            $randoms_vec_expr
        )
        )
    else
        lhs_vec_expr = Expr(:vect, [QuoteNode(c) for c in lhs]...)
        mapping_expr = :(
            $LGMFormula._build_multi_likelihood_mapping(
            $data_expr,
            $lhs_vec_expr,
            $has_intercept,
            $covariates_vec_expr,
            $randoms_vec_expr
        )
        )
    end

    return :(
        $LatentGaussianModels.LatentGaussianModel(
        $family_expr,
        $components_expr,
        $mapping_expr
    )
    )
end

function _components_tuple_expr(has_intercept::Bool,
        covariates::Vector{Symbol}, randoms::AbstractVector, data_expr)
    parts = Any[]
    has_intercept &&
        push!(parts, :($LatentGaussianModels.Intercept()))
    !isempty(covariates) &&
        push!(parts, :($LatentGaussianModels.FixedEffects($(length(covariates)))))
    for term in randoms
        push!(parts, _component_expr_for_term(term, data_expr))
    end
    isempty(parts) &&
        error("@lgm: model has no components — formula must include at least `1`, a covariate, or `f(...)`")
    return Expr(:tuple, parts...)
end

function _component_expr_for_term(term::NamedTuple, data_expr)
    if term.replicate === nothing && term.group === nothing
        return term.comp_expr
    end
    rep_arg = term.replicate === nothing ? :nothing : QuoteNode(term.replicate)
    grp_arg = term.group === nothing ? :nothing : QuoteNode(term.group)
    return :(
        $LGMFormula._wrap_term(
        $(term.comp_expr),
        $data_expr,
        $rep_arg,
        $grp_arg
    )
    )
end

function _randoms_vec_expr(randoms::AbstractVector)
    tuple_exprs = [Expr(:tuple,
                       QuoteNode(term.col),
                       term.comp_expr,
                       term.replicate === nothing ? :nothing : QuoteNode(term.replicate),
                       term.group === nothing ? :nothing : QuoteNode(term.group)
                   ) for term in randoms]
    return Expr(:vect, tuple_exprs...)
end

function _covariates_vec_expr(covariates::Vector{Symbol})
    return Expr(:vect, [QuoteNode(c) for c in covariates]...)
end
