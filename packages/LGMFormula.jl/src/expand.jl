# Build the AST returned by `@lgm`. The expansion is an explicit
# `LatentGaussianModel(...)` constructor call: components and
# likelihood appear literally; only the design matrix / mapping
# construction is deferred to runtime (it depends on the data).

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
"""
function _build_expansion(lhs::Vector{Symbol}, has_intercept::Bool,
        covariates::Vector{Symbol}, randoms::Vector{Tuple{Symbol, Any}},
        data_expr, family_expr)
    components_expr = _components_tuple_expr(has_intercept, covariates, randoms)
    randoms_vec_expr = _randoms_vec_expr(randoms)
    covariates_vec_expr = _covariates_vec_expr(covariates)

    if length(lhs) == 1
        mapping_expr = :(
            $LGMFormula._build_design_matrix(
                $data_expr,
                $(QuoteNode(first(lhs))),
                $has_intercept,
                $covariates_vec_expr,
                $randoms_vec_expr,
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
                $randoms_vec_expr,
            )
        )
    end

    return :(
        $LatentGaussianModels.LatentGaussianModel(
            $family_expr,
            $components_expr,
            $mapping_expr,
        )
    )
end

function _components_tuple_expr(has_intercept::Bool,
        covariates::Vector{Symbol},
        randoms::Vector{Tuple{Symbol, Any}})
    parts = Any[]
    has_intercept &&
        push!(parts, :($LatentGaussianModels.Intercept()))
    !isempty(covariates) &&
        push!(parts, :($LatentGaussianModels.FixedEffects($(length(covariates)))))
    for (_, comp_expr) in randoms
        push!(parts, comp_expr)
    end
    isempty(parts) &&
        error("@lgm: model has no components — formula must include at least `1`, a covariate, or `f(...)`")
    return Expr(:tuple, parts...)
end

function _randoms_vec_expr(randoms::Vector{Tuple{Symbol, Any}})
    pair_exprs = [Expr(:tuple, QuoteNode(col), comp_expr)
                  for (col, comp_expr) in randoms]
    return Expr(:vect, pair_exprs...)
end

function _covariates_vec_expr(covariates::Vector{Symbol})
    return Expr(:vect, [QuoteNode(c) for c in covariates]...)
end
