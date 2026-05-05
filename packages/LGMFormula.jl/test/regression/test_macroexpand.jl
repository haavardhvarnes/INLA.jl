# Structural assertions on `@macroexpand @lgm(...)` output. Verify that
# the expansion is an explicit `LatentGaussianModel(...)` constructor
# call (not a wrapped `lgmformula(...)`), with the components tuple and
# likelihood appearing literally — only the design matrix construction
# is deferred to runtime.

# Walk the expanded AST and pull out the `LatentGaussianModel(...)`
# call. Robust to surrounding `let`-blocks, line-number nodes, and
# module-qualified paths.
function _find_lgm_call(expr)
    if expr isa Expr && expr.head === :call && _is_lgm_constructor(expr.args[1])
        return expr
    elseif expr isa Expr
        for a in expr.args
            r = _find_lgm_call(a)
            r === nothing || return r
        end
    end
    return nothing
end

function _is_lgm_constructor(callee)
    if callee isa Symbol
        return callee === :LatentGaussianModel
    elseif callee isa GlobalRef
        return callee.name === :LatentGaussianModel
    elseif callee isa Expr && callee.head === :.
        return _qualified_name_tail(callee) === :LatentGaussianModel
    else
        return false
    end
end

function _qualified_name_tail(expr::Expr)
    # `Mod.Sub.Name` -> :Name
    if expr.head === :.
        last = expr.args[end]
        if last isa QuoteNode
            return last.value
        elseif last isa Symbol
            return last
        end
    end
    return nothing
end

function _is_call_to(expr, name::Symbol)
    expr isa Expr && expr.head === :call || return false
    callee = expr.args[1]
    if callee isa Symbol
        return callee === name
    elseif callee isa Expr && callee.head === :.
        return _qualified_name_tail(callee) === name
    elseif callee isa GlobalRef
        return callee.name === name
    end
    return false
end

@testset "PR-1 @macroexpand structure" begin
    df_dummy = (y = Float64[], x = Float64[], x1 = Float64[],
        x2 = Float64[], idx = Int[])

    @testset "expansion is a LatentGaussianModel(...) call" begin
        ex = @macroexpand @lgm y ~ 1 data=df_dummy family=GaussianLikelihood()
        call = _find_lgm_call(ex)
        @test call !== nothing
        @test length(call.args) == 4   # callee, family, components, mapping
    end

    @testset "components tuple appears literally — Intercept only" begin
        ex = @macroexpand @lgm y ~ 1 data=df_dummy family=GaussianLikelihood()
        call = _find_lgm_call(ex)
        components = call.args[3]
        @test components isa Expr && components.head === :tuple
        @test length(components.args) == 1
        @test _is_call_to(components.args[1], :Intercept)
    end

    @testset "components tuple — Intercept + FixedEffects(2)" begin
        ex = @macroexpand @lgm y ~ 1 + x1 + x2 data=df_dummy family=GaussianLikelihood()
        call = _find_lgm_call(ex)
        components = call.args[3]
        @test components.head === :tuple
        @test length(components.args) == 2
        @test _is_call_to(components.args[1], :Intercept)
        @test _is_call_to(components.args[2], :FixedEffects)
        # FixedEffects argument is the literal covariate count (2).
        @test components.args[2].args[2] == 2
    end

    @testset "components tuple — Intercept + IID(5)" begin
        ex = @macroexpand @lgm y ~ 1 + f(idx, IID(5)) data=df_dummy family=PoissonLikelihood()
        call = _find_lgm_call(ex)
        components = call.args[3]
        @test components.head === :tuple
        @test length(components.args) == 2
        @test _is_call_to(components.args[1], :Intercept)
        @test _is_call_to(components.args[2], :IID)
        @test components.args[2].args[2] == 5
    end

    @testset "no intercept — components tuple — FixedEffects only" begin
        ex = @macroexpand @lgm y ~ 0 + x1 data=df_dummy family=GaussianLikelihood()
        call = _find_lgm_call(ex)
        components = call.args[3]
        @test components.head === :tuple
        @test length(components.args) == 1
        @test _is_call_to(components.args[1], :FixedEffects)
    end

    @testset "design-matrix call is _build_design_matrix" begin
        ex = @macroexpand @lgm y ~ 1 + x1 + f(idx, IID(5)) data=df_dummy family=PoissonLikelihood()
        call = _find_lgm_call(ex)
        mapping = call.args[4]
        @test _is_call_to(mapping, :_build_design_matrix)
    end

    @testset "family appears literally as the second argument" begin
        ex = @macroexpand @lgm y ~ 1 data=df_dummy family=PoissonLikelihood()
        call = _find_lgm_call(ex)
        family = call.args[2]
        @test _is_call_to(family, :PoissonLikelihood)
    end
end
