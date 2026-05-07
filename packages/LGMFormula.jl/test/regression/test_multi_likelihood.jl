# PR-4 multi-likelihood tuple-LHS coverage. `(y1, y2) ~ rhs` lowers to a
# `LatentGaussianModel((ℓ1, ℓ2), components, StackedMapping(...))` whose
# k blocks each wrap the *same* `LinearProjector(A)` and partition the
# observation rows contiguously: block k owns rows `((k-1)·n+1):(k·n)`.
# Asserts macro / function / hand-written forms produce identical models
# via `_struct_isequal`. Long-format with a `type` column (per-row
# likelihood) is left for a follow-up PR.

@testset "PR-4 multi-likelihood tuple-LHS roundtrip" begin
    rng = Random.Xoshiro(20260507)
    n = 40
    df = (
        y_real=randn(rng, n),
        y_count=rand(rng, 0:5, n),
        y_real2=randn(rng, n),
        x=randn(rng, n),
        idx=rand(rng, 1:5, n)
    )

    @testset "bivariate Gaussian + Intercept" begin
        model_macro = @lgm (y_real, y_real2)~1 data=df family=(
            GaussianLikelihood(), GaussianLikelihood())
        A = sparse(ones(n, 1))
        proj = LinearProjector(A)
        mapping = StackedMapping((proj, proj), [1:n, (n + 1):(2n)])
        model_hand = LatentGaussianModel(
            (GaussianLikelihood(), GaussianLikelihood()),
            (Intercept(),),
            mapping
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Gaussian + Poisson + Intercept + IID" begin
        model_macro = @lgm (y_real, y_count)~1 + f(idx, IID(5)) data=df family=(
            GaussianLikelihood(), PoissonLikelihood())
        A_iid = sparse(1:n, df.idx, 1.0, n, 5)
        A = hcat(sparse(ones(n, 1)), A_iid)
        proj = LinearProjector(A)
        mapping = StackedMapping((proj, proj), [1:n, (n + 1):(2n)])
        model_hand = LatentGaussianModel(
            (GaussianLikelihood(), PoissonLikelihood()),
            (Intercept(), IID(5)),
            mapping
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "three-likelihood: Gaussian + Poisson + Gaussian" begin
        model_macro = @lgm (y_real, y_count, y_real2)~1 + x data=df family=(
            GaussianLikelihood(), PoissonLikelihood(), GaussianLikelihood())
        X = sparse(reshape(Float64.(df.x), n, 1))
        A = hcat(sparse(ones(n, 1)), X)
        proj = LinearProjector(A)
        mapping = StackedMapping(
            (proj, proj, proj),
            [1:n, (n + 1):(2n), (2n + 1):(3n)]
        )
        model_hand = LatentGaussianModel(
            (GaussianLikelihood(), PoissonLikelihood(), GaussianLikelihood()),
            (Intercept(), FixedEffects(1)),
            mapping
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "no intercept + 1 covariate + IID" begin
        model_macro = @lgm (y_real, y_count)~0 + x + f(idx, IID(5)) data=df family=(
            GaussianLikelihood(), PoissonLikelihood())
        X = sparse(reshape(Float64.(df.x), n, 1))
        A_iid = sparse(1:n, df.idx, 1.0, n, 5)
        A = hcat(X, A_iid)
        proj = LinearProjector(A)
        mapping = StackedMapping((proj, proj), [1:n, (n + 1):(2n)])
        model_hand = LatentGaussianModel(
            (GaussianLikelihood(), PoissonLikelihood()),
            (FixedEffects(1), IID(5)),
            mapping
        )
        @test _struct_isequal(model_macro, model_hand)
    end
end

@testset "PR-4 lgmformula function form (tuple-LHS)" begin
    rng = Random.Xoshiro(20260507)
    n = 30
    df = (
        y1=randn(rng, n),
        y2=rand(rng, 0:5, n),
        x=randn(rng, n),
        idx=rand(rng, 1:4, n)
    )

    @testset "function form matches macro form" begin
        model_macro = @lgm (y1, y2)~1 + x + f(idx, IID(4)) data=df family=(
            GaussianLikelihood(), PoissonLikelihood())
        model_func = lgmformula(df;
            lhs=[:y1, :y2],
            intercept=true,
            covariates=[:x],
            randoms=[(:idx, IID(4))],
            family=(GaussianLikelihood(), PoissonLikelihood()))
        @test _struct_isequal(model_macro, model_func)
    end

    @testset "function form: bivariate Gaussian intercept-only" begin
        model_func = lgmformula(df;
            lhs=[:y1, :y2],
            intercept=true,
            family=(GaussianLikelihood(), GaussianLikelihood()))
        A = sparse(ones(n, 1))
        proj = LinearProjector(A)
        mapping = StackedMapping((proj, proj), [1:n, (n + 1):(2n)])
        model_hand = LatentGaussianModel(
            (GaussianLikelihood(), GaussianLikelihood()),
            (Intercept(),),
            mapping
        )
        @test _struct_isequal(model_func, model_hand)
    end
end

@testset "PR-4 @macroexpand structure (tuple-LHS)" begin
    df_dummy = (y1=Float64[], y2=Int[], x=Float64[], idx=Int[])

    @testset "expansion uses _build_multi_likelihood_mapping" begin
        ex = @macroexpand @lgm (y1, y2)~1 data=df_dummy family=(
            GaussianLikelihood(), PoissonLikelihood())
        call = _find_lgm_call(ex)
        @test call !== nothing
        @test length(call.args) == 4
        mapping = call.args[4]
        @test _is_call_to(mapping, :_build_multi_likelihood_mapping)
    end

    @testset "family appears literally as a tuple" begin
        ex = @macroexpand @lgm (y1, y2)~1 data=df_dummy family=(
            GaussianLikelihood(), PoissonLikelihood())
        call = _find_lgm_call(ex)
        family = call.args[2]
        @test family isa Expr && family.head === :tuple
        @test length(family.args) == 2
        @test _is_call_to(family.args[1], :GaussianLikelihood)
        @test _is_call_to(family.args[2], :PoissonLikelihood)
    end

    @testset "components tuple is shared (not duplicated per likelihood)" begin
        ex = @macroexpand @lgm (y1, y2)~1 + x + f(idx, IID(4)) data=df_dummy family=(
            GaussianLikelihood(), PoissonLikelihood())
        call = _find_lgm_call(ex)
        components = call.args[3]
        @test components.head === :tuple
        @test length(components.args) == 3
        @test _is_call_to(components.args[1], :Intercept)
        @test _is_call_to(components.args[2], :FixedEffects)
        @test _is_call_to(components.args[3], :IID)
    end
end

@testset "PR-4 error messages refer to user concepts (tuple-LHS)" begin
    rng = Random.Xoshiro(20260507)
    n = 20
    df = (
        y1=randn(rng, n),
        y2=rand(rng, 0:5, n),
        x=randn(rng, n),
        idx=rand(rng, 1:3, n)
    )

    @testset "unknown LHS column in tuple" begin
        @test_throws ArgumentError begin
            @lgm (y1, missing_y)~1 + x data=df family=(
                GaussianLikelihood(), GaussianLikelihood())
        end
        try
            @lgm (y1, missing_y)~1 + x data=df family=(
                GaussianLikelihood(), GaussianLikelihood())
        catch e
            msg = sprint(showerror, e)
            @test occursin("missing_y", msg)
        end
    end

    @testset "tuple-LHS with single likelihood (length mismatch)" begin
        # Macro accepts the formula; runtime catches the family-length mismatch.
        @test_throws ArgumentError begin
            @lgm (y1, y2)~1 + x data=df family=(GaussianLikelihood(),)
        end
        try
            @lgm (y1, y2)~1 + x data=df family=(GaussianLikelihood(),)
        catch e
            msg = sprint(showerror, e)
            @test occursin("family", msg) || occursin("likelihood", msg)
        end
    end

    @testset "tuple-LHS with non-tuple family" begin
        # `lgmformula` rejects single likelihood for tuple-LHS.
        @test_throws ArgumentError lgmformula(df;
            lhs=[:y1, :y2],
            intercept=true,
            family=GaussianLikelihood())
    end

    @testset "single-LHS with tuple family (lgmformula)" begin
        @test_throws ArgumentError lgmformula(df;
            lhs=:y1,
            intercept=true,
            family=(GaussianLikelihood(), GaussianLikelihood()))
    end

    @testset "tuple-LHS columns of unequal length" begin
        # Construct a NamedTuple by hand with mismatched column lengths.
        bad_df = (y1=randn(rng, 10), y2=randn(rng, 12), x=randn(rng, 10))
        @test_throws DimensionMismatch begin
            @lgm (y1, y2)~1 data=bad_df family=(GaussianLikelihood(), GaussianLikelihood())
        end
    end

    @testset "tuple-LHS with zero columns" begin
        # `()` parses as Expr(:tuple) with empty args; rejected by the parser.
        ex = try
            @eval @lgm ()~1 data=$df family=(GaussianLikelihood(),)
            nothing
        catch e
            e
        end
        @test ex isa LoadError || ex isa ErrorException
    end

    @testset "tuple-LHS entry is non-symbol" begin
        ex = try
            @eval @lgm ($("y1"), y2)~1 data=$df family=(
                GaussianLikelihood(), GaussianLikelihood())
            nothing
        catch e
            e
        end
        @test ex isa LoadError || ex isa ErrorException
    end
end
