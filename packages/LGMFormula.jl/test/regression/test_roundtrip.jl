# Each `@lgm` formula must produce a `LatentGaussianModel` identical
# to the hand-written Tier-1 form. Equality is checked field-by-field
# via `_struct_isequal` (see `test_utils.jl`) since core LGM doesn't
# yet define `==` on `LatentGaussianModel`.

@testset "PR-1 roundtrip equality" begin
    rng = Random.Xoshiro(20260505)
    n = 50
    df = (
        y_real=randn(rng, n),
        y_count=rand(rng, 0:5, n),
        x1=randn(rng, n),
        x2=randn(rng, n),
        idx=rand(rng, 1:5, n)
    )

    @testset "Gaussian + Intercept only" begin
        model_macro = @lgm y_real~1 data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(
            GaussianLikelihood(),
            (Intercept(),),
            sparse(ones(n, 1))
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Gaussian + Intercept + 1 covariate" begin
        model_macro = @lgm y_real~1 + x1 data=df family=GaussianLikelihood()
        X = sparse(reshape(Float64.(df.x1), n, 1))
        A = hcat(sparse(ones(n, 1)), X)
        model_hand = LatentGaussianModel(
            GaussianLikelihood(),
            (Intercept(), FixedEffects(1)),
            A
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Gaussian + Intercept + 2 covariates" begin
        model_macro = @lgm y_real~1 + x1 + x2 data=df family=GaussianLikelihood()
        X = sparse(hcat(Float64.(df.x1), Float64.(df.x2)))
        A = hcat(sparse(ones(n, 1)), X)
        model_hand = LatentGaussianModel(
            GaussianLikelihood(),
            (Intercept(), FixedEffects(2)),
            A
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Poisson + Intercept + IID" begin
        model_macro = @lgm y_count~1 + f(idx, IID(5)) data=df family=PoissonLikelihood()
        A_iid = sparse(1:n, df.idx, 1.0, n, 5)
        A = hcat(sparse(ones(n, 1)), A_iid)
        model_hand = LatentGaussianModel(
            PoissonLikelihood(),
            (Intercept(), IID(5)),
            A
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Poisson + Intercept + covariate + IID" begin
        model_macro = @lgm y_count~1 + x1 + f(idx, IID(5)) data=df family=PoissonLikelihood()
        X = sparse(reshape(Float64.(df.x1), n, 1))
        A_iid = sparse(1:n, df.idx, 1.0, n, 5)
        A = hcat(sparse(ones(n, 1)), X, A_iid)
        model_hand = LatentGaussianModel(
            PoissonLikelihood(),
            (Intercept(), FixedEffects(1), IID(5)),
            A
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Gaussian, no intercept, 1 covariate (`0 +`)" begin
        model_macro = @lgm y_real~0 + x1 data=df family=GaussianLikelihood()
        X = sparse(reshape(Float64.(df.x1), n, 1))
        model_hand = LatentGaussianModel(
            GaussianLikelihood(),
            (FixedEffects(1),),
            X
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Gaussian, no intercept, 1 covariate (`-1 +`)" begin
        model_macro = @lgm y_real~-1 + x1 data=df family=GaussianLikelihood()
        X = sparse(reshape(Float64.(df.x1), n, 1))
        model_hand = LatentGaussianModel(
            GaussianLikelihood(),
            (FixedEffects(1),),
            X
        )
        @test _struct_isequal(model_macro, model_hand)
    end
end

@testset "PR-1 lgmformula function form" begin
    rng = Random.Xoshiro(20260505)
    n = 30
    df = (y=randn(rng, n), x=randn(rng, n), idx=rand(rng, 1:4, n))

    @testset "function form matches macro form" begin
        model_macro = @lgm y~1 + x + f(idx, IID(4)) data=df family=GaussianLikelihood()
        model_func = lgmformula(df;
            lhs=:y,
            intercept=true,
            covariates=[:x],
            randoms=[(:idx, IID(4))],
            family=GaussianLikelihood())
        @test _struct_isequal(model_macro, model_func)
    end

    @testset "function form: intercept-only" begin
        model_func = lgmformula(df;
            lhs=:y, intercept=true, family=GaussianLikelihood())
        model_hand = LatentGaussianModel(
            GaussianLikelihood(), (Intercept(),), sparse(ones(n, 1)))
        @test _struct_isequal(model_func, model_hand)
    end
end
