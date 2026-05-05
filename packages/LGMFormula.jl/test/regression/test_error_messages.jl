# Per `packages/LGMFormula.jl/CLAUDE.md` "Errors refer to user
# concepts": every malformed `@lgm` input must raise an error
# referring to columns, components, or the formula — never
# `FunctionTerm`/`Expr` internals.

@testset "PR-1 error messages refer to user concepts" begin
    rng = Random.Xoshiro(20260505)
    n = 20
    df = (y = randn(rng, n), x = randn(rng, n), idx = rand(rng, 1:3, n))

    @testset "missing data kwarg" begin
        ex = try
            @eval @lgm y ~ 1 family=GaussianLikelihood()
            nothing
        catch e
            e
        end
        @test ex isa LoadError || ex isa ErrorException
        msg = sprint(showerror, ex)
        @test occursin("data", msg)
    end

    @testset "missing family kwarg" begin
        ex = try
            @eval @lgm y ~ 1 data=$df
            nothing
        catch e
            e
        end
        @test ex isa LoadError || ex isa ErrorException
        msg = sprint(showerror, ex)
        @test occursin("family", msg)
    end

    @testset "unknown outcome column" begin
        @test_throws ArgumentError begin
            @lgm not_a_col ~ 1 + x data=df family=GaussianLikelihood()
        end
        try
            @lgm not_a_col ~ 1 + x data=df family=GaussianLikelihood()
        catch e
            msg = sprint(showerror, e)
            @test occursin("not_a_col", msg)
            @test occursin("data", msg)
        end
    end

    @testset "unknown covariate column" begin
        @test_throws ArgumentError begin
            @lgm y ~ 1 + missing_col data=df family=GaussianLikelihood()
        end
        try
            @lgm y ~ 1 + missing_col data=df family=GaussianLikelihood()
        catch e
            msg = sprint(showerror, e)
            @test occursin("missing_col", msg)
        end
    end

    @testset "unknown f-term column" begin
        @test_throws ArgumentError begin
            @lgm y ~ 1 + f(no_idx, IID(3)) data=df family=GaussianLikelihood()
        end
        try
            @lgm y ~ 1 + f(no_idx, IID(3)) data=df family=GaussianLikelihood()
        catch e
            msg = sprint(showerror, e)
            @test occursin("no_idx", msg)
        end
    end

    @testset "f-term level out of range" begin
        bad_df = (y = randn(rng, n), idx = fill(99, n))
        @test_throws ArgumentError begin
            @lgm y ~ 1 + f(idx, IID(3)) data=bad_df family=GaussianLikelihood()
        end
        try
            @lgm y ~ 1 + f(idx, IID(3)) data=bad_df family=GaussianLikelihood()
        catch e
            msg = sprint(showerror, e)
            @test occursin("idx", msg)
            @test occursin("99", msg) || occursin("3", msg)
        end
    end
end
