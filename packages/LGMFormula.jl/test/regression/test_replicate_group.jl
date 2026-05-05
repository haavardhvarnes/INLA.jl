# PR-5 `replicate` / `group` routing.
#
# `f(col, Comp; replicate = id_col)` lowers to runtime
# `Replicate(Comp, R)` with `R = maximum(data.id_col)`. The projector
# block has `R · length(Comp)` columns; row i is 1 in column
# `(id[i] - 1) · length(Comp) + col[i]`. Latent layout matches R-INLA's
# `f(idx, model="ar1", replicate=id)` panel-stacking convention used by
# `LatentGaussianModels.jl/test/oracle/synthetic_replicate_ar1` (Phase
# I-C).
#
# `f(col, Factory; group = grp_col)` lowers to `Group(Factory,
# data.grp_col)` (LGM core factory-form). Per-group sizes come from
# counting the labels in `grp_col`; the projector block has
# `Σ_g s_g` columns laid out group-by-group.

@testset "PR-5 replicate roundtrip" begin
    rng = Random.Xoshiro(20260507)
    n_per = 6
    R = 4
    N = n_per * R
    df = (
        y = randn(rng, N),
        t = repeat(1:n_per, R),
        id = repeat(1:R, inner = n_per),
        x = randn(rng, N),
    )

    @testset "Gaussian + Intercept + Replicate(IID, R)" begin
        model_macro = @lgm y ~ 1 + f(t, IID(n_per); replicate = id) data=df family=GaussianLikelihood()
        # Expected: A_block = sparse(1:N, (id-1)*n_per + t, 1.0, N, R*n_per).
        block_cols = (df.id .- 1) .* n_per .+ df.t
        A_rep = sparse(1:N, block_cols, 1.0, N, R * n_per)
        A = hcat(sparse(ones(N, 1)), A_rep)
        model_hand = LatentGaussianModel(
            GaussianLikelihood(),
            (Intercept(), Replicate(IID(n_per), R)),
            A,
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Gaussian + Intercept + Replicate(AR1, R)" begin
        model_macro = @lgm y ~ 1 + f(t, AR1(n_per); replicate = id) data=df family=GaussianLikelihood()
        block_cols = (df.id .- 1) .* n_per .+ df.t
        A_rep = sparse(1:N, block_cols, 1.0, N, R * n_per)
        A = hcat(sparse(ones(N, 1)), A_rep)
        model_hand = LatentGaussianModel(
            GaussianLikelihood(),
            (Intercept(), Replicate(AR1(n_per), R)),
            A,
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "intercept + covariate + Replicate(IID)" begin
        model_macro = @lgm y ~ 1 + x + f(t, IID(n_per); replicate = id) data=df family=GaussianLikelihood()
        X = sparse(reshape(Float64.(df.x), N, 1))
        block_cols = (df.id .- 1) .* n_per .+ df.t
        A_rep = sparse(1:N, block_cols, 1.0, N, R * n_per)
        A = hcat(sparse(ones(N, 1)), X, A_rep)
        model_hand = LatentGaussianModel(
            GaussianLikelihood(),
            (Intercept(), FixedEffects(1), Replicate(IID(n_per), R)),
            A,
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Replicate + plain f mixed" begin
        df2 = (
            y = randn(rng, N),
            t = repeat(1:n_per, R),
            id = repeat(1:R, inner = n_per),
            s = rand(rng, 1:3, N),
        )
        model_macro = @lgm y ~ 1 + f(t, IID(n_per); replicate = id) + f(s, IID(3)) data=df2 family=GaussianLikelihood()
        block_cols = (df2.id .- 1) .* n_per .+ df2.t
        A_rep = sparse(1:N, block_cols, 1.0, N, R * n_per)
        A_s = sparse(1:N, df2.s, 1.0, N, 3)
        A = hcat(sparse(ones(N, 1)), A_rep, A_s)
        model_hand = LatentGaussianModel(
            GaussianLikelihood(),
            (Intercept(), Replicate(IID(n_per), R), IID(3)),
            A,
        )
        @test _struct_isequal(model_macro, model_hand)
    end
end

@testset "PR-5 group roundtrip (factory form)" begin
    rng = Random.Xoshiro(20260507)
    # Three groups of varying size: G=3, sizes (3, 2, 4) → total 9 obs.
    df = (
        y = randn(rng, 9),
        t = [1, 2, 3, 1, 2, 1, 2, 3, 4],
        grp = [1, 1, 1, 2, 2, 3, 3, 3, 3],
    )
    G = 3
    sizes = [3, 2, 4]
    offsets = [0; cumsum(sizes)[1:(end - 1)]]
    block_cols = offsets[df.grp] .+ df.t

    @testset "Gaussian + Intercept + Group(IID, grp)" begin
        model_macro = @lgm y ~ 1 + f(t, IID; group = grp) data=df family=GaussianLikelihood()
        A_grp = sparse(1:length(df.y), block_cols, 1.0, length(df.y), sum(sizes))
        A = hcat(sparse(ones(length(df.y), 1)), A_grp)
        model_hand = LatentGaussianModel(
            GaussianLikelihood(),
            (Intercept(), Group(IID, df.grp)),
            A,
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Gaussian + Intercept + Group(AR1, grp)" begin
        model_macro = @lgm y ~ 1 + f(t, AR1; group = grp) data=df family=GaussianLikelihood()
        A_grp = sparse(1:length(df.y), block_cols, 1.0, length(df.y), sum(sizes))
        A = hcat(sparse(ones(length(df.y), 1)), A_grp)
        model_hand = LatentGaussianModel(
            GaussianLikelihood(),
            (Intercept(), Group(AR1, df.grp)),
            A,
        )
        @test _struct_isequal(model_macro, model_hand)
    end
end

@testset "PR-5 lgmformula function form" begin
    rng = Random.Xoshiro(20260507)
    n_per = 5
    R = 3
    N = n_per * R
    df = (
        y = randn(rng, N),
        t = repeat(1:n_per, R),
        id = repeat(1:R, inner = n_per),
    )

    @testset "function form: Replicate" begin
        model_macro = @lgm y ~ 1 + f(t, IID(n_per); replicate = id) data=df family=GaussianLikelihood()
        model_func = lgmformula(df;
            lhs = :y,
            intercept = true,
            randoms = [(:t, IID(n_per), :id, nothing)],
            family = GaussianLikelihood())
        @test _struct_isequal(model_macro, model_func)
    end

    @testset "function form: Group" begin
        df2 = (y = randn(rng, 6), t = [1, 2, 1, 2, 3, 1], grp = [1, 1, 2, 2, 2, 3])
        model_macro = @lgm y ~ 1 + f(t, IID; group = grp) data=df2 family=GaussianLikelihood()
        model_func = lgmformula(df2;
            lhs = :y,
            intercept = true,
            randoms = [(:t, IID, nothing, :grp)],
            family = GaussianLikelihood())
        @test _struct_isequal(model_macro, model_func)
    end
end

@testset "PR-5 @macroexpand structure" begin
    df_dummy = (y = Float64[], t = Int[], id = Int[], grp = Int[])

    @testset "replicate term emits _wrap_term in components" begin
        ex = @macroexpand @lgm y ~ 1 + f(t, IID(5); replicate = id) data=df_dummy family=GaussianLikelihood()
        call = _find_lgm_call(ex)
        components = call.args[3]
        @test components.head === :tuple
        @test length(components.args) == 2
        @test _is_call_to(components.args[1], :Intercept)
        # Second component is a `_wrap_term(...)` call, not a bare `IID(5)`.
        @test _is_call_to(components.args[2], :_wrap_term)
    end

    @testset "group term emits _wrap_term in components" begin
        ex = @macroexpand @lgm y ~ 1 + f(t, IID; group = grp) data=df_dummy family=GaussianLikelihood()
        call = _find_lgm_call(ex)
        components = call.args[3]
        @test _is_call_to(components.args[2], :_wrap_term)
    end

    @testset "plain f-term unchanged (no _wrap_term)" begin
        ex = @macroexpand @lgm y ~ 1 + f(t, IID(5)) data=df_dummy family=GaussianLikelihood()
        call = _find_lgm_call(ex)
        components = call.args[3]
        # Bare `IID(5)` call — no _wrap_term.
        @test _is_call_to(components.args[2], :IID)
        @test !_is_call_to(components.args[2], :_wrap_term)
    end
end

@testset "PR-5 error messages refer to user concepts" begin
    rng = Random.Xoshiro(20260507)
    N = 12
    df = (
        y = randn(rng, N),
        t = repeat(1:4, 3),
        id = repeat(1:3, inner = 4),
        grp = repeat(1:3, inner = 4),
    )

    @testset "replicate column missing from data" begin
        @test_throws ArgumentError begin
            @lgm y ~ 1 + f(t, IID(4); replicate = no_id) data=df family=GaussianLikelihood()
        end
        try
            @lgm y ~ 1 + f(t, IID(4); replicate = no_id) data=df family=GaussianLikelihood()
        catch e
            msg = sprint(showerror, e)
            @test occursin("no_id", msg)
        end
    end

    @testset "group column missing from data" begin
        @test_throws ArgumentError begin
            @lgm y ~ 1 + f(t, IID; group = no_grp) data=df family=GaussianLikelihood()
        end
        try
            @lgm y ~ 1 + f(t, IID; group = no_grp) data=df family=GaussianLikelihood()
        catch e
            msg = sprint(showerror, e)
            @test occursin("no_grp", msg)
        end
    end

    @testset "both replicate and group rejected" begin
        ex = try
            @eval @lgm y ~ 1 + f(t, IID(4); replicate = id, group = grp) data=$df family=GaussianLikelihood()
            nothing
        catch e
            e
        end
        @test ex isa LoadError || ex isa ErrorException
        msg = sprint(showerror, ex)
        @test occursin("replicate", msg) && occursin("group", msg)
    end

    @testset "unsupported f-term kwarg" begin
        ex = try
            @eval @lgm y ~ 1 + f(t, IID(4); foo = id) data=$df family=GaussianLikelihood()
            nothing
        catch e
            e
        end
        @test ex isa LoadError || ex isa ErrorException
        msg = sprint(showerror, ex)
        @test occursin("foo", msg) || occursin("kwarg", msg) || occursin("keyword", msg)
    end

    @testset "non-symbol replicate value" begin
        ex = try
            @eval @lgm y ~ 1 + f(t, IID(4); replicate = 7) data=$df family=GaussianLikelihood()
            nothing
        catch e
            e
        end
        @test ex isa LoadError || ex isa ErrorException
    end
end
