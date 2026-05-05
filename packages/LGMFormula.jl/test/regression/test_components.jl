# PR-2 component-coverage tests. For each concrete `AbstractLatentComponent`
# whose projector follows the `f(idx_col, Component)` integer-index
# convention, assert that `@lgm y ~ 1 + f(idx, Component(...))` produces a
# model `_struct_isequal` to the hand-written Tier-1 form.
#
# Components whose projector does NOT fit the integer-index pattern (the
# multi-block components are testable here because their length is `kn`
# and PR-1's projector `sparse(1:n_obs, idx, 1.0, n_obs, kn)` correctly
# zeros the secondary blocks; the truly-deferred ones need richer parser
# machinery) are tracked at the bottom of this file with the PR that will
# pick them up.

@testset "PR-2 component coverage roundtrip" begin
    rng = Random.Xoshiro(20260506)
    n = 30

    # --- single-block components (length == n_levels) -----------------

    @testset "IID(5)" begin
        df = (y = randn(rng, n), idx = rand(rng, 1:5, n))
        comp = IID(5)
        A_block = sparse(1:n, df.idx, 1.0, n, length(comp))
        A = hcat(sparse(ones(n, 1)), A_block)
        model_macro = @lgm y ~ 1 + f(idx, IID(5)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(GaussianLikelihood(), (Intercept(), comp), A)
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "RW1(8)" begin
        df = (y = randn(rng, n), t = rand(rng, 1:8, n))
        comp = RW1(8)
        A_block = sparse(1:n, df.t, 1.0, n, length(comp))
        A = hcat(sparse(ones(n, 1)), A_block)
        model_macro = @lgm y ~ 1 + f(t, RW1(8)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(GaussianLikelihood(), (Intercept(), comp), A)
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "RW2(10)" begin
        df = (y = randn(rng, n), t = rand(rng, 1:10, n))
        comp = RW2(10)
        A_block = sparse(1:n, df.t, 1.0, n, length(comp))
        A = hcat(sparse(ones(n, 1)), A_block)
        model_macro = @lgm y ~ 1 + f(t, RW2(10)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(GaussianLikelihood(), (Intercept(), comp), A)
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "AR1(6)" begin
        df = (y = randn(rng, n), t = rand(rng, 1:6, n))
        comp = AR1(6)
        A_block = sparse(1:n, df.t, 1.0, n, length(comp))
        A = hcat(sparse(ones(n, 1)), A_block)
        model_macro = @lgm y ~ 1 + f(t, AR1(6)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(GaussianLikelihood(), (Intercept(), comp), A)
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Seasonal(12; period=4)" begin
        df = (y = randn(rng, n), m = rand(rng, 1:12, n))
        comp = Seasonal(12; period=4)
        A_block = sparse(1:n, df.m, 1.0, n, length(comp))
        A = hcat(sparse(ones(n, 1)), A_block)
        model_macro = @lgm y ~ 1 + f(m, Seasonal(12; period=4)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(GaussianLikelihood(), (Intercept(), comp), A)
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Besag(W) on 4-node graph" begin
        W = sparse(Float64[0 1 0 0; 1 0 1 0; 0 1 0 1; 0 0 1 0])
        df = (y = randn(rng, n), region = rand(rng, 1:4, n))
        comp = Besag(W)
        A_block = sparse(1:n, df.region, 1.0, n, length(comp))
        A = hcat(sparse(ones(n, 1)), A_block)
        model_macro = @lgm y ~ 1 + f(region, Besag(W)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(GaussianLikelihood(), (Intercept(), comp), A)
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Leroux(W) on 4-node graph" begin
        W = sparse(Float64[0 1 0 0; 1 0 1 0; 0 1 0 1; 0 0 1 0])
        df = (y = randn(rng, n), region = rand(rng, 1:4, n))
        comp = Leroux(W)
        A_block = sparse(1:n, df.region, 1.0, n, length(comp))
        A = hcat(sparse(ones(n, 1)), A_block)
        model_macro = @lgm y ~ 1 + f(region, Leroux(W)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(GaussianLikelihood(), (Intercept(), comp), A)
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Generic0(R)" begin
        R = sparse(Float64[2 -1 0; -1 2 -1; 0 -1 2])
        df = (y = randn(rng, n), idx = rand(rng, 1:3, n))
        comp = Generic0(R)
        A_block = sparse(1:n, df.idx, 1.0, n, length(comp))
        A = hcat(sparse(ones(n, 1)), A_block)
        model_macro = @lgm y ~ 1 + f(idx, Generic0(R)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(GaussianLikelihood(), (Intercept(), comp), A)
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Generic1(R)" begin
        R = sparse(Float64[2 -1 0; -1 2 -1; 0 -1 2])
        df = (y = randn(rng, n), idx = rand(rng, 1:3, n))
        comp = Generic1(R)
        A_block = sparse(1:n, df.idx, 1.0, n, length(comp))
        A = hcat(sparse(ones(n, 1)), A_block)
        model_macro = @lgm y ~ 1 + f(idx, Generic1(R)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(GaussianLikelihood(), (Intercept(), comp), A)
        @test _struct_isequal(model_macro, model_hand)
    end

    # --- multi-block components (length == kn; observations project to ---
    # --- columns 1:n with 0s in the secondary blocks)                  ---

    @testset "BYM(W) — length 2n, idx in 1:n" begin
        W = sparse(Float64[0 1 0 0; 1 0 1 0; 0 1 0 1; 0 0 1 0])
        df = (y = randn(rng, n), region = rand(rng, 1:4, n))
        comp = BYM(W)
        A_block = sparse(1:n, df.region, 1.0, n, length(comp))
        A = hcat(sparse(ones(n, 1)), A_block)
        model_macro = @lgm y ~ 1 + f(region, BYM(W)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(GaussianLikelihood(), (Intercept(), comp), A)
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "BYM2(W) — length 2n, idx in 1:n" begin
        W = sparse(Float64[0 1 0 0; 1 0 1 0; 0 1 0 1; 0 0 1 0])
        df = (y = randn(rng, n), region = rand(rng, 1:4, n))
        comp = BYM2(W)
        A_block = sparse(1:n, df.region, 1.0, n, length(comp))
        A = hcat(sparse(ones(n, 1)), A_block)
        model_macro = @lgm y ~ 1 + f(region, BYM2(W)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(GaussianLikelihood(), (Intercept(), comp), A)
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Generic2(C) — length 2n, idx in 1:n" begin
        C = sparse(Float64[2 -1 0; -1 2 -1; 0 -1 2])
        df = (y = randn(rng, n), idx = rand(rng, 1:3, n))
        comp = Generic2(C)
        A_block = sparse(1:n, df.idx, 1.0, n, length(comp))
        A = hcat(sparse(ones(n, 1)), A_block)
        model_macro = @lgm y ~ 1 + f(idx, Generic2(C)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(GaussianLikelihood(), (Intercept(), comp), A)
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "IID2D(5) — length 2n, idx in 1:n" begin
        df = (y = randn(rng, n), pair = rand(rng, 1:5, n))
        comp = IID2D(5)
        A_block = sparse(1:n, df.pair, 1.0, n, length(comp))
        A = hcat(sparse(ones(n, 1)), A_block)
        model_macro = @lgm y ~ 1 + f(pair, IID2D(5)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(GaussianLikelihood(), (Intercept(), comp), A)
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "IID3D(4) — length 3n, idx in 1:n" begin
        df = (y = randn(rng, n), triple = rand(rng, 1:4, n))
        comp = IID3D(4)
        A_block = sparse(1:n, df.triple, 1.0, n, length(comp))
        A = hcat(sparse(ones(n, 1)), A_block)
        model_macro = @lgm y ~ 1 + f(triple, IID3D(4)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(GaussianLikelihood(), (Intercept(), comp), A)
        @test _struct_isequal(model_macro, model_hand)
    end

    # --- intercept-suppressed combinations ----------------------------

    @testset "no intercept + RW1" begin
        df = (y = randn(rng, n), t = rand(rng, 1:5, n))
        comp = RW1(5)
        A = sparse(1:n, df.t, 1.0, n, length(comp))
        model_macro = @lgm y ~ 0 + f(t, RW1(5)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(GaussianLikelihood(), (comp,), A)
        @test _struct_isequal(model_macro, model_hand)
    end
end

# --- Components NOT covered by PR-1's `f(idx, Component)` parser. -----
# Tracked here so PR-3+ has a checklist:
#
# - `Replicate(component, n_replicates)`: wrapper, requires the
#   `replicate = rep_col` kwarg routing planned for PR-5.
# - `Group(component, group_id)` / `Group(components)`: wrapper,
#   requires the `group = group_col` kwarg routing planned for PR-5.
# - `KroneckerComponent(spatial, temporal)`: composes two children;
#   PR-7 (stretch) covers the `f((s, t), KroneckerComponent(...))`
#   tuple-coordinate syntax.
# - `UserComponent(callable; n, θ0)`: pass-through works at the macro
#   level once a user constructs one, but the `f(col, UserComponent(...))`
#   surface is undertested without a representative end-to-end fixture.
#   Defer until a user-facing example is on hand.
# - `MEB(values)` / `MEC(values)`: per-observation measurement-error
#   components — `length` matches the data vector, not a level set, so
#   the `idx` column is conceptually `1:n_obs`. The natural surface is
#   probably `f(:obs, MEB(measured_x))` or a `mec`/`meb`-specific
#   keyword; defer to PR-5 alongside `replicate`/`group`.
# - `SPDE1D(mesh)` / `SPDE2(mesh)` / `SPDE2NonStationary(...)`: the
#   projector is barycentric over coordinates, not a 0/1 indicator on
#   integer levels. PR-7 (stretch) covers the SPDE-friendly forms; an
#   intermediate ADR is likely needed for the coordinate-column
#   convention.
