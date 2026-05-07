# PR-3 multi-`f` roundtrip coverage. PR-1 already synthesises multiple
# `f(...)` terms by looping over `randoms` and `hcat`-ing per-term
# blocks; PR-3's contract is that `_struct_isequal(model_macro,
# model_hand)` holds for any number of `f(...)` terms (with or without
# covariates / intercept), and for the canonical disease-mapping
# *shapes* (Scotland BYM2 + Poisson, Tokyo RW2 + Binomial).
#
# See `plans/phase-n.md` PR-3 scope correction (2026-05-05) for why
# this lands as test-only — LGM core's `StackedMapping` is a row-
# partition structure for multi-likelihood, not a per-component
# column-partition, so the plan's original `Component => Mapping`
# pseudocode would have required a (last-resort) LGM core change.

@testset "PR-3 multi-`f` roundtrip" begin
    rng = Random.Xoshiro(20260507)
    n = 60

    @testset "Gaussian + BYM2 + AR1 (two effects)" begin
        W = sparse(Float64[0 1 0 0; 1 0 1 0; 0 1 0 1; 0 0 1 0])
        df = (
            y=randn(rng, n),
            region=rand(rng, 1:4, n),
            time=rand(rng, 1:6, n)
        )
        bym2 = BYM2(W)
        ar1 = AR1(6)
        A = hcat(
            sparse(ones(n, 1)),
            sparse(1:n, df.region, 1.0, n, length(bym2)),
            sparse(1:n, df.time, 1.0, n, length(ar1))
        )
        model_macro = @lgm y~1 + f(region, BYM2(W)) + f(time, AR1(6)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(
            GaussianLikelihood(),
            (Intercept(), bym2, ar1),
            A
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Poisson + Besag + IID (two effects)" begin
        W = sparse(Float64[0 1 1; 1 0 1; 1 1 0])
        df = (
            y=rand(rng, 0:5, n),
            area=rand(rng, 1:3, n),
            period=rand(rng, 1:5, n)
        )
        besag = Besag(W)
        iid = IID(5)
        A = hcat(
            sparse(ones(n, 1)),
            sparse(1:n, df.area, 1.0, n, length(besag)),
            sparse(1:n, df.period, 1.0, n, length(iid))
        )
        model_macro = @lgm y~1 + f(area, Besag(W)) + f(period, IID(5)) data=df family=PoissonLikelihood()
        model_hand = LatentGaussianModel(
            PoissonLikelihood(),
            (Intercept(), besag, iid),
            A
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Gaussian + IID + RW1 + Seasonal (three effects)" begin
        df = (
            y=randn(rng, n),
            g=rand(rng, 1:4, n),
            t=rand(rng, 1:8, n),
            m=rand(rng, 1:12, n)
        )
        iid = IID(4)
        rw1 = RW1(8)
        seasonal = Seasonal(12; period=4)
        A = hcat(
            sparse(ones(n, 1)),
            sparse(1:n, df.g, 1.0, n, length(iid)),
            sparse(1:n, df.t, 1.0, n, length(rw1)),
            sparse(1:n, df.m, 1.0, n, length(seasonal))
        )
        model_macro = @lgm y~1 + f(g, IID(4)) + f(t, RW1(8)) + f(m, Seasonal(12; period=4)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(
            GaussianLikelihood(),
            (Intercept(), iid, rw1, seasonal),
            A
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Gaussian + 2 covariates + IID + RW1 (mixed)" begin
        df = (
            y=randn(rng, n),
            x1=randn(rng, n),
            x2=randn(rng, n),
            g=rand(rng, 1:4, n),
            t=rand(rng, 1:6, n)
        )
        iid = IID(4)
        rw1 = RW1(6)
        X = sparse(hcat(Float64.(df.x1), Float64.(df.x2)))
        A = hcat(
            sparse(ones(n, 1)),
            X,
            sparse(1:n, df.g, 1.0, n, length(iid)),
            sparse(1:n, df.t, 1.0, n, length(rw1))
        )
        model_macro = @lgm y~1 + x1 + x2 + f(g, IID(4)) + f(t, RW1(6)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(
            GaussianLikelihood(),
            (Intercept(), FixedEffects(2), iid, rw1),
            A
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "no intercept + IID + RW1" begin
        df = (
            y=randn(rng, n),
            g=rand(rng, 1:4, n),
            t=rand(rng, 1:6, n)
        )
        iid = IID(4)
        rw1 = RW1(6)
        A = hcat(
            sparse(1:n, df.g, 1.0, n, length(iid)),
            sparse(1:n, df.t, 1.0, n, length(rw1))
        )
        model_macro = @lgm y~0 + f(g, IID(4)) + f(t, RW1(6)) data=df family=GaussianLikelihood()
        model_hand = LatentGaussianModel(
            GaussianLikelihood(),
            (iid, rw1),
            A
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    # --- Canonical disease-mapping shape roundtrips -------------------

    @testset "Scotland BYM2 shape: Poisson + covariate + f(area, BYM2)" begin
        # Five-area synthetic adjacency mimicking the Scotland-lip-
        # cancer disease-mapping shape (the actual fit is in
        # LatentGaussianModels.jl/test/oracle; here we only verify the
        # macro produces the right model object).
        W = sparse(Float64[0 1 1 0 0
                           1 0 1 1 0
                           1 1 0 1 1
                           0 1 1 0 1
                           0 0 1 1 0])
        df = (
            cases=rand(rng, 0:10, n),
            x=randn(rng, n),
            area=rand(rng, 1:5, n)
        )
        bym2 = BYM2(W)
        X = sparse(reshape(Float64.(df.x), n, 1))
        A = hcat(
            sparse(ones(n, 1)),
            X,
            sparse(1:n, df.area, 1.0, n, length(bym2))
        )
        model_macro = @lgm cases~1 + x + f(area, BYM2(W)) data=df family=PoissonLikelihood()
        model_hand = LatentGaussianModel(
            PoissonLikelihood(),
            (Intercept(), FixedEffects(1), bym2),
            A
        )
        @test _struct_isequal(model_macro, model_hand)
    end

    @testset "Tokyo RW2 shape: Binomial + f(day, RW2(N; cyclic=true))" begin
        # Tokyo-rainfall-shape: 366 days of the year (cyclic), Binomial
        # over n_trials per day. Synthetic data; oracle is in
        # LatentGaussianModels.jl/test/oracle.
        n_days = 24      # synthetic — oracle uses 366
        n_obs = n
        df = (
            y=rand(rng, 0:5, n_obs),
            day=rand(rng, 1:n_days, n_obs)
        )
        n_trials = fill(5, n_obs)
        rw2 = RW2(n_days; cyclic=true)
        A = hcat(
            sparse(ones(n_obs, 1)),
            sparse(1:n_obs, df.day, 1.0, n_obs, length(rw2))
        )
        model_macro = @lgm y~1 + f(day, RW2(n_days; cyclic=true)) data=df family=BinomialLikelihood(n_trials)
        model_hand = LatentGaussianModel(
            BinomialLikelihood(n_trials),
            (Intercept(), rw2),
            A
        )
        @test _struct_isequal(model_macro, model_hand)
    end
end
