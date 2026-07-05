# Coming from R-INLA

A side-by-side reference for users coming from
[R-INLA](https://www.r-inla.org/). Every R-INLA `family`, `f(model =
...)`, prior, integration strategy, and post-fit accessor that the
Julia ecosystem currently supports has a row here, with the Julia
constructor a user would type and (where applicable) the
[`@lgm`](lgmformula-tutorial.md) formula-DSL analogue.

If you'd rather start from a worked example, the formula-DSL tutorial
walks through Scotland BYM2 and Meuse SPDE end-to-end. This page is
the lookup reference you keep open next to a translation in progress.

!!! note "Mental model — three differences"
    1. **Components and priors are constructors, not strings.**
       R-INLA's `f(idx, model = "bym2", graph = W, hyper = list(prec
       = list(prior = "pc.prec", param = c(1, 0.01))))` becomes
       `BYM2(GMRFGraph(W); hyperprior_prec = PCPrecision(1.0,
       0.01))`. The component is just a Julia struct.
    2. **Models are values, not formulas.** A
       [`LatentGaussianModel`](packages/lgm.md) bundles likelihood,
       components, and observation mapping into a single object you
       can pass to [`inla`](packages/lgm.md), to
       [`LGMTuring.nuts_sample`](https://github.com/HaavardHvarnes/INLA.jl/tree/main/packages/LGMTuring.jl),
       or to
       [`predict_raster`](packages/inlaspderasters.md). The
       [`@lgm`](lgmformula-tutorial.md) macro is sugar that builds
       this object from a formula.
    3. **`inla(model, y)` returns an `INLAResult`.** All summaries
       come through accessor functions —
       [`fixed_effects`](packages/lgm.md),
       [`random_effects`](packages/lgm.md),
       [`hyperparameters`](packages/lgm.md),
       [`log_marginal_likelihood`](packages/lgm.md) — instead of
       `res$summary.fixed`, `res$summary.random`, etc.

## Likelihoods

R-INLA's `family = "..."` string maps to a Julia `<Likelihood>()`
constructor. The `@lgm` column is what you'd write inside
`@lgm y ~ ... family = ...`.

| R-INLA `family` | Julia constructor | `@lgm family =` |
|---|---|---|
| `"gaussian"` | `GaussianLikelihood(; hyperprior = PCPrecision(1.0, 0.01))` | `GaussianLikelihood()` |
| `"poisson"` | `PoissonLikelihood(; E = expected)` | `PoissonLikelihood()` |
| `"binomial"` | `BinomialLikelihood(n_trials)` | `BinomialLikelihood(n)` |
| `"nbinomial"` | `NegativeBinomialLikelihood(; E = expected)` | `NegativeBinomialLikelihood()` |
| `"gamma"` | `GammaLikelihood()` | `GammaLikelihood()` |
| `"beta"` | `BetaLikelihood()` | `BetaLikelihood()` |
| `"betabinomial"` | `BetaBinomialLikelihood(n_trials)` | — |
| `"T"` | `StudentTLikelihood()` | `StudentTLikelihood()` |
| `"skewnormal"` | `SkewNormalLikelihood()` | — |
| `"gev"` | `GEVLikelihood()` | — |
| `"pom"` | `POMLikelihood(n_classes)` | — |
| `"exponentialsurv"` | `ExponentialLikelihood(; censoring)` | — |
| `"weibullsurv"` | `WeibullLikelihood(; censoring)` | `WeibullLikelihood()` |
| `"lognormalsurv"` | `LognormalSurvLikelihood(; censoring)` | — |
| `"gammasurv"` | `GammaSurvLikelihood(; censoring)` | — |
| `"weibullcuresurv"` | `WeibullCureLikelihood(; censoring)` | — |
| `"zeroinflated*0"` | `ZeroInflated*Likelihood0(...)` | — |
| `"zeroinflated*1"` | `ZeroInflated*Likelihood1(...)` | — |
| `"zeroinflated*2"` | `ZeroInflated*Likelihood2(...)` | — |
| `inla.coxph(...)` | `inla_coxph(time, event; nbreakpoints = 15)` | — |

The zero-inflated `*` stands for `Poisson`, `Binomial`, or
`NegativeBinomial`. The three suffixes match R-INLA's three
zero-inflation parameterisations.

`censoring` is an `(int_left, int_right)` pair per observation; pass
`nothing` for fully observed data.

## Latent components

R-INLA's `f(idx, model = "...")` becomes a struct constructor in
Julia. Inside `@lgm`, the same struct is passed positionally to `f()`.

| R-INLA `f(model = ...)` | Julia constructor | `@lgm f(...)` |
|---|---|---|
| `"linear"` | `Intercept(; prec = 1.0e-3, improper = true)` | `1` (in formula) |
| — | `FixedEffects(p; prec = 1.0e-3)` | `x1 + x2 + ...` (in formula) |
| `"iid"` | `IID(n; hyperprior = PCPrecision(1.0, 0.01))` | `f(idx, IID(n))` |
| `"2diid"` | `IID2D(n)` | — |
| `"3diid"` | `IID3D(n)` | — |
| `"rw1"` | `RW1(n; cyclic = false)` | `f(t, RW1(T))` |
| `"rw2"` | `RW2(n; cyclic = false)` | `f(t, RW2(T))` |
| `"ar1"` | `AR1(n)` | `f(t, AR1(T))` |
| `"seasonal"` | `Seasonal(n; period = s)` | `f(t, Seasonal(T; period = s))` |
| `"besag"` | `Besag(graph; scale_model = true)` | `f(area, Besag(GMRFGraph(W)))` |
| `"bym"` | `BYM(graph)` | `f(area, BYM(GMRFGraph(W)))` |
| `"bym2"` | `BYM2(graph; hyperprior_prec = PCPrecision(1.0, 0.01))` | `f(area, BYM2(GMRFGraph(W)))` |
| `"leroux"` | `Leroux(graph)` | `f(area, Leroux(GMRFGraph(W)))` |
| `"generic0"` | `Generic0(R; rankdef = 0)` | `f(idx, Generic0(R))` |
| `"generic1"` | `Generic1(C, G)` | — |
| `"generic2"` | `Generic2(C, G1, G2)` | — |
| `"meb"` | `MEB(values)` | — |
| `"mec"` | `MEC(values)` | — |
| `..., replicate = id` | `Replicate(component, n_replicates)` | `f(idx, IID(n); replicate = id)` |
| `..., group = grp` | `Group(factory, group_id)` | `f(idx, AR1; group = grp)` |
| `"spde2"` (stationary) | `SPDE2(mesh)` (in `INLASPDE.jl`) | `f(loc_x, loc_y, SPDE2(mesh))` |
| `"spde2"` (non-stationary) | `SPDE2NonStationary(mesh; B_τ, B_κ)` | — |
| `rgeneric(...)` | `UserComponent(callable; n)` | — |
| Kronecker space-time | `KroneckerComponent(spatial, temporal)` | — |

For `BYM2`, the precision graph wraps R-INLA's `W` matrix:
`GMRFGraph(W)` from [`GMRFs.jl`](packages/gmrfs.md) builds the
adjacency structure once and the component reuses it. The `α`
keyword on `BYM2` controls the PC prior on φ; the default matches
R-INLA's recommended setting.

## Hyperpriors

R-INLA expresses priors via `hyper = list(prec = list(prior = "...",
param = ...))`. In Julia each prior is a struct passed to the
component or likelihood that owns the hyperparameter.

| R-INLA prior string | Julia constructor | Notes |
|---|---|---|
| `"pc.prec"` | `PCPrecision(U = 1.0, α = 0.01)` | Default for most components |
| `"pc.matern"` | `PCMatern(range_U, range_α, sigma_U, sigma_α)` | Joint PC on (range, σ) for SPDE |
| `"loggamma"` | `GammaPrecision(shape, rate)` | On precision scale (R-INLA's internal) |
| `"lognormal"` | `LogNormalPrecision(μ_log, σ_log)` | |
| `"gaussian"` | `GaussianPrior(μ = 0.0, σ = 1.0)` | On internal unconstrained scale |
| `"beta"` | `BetaPrior(a, b)` | |
| flat | `WeakPrior()` | Improper; rejected by some components |

### Worked priors

#### BYM2 — PC prior on precision and φ

R-INLA:

```r
f(area, model = "bym2", graph = W,
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
               phi  = list(prior = "pc",      param = c(0.5, 2/3))))
```

Julia:

```julia
BYM2(GMRFGraph(W);
     hyperprior_prec = PCPrecision(1.0, 0.01),
     U = 0.5, α = 2/3)
```

The φ prior is parameterised via the PC-mixing-weight construction
(Riebler et al. 2016): `P(φ < U) = α`. R-INLA exposes this through the
`phi` hyper entry; Julia exposes the same through the `U` and `α`
keyword arguments on `BYM2(...)`.

#### RW2 — PC prior on smoothing precision

R-INLA:

```r
f(t, model = "rw2",
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))
```

Julia:

```julia
RW2(T; hyperprior = PCPrecision(1.0, 0.01))
```

The PC prior is on the precision; the user-scale interpretation is
"a marginal SD greater than 1 has 1% prior probability."

#### SPDE2 — PC-Matérn joint prior

R-INLA:

```r
spde <- inla.spde2.pcmatern(mesh,
                            prior.range = c(0.5, 0.5),
                            prior.sigma = c(1.0, 0.01))
f(field, model = spde)
```

Julia:

```julia
SPDE2(mesh; pc = PCMatern(range_U = 0.5, range_α = 0.5,
                          sigma_U = 1.0, sigma_α = 0.01))
```

`PCMatern` is the joint PC prior on the Matérn range ρ and marginal
SD σ (Fuglstad et al. 2019). The four arguments mean: `P(ρ < range_U)
= range_α` and `P(σ > sigma_U) = sigma_α`.

!!! warning "Precision vs SD scale"
    R-INLA documents some priors on the **SD** scale (`pc.prec`'s `U`
    is "P(σ > U) = α", σ being the marginal SD). Julia follows that
    convention for `PCPrecision` and `PCMatern`. But R-INLA's `loggamma`
    is on the **precision** scale (`τ = 1/σ²`). Julia's
    `GammaPrecision(shape, rate)` follows the same precision convention,
    so a literal `param = c(1, 5e-5)` in R-INLA translates to
    `GammaPrecision(1.0, 5.0e-5)` in Julia. See
    [`plans/defaults-parity.md`](https://github.com/HaavardHvarnes/INLA.jl/blob/main/plans/defaults-parity.md)
    for the complete table of where R-INLA and Julia agree, where they
    differ, and why.

## Inference flags

`inla(model, y; kwargs...)` is the Julia analogue of R-INLA's
`inla(formula, family, data, control.inla = ...)`.

| R-INLA flag | Julia kwarg | Default |
|---|---|---|
| `int.strategy = "auto"` | `int_strategy = :auto` | `:auto` (= `:ccd` for dim(θ) > 2, else `:grid`) |
| `int.strategy = "ccd"` | `int_strategy = :ccd` | |
| `int.strategy = "grid"` | `int_strategy = :grid` | |
| `int.strategy = "eb"` | pass `EmpiricalBayes()` to `fit(...)` (or `empirical_bayes(model, y)`) | |
| `control.inla = list(strategy = "gaussian")` | `latent_strategy = Gaussian()` | `Gaussian()` |
| `control.inla = list(strategy = "laplace")` | pass `Laplace()` to `fit(...)` | |
| `control.compute = list(dic = TRUE)` | `dic(res, model, y)` (post-fit) | always available |
| `control.compute = list(waic = TRUE)` | `waic(rng, res, model, y)` | |
| `control.compute = list(cpo = TRUE)` | `cpo(rng, res, model, y)` | |
| `control.compute = list(po = TRUE)` | not yet shipped | |
| `inla.posterior.sample(...)` | `posterior_sample(rng, res, model)` | |
| `control.predictor = list(compute = TRUE)` | `posterior_predictive(...)` | |

The Julia API splits R-INLA's monolithic `control.compute` into
post-fit accessor functions: nothing is computed at fit time beyond
the marginals, and you call exactly the diagnostics you need. This
also means cross-validation diagnostics that need posterior samples
take an explicit `rng::AbstractRNG` argument — [Julia's RNG
state](https://docs.julialang.org/en/v1/stdlib/Random/#Reproducibility)
is never global in this ecosystem.

## Accessors

R-INLA returns a list with `$summary.fixed`, `$summary.random`,
`$summary.hyperpar`, `$mlik`, `$dic`, `$waic`, `$cpo`. Julia returns
an `INLAResult` and gives one accessor per concern.

| R-INLA | Julia |
|---|---|
| `res$summary.fixed` | `fixed_effects(model, res; level = 0.95)` |
| `res$summary.random` | `random_effects(model, res; level = 0.95)` |
| `res$summary.hyperpar` | `hyperparameters(model, res; level = 0.95)` |
| `res$mlik[1]` | `log_marginal_likelihood(res)` |
| `res$dic$dic` | `dic(res, model, y)` |
| `res$waic$waic` | `waic(rng, res, model, y)` |
| `res$cpo$cpo` | `cpo(rng, res, model, y)` |
| `res$cpo$pit` | `pit(rng, res, model, y)` |
| `res$loo$elpd_loo` (via `loo` package) | `psis_loo(res, model, y)` |
| `res$marginals.fixed[[i]]` / `res$marginals.random` | `posterior_marginal_x(res, i)` |
| `res$marginals.hyperpar[[j]]` | `posterior_marginal_θ(res, j; model = model, y = y)` |
| (no equivalent) | `pp_check(rng, res, model, y_obs)` |

`fixed_effects`, `random_effects`, and `hyperparameters` return a
vector of named tuples with `.name`, `.mean`, `.sd`, `.lower`,
`.upper`, where `level` controls the credible-interval coverage.

One scale caveat on the marginal densities: R-INLA's
`marginals.hyperpar` are on the *user* scale (precision, correlation),
while `posterior_marginal_θ` returns the density on the *internal*
scale (log precision, logit correlation). Transform with the usual
change-of-variables Jacobian, e.g. `p(τ) = p_θ(log τ) / τ` for a
precision. By default the density is integrated over the INLA design
(R-INLA's integrated hyperparameter marginals); for models with two or
more hyperparameters pass `model` and `y` to enable it, or it falls
back to the Gaussian approximation at the mode.

## Sampling

R-INLA's `inla.posterior.sample(n, res, ...)` becomes
`posterior_sample(rng, res, model; n_samples = n)`. The Julia version
returns a NamedTuple `(x, θ)` of `Matrix{Float64}`s — `x` is the
latent field samples (rows = latent dim, cols = samples), `θ` is the
hyperparameter samples.

```julia
using Random
rng = Xoshiro(20260506)
res = inla(model, y)
draws = posterior_sample(rng, res, model; n_samples = 1000)

# Posterior on a contrast of two random-effect levels:
contrast = draws.x[i, :] .- draws.x[j, :]
quantile(contrast, [0.025, 0.5, 0.975])
```

To sample at *new* observation points (R-INLA's `predict.inla`-style
workflow), use `posterior_predictive(rng, res, model, A_new)` with a
new observation mapping or matrix `A_new`. To draw `y_new` from the
likelihood given posterior `η`, use
`posterior_predictive_y(rng, res, model)`.

## Prediction

R-INLA exposes prediction through stacking `NA` rows into the design
matrix and reading them back from `res$summary.linear.predictor`.
Julia keeps the original fit clean and provides explicit prediction
endpoints.

| R-INLA workflow | Julia |
|---|---|
| `inla.stack(... NA rows ...)`, then read predictions | `posterior_predictive(rng, res, model, A_new)` |
| Project SPDE field to grid for plotting (`inla.mesh.project`) | `predict_raster(values, mesh, template)` (`INLASPDERasters.jl`) |
| Sample exceedance probability `P(η > θ)` | `predict_raster(...; samples)` with `Exceedance(threshold)` |
| Posterior quantile surfaces | `quantile_rasters(rng, res, model, mesh, template)` |
| Field at observation locations | `extract_at_mesh(raster, mesh)` |

Raster prediction has its own vignette
([Meuse SPDE](vignettes/meuse-spde.md)) walking through the full
INLA fit → posterior raster pipeline.

## SPDE-specific glue

R-INLA wraps the SPDE setup in `inla.spde2.matern` /
`inla.spde2.pcmatern` and the projector in `inla.spde.make.A`; Julia
exposes each step explicitly so users can mix mesh, projector, and
observation mapping freely.

| R-INLA step | Julia |
|---|---|
| `inla.mesh.2d(loc, max.edge = ...)` | `inla_mesh_2d(loc; max_edge = ..., cutoff = ...)` |
| `inla.mesh.1d(...)` | `inla_mesh_1d(...)` (post-Phase M) |
| `inla.spde2.matern(mesh, ...)` | `SPDE2(mesh; pc = PCMatern(...))` |
| `inla.spde.make.A(mesh, loc)` | `MeshProjector(mesh, loc)` |
| `inla.stack(data, A, effects)` | Build observation-mapping matrix `A = hcat(...)` and pass to `LatentGaussianModel(ℓ, components, A)` |

The mesh + SPDE component pieces compose into a `LatentGaussianModel`
the same way as any other component:

```julia
mesh = inla_mesh_2d(loc; max_edge = 0.1, cutoff = 0.05)
spde = SPDE2(mesh; pc = PCMatern(range_U = 0.5, range_α = 0.5,
                                  sigma_U = 1.0, sigma_α = 0.01))
A_field = MeshProjector(mesh, loc) |> as_matrix
A = hcat(ones(n_obs), x_cov, A_field)
model = LatentGaussianModel(GaussianLikelihood(),
                             (Intercept(), FixedEffects(1), spde), A)
res = inla(model, y)
```

## NUTS triangulation

R-INLA does not provide a built-in HMC sampler; users typically
re-fit in Stan or NIMBLE for a sanity check. Julia ships
[`LGMTuring.jl`](https://github.com/HaavardHvarnes/INLA.jl/tree/main/packages/LGMTuring.jl), which evaluates
[`LogDensityProblems`](https://github.com/tpapp/LogDensityProblems.jl)
on the same `LatentGaussianModel` and runs NUTS via AdvancedHMC over
the same hyperparameter posterior:

```julia
using LGMTuring
chain = nuts_sample(model, y, 1000;
                    n_adapts = 200, init_from_inla = inla(model, y))
diff = compare_posteriors(inla(model, y), chain; model = model)
```

`compare_posteriors` returns a side-by-side row per hyperparameter
with `flagged = true` when the posterior summaries disagree beyond
configurable tolerances. This is the load-bearing tier-3 validation
that ships nightly in CI.

## Things R-INLA does that we don't (yet)

The Julia ecosystem is a strict subset of R-INLA's surface. The
following land in v1.x or later:

- **Spatial models on the sphere.** R-INLA's
  `inla.spde2.matern(mesh, B.tau, B.kappa, alpha = 2)` on a spherical
  mesh ships through ADR-031; Julia tracks this in
  [`packages/INLASPDE.jl/plans/plan.md`](https://github.com/HaavardHvarnes/INLA.jl/blob/main/packages/INLASPDE.jl/plans/plan.md)
  as a v0.3+ item.
- **3D SPDE.** Brain-connectivity / subsurface workflows. Tracked in
  the same plan as spherical SPDE.
- **Non-separable space-time SPDE.** The separable case
  (`KroneckerComponent`) ships in Phase M; Lindgren et al.'s 2024
  non-separable generalisation is M+1 territory.
- **Fractional-α SPDE.** Bolin-Kirchner 2020 rational approximation
  for non-integer α. Deferred to v0.2.1+ (ADR-030).
- **`inla.qsample` raw GMRF sampling without conditioning on data.**
  Use [`GMRFs.jl`](packages/gmrfs.md) directly — sample from a GMRF
  via `rand(rng, gmrf)` after constructing the precision matrix.
- **`po` (probability ordinate).** A future Phase K addition.

If you hit a gap not in this list, file an issue at
[`HaavardHvarnes/INLA.jl/issues`](https://github.com/HaavardHvarnes/INLA.jl/issues)
— the v1.x backlog is shaped by user reports.
