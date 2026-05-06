# R-INLA → `@lgm` migration

`LGMFormula.jl` ships a Tier-2 formula macro `@lgm` that expands to the
explicit `LatentGaussianModel(...)` constructor of
[`LatentGaussianModels.jl`](packages/lgm.md). The macro is sugar — it
does no numerical work. Run `@macroexpand @lgm(...)` and you'll see
a constructor call that you could have written by hand.

This page is a side-by-side migration guide for users coming from
R-INLA's `inla(formula, family = ..., data = df)`.

## Anatomy of a formula

```julia
@lgm  y ~ 1 + x + f(idx, IID(n))   data = df   family = PoissonLikelihood()
#     ─┬─ ───┬───   ──────┬──────  ─────┬─────  ──────────┬──────────────
#      │     │            │              │                │
#      │     │            │              │                └ likelihood
#      │     │            │              └ Tables.jl source
#      │     │            └ random effect: column `idx` indexes into `IID(n)`
#      │     └ fixed effects: intercept + scalar covariate `x`
#      └ left-hand-side response column
```

R-INLA equivalent:

```r
inla(y ~ 1 + x + f(idx, model = "iid", hyper = list(prec = list(...))),
     family = "poisson", data = df)
```

The Julia macro mirrors R-INLA's `f(...)` random-effect syntax. The
component class moves from the `model = "iid"` keyword in R to a
constructor `IID(n)` in Julia — every component (`IID`, `AR1`, `RW1`,
`RW2`, `Besag`, `BYM2`, `Leroux`, `Seasonal`, …) is just a struct, so
the formula stays cleanly source-to-source.

## Installing

`LGMFormula.jl` is an optional sub-package; it is not pulled in by
`using INLA`. From a Julia REPL:

```julia
using Pkg
Pkg.develop(url = "https://github.com/HaavardHvarnes/INLA.jl",
            subdir = "packages/LGMFormula.jl")
```

Then:

```julia
using LatentGaussianModels, LGMFormula
```

## Side-by-side: Scotland lip cancer (BYM2)

The full vignette is at [Scotland BYM2](vignettes/scotland-bym2.md).
The R-INLA model:

```r
formula <- cases ~ 1 + x +
  f(area, model = "bym2", graph = W,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))
inla(formula, family = "poisson", E = expected, data = df)
```

In `@lgm`:

```julia
using GMRFs, LatentGaussianModels, LGMFormula

df = (cases = y, x = x, area = collect(1:n))
model = @lgm cases ~ 1 + x + f(area, BYM2(GMRFGraph(W); hyperprior_prec = PCPrecision(1.0, 0.01))) data=df family=PoissonLikelihood(; E = E)
res = inla(model, df.cases)
```

The macro builds the same projector the explicit-constructor vignette
shows — `hcat(ones(n, 1), reshape(x, n, 1), sparse(1:n, area, 1.0, n, 2n))`
— and threads the `BYM2` instance into the components tuple.

## Side-by-side: Tokyo rainfall (cyclic RW2)

The full vignette is at [Tokyo rainfall](vignettes/tokyo-rainfall.md).
The R-INLA model:

```r
formula <- y ~ -1 + f(day, model = "rw2", cyclic = TRUE,
                      hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))
inla(formula, family = "binomial", Ntrials = n_trials, data = df)
```

In `@lgm`:

```julia
df = (y = y, day = collect(1:366))
model = @lgm y ~ 0 + f(day, RW2(366; cyclic = true, hyperprior = PCPrecision(1.0, 0.01))) data=df family=BinomialLikelihood(n_trials)
res = inla(model, df.y)
```

`-1` and `0` both suppress the intercept, matching R-INLA's convention.

## Side-by-side: Meuse SPDE

The full vignette is at [Meuse SPDE](vignettes/meuse-spde.md). R-INLA's
SPDE formula uses `inla.spde.make.A` + `inla.stack` to bind point
coordinates to the mesh:

```r
spde   <- inla.spde2.pcmatern(mesh, prior.range = c(0.5, 0.5), prior.sigma = c(1, 0.5))
A      <- inla.spde.make.A(mesh = mesh, loc = coords)
stack  <- inla.stack(data = list(y = log(zinc)),
                     A = list(A, 1),
                     effects = list(field = 1:spde$n.spde,
                                    list(intercept = 1, dist = dist)))
formula <- y ~ -1 + intercept + dist + f(field, model = spde)
inla(formula, family = "gaussian", data = inla.stack.data(stack), ...)
```

`@lgm` does not yet support coordinate-indexed `f((s, t), SPDE2(mesh))`
forms — that is PR-7 of [Phase N](https://github.com/HaavardHvarnes/INLA.jl/blob/main/plans/phase-n.md)
and lands as `v0.2.2`. For now, build SPDE models with the explicit
constructor as shown in the [Meuse SPDE vignette](vignettes/meuse-spde.md).

The migration when PR-7 ships will look like:

```julia
# Future API — not yet available
model = @lgm logzinc ~ 1 + dist + f((east, north), SPDE2(mesh; pc = PCMatern(...))) data=df family=GaussianLikelihood()
```

## What `@lgm` covers today

| Feature | Syntax | Status |
|---|---|---|
| Intercept-only | `@lgm y ~ 1 …` | ✓ |
| Intercept + scalar covariates | `@lgm y ~ 1 + x1 + x2 …` | ✓ |
| No-intercept | `@lgm y ~ 0 + x …` (`-1` also accepted) | ✓ |
| Single random effect | `@lgm y ~ 1 + f(idx, IID(n)) …` | ✓ |
| Multiple random effects | `@lgm y ~ 1 + f(idx, IID(n)) + f(t, RW1(T)) …` | ✓ |
| Multi-likelihood (tuple-LHS) | `@lgm (y1, y2) ~ rhs … family=(L1, L2)` | ✓ |
| Replicated component | `f(t, AR1(n); replicate = id)` | ✓ |
| Grouped component | `f(t, AR1; group = grp)` | ✓ |
| `Copy(...)` augmentation | `f(...; copy = :name)` | PR-4b |
| Coordinate-indexed SPDE | `f((s, t), SPDE2(mesh))` | PR-7 |

## Restrictions

- Fixed-effects must be bare column symbols. Transformations
  (`log(x)`, interactions `x1*x2`, factor expansions) are not
  supported.
- The `col` of an `f(col, Component)` term must be a column of
  integers in `1:length(Component)` (R-INLA's "index column"
  convention).
- Multi-likelihood `(y1, y2) ~ rhs` is wide-format only — long-format
  with a `type` column is not yet supported.
- `replicate` and `group` are mutually exclusive within a single
  `f(...)` term.

If your model needs anything outside this list, drop back to the
explicit constructor — the formula sugar is strictly opt-in.

## When to use the explicit constructor instead

- Coordinate-indexed projectors (SPDE, mesh-barycentric).
- Custom `StackedMapping` row-partitions for long-format
  multi-likelihood data.
- Any model that needs `Copy(...)` source-component referencing today
  (PR-4b ships this).
- Models with non-trivial transformations on covariates.

The explicit constructor is and remains canonical; `@lgm` is a sugar
layer over it. See the [LatentGaussianModels.jl package
overview](packages/lgm.md) for the constructor API.

## See also

- [`LGMFormula.jl` package page](packages/lgmformula.md) — exported
  API.
- [Phase N plan](https://github.com/HaavardHvarnes/INLA.jl/blob/main/plans/phase-n.md) — the
  PR sequence behind this package.
- [`plans/macro-policy.md`](https://github.com/HaavardHvarnes/INLA.jl/blob/main/plans/macro-policy.md)
  — design constraints on every macro added to this ecosystem.
