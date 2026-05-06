# LGMFormula.jl

Tier-2 formula sugar for [`LatentGaussianModels.jl`](lgm.md). Exposes a
single macro `@lgm` and its function form `lgmformula` that lower a
formula expression bound to a `Tables.jl`-compatible source into an
explicit `LatentGaussianModel(...)` constructor call.

The macro is **strictly source-to-source** per
[ADR-008](https://github.com/HaavardHvarnes/INLA.jl/blob/main/plans/decisions.md)
and [`plans/macro-policy.md`](https://github.com/HaavardHvarnes/INLA.jl/blob/main/plans/macro-policy.md).
Every expansion produces a constructor call that the user could have
written by hand. Run `@macroexpand @lgm(...)` to inspect.

## When to use this

- Migrating models written against R-INLA's `inla(formula, …)` API.
- Wanting concise notation for standard R-INLA-style models with
  named index columns.

For a guided introduction, see the
[migration guide](../lgmformula-tutorial.md).

## Quick example

```julia
using GMRFs, LatentGaussianModels, LGMFormula

df = (y = y, x = x, region = collect(1:n))
model = @lgm y ~ 1 + x + f(region, IID(n)) data=df family=PoissonLikelihood()
res = inla(model, df.y)
```

## Status

`v0.2.0`. Phase N PRs 1–6 closed:

- PR-1: core parser + single-likelihood single-`f` expansion.
- PR-2: component coverage roundtrips.
- PR-3: multi-`f` roundtrip coverage.
- PR-4: tuple-LHS multi-likelihood (wide-format).
- PR-5: `replicate` / `group` term routing.
- PR-6: migration guide + vignette parity (this page).

PR-4b (`Copy` augmentation) and PR-7 (SPDE-friendly coordinate forms)
ship in a follow-up.

## API

```@autodocs
Modules = [LGMFormula]
```
