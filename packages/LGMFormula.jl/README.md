# LGMFormula.jl

Tier-2 formula sugar (`@lgm`) for [`LatentGaussianModels.jl`](../LatentGaussianModels.jl/).

This package is **optional**. Core `LatentGaussianModels.jl` ships a
complete explicit-constructor API (Tier 1); `LGMFormula.jl` adds a
macro that expands to those constructor calls (Tier 2). Run
`@macroexpand @lgm(...)` and you'll see a `LatentGaussianModel(...)`
constructor call you could have written by hand.

See [`plans/macro-policy.md`](../../plans/macro-policy.md) and
[ADR-008](../../plans/decisions.md#adr-008-lgm-macro-lives-in-a-separate-lgmformulajl-package)
for the design constraints.

## Status

`v0.2.0`. Phase N PRs 1–6 closed:

- PR-1: core parser + single-likelihood single-`f` expansion.
- PR-2: component coverage roundtrips (~20 components).
- PR-3: multi-`f` roundtrip coverage.
- PR-4: tuple-LHS multi-likelihood (wide-format).
- PR-5: `replicate` / `group` term routing.
- PR-6: migration guide + vignette parity.

PR-4b (`Copy` augmentation) and PR-7 (SPDE-friendly coordinate forms)
ship in a follow-up.

## Quick example

```julia
using GMRFs, LatentGaussianModels, LGMFormula

df = (y = y, x = x, region = collect(1:n))
model = @lgm y ~ 1 + x + f(region, IID(n)) data=df family=PoissonLikelihood()
res   = inla(model, df.y)
```

The macro builds the projector `A` from the data columns and threads
the components tuple and likelihood through to the constructor. See
the [migration guide](https://haavardhvarnes.github.io/INLA.jl/stable/lgmformula-tutorial/)
for the full R-INLA → `@lgm` correspondence on Scotland BYM2 and
Tokyo rainfall.

## Installation

Not yet on the personal Julia registry. From a fresh Julia REPL:

```julia
using Pkg
Pkg.develop(url = "https://github.com/HaavardHvarnes/INLA.jl",
            subdir = "packages/LGMFormula.jl")
```

## Why a separate package?

- Julia 1.9+ extensions cannot export new symbols. `@lgm` needs to be
  an exported macro, so it cannot live in a weakdep of
  `LatentGaussianModels`.
- StatsModels and its dep tree (StatsBase / StatsFuns / DataAPI /
  ShiftedArrays) is non-trivial. Users who prefer the explicit Tier-1
  API should not pay for it.

## See also

- [`LatentGaussianModels.jl`](../LatentGaussianModels.jl/) — the core
  LGM package. This sub-package is useless without it.
- [Phase N plan](../../plans/phase-n.md) — the PR sequence behind
  this package.
