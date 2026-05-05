"""
    LGMFormula

Tier-2 formula sugar (`@lgm`) for `LatentGaussianModels.jl`.

The macro is a source-to-source transform: every expansion of
`@lgm(...)` produces a `LatentGaussianModel(...)` constructor call
that the user could have written by hand. Run `@macroexpand
@lgm(...)` to inspect.

# Public API

- [`@lgm`](@ref) — formula macro.
- [`lgmformula`](@ref) — function form; the macro lowers to this.

See [`plans/macro-policy.md`](https://github.com/haavardhvarnes/INLA.jl/blob/main/plans/macro-policy.md)
and [ADR-008](https://github.com/haavardhvarnes/INLA.jl/blob/main/plans/decisions.md).
"""
module LGMFormula

using LatentGaussianModels
using SparseArrays
using Tables

export @lgm, lgmformula

include("parse.jl")
include("schema.jl")
include("expand.jl")

"""
    @lgm formula data=df family=Likelihood()

Build a [`LatentGaussianModel`](@ref) from a formula expression bound
to a `Tables.jl`-compatible source.

# Supported (PR-1)

- `@lgm y ~ 1 data=df family=GaussianLikelihood()` — intercept only.
- `@lgm y ~ 1 + x1 + x2 data=df family=GaussianLikelihood()` —
  intercept + scalar covariates.
- `@lgm y ~ 0 + x data=df family=GaussianLikelihood()` —
  no intercept (`-1` also accepted).
- `@lgm y ~ 1 + f(idx, IID(n)) data=df family=PoissonLikelihood()` —
  intercept + one random effect.
- `@lgm y ~ 1 + x + f(idx, IID(n)) data=df family=PoissonLikelihood()` —
  intercept + covariate + random effect.

# Restrictions (PR-1)

- The LHS must be a column name (`Symbol`).
- Fixed-effects terms must be bare column symbols. Transformations
  (`log(x)`, `x1*x2`, factor expansions) are not yet supported.
- A single `f(col, Component)` per formula. The `col` must be a
  column of integers in `1:length(Component)`. Multiple `f(...)` terms
  ship in PR-3.
- Single likelihood. Multi-likelihood (`(y1, y2) ~ ...`) ships in PR-4.

# Expansion

The macro expands to an explicit `LatentGaussianModel(...)` call with
a [`lgmformula`](@ref)-built design matrix; run `@macroexpand` to
inspect. The components tuple and likelihood appear literally in the
expansion; only the design matrix construction is deferred to runtime
(it depends on the data).
"""
macro lgm(args...)
    formula_expr, opts = _parse_args(args)
    haskey(opts, :data) ||
        error("@lgm: missing required keyword `data = ...`")
    haskey(opts, :family) ||
        error("@lgm: missing required keyword `family = ...`")
    lhs, rhs = _parse_formula(formula_expr)
    has_intercept, covariates, randoms = _split_rhs(rhs)
    return esc(_build_expansion(lhs, has_intercept, covariates, randoms,
        opts[:data], opts[:family]))
end

end # module
