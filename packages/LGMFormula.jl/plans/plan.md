# LGMFormula.jl — package plan

> **Status (2026-05-06):** v0.2.0. Phase N PRs 1–6 shipped — the macro
> is feature-complete for the documented Tier-1 surface. PR-4b
> (`Copy` augmentation) and PR-7 (SPDE-friendly coordinate forms)
> remain as follow-ups. See [`plans/phase-n.md`](../../../plans/phase-n.md)
> for the PR ledger.

## Goal

Provide the Tier-2 `@lgm` formula sugar on top of
[`LatentGaussianModels.jl`](../../LatentGaussianModels.jl/). A
source-to-source transform that expands to an explicit
`LatentGaussianModel(...)` constructor call. No numerical code.

## Module layout

```
src/
├── LGMFormula.jl   # main module, exports @lgm and lgmformula
├── parse.jl        # Expr-tree walker: LHS/RHS, fixed effects, f(...) terms
├── schema.jl       # Tables.jl binding; per-term projector blocks; runtime wrap helper
└── expand.jl       # build the LatentGaussianModel(...) AST

test/
├── runtests.jl
├── test_utils.jl
└── regression/
    ├── test_macroexpand.jl       # @macroexpand vs frozen AST snapshots
    ├── test_roundtrip.jl         # model == handwritten Tier-1 (struct equality)
    ├── test_error_messages.jl    # malformed input errors are user-readable
    ├── test_components.jl        # per-component coverage roundtrips (PR-2)
    ├── test_multi_f.jl           # multi-`f` shapes (PR-3)
    ├── test_multi_likelihood.jl  # tuple-LHS multi-likelihood (PR-4)
    └── test_replicate_group.jl   # replicate / group routing (PR-5)
```

## Shipped (Phase N PR-1..PR-6)

- `@lgm y ~ rhs data=df family=L()` parses LHS, fixed-effect terms,
  and `f(col, Component(...))` terms.
- `lgmformula(data; lhs, intercept, covariates, randoms, family)`
  is the function form the macro lowers to.
- Multi-likelihood tuple-LHS: `(y1, y2) ~ rhs ... family=(L1, L2)`
  (wide-format).
- `f(col, Component(n); replicate = id_col)` → runtime
  `Replicate(component, R)` with `R = maximum(data.id_col)`.
- `f(col, Factory; group = grp_col)` → runtime
  `Group(Factory, data.grp_col)`.
- Errors refer to user concepts: column names, kwarg names,
  component class names. No `FunctionTerm`/`Expr` internals leak.
- Migration guide + Scotland/Tokyo `@lgm` vignette sections + a Meuse
  pointer to PR-7.

## Open work (post-PR-6)

- **PR-4b** — `Copy` augmentation: `f(...; copy = :name)` resolves to
  `Copy(target_indices)` plus `CopyTargetLikelihood(base, copies)`.
  Needs ADR-034 on component-naming syntax.
- **PR-7 (stretch)** — SPDE-friendly coordinate forms:
  `f((s, t), KroneckerComponent(SPDE2(mesh), AR1()))` and the
  tuple-coordinate parser path. Mesh threaded through unchanged via
  expression interpolation.
- **Long-format multi-likelihood** — `(y, type)` style LHS where rows
  are stacked and `type` selects the likelihood per row. Deferred to
  v0.3+; wide-format covers the canonical R-INLA `inla.stack` flows.

## Out of scope

- Any Bayesian inference. `LGMFormula.jl` returns a model; fitting
  uses `LatentGaussianModels.jl`.
- DataFrames-specific fast paths. Tables.jl is the contract.
- Interaction terms `x1:x2` and covariate transformations
  (`log(x)`, factor expansions) — defer until a user asks.
- `cgeneric` C-callable bindings — explicit out of scope per
  [`plans/macro-policy.md`](../../../plans/macro-policy.md).

## Validation

- `test/regression/` covers `@macroexpand` AST shape, model roundtrip
  equality against hand-written Tier-1 forms, and user-readable error
  messages.
- No oracle / triangulation tier — parsing is deterministic and the
  numerical content is whatever `LatentGaussianModels.jl` already
  validates.
