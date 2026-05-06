# Guidance for Claude Code in LGMFormula.jl

Extends [`/CLAUDE.md`](../../CLAUDE.md). Scoped to the formula macro.

## Scope

This package owns:
- The `@lgm` macro and its supporting parsing logic.
- The `lgmformula(data; lhs, intercept, covariates, randoms, family)`
  function form the macro lowers to.
- Schema handling: binding `region` in `f(region, Besag(W))` to a
  column of a `Tables.jl`-compatible source via runtime helpers in
  `schema.jl`.
- Per-term sparse projector block construction — pure index-shuffling
  on top of `Tables.getcolumn`; no GMRF or LGM arithmetic.

Out of scope:
- Anything numerical. This package never does arithmetic beyond what
  `sparse(rows, cols, vals, m, n)` does for the design matrix.
- Component construction. `Besag(W)` etc. are constructed by
  `LatentGaussianModels`; the macro just threads arguments through.

## Design constraints

- **Tier 1 completeness.** Every expansion of `@lgm(...)` must produce
  a `LatentGaussianModel(...)` constructor call that a user could have
  written by hand. Run `@macroexpand` and it must be readable.
- **No runtime semantics.** The macro is purely a source-to-source
  rewrite. It does not introduce any state or behavior that the
  explicit constructor lacks.
- **Errors refer to user concepts.** An error from a malformed
  `f(region, Besag(W))` term should say "unknown column `region` in
  data" or "component `Besag(W)` is not a valid component", not
  something about `FunctionTerm` internals.

## Dependencies allowed

Core:
- `LatentGaussianModels` — the host.
- `StatsModels` — `@formula` parsing helpers.
- `Tables` — data-source interface.
- `SparseArrays` (stdlib) — design-matrix construction.

Nothing else without an ADR.

## Testing

All tests live under `test/regression/` (parsing is deterministic, so
no oracle / triangulation tiers):
- `test_macroexpand.jl` — `@macroexpand @lgm(...)` matches frozen AST
  shapes.
- `test_roundtrip.jl` — `model == handwritten_constructor` via struct
  equality (`_struct_isequal` helper in `test_utils.jl`).
- `test_error_messages.jl` — malformed input errors refer to user
  concepts.
- `test_components.jl` — per-component coverage (PR-2).
- `test_multi_f.jl` — multi-`f` shapes (PR-3).
- `test_multi_likelihood.jl` — tuple-LHS multi-likelihood (PR-4).
- `test_replicate_group.jl` — `replicate` / `group` routing (PR-5).

## Review criteria for a new macro

See [`plans/macro-policy.md`](../../plans/macro-policy.md) "Review
criteria." All of those apply to every macro added to this package.
