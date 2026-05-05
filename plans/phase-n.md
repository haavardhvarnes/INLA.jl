# Phase N — `LGMFormula.jl` maturity

## Context

`LGMFormula.jl` is the Tier-2 formula sugar deferred from v0.1 per
[ADR-015](decisions.md). The package was sized at 4–6 weeks
when [`replan-2026-04-28.md:408-431`](replan-2026-04-28.md) was
written (only Phases F–H of LGM had shipped). With Phases F–M now in
v0.2.0, the formula DSL must cover a substantially larger Tier-1
surface: multi-likelihood + `Copy` + `StackedMapping` (Phase G), seven
new likelihoods (Phase H + J), `IIDND`/`MEB`/`MEC`/`Replicate`/`Group`
(Phase I), six new posterior accessors that the formula never sees but
the docs must explain (Phase K), `UserComponent` (Phase L), and 1D /
non-stationary / space-time SPDE components (Phase M).

The macro itself stays source-to-source per
[`plans/macro-policy.md`](macro-policy.md). Tier 1 remains
canonical; `@lgm` is sugar that expands to a `LatentGaussianModel(...)`
constructor call with `@macroexpand`-readable output.

## State of the world (read-only audit)

`packages/LGMFormula.jl/` today: `Project.toml` (0.1.0-DEV,
`LatentGaussianModels = "0.1"` compat), `CLAUDE.md`, `README.md`,
`plans/plan.md`, `LICENSE`, `.gitignore`. **No `src/`, no `test/`.**
The package currently:
- Cannot be `using`'d (no module file).
- Has compat pinned to LGM 0.1 — needs bump to 0.2 alongside any
  implementation.
- Ships an M1 / M2 / M3 milestone breakdown in `plans/plan.md` that
  predates Phases G–M and is now scope-shy.

The Tier-1 surface the macro must cover, after Phase M close at v0.2.0:

| Tier-1 feature | Phase | `@lgm` requirement |
|---|---|---|
| `LatentGaussianModel(ℓ, components, mapping)` 4-arg | F | Always-present |
| `LikelihoodMap` + multi-likelihood tuple | G | `@lgm (y1, y2) ~ ...` joint syntax |
| `StackedMapping` per-likelihood blocks | G | Built from `+ f(...)` term grouping |
| `Copy(source_indices; β, hyper_β)` | G | `f(...; copy = :name)` term |
| `Replicate(component, n)` | I | `f(col; replicate = rep_col)` kwarg |
| `Group(component, group_id)` | I | `f(col; group = group_col)` kwarg |
| Twenty likelihoods (Gaussian / Poisson / Binomial / NegBin / Gamma / Beta / BetaBinomial / StudentT / SkewNormal / GEV / POM / Weibull{PH,AFT} / Exponential / LognormalSurv / GammaSurv / Coxph / WeibullCure / ZIP{0,1,2} / ZINegBin{0,1,2}) | F+H+J | `family = ...` kwarg |
| Twenty components (Intercept / FixedEffects / IID / IIDND / RW1 / RW2 / AR1 / Seasonal / Besag / BYM / BYM2 / Leroux / Generic0 / Generic1 / Generic2 / MEB / MEC / SPDE1D / SPDE2 / SPDE2NonStationary / KroneckerComponent / UserComponent) | F-M | Pass-through `f(col, Component(...))` |
| Spatial `f(coords, SPDE2(mesh))` | M | Coords as a 2-column term |

R-INLA's `inla(formula, family = ..., data = ...)` is the calling
convention to match where it doesn't conflict with Julia idioms. The
predecessor [`IntegratedNestedLaplace.jl`](https://github.com/haavardhvarnes/IntegratedNestedLaplace.jl)
has a `_build_latent_effect` reference that documents the field-by-
field translation; port the pattern, not the file.

## PR sequence

Six PRs + one stretch tail, similar shape to Phase M but lighter on
numerical content. No new oracle fixtures — the contract is "model
`isequal` to the hand-written Tier-1 form."

### PR-1 — Core parser + single-likelihood single-`f` expansion

Stand up `LGMFormula.jl/src/`:
- `LGMFormula.jl` — module entry; exports `@lgm`, `lgmformula`.
- `parse.jl` — walk the `Expr` tree, separate the LHS, the RHS
  fixed-effects part, and the `f(...)` calls.
- `schema.jl` — bind `f(col, Component(...))` to a `Tables.jl` column;
  emit a clear error for unknown columns.
- `expand.jl` — rebuild as `LatentGaussianModel(ℓ, (...), mapping)`.

Coverage:
- `@lgm y ~ 1 + x1 + x2 data=df family=Gaussian()` → fixed effects only.
- `@lgm y ~ 1 + x + f(idx, IID(n)) data=df family=Poisson()` → one
  random effect.
- `lgmformula(formula::FormulaTerm, df; family, ...)` function form.

**Tests** (`test/regression/test_macroexpand.jl`,
`test/regression/test_roundtrip.jl`): `@macroexpand` against a frozen
AST snapshot for ~6 minimal cases; `model_lgm == model_explicit` via
`isequal` for those 6.

### PR-2 — Component coverage tests

No source change. For each of the ~20 concrete components, add a
roundtrip test asserting that an `@lgm` formula expressing that
component produces an `isequal` model to the hand-written form. The
test surface is the implicit acceptance gate on PR-1's expansion
machinery: if a component constructor needs a kwarg the parser doesn't
forward, this PR catches it.

Half a day of work; ~150 LOC of test scaffolding. Lands as the
quality gate before PR-3 enlarges the parser.

### PR-3 — Multi-`f` roundtrip coverage

**Scope correction (recorded 2026-05-05).** This PR's original scope
above envisaged a per-component `StackedMapping` synthesis with
`Component => Mapping` pairs. That pseudocode does not fit LGM core's
actual `StackedMapping` API: `StackedMapping(blocks, rows)` is a
**row-partition** structure for multi-likelihood (`inla.stack`-style
disjoint observation slices sharing a column-stacked latent), not a
per-component column-partition. There is no LGM core mapping type that
keys per-component mappings — and per
[`LatentGaussianModels.jl/CLAUDE.md`](../packages/LatentGaussianModels.jl/CLAUDE.md)
LGMFormula must consume the existing API, not extend it.

**General guidance** (applies to all Phase N PRs): LGM core changes
are last-resort. If a Phase N PR's plan seems to need an LGM core API
addition, that's a replan signal — re-scope the PR around what
LGMFormula can deliver alone, or escalate as a separate Tier-1
follow-up.

PR-1's `_build_design_matrix` already handles multi-`f` correctly: it
loops over every `f(...)` term, pushes a per-term `sparse(1:n_obs,
idx, 1.0, n_obs, length(comp))` block, and returns
`reduce(hcat, blocks)`. The LGM constructor wraps that single
`SparseMatrixCSC` as `LinearProjector(A)` — the v0.2 default per
ADR-017. A multi-`f` formula like `y ~ 1 + f(region, BYM2(W)) +
f(time, AR1(T))` therefore lowers to:

```julia
LatentGaussianModel(
    family,
    (Intercept(), BYM2(W), AR1(T)),
    A,        # = hcat(intercept_block, A_region, A_time) :: SparseMatrixCSC
)
```

**Revised PR-3 scope: test-only coverage.** No source change. New
file `test/regression/test_multi_f.jl` adds:
- Synthetic two-effect roundtrips (BYM2 + AR1, Besag + IID, IID +
  RW1).
- Synthetic three-effect roundtrip (IID + RW1 + Seasonal).
- Mixed covariate + multi-`f` (Gaussian, several covariates plus two
  random effects).
- No-intercept multi-`f` (`0 + f(...) + f(...)`).
- Scotland BYM2 *shape* roundtrip (Poisson + intercept +
  `f(area, BYM2(W))` on a small synthetic adjacency).
- Tokyo RW2 *shape* roundtrip (Binomial + `f(day, RW2(...; cyclic=true))`
  on a tiny synthetic dataset).

These are roundtrip equality tests against the hand-written Tier-1
form — no full inference. Full disease-mapping oracle parity is
already covered in `LatentGaussianModels.jl/test/oracle/`; PR-3's job
is to prove the macro produces an `_struct_isequal` model for these
canonical shapes.

### PR-4 — Multi-likelihood (tuple-LHS)

**Scope split (recorded 2026-05-05).** PR-4's original scope bundled
multi-likelihood and `Copy` into one PR. Both are real expansion-
machinery work — multi-likelihood needs parser + schema + expand
changes for tuple-LHS, and `Copy` needs a syntax-level decision on
naming source components in the formula plus `CopyTargetLikelihood`
plumbing. Splitting keeps each PR scoped to one concept and gives
`Copy` its own ADR. `Copy` ships as PR-4b below.

`@lgm (y_long, y_surv) ~ rhs data=df family=(Gaussian(), Weibull())`
lowers to:

```julia
LatentGaussianModel(
    (Gaussian(), Weibull()),
    (Intercept(), ...),
    StackedMapping(
        (LinearProjector(A), LinearProjector(A)),
        [1:n, (n+1):(2n)],
    ),
)
```

The shared RHS is built once into a single sparse `A` (length-n rows);
each likelihood block wraps the same projector. Observations are
expected in wide-format — each LHS column has length `n`, the stacked
observation vector is `vcat(y1, y2, ...)`. Long-format with a `type`
column is left for a follow-up.

**ADR-033** (PR-4): multi-likelihood formula syntax. Candidates were
tuple-LHS `(y1, y2) ~ ...` vs. R-INLA's per-row `family` plus
`Y ~ stack(...)`. Recommendation: ship tuple-LHS as Julia-idiomatic;
document R-INLA's stacked-Y form as a fallback users can write
explicitly.

**Per-likelihood RHS variation** (e.g. Baghfalaki joint
longitudinal-survival, where the survival predictor differs from the
longitudinal predictor) is not in PR-4 — it requires PR-4b's `Copy`
plus a syntax for "this f-term applies only to likelihood k", which
is its own design call.

**Tests**: tuple-LHS roundtrip on `(y1, y2) ~ 1 + f(idx, IID(n))`,
`(y_g, y_p) ~ 1 + f(idx, IID(n))` matching the
`test_synthetic_joint_gauss_pois` form, and a 3-likelihood case.
Macroexpand structural assertions on the tuple-LHS expansion. Error-
message tests for mismatched lengths and missing LHS columns.

### PR-4b — `Copy` augmentation

`f(...; copy = :name)` resolves `:name` against named f-terms in the
formula and emits a `Copy(target_indices)` augmentation, wrapping the
target likelihood with `CopyTargetLikelihood(base, copy_spec)`. The
LGM core API (`packages/LatentGaussianModels.jl/src/likelihoods/copy.jl`)
ships in v0.2: `Copy(source_indices; β_prior, β_init, fixed)` and
`CopyTargetLikelihood(base, copies::Tuple)`.

**ADR-034** (PR-4b): component naming syntax in the formula. Two
candidates: (a) implicit naming via the column symbol (`f(subject,
IID(n))` is named `:subject`), (b) explicit `name = :foo` kwarg. Lean
(a); the column symbol is already unique per f-term in the cases we
care about (R-INLA's `inla.stack` makes the same assumption). Per-
likelihood targeting (which Copy goes into which likelihood's
predictor) needs its own kwarg — `f(subject; copy = :u, target = 2)`
or similar.

**Tests**: roundtrip on the Phase G Baghfalaki joint longitudinal-
survival fit shape. R-INLA's `inla.stack` for that fit is the ground
truth.

### PR-5 — `replicate` / `group` routing

`f(col; replicate = rep_col)` → `Replicate(component, n_levels(rep_col))`.
`f(col; group = group_col)` → `Group(component, group_col_vec)`.

This is the one place where the macro looks at `data` columns at
expansion time (to compute level counts and group ids); document why
in the ADR draft.

**Tests**: roundtrip on the Phase I synthetic_replicate_ar1 fixture
(longitudinal panel, replicated AR1).

### PR-6 — Migration guide + vignette parity

`docs/src/lgmformula-tutorial.md` — R-INLA → `@lgm` migration with
side-by-side examples for the three flagship vignettes (Scotland BYM2,
Tokyo rainfall, Meuse SPDE). Add the `@lgm` versions of these
vignettes as a second tab in their existing pages — the explicit
constructor stays canonical; the formula version is the migrant
on-ramp.

Add `getting-started.md` callout: "Coming from R-INLA?" → links to
`@lgm`. Update `README.md` to mention LGMFormula in the "What ships"
list.

**No tests** beyond the docs build.

### PR-7 (stretch) — `KroneckerComponent` and SPDE-friendly forms

`@lgm y ~ 1 + f((s, t), KroneckerComponent(SPDE2(mesh), AR1()))` —
the comma-separated coordinate form. Two design calls, surfaced in
PR-7's ADR draft:
- **Tuple-coordinate syntax** `f((s, t), ...)` vs. R-INLA's
  `f(s, ..., group = t, control.group = list(model = "ar1"))`.
- **Mesh threading**: the `SPDE2(mesh)` constructor takes a mesh
  argument that the formula has to thread through unchanged. Pure
  pass-through via Julia's expression interpolation `$(mesh)`.

**Stretch criterion**: ship if PRs 1–6 close inside week 5; defer to
Phase N+1 if time-pressed (mirrors Phase M's PR-7 deferral pattern,
this time on grounds of bandwidth — KroneckerComponent's mesh-binding
ergonomics need a separate ADR, and the explicit constructor is a
viable workaround for users who need it).

## Out of scope for Phase N

- `cgeneric` C-callable bindings — explicit out of scope per
  [`packages/LGMFormula.jl/CLAUDE.md`](../packages/LGMFormula.jl/CLAUDE.md)
  and [`plans/macro-policy.md`](macro-policy.md).
- Interaction terms `x1:x2` — defer until a user asks.
- DataFrames-specific fast paths. Tables.jl is the contract.
- Any numerical code. The macro is purely a source-to-source rewrite
  per ADR-008.

## Critical files

| Concern | Path |
|---|---|
| Module entry (new) | `packages/LGMFormula.jl/src/LGMFormula.jl` |
| Parser (new) | `packages/LGMFormula.jl/src/parse.jl` |
| Schema binding (new) | `packages/LGMFormula.jl/src/schema.jl` |
| AST expansion (new) | `packages/LGMFormula.jl/src/expand.jl` |
| StatsModels integration (new) | `packages/LGMFormula.jl/src/terms.jl` |
| Tier-1 constructor reference | [`packages/LatentGaussianModels.jl/src/model.jl:54-150`](../packages/LatentGaussianModels.jl/src/model.jl) |
| Observation mapping types | [`packages/LatentGaussianModels.jl/src/observation_mapping.jl`](../packages/LatentGaussianModels.jl/src/observation_mapping.jl) |
| `Copy` source | [`packages/LatentGaussianModels.jl/src/likelihoods/copy.jl`](../packages/LatentGaussianModels.jl/src/likelihoods/copy.jl) |
| Replicate / Group source | [`packages/LatentGaussianModels.jl/src/components/{replicate,group}.jl`](../packages/LatentGaussianModels.jl/src/components/replicate.jl) |
| Predecessor reference | `_build_latent_effect` in `IntegratedNestedLaplace.jl` (port pattern, not file) |
| LGMFormula plan ledger | [`packages/LGMFormula.jl/plans/plan.md`](../packages/LGMFormula.jl/plans/plan.md) |
| Macro policy | [`plans/macro-policy.md`](macro-policy.md) |
| ADR registry | [`plans/decisions.md`](decisions.md) — append ADR-033 (multi-likelihood syntax), optional ADR-034 (KroneckerComponent ergonomics) |

## Verification

Phase close gates:
- `@macroexpand @lgm(...)` against frozen AST snapshots for ~12
  representative cases (one per component family).
- Roundtrip `isequal` on the three flagship vignettes
  (Scotland BYM2, Tokyo rainfall, Meuse SPDE) plus the Baghfalaki
  joint and the synthetic replicate fit.
- All Phase F–M oracle fixtures still pass (no Tier-1 source changes —
  this is just a sanity check that the LGMFormula compat bumps didn't
  pull in something incompatible).
- `JET.@report_opt` clean on the parser hot path.
- Aqua check clean.
- Migration guide published; three vignettes have a `@lgm` section
  alongside the explicit-constructor section.

## Release target

Phase N closes as **`v0.2.1`** — `LGMFormula.jl` ships at `v0.2.0`
(first release of the package), umbrella `INLA.jl` bumps to `v0.2.1`
and gains LGMFormula as a `[deps]` (or as a `[weakdeps]` re-export —
design call in PR-1's ADR draft, leaning toward `[deps]` since
LGMFormula is small and the umbrella is the curated convenience
surface).

PR-7 (KroneckerComponent stretch) is the only candidate for promoting
to `v0.2.2` if shipped post-PR-6.

Tag at phase close: `v0.2.1` on `main`, with a CHANGELOG entry
covering all six (or seven) PRs, and a release commit
`chore(release): v0.2.1 — Phase N close (LGMFormula maturity)`
mirroring the v0.2.0 close pattern.

## Estimated cadence

At observed solo-developer pace (~1 PR/week through Phases L and M):

- Week 1: PR-1 core parser
- Week 2: PR-2 component coverage tests + PR-3 multi-`f` (combined,
  PR-2 is a half-day)
- Week 3: PR-4 multi-likelihood + Copy
- Week 4: PR-5 replicate/group
- Week 5: PR-6 migration guide
- Stretch / Phase N+1: PR-7 KroneckerComponent ergonomics

Replan estimate was 4–6 weeks; this plan targets 5 weeks with PR-7 as
slack. Phase N is lighter than Phase M (no FEM, no oracle generation,
no R cross-runs) — calendar risk is mostly in PR-3's mapping-synthesis
correctness.
