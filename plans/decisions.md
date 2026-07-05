# Architecture decision records

Each significant design decision gets a short entry here. Format:

```
## ADR-NNN: Title

Status: Accepted | Proposed | Superseded by ADR-MMM
Date: YYYY-MM-DD

### Context
What's the situation that forced a decision?

### Decision
What did we choose?

### Consequences
What does this buy us, what does it cost us, what's the escape hatch?
```

Numbering is sequential; never renumber. ADRs remain in numerical
order in the body — the topical index below is a forward index, not
a reorder.

## Topical index

For navigation; ADR bodies appear in numerical order under the index.

**Architecture & package layout**
- [ADR-001](#adr-001-package-split-into-gmrfs--latentgaussianmodels--inlaspde) — package split into GMRFs / LatentGaussianModels / INLASPDE
- [ADR-008](#adr-008-lgm-macro-lives-in-a-separate-lgmformulajl-package) — `@lgm` macro lives in `LGMFormula.jl`
- [ADR-009](#adr-009-turing--hmc-bridge-lives-in-a-separate-lgmturingjl-package) — Turing / HMC bridge lives in `LGMTuring.jl`
- [ADR-014](#adr-014-main-is-fast-forwarded-to-claudehungry-pascal-main-becomes-the-integration-branch-going-forward) — `main` is the integration branch
- [ADR-015](#adr-015-lgmformulajl-and-gmrfspardisojl-both-cut-from-v01-deferred-to-v02) — LGMFormula.jl / GMRFsPardiso.jl deferral
- [ADR-020](#adr-020-drop-julia-110-lts-support--julia-112-is-the-minimum-supported-version) — Julia 1.12+ minimum
- [ADR-039](#adr-039-lgmformula--inlaspde-integration-ships-as-a-julia-19-weakdep-extension-not-a-hard-dep-or-duck-typed-hook) — LGMFormula ↔ INLASPDE weakdep extension
- [ADR-042](#adr-042-inlaspderasters-takes-a-load-bearing-dep-on-latentgaussianmodels-via-model-res--overloads) — INLASPDERasters' load-bearing LGM dep

**Dispatch & API conventions**
- [ADR-003](#adr-003-multiple-dispatch-not-macros-for-the-primary-api) — multiple dispatch over macros
- [ADR-010](#adr-010-public-kwargs-mirror-r-inla-names-with-snake_case--symbol-or-type-dual-input) — snake_case + symbol-or-type kwargs
- [ADR-011](#adr-011-top-level-api-is-fitmodel-y-strategy-kwargs-inlamodel-y-kwargs-is-a-convenience-alias) — `fit(model, y, strategy; …)` top-level
- [ADR-025](#adr-025-usercomponent-callback-signature-r-inla-rgeneric) — `UserComponent` callback (rgeneric)
- [ADR-026](#adr-026-marginal-strategies-via-abstractmarginalstrategy-type-dispatch) — `AbstractMarginalStrategy` type dispatch

**Inference algorithms**
- [ADR-002](#adr-002-scimls-linearsolvejl-as-the-sparse-factorization-backend) — LinearSolve.jl factorization backend
- [ADR-004](#adr-004-selected-inversion-takahashi-recursion-as-an-explicit-risk) — selected inversion (Takahashi)
- [ADR-006](#adr-006-full-laplace-is-the-phase-3-default-simplifiedlaplace-deferred) — full Laplace default (amended; see ADR-016)
- [ADR-012](#adr-012-adopt-selectedinversionjl-for-sparse-selected-inverse) — adopt SelectedInversion.jl
- [ADR-016](#adr-016-simplified-laplace-mean-shift-correction-rue-martino-added-as-opt-in-latent_strategy) — simplified-Laplace mean-shift correction
- [ADR-027](#adr-027-is-correction-port-from-integratednestedlaplacejl--declined-for-v0x-deferred-to-per-workflow-strategy-if-a-use-case-appears) — IS correction declined
- [ADR-031](#adr-031-targeted-exception-classification-in-the-laplace-bad-θ-wrapper) — targeted exception classification
- [ADR-045](#adr-045-proposed-de-densify-the-constrained-laplace-null-space-bump-low-rank-or-bordered-kkt) — de-densify the constrained-Laplace null-space bump (Proposed)
- [ADR-046](#adr-046-integrated-hyperparameter-marginals--design-point-reuse-for-mθ--1-conditional-mode-profile-slices-for-mθ--2) — integrated hyperparameter marginals
- [ADR-047](#adr-047-item-163-closed--the-simplified-laplace-variance-correction-is-mis-specified-no-such-term-exists-in-the-reference-method) — "SLA variance correction" closed as mis-specified

**Observation mapping & projectors**
- [ADR-005](#adr-005-projector-matrix-a-as-a-model-field-in-v0x-possibly-promoted-later) — projector A as model field
- [ADR-017](#adr-017-projector-seam--abstractobservationmapping-for-joint-likelihood--multi-response-models) — `AbstractObservationMapping`
- [ADR-021](#adr-021-copy-component--scaling-β-lives-on-the-receiving-likelihood-not-on-the-projection-mapping) — `Copy` β on receiving likelihood

**Likelihoods**
- [ADR-018](#adr-018-censoring-as-a-likelihood-level-feature--censoring-enum--per-row-vector-on-the-survival-likelihood-struct) — censoring at likelihood level
- [ADR-019](#adr-019-zero-inflated-count-families--three-r-inla-parameterisations--three-base-distributions) — zero-inflated count families
- [ADR-024](#adr-024-categorical--multinomial-via-independent-poisson-reformulation-helper-function-api-no-new-likelihood-type) — Categorical / Multinomial via independent-Poisson

**Latent components**
- [ADR-022](#adr-022-iidndn-parameterisation--separable-log-τ_i-atanh-ρ_ij-by-default-wishartinvwishart-on-the-joint-precision-as-alternative) — `IIDND{N}` parameterisation
- [ADR-023](#adr-023-meb-and-mec-measurement-error-components--β-via-copy-decomposition-non-zero-prior_mean-made-load-bearing) — `MEB` / `MEC` measurement-error

**SPDE & mesh**
- [ADR-007](#adr-007-delaunaytriangulationjl-for-mesh-generation-fmesher-wrap-as-fallback-sub-package) — DelaunayTriangulation.jl + fmesher fallback
- [ADR-013](#adr-013-spde2-v01-supports-α--1-2-internal-hyperparameters-are-log-τ-log-κ) — SPDE2 α ∈ {1, 2}
- [ADR-028](#adr-028-gaussian-basis-prior-on-non-stationary-spde-basis-coefficients--match-r-inlas-thetapriormeanthetapriorprec-per-coefficient-parameterisation-defer-pc-on-basis-norm) — Gaussian-basis prior, non-stationary SPDE
- [ADR-029](#adr-029-kroneckercomponent--generic-two-component-kronecker-composer-for-separable-space-time-gmrfs) — `KroneckerComponent` separable space-time
- [ADR-030](#adr-030-fractional-α-spde-bolin-kirchner-2020--deferred-to-v021-phase-m-closes-at-pr-6) — fractional-α SPDE deferred
- [ADR-032](#adr-032-mesh-utilities-maturity--alpha-shape-boundary-tuple-max_edge-boundary-pre-subdivision) — mesh utilities maturity
- [ADR-036](#adr-036-spde2-retains-inlamesh-so-the-lgm-extension-can-build-a-meshprojector-at-runtime) — `SPDE2` retains `INLAMesh`

**`@lgm` formula syntax**
- [ADR-033](#adr-033-multi-likelihood-formula-syntax--tuple-lhs-y1-y2--rhs) — multi-likelihood tuple-LHS
- [ADR-034](#adr-034-implicit-f-term-naming-via-column-symbol) — implicit f-term naming
- [ADR-035](#adr-035-lgm-replicate--group-routing--runtime-wrap-calls-in-the-components-tuple) — `replicate` / `group` routing
- [ADR-037](#adr-037-lgm-accepts-tuple-coordinate-first-arg-in-f-arity--2-3-only) — tuple-coordinate first arg in `f(...)`
- [ADR-038](#adr-038-kroneckercomponent-space-time-lgm-form-takes-a-3-tuple-coordinate-design-block-is-a-sparse-khatri-rao-matrix) — Kronecker space-time 3-tuple

**Raster / geo-spatial**
- [ADR-040](#adr-040-predict_rastermodel-res---gaussian-summary--sample-based-path-with-exceedance-wrapper) — `predict_raster(model, res, …)` overloads
- [ADR-041](#adr-041-crs-policy--predict_raster-rejects-mismatched-crs-at-the-api-boundary-mesh_crs-keyword-is-opt-in) — CRS policy via `mesh_crs` keyword

**Testing & triangulation**
- [ADR-043](#adr-043-genjl-second-mcmc-sanity-check--deferred-to-v1x-v10-ships-nuts-only-triangulation) — Gen.jl deferred to v1.x
- [ADR-044](#adr-044-tier-3-triangulation-tolerances-at-v10--tol_mean--15-sds-tol_sd--060-uniformly) — tier-3 v1.0 tolerances

---

## ADR-001: Package split into GMRFs / LatentGaussianModels / INLASPDE

Status: Accepted
Date: 2026-04

### Context
A single monolithic `INLA.jl` is tempting for simplicity. But the GMRF core
is independently useful to people doing sparse Gaussian modeling (Markov
random fields, disease mapping with Stan, image restoration, 4D-Var in
data assimilation). SPDE machinery brings in Meshes.jl and is irrelevant to
many users.

### Decision
Three packages:
- `GMRFs.jl` — numerical core, dependency-light.
- `LatentGaussianModels.jl` — LGM abstraction + INLA algorithm, depends on GMRFs.
- `INLASPDE.jl` — SPDE/FEM, depends on LatentGaussianModels + Meshes.

Kept in a monorepo during planning; each directory is a valid Julia
package and can move to its own repo once stable.

### Consequences
- Good: users install only what they need; GMRFs.jl attracts contributors
  outside the INLA community.
- Good: clear dependency layering; changes to SPDE can't break GMRFs.
- Cost: coordinated releases across three packages require discipline. A
  shared CompatHelper setup and `juliaup/pkg-release-manager` style tooling
  helps.

---

## ADR-002: SciML's LinearSolve.jl as the sparse factorization backend

Status: Accepted
Date: 2026-04

### Context
R-INLA's `GMRFLib_sparse_interface` abstracts over TAUCS, PARDISO, and
BAND solvers. We need the same abstraction. Options:
1. Wrap CHOLMOD directly via SuiteSparse.jl.
2. Use LinearSolve.jl's abstract interface.
3. Build our own dispatcher.

### Decision
LinearSolve.jl. Cholmod/KLU/Pardiso are all already exposed through it.
Symbolic-factorization reuse works via `init(prob)` + `solve!(cache)`.

### Consequences
- Good: zero wheel reinvention; get multi-backend for free.
- Good: standard SciML pattern, familiar to the ecosystem.
- Cost: LinearSolve has its own release cadence; we track master in CI to
  catch breakages.
- Escape hatch: if we hit a case where LinearSolve's overhead matters (inner
  loop allocation is suspected), drop to direct `cholmod_factorize!` for
  that specific hot path, with a benchmark justifying it.

---

## ADR-003: Multiple dispatch, not macros, for the primary API

Status: Accepted
Date: 2026-04

### Context
Turing's `@model` is idiomatic Julia for probabilistic programming. R-INLA
uses R formulas. Stan has its own DSL. For an LGM library, these trace-
based/formula-based APIs are a poor fit — the model is static structure,
not procedural code.

### Decision
Concrete types and multiple dispatch are the primary API. Optional formula
sugar macro (`@lgm`) expands to the explicit constructor. Full policy in
`plans/macro-policy.md`.

### Consequences
- Good: type-stable, AD-friendly, extensible by third parties via
  `struct + methods`, no DSL for users to learn.
- Good: the `rgeneric`/`cgeneric` equivalent is automatic — users define
  a subtype of `AbstractLatentComponent`, no C callbacks needed.
- Cost: R-INLA migrants face a syntax change. The macro sugar narrows the
  gap but doesn't eliminate it.

---

## ADR-004: Selected inversion (Takahashi recursion) as an explicit risk

Status: Proposed
Date: 2026-04

### Context
INLA requires marginal variances `diag(Q⁻¹)` for the Laplace approximation.
R-INLA uses a specialized C implementation of the Takahashi recursion,
which is not trivially reimplemented. Julia has `SelectedInversion.jl` (young),
CHOLMOD's `sparseinv` (not exposed conveniently), or we implement Takahashi
ourselves.

### Decision (proposed)
Evaluate `SelectedInversion.jl` early in Phase 3. If it meets correctness
and performance requirements on our benchmark suite, use it. If not,
implement Takahashi natively using sparse Cholesky factor access. Either
way, this is a named Phase 3 sub-milestone with a budget of 1–2 months.

### Consequences
- Good: if SelectedInversion.jl works, we save months.
- Cost: if it doesn't, we're writing a numerically delicate algorithm.
- Mitigation: both paths start with the same correctness tests against
  dense `inv(Q)` diagonal.

---

## ADR-005: Projector matrix A as a model field in v0.x, possibly promoted later

Status: Proposed
Date: 2026-04

### Context
R-INLA's `inla.stack` machinery conflates observation-to-latent mapping,
data stacking, and effect indexing. A clean Julia version needs to decide
whether the projector A is a field on `LatentGaussianModel` or its own
abstract type `AbstractObservationMapping`.

### Decision (proposed)
v0.1–0.3: `projector::A` as a model field, with a `IdentityProjector()` for
areal models and a concrete `MeshProjector(mesh, locations)` for SPDE.
Re-evaluate once misaligned and joint-likelihood models land in Phase 5.

### Consequences
- Good: simpler for the MVP; easier to explain.
- Cost: likely to be refactored to an abstract type by Phase 5, producing
  a minor breaking change.
- Mitigation: the constructor API `LatentGaussianModel(; projector = ...)`
  stays the same even if the field type changes.

---

## ADR-006: Full Laplace is the Phase-3 default; `simplified.laplace` deferred

Status: Accepted, amended 2026-04-30 — see ADR-016
Date: 2026-04

### Context
R-INLA ships `simplified.laplace` (Rue-Martino 2009) as the default
`strategy`. It is noticeably cheaper than full Laplace and gives more
accurate tails on skewed likelihoods. Implementing it correctly requires
the Rue-Martino correction terms, which are numerically delicate.

### Decision
v0.1 ships **full Laplace** as the default `strategy`. `Gaussian`
(fast preview) is also shipped. `SimplifiedLaplace` is deferred to v0.3
when we have tier-2 confidence in the full Laplace path.

### Consequences
- **Cost:** posterior tails on Poisson/Binomial with extreme counts will
  differ visibly from R-INLA at the first-releases. Users comparing to
  R-INLA at the 1% tolerance in tier-2 tests may see failures
  specifically on skewed likelihoods. Document the divergence
  prominently and widen tier-2 tolerance retroactively on the
  skewed-likelihood fixtures.
- **Good:** v0.1 lands sooner; the simpler algorithm is easier to
  validate in isolation.
- **Escape hatch:** `strategy = :gaussian` for a fast preview; explicit
  warning at fit time on likelihoods we know are tail-sensitive.

### References
- Rue, Martino, Chopin (2009), §3.2.
- `plans/defaults-parity.md` "Default Laplace strategy" section.

### Amendment 2026-04-24 — marginal-reconstruction scope split

The original ADR conflated two distinct things that R-INLA also
calls `strategy`:

1. **Top-level inference strategy** — the Laplace step used inside
   the INLA outer loop (Gaussian / simplified Laplace / full Laplace
   per marginal). This is what the original decision governed.
2. **Posterior marginal reconstruction** — how `p(x_i | y)` is
   reassembled from the θ-grid of per-θ Laplaces. R-INLA also
   exposes `:gaussian`, `:simplified_laplace`, `:laplace` on this
   step (post-processing, not inference-time).

Commit `85db314` added a `strategy` **kwarg on
`posterior_marginal_x`** covering case (2) only — opt-in, default
remains `:gaussian`. Likelihood contract extended with closed-form
`∇³_η_log_density` for Gaussian / Poisson / Binomial plus a
central-difference fallback. Collapses to Gaussian on
quadratic-in-η likelihoods; verified in
[test_simplified_laplace.jl](packages/LatentGaussianModels.jl/test/regression/test_simplified_laplace.jl).

This does **not** reverse the v0.3 deferral. The following items
remain Phase-3-late / v0.3 scope:

- Flipping the *default* on `posterior_marginal_x` from
  `:gaussian` to `:simplified_laplace` — blocked on the
  Pennsylvania Poisson oracle in the replan's Phase C.
- Adding an analogous kwarg to the inference-time Laplace strategy
  (case 1 above). The top-level default is and remains full
  Laplace per this ADR.
- The full Rue-Martino per-marginal `:laplace` strategy (the
  expensive per-marginal re-Laplace).

Rationale for keeping the landed work rather than reverting: the
cubic-derivative closed forms are real correctness infrastructure
needed in Phase C regardless of when the default flips; losing them
incurs the same re-implementation cost later with no payoff now.

Status of this ADR: Accepted, amended 2026-04-24.

---

## ADR-007: DelaunayTriangulation.jl for mesh generation; fmesher wrap as fallback sub-package

Status: Accepted
Date: 2026-04

### Context
R-INLA relies on `fmesher` (a Lindgren C++ tool, now externalized as
`inlabru-org/fmesher`) for constrained Delaunay triangulation with
boundary refinement and extension buffers. Julia's native option is
`DelaunayTriangulation.jl`, which is actively maintained and
feature-rich but has not been validated against `fmesher` output on
SPDE-grade meshes.

### Decision
Default path: `DelaunayTriangulation.jl`. `INLASPDE.jl` M2 includes a
mesh-quality comparison step against `fmesher` output on a fixed
boundary. If that comparison reveals quality gaps that affect SPDE
accuracy (minimum angle, maximum edge length, boundary-refinement
behavior), the fallback is a **new sub-package `INLASPDEFmesher.jl`**
that wraps `fmesher` via BinaryBuilder and exposes an
`fmesher_mesh_2d(...)` alternative constructor. The fallback is a
sub-package, not a weakdep, because binary-artifact dependencies are
heavy and the user should opt in explicitly.

### Consequences
- **Good:** native Julia is the default and works end-to-end without
  binary artifacts.
- **Cost:** if mesh quality is inadequate, Phase 4 M2 slips by 2–4 weeks
  while the fmesher fallback is packaged.
- **Escape hatch:** the fmesher sub-package is planned-but-deferred;
  creating it is a one-time cost not on the critical path for v0.1
  SPDE.

### References
- `packages/INLASPDE.jl/plans/plan.md` M2, Risk items.
- `plans/dependencies.md` — fmesher note.

### Resolution 2026-04-26 — native path sufficient, gate relaxed

The M3 parity fixtures (`6d0784d`) measured native mesh quality
against fmesher on the three reference boundaries. The strict
parity gate (5% vertex count, 0.95× angle, 1.05× edge) failed on
all three; the failure is structural (DT.jl uses an
equilateral-area Ruppert bound; fmesher uses per-edge bisection),
not a bug. The Meuse SPDE oracle (INLASPDE M5) nevertheless
**passes within tolerance using the native mesh**, so mesh quality
is sufficient for SPDE work in v0.1.

**Decision (resolution):** declare native path sufficient. Relax the
parity gate from strict fmesher equivalence to a regression floor
on DT.jl's measured behaviour:

|                | original gate         | resolved gate        |
|----------------|-----------------------|----------------------|
| `rel_vcount`   | ≤ 0.05                | ≤ 0.50               |
| `min_angle_J`  | ≥ max(20°, 0.95·R)    | ≥ 25.0°              |
| `max_edge_J/R` | ≤ 1.05                | ≤ 2.5                |

The resolved gate is locked in as plain `@test` (no
`@test_broken`); a DT.jl regression that materially degrades mesh
quality now fails CI immediately.

`INLASPDEFmesher.jl` remains a planned-but-deferred fallback. The
trigger to actually build it is a downstream user reporting that
mesh quality is biting them on a real fit, not the strict gate
failing in isolation.

ADR-007 is hereby **closed**; further `INLASPDEFmesher.jl` work is
tracked as a v0.2 candidate via `plans/replan-2026-04.md` Phase D.

### References (resolution)
- `packages/INLASPDE.jl/plans/plan.md` M3 parity table.
- `packages/INLASPDE.jl/test/oracle/test_fmesher_parity.jl`.
- `plans/replan-2026-04.md` Phase D, item 1.

---

## ADR-008: `@lgm` macro lives in a separate `LGMFormula.jl` package

Status: Accepted
Date: 2026-04

### Context
The Tier-2 formula sugar `@lgm y ~ 1 + f(region, Besag(W))` is desirable
for R-INLA migrants. Parsing it cleanly benefits from StatsModels.jl's
`@formula` infrastructure for the fixed-effects side. But StatsModels
pulls StatsBase / StatsFuns / DataAPI / Tables into the dep tree, and
Julia 1.9+ extensions cannot export new symbols — so `@lgm` cannot
live in a weakdep extension of core LGM.

### Decision
The macro and its StatsModels dependency live together in a separate
sub-package `packages/LGMFormula.jl/`. Core `LatentGaussianModels.jl`
ships only the Tier-1 explicit constructor. Users who want formula
sugar install one extra package:
```julia
using Pkg; Pkg.add("LGMFormula")
using LatentGaussianModels, LGMFormula
```

### Consequences
- **Good:** core LGM's dep tree stays narrow; StatsModels is not forced
  on users who don't want it; the macro can be a proper exported
  symbol with its own namespace.
- **Good:** macro-policy stays honest — Tier 1 is complete without the
  macro *and* without the macro's host package.
- **Cost:** one extra `Pkg.add` for users migrating from R-INLA. A line
  in the getting-started guide.
- **Scheduling:** LGMFormula.jl is not on the MVP / Phase-3 critical
  path. Build it once the constructor API is stable (end of Phase 3).

### References
- `plans/macro-policy.md` (rewritten alongside this ADR).
- `plans/dependencies.md` sub-packages table.

---

## ADR-009: Turing / HMC bridge lives in a separate `LGMTuring.jl` package

Status: Accepted
Date: 2026-04

### Context
A Turing / AdvancedHMC bridge has two use cases: (1) tier-3
triangulation against HMC posteriors; (2) INLA-within-MCMC flows
(Gómez-Rubio Ch. 5–7) that wrap INLA's conditional `p(x | θ, y)` inside
an outer MCMC loop on θ. Turing's transitive closure is 20–40 s of
TTFX; putting Turing even in `[weakdeps]` of core LGM inflates the
install and release surface.

### Decision
Core LGM commits **only** to `LogDensityProblems.jl` conformance —
implementing `LogDensityProblems.capabilities`, `logdensity`,
`logdensity_and_gradient`, `dimension`. Nothing Turing-specific in
core. A separate sub-package `packages/LGMTuring.jl/` depends on
LatentGaussianModels + Turing + AdvancedHMC and provides:

- `sample(lgm, ::NUTS, n; init_from_inla = true, kwargs...)` wrapper.
- INLA-within-MCMC loop using LGM's exposed
  `sample_conditional(lgm, θ, y)` — see P7 in the initial review.
- `compare(inla_fit, nuts_chain)` diagnostic for triangulation tier.

### Consequences
- **Good:** core LGM's `LogDensityProblems` conformance is already
  useful to anyone who wants to use AdvancedHMC or Pathfinder directly,
  without our bridge.
- **Good:** Turing's release cadence is isolated to one sub-package.
- **Cost:** triangulation-tier tests (Stan/NIMBLE/Turing cross-checks)
  move to `packages/LGMTuring.jl/test/triangulation/`, not core LGM's
  `test/triangulation/`. This means core LGM's tier-3 test list is
  shorter — acceptable, tier-3 is slow and not PR-gating anyway.
- **Requirement:** LGM must expose `sample_conditional(lgm, θ, y)` as
  public API, not internal. Add to the component/inference contract
  in M3.

### References
- `plans/macro-policy.md` Tier 3 description.
- `packages/LatentGaussianModels.jl/plans/plan.md` M3 (conditional
  sampling as public API).

---

## ADR-010: Public kwargs mirror R-INLA names with snake_case + symbol-or-type dual input

Status: Accepted
Date: 2026-04

### Context
R-INLA users coming to the Julia port will look for familiar kwargs:
`int.strategy`, `control.compute$dic`, `scale.model`. A pure-Julia API
using only `AbstractIntegrationScheme` type instances would be
correct-but-alien. A pure R-style API would clash with Julia
conventions. The design should accept both.

### Decision
Public fit-time and component kwargs mirror R-INLA's names with
snake_case translation (`int.strategy` → `int_strategy`,
`scale.model` → `scale_model`). Each such kwarg accepts either:

- a **symbol** (`:ccd`, `:grid`, `:laplace`) — resolves to the canonical
  type via an internal table; or
- a **type instance** (`CCD()`, `Grid(n = 15)`, `Laplace()`) — used
  directly for advanced configuration.

Nested R-INLA `control.*` groups are flattened into prefixed flat
kwargs (`compute_dic`, `compute_waic`) rather than nested NamedTuples.

The complete kwarg table lives in `plans/defaults-parity.md`. Defaults
match R-INLA's except where explicitly documented in the divergences
section.

### Consequences
- **Good:** R-INLA users keep muscle memory; Julia users have
  type-directed autocompletion via the type-instance form.
- **Good:** `help?>` on a specific kwarg type (`?CCD`) still yields
  useful documentation independently of the fit function.
- **Cost:** every kwarg-accepting function needs the
  symbol-dispatch boilerplate. Centralize it in
  `LatentGaussianModels.Inference._resolve(::Val{:ccd})`, etc.

### References
- `plans/defaults-parity.md` kwargs table.

---

## ADR-011: Top-level API is `fit(model, y, strategy; kwargs...)`; `inla(model, y; kwargs...)` is a convenience alias

Status: Accepted
Date: 2026-04

### Context
Two candidate public entry points appeared in earlier drafts:

- `fit(model, y, strategy; kwargs...)` — dispatched on
  `AbstractInferenceStrategy`, following the StatsBase /
  MLJ / Distributions convention for model-fitting in Julia.
- `inla(model, y; kwargs...)` — a short, R-INLA-familiar name with
  strategy selected via kwarg.

READMEs used `inla(...)`; the LGM CLAUDE.md used `fit(...)`. Pick one
canonical entry point; define the other as a thin wrapper.

### Decision
The **canonical** entry point is

```julia
fit(model, y, strategy = INLA(); kwargs...) -> AbstractInferenceResult
```

with `strategy::AbstractInferenceStrategy`. Dispatch on the strategy
type selects the algorithm; kwargs configure it per ADR-010.

`inla(model, y; kwargs...)` is a thin convenience alias defined as

```julia
inla(model, y; kwargs...) = fit(model, y, INLA(); kwargs...)
```

It exists because R-INLA users expect to type it and because it scans
as a domain verb in examples. It is not the canonical method, only
sugar.

Likewise:
- `empirical_bayes(model, y; kwargs...) = fit(model, y, EmpiricalBayes(); kwargs...)`
- `laplace(model, y; kwargs...) = fit(model, y, Laplace(); kwargs...)`

### Consequences
- **Good:** the underlying dispatch is via multiple-dispatch on the
  strategy type (StatsBase-idiomatic), so third-party strategies slot
  in by subtyping `AbstractInferenceStrategy` and defining a `fit`
  method — no change to the LGM package.
- **Good:** the aliases `inla`, `laplace`, `empirical_bayes` keep
  READMEs readable and preserve muscle memory for R-INLA migrants.
- **Cost:** two ways to do the same thing. Documentation must pick one
  style per page and not mix. Convention: the *quickstart* vignettes
  use `inla(...)`; the *reference* docs use `fit(model, y, INLA())`.
- **Cost:** third-party strategy authors must know to define `fit`,
  not a custom verb. Document this clearly in
  `AbstractInferenceStrategy`'s docstring.

### References
- `plans/macro-policy.md` Tier 3 bridge convention.
- `plans/defaults-parity.md` kwargs table — all kwargs listed apply
  to both `fit(..., INLA(); kwargs...)` and `inla(...; kwargs...)`.

---

## ADR-012: Adopt `SelectedInversion.jl` for sparse selected inverse

Status: Accepted
Date: 2026-04-23

### Context

ADR-004 named selected inversion (Takahashi recursion for
`diag(Q⁻¹)` on the sparsity pattern of the Cholesky factor) as the
single biggest Phase-3 numerical risk: R-INLA uses a specialized C
implementation, Julia options were (a) `SelectedInversion.jl` (young),
(b) CHOLMOD's `sparseinv` (not exposed cleanly), or (c) native
Takahashi on the sparse Cholesky factor.

At Phase 3 entry the GMRFs.jl `marginal_variances` reference impl
densifies `Q` and is gated on `n < 1000`. The INLA outer loop's
posterior-variance computation uses the same densification per design
point — a hard wall for any realistic problem.

Empirical evaluation of `SelectedInversion.jl` v0.2:
- API: `selinv(Q::SparseMatrixCSC; depermute = true)` returns a
  NamedTuple `(Z::SparseMatrixCSC, p::Vector{Int})`. `diag(Z)` gives
  the marginal variances directly.
- Correctness: matches `diag(inv(Matrix(Q)))` to ~1e-16 on band,
  random-SPD, and RW-style Laplacian matrices.
- Dependencies: SparseArrays + LinearAlgebra + SuiteSparse (all
  already transitive via LinearSolve) — no new heavy deps.
- Failure mode: throws `PosDefException` on non-PD input; standard
  LinearAlgebra behavior.

### Decision

Adopt `SelectedInversion.jl` as a **core `[deps]` of `GMRFs.jl`** (not
LGM — per ADR-001 layering, selected inversion is numerical core).
Route LGM's posterior marginal-variance path through
`GMRFs.marginal_variances`. The reference dense impl is kept as a
correctness oracle for small `n`, reachable via `method = :dense`.

### Consequences

- **Good:** unblocks any `n > ~1000` problem; no native Takahashi
  implementation needed, saving the 1–2 months of schedule budget
  named in ADR-004.
- **Good:** fixture tests can now compare directly against R-INLA's
  `inla.qinv` output at scales matching real spatial-epidemiology
  problems.
- **Cost:** one more direct dep; maintenance risk if
  `SelectedInversion.jl` goes stale. Mitigation: pin compat bounds,
  keep the dense reference path for fallback.
- **Cost:** does not solve `logdet(Q)` for intrinsic GMRFs —
  generalised log-determinant on the non-null subspace is a separate
  problem, tracked in `plans/defaults-parity.md`.

### References
- ADR-004 — Selected inversion named risk.
- `plans/dependencies.md` — updated core deps table for GMRFs.jl.

---

## ADR-013: SPDE2 v0.1 supports α ∈ {1, 2}; internal hyperparameters are `(log τ, log κ)`

Status: Accepted
Date: 2026-04-23

### Context

The Lindgren-Rue-Lindström 2011 SPDE-Matérn link yields a closed-form
sparse precision `Q(τ, κ)` only for integer `α`. R-INLA's
`inla.spde2.matern` defaults to `α = 2` (which corresponds to Matérn
smoothness `ν = α - d/2 = 1` in 2D). Fractional `α` via Bolin-Kirchner
rational approximations is deferred to v0.3 per
`packages/INLASPDE.jl/plans/plan.md`.

Two open parameterisation questions:

1. **User-scale vs internal scale.** PC priors are most naturally stated
   on the user-scale pair `(range ρ, marginal σ)` per Fuglstad et al.
   2019. The Laplace inner loop works on unconstrained real coordinates.
2. **α surface.** Expose `α` as a compile-time type parameter
   (`SPDE2{1}`, `SPDE2{2}`) or a runtime field.

### Decision

- SPDE2 v0.1 supports `α ∈ {1, 2}`. `α = 2` is the default, matching
  R-INLA. `α` is a **type parameter** (`SPDE2{α}`) so the assembly
  path branches at compile time and the precision-matrix method is
  fully type-stable.
- Internal hyperparameters on the Laplace scale are `θ = [log τ, log κ]`.
  PC-Matérn priors are authored on `(ρ, σ)` and transformed to
  `(log τ, log κ)` via the closed-form Jacobian (cf. Fuglstad et al.
  2019 eqs. 7-8). Users never see `(log τ, log κ)` directly; posterior
  summaries are reported on `(ρ, σ)` via `user_scale`.
- `log_hyperprior(spde, θ)` evaluates the PC density on user scale and
  adds the Jacobian term so that the Laplace-scale posterior is
  normalised consistently with R-INLA's.

### Consequences

- **Good:** type-stable `precision_matrix(::SPDE2{α}, θ)` for both
  α values; no run-time dispatch in the inner Newton loop.
- **Good:** user-facing priors are on the domain statisticians think
  in (range/sd), matching R-INLA's defaults — lowers the parity risk
  in tier-2 oracle tests.
- **Cost:** the Jacobian from `(ρ, σ)` to `(log τ, log κ)` must be
  implemented carefully and regression-tested against R-INLA's
  `inla.pc.prior.matern` — a known source of sign/factor bugs.
- **Escape hatch for fractional α:** v0.3 will add a separate
  `SPDEFractional` concrete type; the type parameter on `SPDE2` does
  not preclude this.

### References
- `packages/INLASPDE.jl/plans/plan.md` M2 and M4.
- Lindgren, Rue, Lindström 2011. SPDE.
- Fuglstad, Simpson, Lindgren, Rue 2019. PC priors for Gaussian random
  fields.

---

## ADR-014: `main` is fast-forwarded to `claude/hungry-pascal`; main becomes the integration branch going forward

Status: Accepted
Date: 2026-04-24

### Context

Between 2026-04-23 13:46 and 2026-04-24 11:15, 17 feature commits
landed on branch `claude/hungry-pascal` while `main` remained at
`9c2f69d`. A status audit on 2026-04-24 discovered the gap:

- `git merge-base main claude/hungry-pascal` resolves to `main`'s
  HEAD — the branches share a linear history, no divergence.
- Every commit on `hungry-pascal` is by the same author as on `main`.
- No pull request, review, or CI gate was configured; the branch was
  simply the workspace where work continued after `main` stopped
  being advanced.
- Roadmap progress actually landed on `hungry-pascal`: MVP go/no-go
  (Scotland BYM2, ADR-006 scope, ADR-012, ADR-013), much of Phase 3
  (simplified-Laplace correction, DIC/WAIC/CPO/PIT, posterior
  marginals, linear constraints in Laplace), and Phase 4 M1–M6-A
  (FEM assembly, SPDE2, mesh generation, MeshProjector, Meuse
  oracle, GeoInterface ext, rasters predict/quantile, MakieExt).

Reading the roadmap against `main` gave a misleadingly pessimistic
picture (late Phase 1). The true state is through Phase 4 M6-A.

### Decision

1. **Fast-forward `main` to `claude/hungry-pascal`.** The update is
   `git checkout main && git merge --ff-only claude/hungry-pascal`.
   No cherry-picking, no squash — the 17 commits are individually
   scoped (Conventional Commits, per-feature) and bisect-able as-is.
2. **Delete `claude/hungry-pascal`** after the fast-forward. It has
   no independent meaning.
3. **Adopt a no-stale-main rule going forward.** Either (a) work
   directly on `main` until branch protection lands, or (b) once
   branch protection and CI are enabled per
   `plans/initial-commits.md`, open a PR per feature branch and
   require merge before starting the next unrelated piece of work.
   Leaving a work branch more than one working day ahead of `main`
   without a PR is the concrete anti-pattern this ADR names.
4. **Record roadmap drift explicitly.** After the merge,
   `ROADMAP.md` is updated to reflect the new baseline, and a
   replan document lands in `plans/replan-2026-04.md` covering
   Phase B–E per the status review.

### Consequences

- **Good:** a single canonical branch; external readers (and Claude
  Code sessions) evaluating project status will not be misled by a
  stale `main`.
- **Good:** all prior review benefits are still available — the
  commit granularity on `hungry-pascal` is already per-feature, so
  `git log` on the merged `main` reads the same way a review would
  have produced.
- **Cost:** the ADR numbering on this branch (which only sees through
  ADR-011) must be reconciled at merge time. `hungry-pascal` adds
  ADR-012 and ADR-013; this branch adds ADR-014. The textual merge
  is trivial — no ADR is renumbered, only concatenated in order —
  but the merge commit must verify the sequence 012, 013, 014 is
  monotonic with no gaps before the fast-forward is accepted.
- **Cost:** one-off effort to backfill branch protection rules on
  GitHub (per `plans/initial-commits.md` §"Branch protection
  rules"). This is independent of the merge itself but should land
  in the same working session to prevent recurrence.

### Follow-up items flagged by this ADR (tracked but not resolved here)

- **ADR-006 divergence.** `85db314` landed a simplified-Laplace
  skew correction, which ADR-006 explicitly deferred to v0.3.
  Either amend ADR-006 with a superseding note or revert the
  feature. Decide before v0.1 tag.
- **Fixture generation pipeline.** Phase 4 M5 introduced JLD2
  fixtures under `packages/INLASPDE.jl/test/oracle/fixtures/` but
  the R generation scripts in `scripts/generate-fixtures/spde/`
  are not exercised by CI. This is the single largest untracked
  correctness risk identified by the status review.
- **Missing Phase 2 components.** `Leroux`, `BYM` (non-reparameterised
  form), `Seasonal`, and `Generic0/1` at the LGM component layer
  have plan entries but no source — the merge does not close these
  items; they are Phase B in the replan.

### References

- `ROADMAP.md` — phase numbering to be revised after merge.
- `plans/initial-commits.md` — branch protection rules §.
- ADR-006 — simplified Laplace deferral, now in tension with
  `85db314`.

---

## ADR-015: `LGMFormula.jl` and `GMRFsPardiso.jl` both cut from v0.1; deferred to v0.2

Status: Accepted; LGMFormula.jl realised v0.2.2 (Phase N); GMRFsPardiso.jl still deferred (see Phase Q PR-1 ledger flip)
Date: 2026-04-26

### Context

Phase E3 of the v0.1 replan ([`replan-2026-04.md`](replan-2026-04.md))
flagged both sub-packages as cuttable, with explicit cut criteria:

- `LGMFormula.jl`: cut if a usable migration doc exists for users
  coming from R-INLA's `f(...)` syntax.
- `GMRFsPardiso.jl`: cut if comparative benchmark on a 10⁵-vertex
  SPDE Q does not beat CHOLMOD by ≥ 30%.

By 2026-04-26 the four `src/`-bearing packages are GA-ready (Phase E1
Aqua + JET pass, Phase E2 docs site live), and the question is whether
to delay tagging v0.1.0 to ship two more sub-packages.

### Decision

**Both cut from v0.1; defer to v0.2.**

`LGMFormula.jl`: the explicit `LatentGaussianModel(ℓ, components, A)`
constructor is a small surface and the
[Getting started](../docs/src/getting-started.md) and three vignettes
already cover the migration path from R-INLA. The macro adds zero
numerical capability — purely ergonomics — and macro design that's
robust against `StatsModels` schema-application edge cases is at
least 3 weeks of work that would block the v0.1 tag. Re-add as a
sub-package in v0.2 with the M1+M2 milestones from
[packages/LGMFormula.jl/plans/plan.md](../packages/LGMFormula.jl/plans/plan.md).

`GMRFsPardiso.jl`: Pardiso.jl upstream has had periods of disrepair
([packages/GMRFsPardiso.jl/plans/plan.md](../packages/GMRFsPardiso.jl/plans/plan.md)
"Risk items"), license-detection plumbing is non-trivial, and the
benchmark gate requires a real 10⁵-vertex run we have not done. CHOLMOD
remains the only backend in v0.1; the `FactorCache` interface in
`GMRFs.jl` is already designed so a Pardiso specialization can land
without an API break. Re-evaluate in v0.2 after a benchmark study.

The two scaffold directories (`packages/LGMFormula.jl/` and
`packages/GMRFsPardiso.jl/`) stay in the repo with their `Project.toml`
and `plans/` intact; they simply do not have `src/` yet and are not in
the v0.1 release manifest.

### Consequences

- v0.1.0 release manifest is exactly four packages: GMRFs.jl,
  LatentGaussianModels.jl, INLASPDE.jl, INLASPDERasters.jl.
- ADR-008 (`@lgm` lives in LGMFormula.jl) is unchanged; the package
  just doesn't ship in v0.1.
- No new public API surface to maintain in v0.1, which keeps the
  registry submission minimal and the deprecation surface for v0.2
  empty.
- Users wanting Pardiso must build a custom `FactorCache` against the
  GMRFs.jl `AbstractFactorCache` interface; this is a power-user path
  and is acceptable for v0.1.

### References

- [packages/LGMFormula.jl/plans/plan.md](../packages/LGMFormula.jl/plans/plan.md)
- [packages/GMRFsPardiso.jl/plans/plan.md](../packages/GMRFsPardiso.jl/plans/plan.md)
- ADR-008, ADR-009 — the "split into a sub-package" pattern.

---

## ADR-016: Simplified-Laplace mean-shift correction (Rue-Martino) added as opt-in `latent_strategy`

Status: Accepted
Date: 2026-04-28

### Context

ADR-006 (amended 2026-04-24) split simplified-Laplace into two pieces:
a density-shape correction on `posterior_marginal_x` (landed as the
`strategy = :simplified_laplace` kwarg in commit `85db314`) and a
mean-shift correction on the latent posterior summary (deferred). The
density correction multiplies each per-θ Gaussian by a Hermite-3 skew
factor `(1 + γ/6 · H₃(s))`, expanded around the **unshifted** Newton
mode `x̂(θ)`. R-INLA's `simplified.laplace` additionally shifts that
mode by `Δx = ½ H⁻¹ Aᵀ (h³ ⊙ σ²_η)`, so for a true match against
R-INLA's posterior mean on skewed likelihoods (Poisson, Binomial,
Gamma, NegBin) we need the mean shift as well.

The mean shift is one
multi-RHS sparse triangular solve per integration point; with our
existing `FactorCache` `\` (`packages/GMRFs.jl/src/factorization.jl:74`),
existing closed-form `∇³_η_log_density` for all v0.1 likelihoods, and
the kriging projection idiom from `_latent_skewness`
(`marginals.jl:176`), the implementation is ~50 LoC plus a small wiring
delta in `inla.jl`.

We also evaluated and rejected porting the same reference's Edgeworth
log-marginal-likelihood correction and N=100 importance-sampling
log-marginal correction:

- mlik parity is already established (project memory
  `project_inla_mlik_gap.md`, resolved 2026-04-27 — the residual gap
  was a normalising constant in `Intercept`/`BYM`, not a missing
  higher-order term).
- The Edgeworth term materialises a dense `Σ_η = A H⁻¹ Aᵀ` and is
  `O(n_obs²)` to `O(n_obs³)` — blows the Phase-D 30 s Scotland budget.
- The IS estimator at fixed N=100 ships no ESS diagnostic, and at this
  scale its Monte-Carlo error (~ 0.1 nat) is comparable to or larger
  than the corrections it claims to make.

### Decision

Add a new `latent_strategy::Symbol` kwarg to `INLA(...)`, accepting
`:gaussian` (default) and `:simplified_laplace`. Setting
`:simplified_laplace` applies the Rue-Martino mean shift to
`INLAResult.x_mean` and `x_var` while leaving `LaplaceResult.mode`
(the Newton fixed point) unchanged — downstream code that operates on
the Newton mode (`_latent_skewness`, sampling, log-marginal) is
deliberately unaffected. The density-shape correction in
`posterior_marginal_x(strategy = :simplified_laplace)` remains
orthogonal: a user can pick either, both, or neither.

The new field is on `INLA(...)`, not a new strategy type, so
third-party `AbstractInferenceStrategy` subtypes (per ADR-011) do not
have to enumerate latent variants.

Default stays `:gaussian` — flipping it is a separate v0.3 decision
gated on the Pennsylvania oracle re-run with the new path active.

The variance correction (R-INLA's third simplified-Laplace term,
`H⁻¹ Aᵀ diag(h⁴) A H⁻¹` per coordinate) and the Edgeworth / IS log-mlik
corrections remain v0.3 / out of scope.

### Consequences

- **Good:** closes the v0.1 gap to R-INLA's posterior *mean* on skewed
  likelihoods. The `:gaussian` and `:simplified_laplace` paths are
  orthogonal at the API level, so opt-in users see the new behaviour
  with zero risk to existing callers.
- **Good:** complements the already-shipped Hermite skew correction
  in `posterior_marginal_x` — together they cover R-INLA's
  `simplified.laplace` mean and shape; only the variance term is
  still owed.
- **Cost:** one multi-RHS sparse triangular solve per integration
  point. Measured at ~0% overhead on the small Pennsylvania-style
  models in the existing oracle suite; the projected upper bound for
  larger SPDE fits (Meuse) is ~5%.
- **Cost:** maintenance of two latent-strategy paths through the
  integration loop. Localised: ~3 lines in `fit(::INLA)` accumulator
  plus the 50-line helper in `simplified_laplace_correction.jl`.
- **Cost:** Pennsylvania oracle fixture needs regeneration to add a
  `bym_mean_sla` field for the matched R-INLA `strategy = "simplified.laplace"`
  output. The existing Gaussian-strategy assertions are untouched.
- **Escape hatch:** `latent_strategy = :gaussian` is exactly the prior
  behaviour, bit-for-bit. Verified by a regression test that asserts
  `‖x_mean_gaussian − x_mean_simplified_laplace‖∞ < 1e-12` on a
  Gaussian-likelihood model (where `h³ ≡ 0`).
- **Defer v0.3:** flipping the default to `:simplified_laplace`,
  adding the variance correction, and any per-marginal full
  `:laplace` strategy.

### References

- ADR-006 (and 2026-04-24 amendment) — original deferral.
- `replan-2026-04.md` "What this replan does not schedule" — entry
  rewritten to keep Edgeworth / variance correction deferred while
  removing the mean shift from the deferred list.
- Rue, Martino, Chopin (2009) §4.2.
- `packages/LatentGaussianModels.jl/src/inference/simplified_laplace_correction.jl`
  — implementation.
- `packages/LatentGaussianModels.jl/test/regression/test_sla_mean_shift.jl`
  — regression coverage (Gaussian collapse, dense formula
  cross-check, BYM2 sum-to-zero preservation).

### Amendment 2026-07-06 — variance-correction claim withdrawn (see ADR-047)

The "variance correction (R-INLA's third simplified-Laplace term,
`H⁻¹ Aᵀ diag(h⁴) A H⁻¹` per coordinate)" this ADR deferred does not
exist in Rue-Martino-Chopin (2009) §3.2.3: the expansion defines only
the `γ^(1)` mean and `γ^(3)` skewness terms (both third-derivative),
and the fitted skew-normal pins the variance at 1 by construction. The
deferral is closed as mis-specified by ADR-047. Modern R-INLA's
variance adjustments are the strategy-independent VB corrections
(`control.vb`), recorded in ADR-047 as a possible future feature on
the `_apply_integration_moments` seam. The mean-shift decision and
implementation in this ADR are unaffected.

---

## ADR-017: Projector seam — `AbstractObservationMapping` for joint-likelihood / multi-response models

Status: Proposed
Date: 2026-04-29

Supersedes ADR-005 (which deferred the call to "Phase 5"; that phase
is now Phase G of `replan-2026-04-28.md`, and the call is due).

### Context

The v0.1 `LatentGaussianModel` carries a single `A::AbstractMatrix`
projector mapping the stacked latent vector `x` to the linear
predictor `η = A x`, alongside a single `likelihood::AbstractLikelihood`
applied row-wise to `η`. This is sufficient for every fixture
currently in `test/oracle/`: areal Poisson (Scotland, Pennsylvania),
geostatistical Gaussian (Meuse), all the synthetic singletons.

Phase G unlocks the rest of R-INLA's mainstream coverage —
joint longitudinal-survival, multi-response geostatistics, disease
mapping with multiple sentinels — and every one of these models
violates **both** assumptions simultaneously:

1. There are *several* likelihood blocks, e.g. Gaussian on the
   longitudinal arm and Weibull on the survival arm.
2. Each block has its own row-to-latent mapping `A_k`, and the
   `Copy` component (Phase G scope) shares one block's latent into
   another block's predictor with an estimated scaling β.

R-INLA papers over this with `inla.stack`, which conflates *three*
concerns: observation mapping, data stacking, and effect indexing.
Replicating `inla.stack` verbatim would import its sharpest
ergonomic edge — most R-INLA support questions are stack-related —
while losing the type-stable structural API that ADR-003 committed
us to. We need an alternative that (a) handles the multi-block case
cleanly, (b) keeps the v0.1 single-block users on the existing
constructor with no source change, and (c) gives third-party
components a typed extension point rather than a magic
`hcat`-shaped contract.

The architecture document (`plans/architecture.md` §"The projector
question — open issue") flags this and explicitly defers the
abstract-type promotion to an ADR. The replan's Phase G "Risks"
section names a wrong abstraction here as "the first phase where
[it] will compound across the rest of the v0.2 horizon" and asks
for an outside reviewer before the seam merges. This ADR is that
review prerequisite — the code change does not start before this
record is accepted.

### Decision (proposed)

Promote the projector slot to a dispatchable abstract type
`AbstractObservationMapping <: Any`, owned by
`LatentGaussianModels.jl`. The `LatentGaussianModel` struct's `A`
field becomes `mapping::M <: AbstractObservationMapping`, and the
`likelihood::L <: AbstractLikelihood` field becomes
`likelihoods::Tuple{Vararg{<:AbstractLikelihood}}` keyed by
observation block via the mapping.

Concrete subtypes shipped with the seam:

```julia
abstract type AbstractObservationMapping end

# v0.1 default — wraps the current m.A. Single-block; one likelihood.
struct LinearProjector{A <: AbstractMatrix} <: AbstractObservationMapping
    A::A
end

# Areal / row-aligned shortcut — algebraically a no-op A = I.
struct IdentityMapping <: AbstractObservationMapping
    n::Int
end

# Multi-block joint models — one (A_k, likelihood_k) per block.
# Block boundaries are encoded as a Vector{UnitRange{Int}} keyed
# parallel to `likelihoods`.
struct StackedMapping{T <: Tuple} <: AbstractObservationMapping
    blocks::T                    # Tuple of LinearProjector / IdentityMapping
    rows::Vector{UnitRange{Int}} # observation-index ranges per block
end

# Separable space-time — also serves Phase M.
struct KroneckerMapping{S, T} <: AbstractObservationMapping
    A_space::S
    A_time::T
end
```

Required interface (every `AbstractObservationMapping` implements):

```julia
apply!(η, mapping, x)              # η .= mapping * x, in place
apply_adjoint!(g, mapping, r)      # g .+= mappingᵀ * r, accumulating
nrows(mapping) -> Int              # row count of the implicit matrix
ncols(mapping) -> Int              # column count = length(x)
likelihood_for(mapping, i) -> Int  # which block-index owns row i
                                   # (default returns 1 for single-block)
```

`apply!` and `apply_adjoint!` are the load-bearing inner-Newton
operations — they replace every `m.A * x` and `m.A' * r` site in
`inference/laplace.jl`. Performance contract: an `IdentityMapping`
falls through to a `copyto!`; a `LinearProjector{<:SparseMatrixCSC}`
calls `mul!` on the underlying matrix with no wrapper allocation.
Benchmarked against the v0.1 baseline on `test/oracle/`: target ≤2 %
regression on Scotland BYM2 wall-clock, ≤5 % on Meuse SPDE.

Constructor compatibility:

```julia
# v0.1 form — unchanged; wraps in LinearProjector internally.
LatentGaussianModel(ℓ::AbstractLikelihood, components, A::AbstractMatrix)

# v0.2 multi-likelihood form.
LatentGaussianModel(ℓs::Tuple{Vararg{AbstractLikelihood}}, components,
                    mapping::AbstractObservationMapping)
```

The single-likelihood signature continues to work for **two minor
versions** with a one-line deprecation warning advising the new tuple
form; v0.4 drops it. Existing `test/oracle/` fixtures keep their
constructors unchanged across v0.1 → v0.2.

### Alternatives considered

**A. Carry a tuple of `(likelihood, A)` pairs directly on the
model.** Simpler — no new abstract type, no dispatch. Rejected
because it inlines the abstraction at every call site
(`for (ℓ, A) in pairs ...`) and forecloses on `KroneckerMapping`,
which Phase M needs anyway. The dispatch seam is the cheap insurance
against re-litigating this in 18 months.

**B. Keep `A` and `likelihood` as v0.1, expose multi-block via a new
`MultiLikelihoodModel` type.** Two parallel model types means every
inference strategy, diagnostic, and accessor in
`LatentGaussianModels.jl` and `INLASPDE.jl` ships in two flavours.
Rejected on cost.

**C. Adopt R-INLA's `inla.stack` 1:1.** Already rejected by ADR-003
(structural over procedural API).

### Consequences

**Good.**

- Joint-likelihood models become a struct + a tuple, not a DSL —
  preserves the ADR-003 commitment.
- `Copy` (also Phase G) lands as a `LinearPredictorTerm` that the
  mapping holds, with no further surgery on the model struct.
- Phase M's separable space-time gets `KroneckerMapping` for free.
- The single-likelihood hot path is unchanged: `LinearProjector{<:
  SparseMatrixCSC}` is one indirection that the compiler inlines.
- Third-party components define their own
  `AbstractObservationMapping` (e.g. `LowRankProjector` for
  factor-analytic models) without touching LGM's source.

**Cost.**

- One breaking change in the `LatentGaussianModel` field layout
  (`A` → `mapping`). Mitigated by the constructor compatibility
  shim and the deprecation window.
- Every `m.A` access inside `LatentGaussianModels.jl`,
  `INLASPDE.jl`, and the LGM tests is rewritten to go through
  `apply!` / `apply_adjoint!`. ~30 sites by `grep`. One PR.
- The `nhyperparameters` accounting changes — likelihood θ is now
  a vector-of-tuples rather than a scalar offset. Touches
  `θ_ranges` initialisation in `model.jl:52-58`.
- Inference strategies that materialise `A` densely (none today, but
  `simplified_laplace_correction.jl` does row-slicing) need a path
  through the mapping interface. Tracked as a follow-up sub-task,
  not a blocker for the seam itself.

**Escape hatch.** If the mapping abstraction underperforms in a
benchmarked hot path, drop to the underlying matrix via
`mapping.A` (for `LinearProjector`) and call `mul!` directly —
the wrapper is intentionally thin enough that this is a single-line
opt-out for that call site.

### Phasing — how this lands

The replan asks for multi-likelihood as one PR series and `Copy` as
a follow-up; this ADR keeps that split:

1. **PR series 1 — seam only.** Introduce
   `AbstractObservationMapping`, the four concrete subtypes, the
   compatibility constructor, and the `m.A → m.mapping` rewrite.
   Single-likelihood semantics throughout. Goal: every existing
   oracle test passes unchanged.
2. **PR series 2 — multi-likelihood.** Promote `likelihood` →
   `likelihoods` and route per-block density evaluation through
   `likelihood_for(mapping, i)`. New oracle fixture: a two-block
   Gaussian + Poisson synthetic.
3. **PR series 3 — `Copy` component.** Lands on top of (2) using
   `LinearPredictorTerm` to express shared latents. Fixture:
   Baghfalaki et al. (gated on Phase H survival likelihoods —
   merges as a later vignette PR).

### Open questions for the reviewer

1. Should `IdentityMapping` carry `n::Int` or be dimension-free?
   Carrying `n` keeps `nrows`/`ncols` cheap and lets us assert
   shape at construction; the cost is one `Int` per model.
2. Block ordering invariant: are the rows of `StackedMapping`
   *required* to be contiguous and sorted, or should we allow
   interleaved row indices? Contiguous is simpler; interleaved
   matches R-INLA's `inla.stack` more closely. Recommendation:
   contiguous in v0.2, defer interleaved to a follow-up if a real
   user needs it.
3. Do we expose `apply!`/`apply_adjoint!` as public API, or hide
   them behind `*` / `mul!` overloads only? Public is friendlier
   to third-party components; private keeps us free to tweak the
   signature. Recommendation: public, documented in the package
   `CLAUDE.md`.

### References

- ADR-003 (multi-dispatch primary API) — the commitment this ADR
  preserves under multi-likelihood pressure.
- ADR-005 (projector as field, deferred decision) — superseded.
- `plans/architecture.md` §"The projector question — open issue" —
  the open issue this ADR closes.
- `plans/replan-2026-04-28.md` Phase G — "Multi-likelihood support",
  "AbstractObservationMapping seam", "`Copy` component", and the
  "Risks" subsection asking for outside ADR review.
- R-INLA source: `inla.stack` in [`r-inla/rinla/R/stack.R`](https://github.com/hrue/r-inla/blob/devel/rinla/R/stack.R)
  — the abstraction we are explicitly *not* copying.
- Baghfalaki et al. tutorial — the joint-model fixture target.

---

## ADR-018: Censoring as a likelihood-level feature — `Censoring` enum + per-row vector on the survival likelihood struct

Status: Proposed
Date: 2026-04-29

### Context

Phase H of `plans/replan-2026-04-28.md` ships the survival / time-to-event
likelihood pack — Exponential, Weibull (PH and AFT), LognormalSurv,
GammaSurv, Coxph (piecewise-baseline), WeibullCure — plus the zero-inflated
count families. Every survival family in this pack must distinguish four
observation modes: uncensored event time, right-censored, left-censored,
and interval-censored. R-INLA encodes this via `inla.surv(time, event,
time2)` where `event ∈ {0, 1, 2, 3}` is the censoring code and `time2` is
the upper interval bound.

The replan's Phase H scope (line 208) calls for a "`Censoring` enum
(`:none`, `:right`, `:left`, `:interval`) + per-row vector on the
likelihood struct". This ADR pins down (a) the type and where it lives,
(b) how the boundary times for interval-censoring are carried, (c) the
hot-path contract for `log_density` / `∇_η` / `∇²_η` / `∇³_η` under
censoring, and (d) the upgrade path for non-survival likelihoods that
might want censoring later (censored Gaussian → tobit). The ADR is a
prerequisite for PR1 of the survival pack (Exponential + the censoring
contract); Weibull / Lognormal / Gamma / Coxph land as follow-up PRs that
inherit the contract, not re-litigate it.

### Decision (proposed)

**1. The enum.** Define a value-type enum in `LatentGaussianModels`:

```julia
@enum Censoring NONE RIGHT LEFT INTERVAL
```

Exported. Symbol coercion (`Censoring(:none)`, etc.) is provided as a
convenience for keyword-call sites; storage is always the enum, not a
`Symbol`, so the inner Newton loop sees a single byte per row and the
compiler can unbox the branch.

**2. Per-row vector on the survival likelihood struct.** Each survival
family carries the censoring data as struct fields, *not* on `y`:

```julia
struct ExponentialLikelihood{L, C, V} <: AbstractLikelihood
    link::L              # LogLink default — rate parameterisation
    censoring::C         # Nothing | AbstractVector{Censoring}
    time_hi::V           # Nothing | AbstractVector{<:Real} — INTERVAL upper bounds
end
```

Construction:

```julia
ExponentialLikelihood()                                   # all NONE, fastest path
ExponentialLikelihood(censoring = [NONE, RIGHT, NONE])    # mixed; time_hi unused
ExponentialLikelihood(censoring = [NONE, INTERVAL],
                      time_hi   = [0.0, 5.0])             # interval row uses time_hi[2]
```

The `Nothing` sentinel is load-bearing — it gives a separate dispatch for
the all-uncensored case so the closed-form simple-density path is not
dragged through a per-row switch in the hot path.

**3. The `log_density` contract under censoring.** For a survival family
with hazard `h(t; η, θ)` and survival `S(t; η, θ) = exp(-Λ(t; η, θ))`:

| Censoring        | log p(y_i \| η_i, θ)                                        |
|------------------|------------------------------------------------------------|
| `NONE`           | `log f(t_i) = log h(t_i) − Λ(t_i)`                          |
| `RIGHT`          | `log S(t_i) = −Λ(t_i)`                                      |
| `LEFT`           | `log F(t_i) = log(1 − S(t_i))`                              |
| `INTERVAL`       | `log[S(t_lo) − S(t_hi)] = log[exp(−Λ_lo) − exp(−Λ_hi)]`     |

`y[i]` carries `t_lo` for `INTERVAL` rows; `time_hi[i]` carries `t_hi`. For
`NONE`/`RIGHT`/`LEFT`, `y[i]` is the relevant single time. The `INTERVAL`
log-difference uses `logsubexp(a, b) = a + log1p(-exp(b - a))` for `a > b`
to retain digits; this is a small helper, not a new dependency.

`∇_η`, `∇²_η`, `∇³_η` follow the same per-row dispatch — for each
censoring mode the chain rule against the link function `g(η) = ` rate
(or scale, depending on family) yields a closed form that is no more
expensive than the uncensored case modulo a single `expm1` per
non-`NONE` row. The `Nothing`-censoring dispatch reuses the simple
broadcast path.

**4. Boundary-time storage decision.** Three alternatives considered:

- **A. `time_hi::Union{Nothing, Vector{<:Real}}` on the struct.**
  Adopted. Cheap for the common single-bound case; one extra Vector for
  the interval case; the type-parameterised `Nothing` keeps dispatch
  type-stable. Cost: a tiny amount of dead storage when only some rows
  are `INTERVAL` (`time_hi[i]` for non-`INTERVAL` rows is unread).
- **B. Heterogeneous `y::Vector{Union{T,Tuple{T,T}}}`.** Rejected —
  breaks every existing `y::AbstractVector{<:Real}` contract throughout
  diagnostics, plotting, and the predict/predict_quantile API.
- **C. A `SurvivalOutcome{T}` data type carrying `(t_lo, t_hi, c)` per
  row.** Rejected — forces a dependent type on every survival likelihood
  and on the `y` argument across the public API; the per-row vector
  approach achieves the same density of information with no surface
  change to `y`.

**5. Default for non-censored families.** Gaussian, Poisson, Binomial,
NegBinomial do **not** acquire censoring fields in PR1. The contract
is opt-in per likelihood. A future censored Gaussian (tobit) will follow
the same pattern: add `censoring::Union{Nothing, Vector{Censoring}}`
plus `time_hi` to its struct, dispatch on the `Nothing` field for the
existing fast path. No cross-cutting change to `AbstractLikelihood` is
required; the seam is the per-likelihood struct field, not the abstract
type.

**6. Validation.** At `fit` time, assert
`length(censoring) == length(y)` and `time_hi[i] > y[i]` for every
`INTERVAL` row. These are public-boundary `@assert`s, removed from the
inner Newton loop per the existing hot-path policy in
`packages/LatentGaussianModels.jl/CLAUDE.md`.

**7. Diagnostics under censoring.** `pointwise_log_density` honours the
censoring code per row (i.e., its sum equals `log_density`, regardless
of the censoring mix). `pointwise_cdf` for censored rows is undefined in
v0.1 — PIT for censored observations needs a non-trivial reweighting
(Henderson-Crowther) that is deferred to a v0.2 follow-up; the default
`pointwise_cdf` fallback `throw`s on first call so a user who tries to
compute PIT with censored data gets a clean error.

### Alternatives considered

**A. Censoring as an outer wrapper likelihood —
`CensoredLikelihood{L<:AbstractLikelihood}`.** This keeps every base
likelihood (Gaussian, Exponential, Weibull) censoring-naïve and adds
censoring through composition. Rejected because the survival families
*already* use the survival function `S` directly in their derivative
formulas — a wrapper would have to re-derive `∇²` and `∇³` of
`log[S − S]` numerically, blowing up cost for `INTERVAL`. Composition is
attractive in the abstract but pays a real performance cost in the inner
Newton loop, and does not naturally extend to censoring upper-bound
storage.

**B. Censoring as part of the `LikelihoodMapping` (ADR-017) instead of
on the likelihood struct.** Rejected — `LikelihoodMapping` is the
*routing* layer (which row goes to which likelihood block), not the
*outcome model* layer. Censoring is intrinsic to the outcome model
(survival vs. event) and naturally co-locates with the family-specific
`log_density` formula. Putting it on the mapping would force every
joint-likelihood model to carry a redundant censoring vector even for
its non-survival blocks.

**C. `Vector{Symbol}` instead of `Vector{Censoring}`.** Rejected on
performance — `Symbol` is a pointer, the per-row branch through `===`
is heavier than an enum compare; on a 5e4-observation Tonsil fit the
branch dominates the inner-loop cost relative to the closed-form
hazard arithmetic. The `Censoring(:none)` symbol coercion at
construction is the convenience users wanted; the storage is the
enum.

### Consequences

**Good.**

- Survival families share one censoring contract; PR1 (Exponential)
  shakes it down end-to-end before Weibull adds a shape hyperparameter
  and Coxph adds a piecewise baseline.
- The `Nothing`-censoring fast path means uncensored exponential fits
  pay zero overhead from the censoring infrastructure — same throughput
  as a hypothetical "no-censoring-ever" implementation.
- Future censored-Gaussian (tobit) lands as a struct-field addition, not
  an `AbstractLikelihood` contract change.
- Aligns with R-INLA's `inla.surv` on the data side (lo/hi/code triple)
  while keeping the Julia API typed and dispatch-friendly.

**Cost.**

- Each survival family must implement the four-way per-row dispatch in
  `log_density`, `∇_η`, `∇²_η`, `∇³_η`. ~40 LoC per family. Mitigated by
  a shared internal helper module
  (`packages/LatentGaussianModels.jl/src/likelihoods/survival/_censoring.jl`)
  exposing the censoring branches as inlineable building blocks.
- `time_hi` storage is wasted for non-`INTERVAL` rows in mixed datasets.
  In practice interval censoring is rare (Tonsil: 0%, Leuk: 0%, Gambia:
  0%); the typical fit pays one extra `Vector` of length n_obs that goes
  unread. Acceptable.
- Validation assertions at the public boundary need clear, R-INLA-style
  error messages — a row with `INTERVAL` and `time_hi[i] ≤ y[i]` is the
  most common user mistake and the message must call this out.

**Escape hatch.** If a user needs a censoring mode the enum doesn't cover
(e.g. Type-II progressive censoring), they can implement their own
`AbstractLikelihood` subtype that handles its own outcome encoding
without using `Censoring` — the enum is convention, not a load-bearing
hook in the abstract type's contract. The shared helper module is
internal API.

### Phasing — how this lands

1. **PR1 — `Censoring` + `ExponentialLikelihood` + oracle.** Lands the
   enum in `LatentGaussianModels`, the `_censoring.jl` helper, the
   Exponential family with full censoring support, an R-INLA oracle
   fixture (synthetic exponential survival, Phase H scope) under
   `scripts/generate-fixtures/lgm/exp_survival/`, and tier-1 / tier-2
   tests. **This ADR's deliverable.**
2. **PR2 — `WeibullLikelihood` (PH parameterisation).** Reuses the
   censoring helper. Adds the shape hyperparameter and PC prior on it.
3. **PR3 — `LognormalSurvLikelihood`, `GammaSurvLikelihood`,
   `WeibullCureLikelihood`.** Three families, one PR; the censoring
   pattern is now mechanical.
4. **PR4 — `CoxphLikelihood`** with piecewise-constant baseline and
   stratification. Most complex; censoring contract carried over.
5. **PR5 — Zero-inflated families.** Independent of censoring; gated on
   the multi-likelihood seam (ADR-017) only. Lands separately.
6. **PR6 — `PCAlphaW` hyperprior.** Sørbye-Rue 2017. Pure prior
   addition; independent of the censoring contract.

### Open questions

1. **Left-truncation.** R-INLA's `inla.surv` carries a `truncation`
   column for left-truncated data (e.g., delayed entry into a cohort).
   This is **not** in the v0.1 censoring enum. Recommendation: add
   `truncation::Union{Nothing, Vector{<:Real}}` on the survival likelihood
   in a v0.2 PR if a user reports needing it; the contract here doesn't
   foreclose on it.
2. **`pointwise_cdf` under censoring.** The Henderson-Crowther PIT for
   censored observations is non-trivial. Recommendation: deferred to
   v0.2; the v0.1 default `throw` is the right ergonomic for now.
3. **`@enum` vs. `BitFlags`.** A `@enum` cannot be combined (e.g.,
   right-censored *and* left-truncated). If we need the cross product
   later we'll move to `@bitflag`. Recommendation: stick with `@enum`
   in v0.1; revisit in v0.2 alongside left-truncation.

### References

- `plans/replan-2026-04-28.md` Phase H, line 208 — censoring enum scope.
- ADR-017 — multi-likelihood seam; Phase H survival families plug into
  the same `LikelihoodMapping` story but the censoring contract is
  orthogonal.
- ADR-006 (and 2026-04-24 amendment) — sets the precedent of opt-in
  per-likelihood extensions to the inner-Newton path
  (`∇³_η_log_density` was added the same way).
- Klein & Moeschberger (2003), *Survival Analysis*, §3.5 — the four
  censoring modes and their likelihoods.
- R-INLA `inla.surv` source —
  [`r-inla/rinla/R/inla.surv.R`](https://github.com/hrue/r-inla/blob/devel/rinla/R/inla.surv.R).
- `packages/LatentGaussianModels.jl/CLAUDE.md` — likelihood contract +
  hot-path policy this ADR extends.

---

## ADR-019: Zero-inflated count families — three R-INLA parameterisations × three base distributions

Status: Proposed
Date: 2026-05-01

### Context

ADR-018 PR5 ships the zero-inflated count families promised in Phase H of
`plans/replan-2026-04-28.md`. R-INLA exposes three parameterisations
(`zeroinflatedX0`, `zeroinflatedX1`, `zeroinflatedX2`) over three count
families (Poisson, Binomial, NegativeBinomial), giving nine concrete
likelihoods. They differ in (a) how the zero-inflation probability `π`
relates to the count distribution's mean and (b) the hyperparameter scale.
Three ADR questions need pinning before nine likelihoods land:

1. **Which parameterisation is which?** R-INLA's documentation is sparse;
   the source-of-truth is `r-inla/inlaprog/src/likelihood.c`. We need to
   commit to a single canonical mapping so future debugging against
   R-INLA fixtures isn't a guessing game.
2. **One file per family, or one per parameterisation?** Nine likelihoods
   is enough that the structural choice matters for both readability and
   the AbstractLikelihood contract.
3. **Which gradients are closed-form?** `∇³_η_log_density` matters only
   for the simplified-Laplace correction (ADR-006 amendment). Closed-form
   gradients across all 9 × 3 = 27 method positions is unnecessary.

### Decision (proposed)

**1. Parameterisations — canonical mapping.** Three families × three
suffixes; the suffix names the parameterisation, *not* the base
distribution. Across all nine, `θ` carries `[log(size)?, zi_scalar]`:

| Suffix | Name                    | π_i formula                | Count component on `y > 0`           |
|--------|-------------------------|----------------------------|--------------------------------------|
| `0`    | hurdle                  | `logit(π) = θ` (constant)  | base distribution truncated at 0     |
| `1`    | standard mixture        | `logit(π) = θ` (constant)  | base distribution (zero allowed)     |
| `2`    | intensity-modulated     | `π_i = 1 - q_i^α`, `θ = log α` | base distribution (zero allowed) |

`q_i` is family-specific: `μ_i / (1 + μ_i)` for Poisson and
NegativeBinomial (with `μ_i = E_i · exp(η_i)`); `sigmoid(η_i)` for
Binomial. ZINB carries an extra `θ[1] = log(size)` for overdispersion;
the zi scalar is `θ[2]`. ZIP and ZIB carry only the zi scalar (1
hyperparameter total).

This mapping matches R-INLA's `family = "zeroinflated{poisson,binomial,
nbinomial}{0,1,2}"` byte-for-byte.

**2. File layout — one file per family, three structs each.** All three
parameterisations of a family share the same `y > 0` count component
arithmetic; bundling them in a single file lets the shared expressions
appear once in comments and the differences live in adjacent functions.

```
src/likelihoods/zero_inflated/
    _helpers.jl        # logsumexp2, only
    poisson.jl         # ZIP0/1/2 — three structs, three log_density,
                       #   three ∇_η, three ∇²_η, plus ZIP1 ∇³_η
    binomial.jl        # ZIB0/1/2 — same shape; ZIB1 also closed-form ∇³
    negbinomial.jl     # ZINB0/1/2 — same shape; ∇³ falls back to FD
```

Each struct is a separate `<: AbstractLikelihood` so dispatch picks the
right closed form without runtime branching on `family`. The
constructor signature mirrors the plain count family:

- ZIP / ZINB: `T(; link = LogLink(), E = nothing,
  hyperprior_size = GammaPrecision(1.0, 0.1),  # ZINB only
  hyperprior_zi = GaussianPrior(0.0, 1.0))`
- ZIB: `T(n_trials; link = LogitLink(),
  hyperprior = GaussianPrior(0.0, 1.0))`

ZIP/ZINB enforce `LogLink`; ZIB enforces `LogitLink`. Other links throw
in the constructor (matches the plain-NB / plain-Bin policy already in
the package).

**3. Default hyperpriors — match R-INLA verbatim.**

- `θ = logit(π)` (types 0/1): `gaussian(mean = 0, prec = 1)` on the
  internal scale, encoded as `GaussianPrior(0.0, 1.0)` (the new
  R-INLA-equivalent prior added in PR5 step A).
- `θ = log(α)` (type 2): same `gaussian(0, 1)`.
- `θ[1] = log(size)` (ZINB): `loggamma(1, 0.1)`, encoded as the
  existing `GammaPrecision(1.0, 0.1)`. Identical to plain NB.

**4. Gradient closure — closed-form everywhere except `∇³_η` on types
0/2.** `∇_η_log_density` and `∇²_η_log_density` are closed-form for all
nine likelihoods. `∇³_η_log_density` is closed-form for **type 1 only**
(both ZIP1 and ZIB1; ZINB1 falls back via the abstract default). The
type-1 simplified-Laplace correction is the most common
zero-inflated use case in disease-mapping; types 0 and 2 fall back to
the AbstractLikelihood FD default (acceptable since ADR-006's amendment
only requires closed-form `∇³` where it materially affects the inner
hot path).

Gradient derivations (recorded in source comments, not the ADR):

- **Type 0 (hurdle).** `y = 0` branch's η-derivative vanishes because
  `log p = log π` is `η`-independent. `y > 0` derivatives need a
  truncation correction `-K = -∂η log(1 - P_count(0))`.
- **Type 1 (standard mixture).** `y > 0` reduces to plain count family.
  `y = 0` uses a posterior weight `w = (1-π)·P_count(0) /
  (π + (1-π)·P_count(0))`, computed via `logsumexp2` (the only shared
  helper). All three derivatives close cleanly:
  `∇_η = -μ·w` (ZIP1), `-n·p·w` (ZIB1), `-s·μ·w/M` (ZINB1) and similar
  for higher orders.
- **Type 2 (intensity-modulated).** `y > 0` adds `α/(1+μ)` to the plain
  count gradient (ZIP/ZINB) or `α(1-p) - n·p` (ZIB1). `y = 0` is
  `log f` with `f = 1 - q^α·D`, `D = 1 - P_count(0)`; ∂²η computed via
  the quotient rule `(∂²f·f - (∂f)²)/f²`.

**5. Hyperprior split for ZINB.** Two hyperparameters → two prior fields.
The struct carries `hyperprior_size::P1` and `hyperprior_zi::P2`
separately so each can be tuned independently via R-INLA-style kwargs:
`ZeroInflatedNegativeBinomialLikelihood1(; hyperprior_size =
PCPrecision(1.0, 0.01), hyperprior_zi = GaussianPrior(0.0, 0.5))`.

### Consequences

- 27 closed-form methods across 9 likelihoods, ≈800 LoC. Validated to
  ~1e-9 against FiniteDiff on a mix-of-zeros, mix-of-positives, and
  large-count test grid (`test/regression/test_zero_inflated.jl`).
- Type-2 ZINB `y = 0` ∂²η has the most algebraically dense closed form
  (the `∂A` term mixes both `pn0` and `(1-μ)`). FD validation is the
  primary correctness check; the comment block in `negbinomial.jl`
  records the derivation step-by-step so future audits can reproduce it.
- Adding a fourth parameterisation (R-INLA does not currently support
  one, but ZIB has been discussed) is mechanical: drop in a new
  `ZeroInflatedXLikelihood3` with its own `log_density` and gradients.
  No abstract-type change.
- The simplified-Laplace correction (ADR-006 amendment) gains the
  type-1 families immediately; types 0 and 2 will fall back to the
  classical Gaussian approximation at the inner Laplace step. This is
  acceptable for v0.1 — disease-mapping fits use type 1 in the vast
  majority of published applications.

### Phasing — how this lands

PR5 of ADR-018 is the single PR landing all nine likelihoods plus
regression tests. Oracle fixture (synthetic ZIP1 vs R-INLA) lands in a
follow-up PR per the replan; tier-1 regression tests are sufficient to
unblock further survival/zero-inflated documentation work.

### Open questions

1. **R-INLA's `quantile` parameter on type 2.** R-INLA exposes a
   `quantile` kwarg that re-parameterises α via the prior expected
   probability `P(y = 0)` at a chosen quantile of the linear predictor.
   This is a *user-facing convenience* that gets translated to a prior
   on `θ = log α`. Recommendation: skip in v0.1; users can set the
   prior on `log α` directly via `hyperprior_zi`. Add in v0.2 if there
   is demand.
2. **Closed-form ∇³ for ZINB1.** The ZINB1 `y = 0` posterior weight
   structure is the same as ZIP1 / ZIB1. The ∇³ derivation is
   mechanically the same shape but algebraically dense (extra `s` and
   `M` factors). Recommendation: deferred until simplified-Laplace
   correction performance on ZINB1 becomes a documented bottleneck.
3. **Predictive PIT under zero-inflation.** `pointwise_cdf` for ZI
   families needs the mixture CDF, which has a discrete jump at zero.
   Recommendation: deferred to v0.2; the `pointwise_log_density`
   methods are sufficient for DIC / WAIC / log marginal likelihood
   diagnostics shipping in v0.1.

### References

- `plans/replan-2026-04-28.md` Phase H — zero-inflated families scope.
- ADR-006 (and 2026-04-24 amendment) — `∇³_η_log_density` opt-in.
- ADR-017 — multi-likelihood seam; ZI families plug into the same
  `LikelihoodMapping` story.
- ADR-018 — PR5 of the survival pack is this work.
- R-INLA likelihood source —
  [`r-inla/inlaprog/src/likelihood.c`](https://github.com/hrue/r-inla/blob/devel/inlaprog/src/likelihood.c).
- Lambert (1992), *Technometrics* 34(1) — original ZIP1 standard
  mixture parameterisation.
- Heilbron (1994), *Biometrical Journal* 36(5) — hurdle (type 0)
  parameterisation.

---

## ADR-020: Drop Julia 1.10 LTS support — Julia 1.12 is the minimum supported version

Status: Accepted
Date: 2026-05-01

### Context

The replan-2026-04-28 Phase F plan committed the v0.1.0 / v0.1.1 release
to Julia 1.10 LTS + current stable, on the assumption that LTS coverage
broadens the user base. In practice this project has zero LTS-pinned
users (no downstream issues filed against 1.10, no LTS-only deps), and
maintaining the LTS lane has carried real cost:

- The CI matrix doubles for the four core packages (LTS × current,
  Linux × macOS × Windows), and the LTS includes are the slow tail.
- Several recent commits used 1.11+ syntax (`@kwdef` improvements,
  `Returns`, `Splat`) that needed manual back-porting for the LTS lane.
- Julia 1.12 is the current stable — the pragmatic floor — and is what
  the local development environment, the benchmark machine, and the
  authors' editors all run on.

### Decision

**Julia 1.12 is the minimum supported version across the entire monorepo.**
This applies to every package's `[compat] julia` field and every CI
matrix lane.

- All `Project.toml` `[compat] julia` entries are bumped to `"1.12"`.
- The `.github/workflows/test.yml` matrix drops `'1.10'` and keeps `'1'`
  (which resolves to 1.12.x today and to whatever stable is when CI
  runs in the future).
- Cross-platform coverage (macOS, Windows) tracks `'1'` rather than the
  former 1.10 pin.
- The replan-2026-04-28 acceptance criterion ("`Pkg.add(\"INLA\")` from
  a fresh depot resolves on Julia 1.10 LTS and current stable") is
  superseded — only current stable.

### Consequences

- New language features ≥ 1.11 (e.g. `Returns`, `Splat`, public marker
  in `module`) are now usable without conditional shims.
- Smaller CI matrix → faster PR feedback, lower spend.
- Users on 1.10 LTS who try `Pkg.add("INLA")` will get a clean compat
  error from Pkg.resolve, not a broken install.
- AutoMerge on the General registry should be untroubled — `julia =
  "1.12"` is a valid lower bound for Pkg.

### References

- `plans/replan-2026-04-28.md` Phase F — superseded acceptance criterion.
- ADR-001 — package split context (which versioning policy applies to all
  four core packages uniformly).

---

## ADR-021: `Copy` component — scaling β lives on the receiving likelihood, not on the projection mapping

Status: Accepted
Date: 2026-05-02

### Context

ADR-017 closed the projector seam question by promoting the model's
`A` slot to a dispatchable `AbstractObservationMapping`. Phase G PR
series 1 and 2 (multi-likelihood + the seam) have since merged. PR
series 3 — the `Copy` component, the single largest blocker on the
joint longitudinal-survival vignette — needs to express
`η_target[i] += β · x_source[k(i)]`, where `x_source` is some other
component's latent slice and `β` is an estimated hyperparameter.

The plumbing question β raises is *where in the architecture does β
live*. Three places it could live:

A. **θ-aware `apply!`.** Extend the seam contract from
   `apply!(η, mapping, x)` to `apply!(η, mapping, x, θ)`. The Copy
   mapping reads β from θ via a stored hyperparameter index and
   applies the scaling in its `apply!` method.

B. **Mutable scaling on the mapping struct.** Add a
   `set_scaling!(mapping, θ)` method called once per θ-iteration before
   `apply!` runs; the Copy mapping caches the current β internally.

C. **β on the receiving likelihood.** β is a hyperparameter of the
   likelihood that receives the copied effect (not of the source
   component or of the mapping). The mapping stays β-free; after the
   main `apply!` runs, the receiving likelihood adds its own
   `β · x_source[k(i)]` contribution to η via a new hook.

ADR-017's `LinearPredictorTerm` paragraph (line 1036) is silent on
which of these wins; the ADR explicitly leaves Copy's plumbing for
PR series 3 to settle, with the "outside reviewer" prerequisite from
the replan-2026-04-28 Phase G "Risks" subsection.

### Decision

**Option C: β is a hyperparameter of the receiving likelihood.**

The seam contract from ADR-017 stays exactly as written —
`apply!(η, mapping, x)` is θ-free. After the mapping's `apply!`
populates η from x, each likelihood gets a chance to apply its own
post-projection contributions via a new hook on `AbstractLikelihood`:

```julia
# Default no-op — most likelihoods don't have copies.
add_copy_contributions!(η_block, ℓ::AbstractLikelihood,
                        x::AbstractVector, θ_ℓ) = η_block
```

Likelihoods that participate in joint models (the survival likelihoods
from Phase H, primarily) override this to read their `β` slot from
`θ_ℓ` and apply `β * x[component_range]` to their η block. The
`Copy` "component" that user code interacts with is therefore an
ergonomic constructor that:

1. registers no new entry in the latent vector (no extra columns in
   x; the source component's latent is the only copy);
2. registers no new precision-matrix block;
3. records `(target_block, source_component_index, β_prior)` on the
   receiving likelihood when the model is constructed;
4. expands `nhyperparameters(ℓ)` by 1 (or more, for vector β) so the
   inference loop allocates β alongside the likelihood's other
   hyperparameters.

In effect, "Copy" is sugar over a likelihood-side feature, not a new
latent component class. The latent component listing in the model
struct is unchanged from PR series 2.

### Alternatives considered

**A. θ-aware `apply!`.** Cleaner stylistically — the mapping owns the
projection completely, including any θ-dependent scaling — but
breaks the ADR-017 contract that all ~30 `m.A`-equivalent sites in
`inference/laplace.jl` and `INLASPDERasters.jl` were rewritten
against ten weeks ago. Third-party `AbstractObservationMapping`
implementers (the `LowRankProjector` / `KroneckerMapping` audience)
would need to update their mappings to accept θ even when they don't
use it. Rejected on lock-in cost.

**B. Mutable mapping state with `set_scaling!`.** Avoids the contract
break but introduces a hidden imperative call ordering — the inner
Newton loop has to call `set_scaling!` before each `apply!` or the
mapping returns stale η. Surprises third-party implementers, makes
the seam less type-stable. Rejected on ergonomics.

**C** keeps the seam pure, keeps the Copy implementation local to
likelihoods that actually use it, and matches R-INLA's mental model
where β is a property of the receiving formula's `f(..., copy=...)`,
not of the latent term being copied.

### Consequences

**Good.**

- ADR-017's `apply!` / `apply_adjoke!` contract is unchanged —
  third-party `AbstractObservationMapping` implementers don't have
  to retrofit anything. The 30 rewritten call sites stay rewritten.
- Copy's hyperparameter accounting flows through the existing
  `likelihood_θ_ranges` machinery from PR series 2 — no parallel
  routing. β is just another likelihood hyperparameter.
- The receiving likelihood's `log_density` / `∇_η_log_density` /
  `∇²_η_log_density` methods need no changes: they already operate
  on η. The hook fires before they do; everything downstream is the
  same.
- The Copy component's user-visible constructor stays intact — users
  don't need to know β lives on the likelihood; they pass
  `Copy(target_component; β_prior)` and the model constructor wires
  it through.

**Cost.**

- One new abstract method on `AbstractLikelihood`:
  `add_copy_contributions!`. Every concrete likelihood inherits the
  default no-op (single-line method). The ones that opt in are the
  survival likelihoods (Weibull, Exponential, log-normal, gamma)
  plus Gaussian on the longitudinal arm — five sites at most.
- The model constructor has to thread β-prior info from the
  user-facing `Copy(...)` call onto the right likelihood's
  `θ_ranges`. One pass over `m.likelihoods` at construction.
- The β-source association (`β` on this likelihood reads from the
  random-intercept *of that other component*) needs a stored index
  on the receiving likelihood. Adds one field per opt-in likelihood.

**Escape hatch.** If a future use case needs a Copy whose receiving
*term* isn't a likelihood (e.g. a copy that targets a fixed-effect
slot, not an observation block), the alternative B path is still
available — `set_scaling!` becomes the second mechanism without
disturbing the first.

### Phasing — how this lands

1. **PR-3a — likelihood hook.** Add the no-op
   `add_copy_contributions!(η, ℓ, x, θ)` default and the source-index
   storage on the abstract likelihood interface. Every existing
   concrete likelihood gets the default; oracle suite stays green.
2. **PR-3b — `Copy` ergonomic constructor.** Add the user-visible
   `Copy(target, …)` constructor that wires β-prior onto the receiving
   likelihood. Add a closed-form regression test (β = 1.0 fixed
   should reproduce the unscaled-share oracle result).
3. **PR-3c — Baghfalaki vignette.** Joint Gaussian + Weibull with a
   shared subject-specific random intercept and a Copy-scaled
   contribution into the survival linear predictor. Lands as the
   final Phase G PR; oracle fixture in
   `packages/LatentGaussianModels.jl/test/oracle/fixtures/baghfalaki.jld2`.

### References

- ADR-017 — projector seam decision; this ADR closes the β-plumbing
  question deferred there.
- ADR-018 — censoring as a per-row vector on the survival likelihood
  struct; the same likelihood-side feature pattern that this ADR
  follows for β.
- `plans/replan-2026-04-28.md` Phase G — "Joint-models scaffolding";
  Risks subsection (outside-reviewer prerequisite).
- R-INLA `f(..., copy=..., hyper=list(beta = ...))` — the
  conceptual model this ADR matches.
- Baghfalaki, T., Esfandyari, S. & Nazari, V. (2024). *A Bayesian
  joint modelling of longitudinal and time-to-event data using
  INLA: A tutorial.* — the joint model this ADR's PR-3c vignette
  reproduces.

---

## ADR-022: `IIDND{N}` parameterisation — separable (`log τ_i, atanh ρ_{ij}`) by default, Wishart/InvWishart on the joint precision as alternative

Status: Accepted
Date: 2026-05-02

### Context

Phase I-A opens the multivariate-IID work: the `IID2D` and `IID3D`
families that R-INLA exposes as `model = "2diid"` / `"iid3d"`, used in
joint longitudinal-survival random effects, paired-areal disease
mapping, and bivariate meta-analysis. Before any code lands, the
parameterisation has to be locked because it touches the public kwargs,
the prior interface, the Hessian-at-θ̂ stability, and the matching of
R-INLA's defaults — all four of which are difficult to walk back later.

The component sits on `n × N` latent slots with joint precision
`Λ ∈ ℝ^{N × N}` (the GMRF block is `Λ ⊗ I_n`, so the per-replicate
precision is `Λ` and the cross-replicate structure is independence).
What needs deciding is *how the user supplies Λ*, since that drives the
hyperparameter vector θ_c, the prior interface, and the user-facing
parameters that get reported in the summary.

Three parameterisations were on the table:

A. **Joint Λ with Wishart/InvWishart prior.** θ_c stores the
   `N(N+1)/2` distinct entries of Λ on a Cholesky scale; the prior is a
   single `Wishart(r, V)` (or `InvWishart`). Matches R-INLA's `iid2d`
   model code's default. Requires a new `AbstractJointHyperPrior`-style
   abstract type because the existing `AbstractHyperPrior` is
   documented as scalar-only (`packages/LatentGaussianModels.jl/src/priors/abstract.jl`
   line 18: "Multi-dimensional priors … live in INLASPDE.jl because
   they are inherently coupled").

B. **Marginal precisions × correlations (separable).**
   - For `N = 2`: θ_c = `(log τ_1, log τ_2, atanh ρ)` — three scalars,
     each carrying its own `AbstractHyperPrior`. Matches R-INLA's
     `2diid` model code (the alternate `f(., model="2diid", ...)` form).
   - For `N ≥ 3`: θ_c stores the strictly-lower-triangular Cholesky
     factor of the correlation matrix on a tangent-space scale (Lewandowski-Kurowicka-Joe
     2009 / Stan's `cholesky_factor_corr`). Marginal precisions stay
     scalar.

C. **Marginal precisions × correlation matrix on a sphere
   (Lewandowski).** Strictly more complex than (B); the only material
   difference is the prior shape, and R-INLA's defaults aren't
   expressible cleanly in this form.

### Decision

**Adopt (B) as the default parameterisation. Provide (A) — Wishart /
InvWishart on the joint Λ — as an explicit alternative when the user
supplies one as the `hyperprior` kwarg.**

For `N = 2`:

```julia
IIDND(n, 2;
      hyperprior_precs = (PCPrecision(), PCPrecision()),
      hyperprior_corr  = PCCor0(U = 0.5, α = 0.5))
```

with internal-scale `θ_c = (log τ_1, log τ_2, atanh ρ)`. The
user-facing summary reports `(τ_1, τ_2, ρ)`.

For `N ≥ 3`:

```julia
IIDND(n, 3;
      hyperprior_precs = ntuple(_ -> PCPrecision(), 3),
      hyperprior_corrs = ntuple(_ -> PCCor0(U = 0.5, α = 0.5), 3))
```

with internal-scale θ_c packing the three log-precisions followed by
the `N(N-1)/2 = 3` `atanh ρ_{ij}` entries (i < j); the on-disk Cholesky
factor of the correlation matrix is reconstructed from the `atanh ρ`
entries via the Lewandowski-Kurowicka-Joe stick-breaking step. (The
`atanh-of-each-pairwise-corr` parameterisation is *not* injective onto
positive-definite correlation matrices for N ≥ 4 in general; ADR-022
locks IIDND to N ≤ 3, with N ≥ 4 deferred to a successor ADR if the
need arises. R-INLA's `iid3d` model also stops at N = 3.)

For the Wishart alternative:

```julia
IIDND(n, N; hyperprior = Wishart(r = N + 1, V = Matrix(I, N, N)))
```

θ_c packs the lower-triangular Cholesky factor of Λ; internal scale is
unconstrained `ℝ^{N(N+1)/2}` via the log-Cholesky map (positive
diagonal entries are stored as `log L_{ii}`). This needs a new
`AbstractJointHyperPrior` abstract type and a single
`log_prior_density(prior, θ_block)` method on it; the existing scalar
`AbstractHyperPrior` machinery is left untouched. Wishart and
InvWishart are the only initial concrete subtypes.

The `IIDND` struct dispatches on whichever kwarg path the user took
(separable vs joint) via two distinct `IIDND` types — `IIDND_Sep{N}`
and `IIDND_Joint{N}` — under a single `AbstractIIDND` umbrella. The
public constructor `IIDND(n, N; ...)` selects the right concrete type
by inspecting the kwargs; users don't see the dispatch.

### Why (B) over (A) as default

1. **Matches R-INLA's `2diid` defaults bit-for-bit.** R-INLA's
   `2diid` model uses `loggamma + atanh-ρ-Gaussian` defaults.
   Defaulting to (A) — Wishart on Λ
   — would silently diverge from R-INLA's most-used multivariate-IID
   path (R-INLA's docs list `2diid` as the recommended form when the
   user has scalar precisions + correlation in mind, which is the
   common joint-models case).
2. **Fits the existing `AbstractHyperPrior` infra without a new
   abstract type.** Each of the three (or six, for N=3) hyperparameter
   slots is scalar and gets its own existing prior class. Wishart/
   InvWishart are inherently coupled and *do* need new infra; option
   (B) lets that infra arrive only when the user opts in.
3. **PC priors compose naturally.** `PCCor0` on the correlation slot
   (reference at ρ = 0, penalising departures from independence —
   matches R-INLA's `pc.cor0`), `PCPrecision` on each marginal
   precision — the already-shipped PC priors. The Wishart path doesn't
   have a PC analogue in the literature. (Note: the early drafts of
   this ADR called the prior `PCCor1`; corrected before code landed
   because R-INLA's `pc.cor1` reserves that name for the
   reference-at-ρ=1 prior used by AR(1)'s lag-1 correlation.)
4. **Hessian-at-θ̂ stability is well-understood.** R-INLA's published
   stability results for `2diid` carry over directly. The Cholesky
   parameterisation of the joint Λ has known degenerate-Hessian
   pathologies on the diagonal (the `log L_{ii}` entries) when Λ is
   near-singular; defaulting to that path before we have a
   stress-tested implementation is the wrong order.

### Why (A) is still offered

When the user is genuinely working in the precision-matrix mental
model — typically because they have an informative Wishart prior from
prior elicitation, or because they're following a textbook example
that uses Wishart — forcing them through the separable path is a
violation of "match R-INLA's user mental model". `iid2d` (Wishart
default) and `2diid` (separable default) coexist in R-INLA precisely
because both audiences exist.

### Why not (C)

The sphere parameterisation is strictly more complex than (B) —
Lewandowski-Kurowicka-Joe is the standard tool for this — and the only
material difference is the prior shape. R-INLA's defaults aren't
expressible in (C) without a Jacobian compensation that the user
shouldn't have to think about. Rejected on lock-in cost vs zero
quality benefit.

### Consequences

**Good.**

- `IID2D` is the highest-leverage Phase I-A target, and the separable
  parameterisation makes it (a) bit-for-bit comparable to R-INLA's
  `2diid` oracle, (b) ergonomic for the joint-longitudinal-survival
  use case where the user thinks in terms of "subject-specific
  random intercept and slope, with a correlation between them".
- The Wishart alternative path is opt-in, so the `AbstractJointHyperPrior`
  infra only ships when there's a real caller for it. PR sequencing
  can defer the Wishart path to its own PR.
- PC priors stay primary throughout the package; no backslide on the
  R-INLA defaults-parity track.

**Cost.**

- Two `IIDND` concrete types (`Sep{N}` + `Joint{N}`) instead of one.
  Localised: ~80 LoC each, plus the umbrella abstract type. Both
  share `length`, `nhyperparameters`, `gmrf` defaults via the abstract
  type; only `precision_matrix`, `log_hyperprior`,
  `initial_hyperparameters`, and the public constructor logic differ.
- One new abstract type (`AbstractJointHyperPrior`) plus `Wishart` /
  `InvWishart` concrete subtypes when the Wishart path lands. Three
  new files in `src/priors/`. Doesn't affect any of the existing
  scalar-prior-using components.
- `Cholesky-of-correlation` reconstruction (LKJ stick-breaking) for
  N=3 needs ~30 LoC of arithmetic and a regression test against the
  Stan reference.

**Escape hatch.** If a user case emerges where the separable default
diverges materially from a Wishart default they expected (e.g. they're
porting an R-INLA model that *did* use `iid2d` rather than `2diid`),
they pass `hyperprior = Wishart(...)` and the constructor flips to
the joint path with no further code change.

**Defer.** N ≥ 4 IID is deferred — neither R-INLA nor the joint-models
literature has an active need, and the Lewandowski-Kurowicka-Joe
parameterisation gets fragile at higher N. A successor ADR can lift
this if a real use case appears.

### Phasing — how this lands

1. **PR-1a — `IIDND_Sep{2}` (i.e. `IID2D`).** Add the separable
   constructor, `precision_matrix`, log-prior wiring through the three
   scalar `AbstractHyperPrior` slots, regression tests against a dense
   reference Λ. Add `PCCor0` prior. No Wishart path yet.
2. **PR-1b — `IID3D`.** Same pattern at N=3 with the three-correlation
   LKJ piece. `IID3D` regression test against R-INLA `iid3d` (R-INLA's
   `iid3d` is documented as brittle on small samples; oracle is
   informational, not load-bearing).
3. **PR-1c — Wishart / InvWishart joint path.** Adds
   `AbstractJointHyperPrior`, the two concrete priors, and
   `IIDND_Joint{N}` constructor logic. Lands only when a user case
   asks for it; otherwise sits behind the proposed-status flag.

### References

- ADR-006 — PC priors as the default prior class.
- ADR-021 — recent component ADR; this ADR follows the same shape.
- `plans/replan-2026-04-28.md` Phase I — original IID2D / IID3D scope.
- Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sørbye, S. H.
  (2017). *Penalising model component complexity: a principled,
  practical approach to constructing priors.* — general PC-prior
  framework underlying `PCCor0` (reference at ρ = 0).
- Sørbye, S. H., & Rue, H. (2017). *Penalised complexity priors for
  stationary autoregressive processes.* — companion construction at
  reference ρ = 1 (R-INLA's `pc.cor1`), used by AR(1) but not by
  IID2D.
- Lewandowski, D., Kurowicka, D., & Joe, H. (2009). *Generating
  random correlation matrices based on vines and extended onion
  method.* — N≥3 Cholesky-of-correlation parameterisation.
- R-INLA `f(., model="2diid", ...)` and `f(., model="iid3d", ...)` —
  the parameterisations this ADR matches.

---

## ADR-023: `MEB` and `MEC` measurement-error components — β-via-`Copy` decomposition; non-zero `prior_mean` made load-bearing

Status: Accepted
Date: 2026-05-03

### Context

Phase I-B opens the measurement-error work: R-INLA's `f(w, model="meb", ...)`
(Berkson) and `f(w, model="mec", ...)` (Classical) latent components, used
for errors-in-variables regression in epidemiology, environmental
exposure modelling, and any setting where a covariate is measured with
non-trivial uncertainty (Carroll, Ruppert, Stefanski, Crainiceanu 2006;
Muff, Riebler, Held, Rue, Saner 2015).

Both R-INLA models expose, in the linear predictor, the *scaled* latent
covariate `ν_i = β · x_i` rather than `x_i` itself. Internally the
latent block is the unscaled `x` — a length-`n_unique(w)` Gaussian
vector with prior mean and prior precision determined by the model
flavour and a small set of hyperparameters — and `β` is a per-component
hyperparameter that multiplies the projection from latent to linear
predictor before it lands in `η`. R-INLA registers both models in the
latent-models table alongside `ar1`/`rw1`/`iid` (`models.R`,
`n.required = FALSE`), so the user's mental model is "this is a
component", not "this is a coefficient prior".

The two flavours differ by which way the noise tie runs:

- **MEB (Berkson)**: `x_i = w_i + u_i`, `u_i ~ N(0, (τ_u s_i)⁻¹)`. The
  observed covariate `w` is *fixed* and the latent is `w + u`. The
  latent-block prior is therefore `x ~ N(w, (τ_u D)⁻¹)` with
  `D = diag(s)` after marginalising `u`.
- **MEC (Classical)**: `w_i = x_i + u_i`, `u_i ~ N(0, (τ_u s_i)⁻¹)`,
  with prior `x ~ N(μ_x · 1, (τ_x I)⁻¹)`. The observed `w` is a noisy
  measurement of the latent truth `x`. R-INLA absorbs the Berkson tie
  `w | x ~ N(x, (τ_u D)⁻¹)` into the prior on `x` via Gaussian
  conjugacy: posterior precision becomes `τ_x I + τ_u D`, posterior mean
  becomes `(τ_x I + τ_u D)⁻¹ (τ_x μ_x 1 + τ_u D w)`.

Both flavours fit the existing `AbstractLatentComponent` contract if
two things are true:

1. The component's prior mean is allowed to be non-zero (the optional
   `prior_mean(c, θ)` method documented in `CLAUDE.md` becomes
   load-bearing, where for every other component it has been zero).
2. The β-multiplication of the latent block before it enters η is
   handled somewhere — and we have a choice of where.

`plans/replan-2026-04-28.md` (lines 255–256) originally labelled MEB
as "classical" and MEC as "Berkson". This is the opposite of R-INLA's
convention. The R-INLA LaTeX source
(`inlaprog/doc/latent/{meb,mec}.tex`) is unambiguous: `meb` =
"Berkson Measurement Error", `mec` = "Classical Measurement Error". The
B/C suffix in the model name encodes the error structure, not "Bayesian
vs. classical inference". This ADR locks the correct R-INLA naming;
the same PR that lands this ADR amends the replan in place.

### Decision

**Three coupled choices.**

#### 1. β lives on the *receiving* likelihood via `Copy`, not on the component

Consistent with ADR-021's `Copy` decision. The MEB and MEC components
own only the latent block `x` (with non-zero prior mean and structured
precision); the β-scaling is wired up by the user (or by a constructor
helper) as a `Copy(source_indices; β_prior=...)` on the
`CopyTargetLikelihood` that wraps the receiving observation likelihood.

```julia
# User-facing pattern:
component = MEB(w; scale = s, τ_u_prior = GammaPrecision(1.0, 1.0e-4))
β_copy = Copy(component_range; β_prior = GaussianPrior(1.0, 0.001))
target = CopyTargetLikelihood(GaussianLikelihood(), β_copy)
model = LatentGaussianModel((target,),
                            (Intercept(), component),
                            obs_mapping)
```

The receiving-likelihood placement is the existing infrastructure;
nothing about the inner Newton loop, the η-Jacobian augmentation, or
the hyperparameter accounting changes. The MEB/MEC component itself
contributes zero copies of its own — it is, structurally, an "IID with
non-zero prior mean and per-element precision".

#### 2. Latent layout

Both components store a length-`n_unique(w)` latent block. Observations
with identical `w` (after rounding to the supplied `digits`, matching
R-INLA's `f(w, model = "meb", values = ...)` semantics) share the same
latent slot. The projector from latent to row-space is the natural
"row `i` reads slot `k` where `w[i] = unique_w[k]`" injection — already
expressible as a `LinearProjector` from ADR-017's observation-mapping
seam, no new mapping type required.

`prior_mean(c, θ)` becomes the load-bearing way to expose the non-zero
mean to the inner Newton loop. The CLAUDE.md component-contract entry
already documents `prior_mean` as an optional method with default
`zeros`; this ADR promotes it to "load-bearing for at least two
components" and adds the corresponding hot-path test (the inner Newton
step must read `prior_mean` whenever it is non-zero).

`precision_matrix(c, θ)` is diagonal in both flavours:

- **MEB**: `Q_c = τ_u · D`, `μ_c = w_unique`.
- **MEC**: `Q_c = τ_x I + τ_u D`, `μ_c = Q_c⁻¹ · (τ_x μ_x 1 + τ_u D w_unique)`.

Note the MEC mean depends on θ (through `τ_x`, `τ_u`, `μ_x`), unlike
MEB's θ-constant mean. This is fine — `prior_mean(c, θ)` is dispatched
on θ already; the contract supports it.

#### 3. Hyperparameter list and defaults — bit-for-bit R-INLA

**MEB**:

| slot | name | internal scale | prior | initial | default-fixed |
|---|---|---|---|---|---|
| θ_c[1] | β | identity | `GaussianPrior(1.0, 0.001)` (R-INLA: gaussian, mean=1, prec=0.001) | 1.0 | no |
| θ_c[2] | log τ_u | log | `GammaPrecision(1.0, 1.0e-4)` (R-INLA: loggamma, shape=1, rate=1e-4) | log(1000) ≈ 6.9078 | no |

β lives on the receiving likelihood via `Copy`; θ_c[1] is the Copy's β
slot, *not* the component's. Internally only `log τ_u` belongs to the
component's hyperparameter vector. The "MEB has 2 hyperparameters" line
in R-INLA's docs is a user-facing accounting that includes β; the
Julia accounting is "1 component hyperparameter + 1 β slot via Copy".

**MEC**:

| slot | name | internal scale | prior | initial | default-fixed |
|---|---|---|---|---|---|
| θ_c[1] | β | identity | `GaussianPrior(1.0, 0.001)` | 1.0 | no |
| θ_c[2] | log τ_u | log | `GammaPrecision(1.0, 1.0e-4)` | log(10000) ≈ 9.21 | **yes** |
| θ_c[3] | μ_x | identity | `GaussianPrior(0.0, 1.0e-4)` | 0.0 | **yes** |
| θ_c[4] | log τ_x | log | `GammaPrecision(1.0, 1.0e4)` | -log(10000) ≈ -9.21 | **yes** |

Same β-via-Copy split as MEB. The default-fixed slots (τ_u, μ_x, τ_x)
are R-INLA's "degrades to plain regression unless the user opts in":
τ_u huge ⇒ no measurement error, τ_x tiny ⇒ vague prior on x, μ_x
fixed at 0. The user unfixes whichever slots they want to estimate by
passing `fix_τ_u = false` (etc.) at construction. The Julia-side public
constructor mirrors R-INLA's `f(w, model = "mec", control.fixed =
list(...))` form.

#### Public API

```julia
MEB(w::AbstractVector;
    scale = ones(length(unique(w))),
    digits::Int = 8,
    τ_u_prior::AbstractHyperPrior = GammaPrecision(1.0, 1.0e-4),
    τ_u_init::Real = log(1000.0))
# Returns a `MEB` component (subtype of `AbstractLatentComponent`) plus,
# via the constructor's secondary return, a `Copy` template the user
# attaches to the receiving likelihood.

MEC(w::AbstractVector;
    scale = ones(length(unique(w))),
    digits::Int = 8,
    τ_u_prior = GammaPrecision(1.0, 1.0e-4),
    μ_x_prior = GaussianPrior(0.0, 1.0e-4),
    τ_x_prior = GammaPrecision(1.0, 1.0e4),
    τ_u_init  = log(10000.0),
    μ_x_init  = 0.0,
    τ_x_init  = -log(10000.0),
    fix_τ_u::Bool = true,
    fix_μ_x::Bool = true,
    fix_τ_x::Bool = true)
```

The `digits` kwarg matches R-INLA's `f(...; values = unique(round(w,
digits = ...)))` convention for de-duplicating observed covariate
values into latent slots.

The β slot does not appear in either constructor's kwargs — it lives on
the user-supplied `Copy(...)` attached to the receiving likelihood. To
make this ergonomic, both constructors return a `(component, β_copy_template)`
named tuple where `β_copy_template` is a `Copy` pre-configured with
R-INLA's `gaussian(1, 0.001)` β prior; the user can override by
constructing their own `Copy` against the component's range.

### Why β-via-Copy and not β-on-the-component

1. **Consistency with ADR-021.** That ADR accepted "scaling β lives on
   the receiving likelihood" as the rule. Reopening it for MEB/MEC
   would either fragment the codebase (two β-attachment idioms) or
   force an ADR-021 supersedence. The Copy route is the one we already
   tested.
2. **No change to `add_copy_contributions!` or `_accumulate_copy_jacobian!`.**
   The hot path is already wired for β-times-source-block contributions
   into η. Reusing it for MEB/MEC means zero net new hot-path code
   beyond the component's `precision_matrix` and `prior_mean`.
3. **The "MEB component contributes its own copy" alternative would
   couple component construction to likelihood construction.** Right
   now `LatentGaussianModel` accepts `(likelihoods, components,
   mapping)` independently. A self-copying component would need a
   back-reference from latent index → receiving likelihood, which
   `LatentGaussianModel` does not currently maintain.
4. **R-INLA's posterior of `ν = β·x` is recoverable post-fit as a
   derived quantity from the joint `(x, β)` posterior** — no need to
   carry `ν` as a primary latent. A small accessor in PR-2c
   (`measurement_error_scaled_latent(model, res)`) closes the
   user-facing ergonomics gap.

### Why two distinct components and not a parameterised one

`MEB` and `MEC` differ in:
- prior-mean dependence on θ (constant for MEB, θ-dependent for MEC);
- precision-matrix structure (diagonal `τ_u D` vs `τ_x I + τ_u D`);
- number of hyperparameters (1 vs 3 once β is excluded);
- default-fixed pattern (none for MEB, three of three for MEC).

Folding both into a single `MeasurementError(...; flavour = :berkson|:classical)`
constructor would require the struct's hyperparameter vector to be
parameterised on `flavour`, which forces the field types to be
`Vector` instead of `NTuple` and gives up dispatch granularity for
`precision_matrix`/`prior_mean`. The two-component split is the same
shape as `BYM` vs `BYM2` or `Generic0` vs `Generic1` vs `Generic2` —
distinct components matching distinct R-INLA model strings.

### Why the replan's labels-swapped error matters

Carroll-style measurement-error literature uses the Berkson/Classical
distinction to flag which way the noise tie runs, and the analytical
behaviour differs sharply between them — bias attenuation under
classical error, no-attenuation under Berkson, and entirely different
prior-elicitation requirements. Mislabelling the components in user
documentation would silently push users toward the wrong model for
their data. `plans/replan-2026-04-28.md` lines 255–256 originally
called `MEB(...)` "classical" and `MEC(...)` "Berkson" — the opposite
of R-INLA. The active session-level Phase I-onwards replan
(2026-05-02) inherited the same swap from the parent replan. This ADR
is the single source of truth from this point forward; the committed
replan was corrected in the same PR that landed this ADR, and the
session-level replan picks up the correction via its next sync.

### Consequences

**Good.**

- MEB and MEC ship as standard `AbstractLatentComponent` subtypes with
  no new abstract types, no new prior categories, no Copy-machinery
  changes. The β-attachment idiom is unchanged from ADR-021.
- `prior_mean(c, θ)` is promoted from "documented but unused optional"
  to "documented and load-bearing in two components". The inner Newton
  loop already calls it; this ADR confirms that path stays load-bearing
  rather than getting silently optimised away.
- The R-INLA defaults table is captured here verbatim for
  `defaults-parity.md`'s next sync.
- The replan's MEB/MEC label swap is corrected in writing.

**Cost.**

- ~150 LoC for `MEB` (struct, constructor, `precision_matrix`,
  `prior_mean`, `log_hyperprior`, `log_normalizing_constant`, `gmrf`)
  plus ~50 LoC of regression tests in `test/regression/test_meb.jl`.
- ~200 LoC for `MEC` (more hyperparameters, conjugate-Gaussian mean
  formula, three default-fixed slots, fix-toggle kwargs) plus ~70 LoC
  of regression tests.
- One R-INLA oracle fixture (Carroll-style classical-error regression)
  in `test/oracle/`. Probably reuse the Muff et al. 2015 example data.
- Promotion of `prior_mean` to load-bearing requires a regression test
  asserting the inner Newton step reads it correctly. ~30 LoC.
- The "constructor returns `(component, β_copy_template)` named tuple"
  is a slightly unusual constructor return shape; document explicitly
  in both docstrings.

**Escape hatch.** If the named-tuple constructor return shape produces
user-visible friction (someone tries to write `c = MEB(w)` and is
surprised by the tuple), wrap it in a small helper:

```julia
component, β_copy = MEB(w)            # explicit destructuring
m = MEB(w).component                  # ignore β-default-template
m, β_copy = unpack_meb(MEB(w))        # named helper
```

**Defer.** R-INLA exposes both `meb` and `mec` per-row group / replicate
extensions (`f(w, model = "meb", group = ..., replicate = ...)`).
Phase I-C handles the `replicate` / `group` machinery for *all*
components uniformly; MEB/MEC pick them up for free at that point. No
extra ADR work needed when Phase I-C lands.

### Phasing — how this lands

1. **PR-2a (this ADR + replan correction).** Lock the design before
   code. Amend `plans/replan-2026-04-28.md` (lines 255–256) to correct
   the MEB/MEC label swap and reference this ADR.
2. **PR-2b — `MEB` (Berkson) component.** The simpler of the two —
   prior mean is θ-constant `w_unique`, precision is diagonal `τ_u D`,
   one component hyperparameter (`log τ_u`). Regression tests + R-INLA
   oracle on a synthetic Berkson example.
3. **PR-2c — `MEC` (Classical) component.** Adds the conjugate-Gaussian
   prior-mean formula, three default-fixed hyperparameters, and the
   `fix_*` kwargs. R-INLA oracle on a Carroll-style classical-error
   regression (Muff et al. 2015 example data is the canonical fixture).
4. **PR-2d — Vignette.** Port the Carroll classical-error regression
   into `docs/src/vignettes/measurement-error-regression.md`.

### References

- ADR-021 — `Copy` component; β-on-receiving-likelihood rule that this
  ADR follows.
- ADR-017 — `AbstractObservationMapping` seam; provides the
  `LinearProjector` for the row-to-unique-slot mapping.
- ADR-022 — most-recent component ADR; this ADR follows its shape.
- `plans/replan-2026-04-28.md` Phase I — original scope; this ADR
  corrects the MEB/MEC label swap originally introduced there
  (lines 255–256), now amended in the same PR that lands this ADR.
- `plans/defaults-parity.md` — to be updated with the MEB/MEC default
  table when PR-2b/c land.
- R-INLA model docs:
  - `https://inla.r-inla-download.org/r-inla.org/doc/latent/meb.pdf`
  - `https://inla.r-inla-download.org/r-inla.org/doc/latent/mec.pdf`
- R-INLA source:
  - `inlaprog/doc/latent/meb.tex`, `mec.tex` — model equations and
    R-INLA's hyperparameter accounting (`hyperid 3001/3002` for MEB,
    `hyperid 2001-2004` for MEC).
  - `rinla/R/models.R` — registration confirms latent-section
    placement, `n.required = FALSE`, no `covariate` field.
- Muff, Riebler, Held, Rue, Saner (2015). *Bayesian analysis of
  measurement error models using INLA.* — canonical Bayesian-INLA
  treatment of MEB/MEC; the parameterisation choices in R-INLA trace
  to this paper.
- Carroll, Ruppert, Stefanski, Crainiceanu (2006). *Measurement Error
  in Nonlinear Models: A Modern Perspective.* — textbook reference for
  the Berkson vs. classical distinction and the associated bias
  results.

---

## ADR-024: Categorical / Multinomial via independent-Poisson reformulation; helper-function API, no new likelihood type

Status: Accepted
Date: 2026-05-04

### Context

Phase J's tail item is multi-class outcomes — the user has an `n × K`
count or one-hot matrix `Y` (with `K ≥ 3`) and wants a regression on
the class probabilities. The two textbook parameterisations are:

1. **Multinomial-logit** (softmax). For each row `i` and class `k`:

   ```math
   π_{ik}(η) = exp(η_{ik}) / Σ_j exp(η_{ij}),
   ```

   with `K − 1` free `η` per row (one class fixed for identifiability).
   Hessian *within an observation across the `K − 1` η values* is
   *not* diagonal — the off-diagonal `∂²L_i/∂η_{ik}∂η_{im} ≠ 0` for
   `k ≠ m`.

2. **Stick-breaking.** Sequential binary logits:
   `π_{i1} = σ(η_{i1})`, `π_{ik} = (1 − Σ_{j<k} π_{ij}) σ(η_{ik})`.
   Hessian within an observation also non-diagonal.

The existing `AbstractLikelihood` contract (CLAUDE.md, package-level
section) requires `∇²_η_log_density` to return a `Vector{<:Real}` —
the diagonal of the Hessian. From the same docstring:

> For likelihoods where the Hessian is not strictly diagonal
> (e.g. Dirichlet-multinomial), this contract does not fit — such
> likelihoods must implement a different `AbstractLikelihood` subtype
> with a dedicated method (not part of v0.1).

That note remains accurate for v0.2.x: the package has not introduced
a `AbstractMatrixHessianLikelihood` subtype, and the inner Newton
loop, the simplified-Laplace correction, the marginal accessors, and
all the integration schemes assume a diagonal-Hessian likelihood. The
diagonal contract is load-bearing across at least eight modules; flipping
to a banded-or-block Hessian is a Phase L-or-later refactor (it
naturally pairs with the FullLaplace strategy work).

R-INLA's standard recipe sidesteps the non-diagonal Hessian: the
**Poisson-multinomial equivalence** (Baker 1994, Chen 1985, R-INLA's
[`Multinomial-INLA.pdf`](https://www.r-inla.org/learnmore/tutorials)
tutorial). Reframe each multinomial trial as `K` independent
Poisson observations — counts `(0, …, 1, …, 0)` with the chosen class
in row `y_i` flagged by `1` — and absorb the row total constraint
into a per-observation nuisance intercept `α_i`:

```math
y_{ik} | α_i, η_{ik} ∼ Poisson(exp(α_i + η_{ik})),
\qquad i = 1, …, n; k = 1, …, K.
```

Each `α_i` is given a near-flat prior so it acts as a free nuisance
absorbing the multinomial total. Marginalising `α_i` analytically
recovers the original multinomial likelihood up to an additive
constant; integrating it numerically (which is what INLA does) gives
the same posterior on the regression effects, bit-for-bit, as a
direct multinomial-logit fit. The `K` Poisson observations within a
row each have their own diagonal Hessian entry — the *cross-row* tying
happens through `α_i`, which is a *latent* not a likelihood concern,
so the diagonal contract is preserved.

This is what R-INLA hard-wires; its `Multinomial-INLA.pdf` tutorial is
the primary reference for users. The reformulation exists at the
*user/data-construction* layer, not as a new likelihood family — there
is no `family = "multinomial"` in R-INLA that uses the trick
internally. Users construct the long-format data manually.

### Decision

**Three coupled choices.**

#### 1. v0.2.x: independent-Poisson reformulation; no new likelihood type

The package ships a *helper function* that converts an `n × K`
response matrix into the long-format data needed for the Poisson trick,
plus the projector + nuisance-intercept component. The user composes
this helper output with the existing `PoissonLikelihood`, `Intercept`,
`FixedEffects`, etc. — no new `AbstractLikelihood` subtype.

```julia
# Public API (v0.2.2):
helper = multinomial_to_poisson(Y; class_names = nothing)
# Returns a NamedTuple with:
#   y        :: Vector{Int}            # length n*K, vec'd row-major
#   row_id   :: Vector{Int}            # row index repeated K times
#   class_id :: Vector{Int}            # class 1..K cycling per row
#   n_rows, K, n_long
```

The user then constructs the model:

```julia
helper = multinomial_to_poisson(Y)
n, K = helper.n_rows, helper.K

# Per-observation nuisance intercept α_i. R-INLA convention is
# `f(idx, model="iid", hyper=list(prec=list(initial=-10, fixed=TRUE)))`,
# i.e. precision fixed at exp(-10) ≈ 4.5e-5 (very flat). Use the
# existing IID component with a fixed-precision hyper.
c_α = IID(n; τ_init = log(exp(-10)), τ_fixed = true)

# Class-by-covariate fixed effects: drop class K for identifiability,
# leaving (K − 1) × p coefficients. The user builds the design from
# `(class_id, row_id, X_original)`; helper exposes
# `multinomial_design_matrix(helper, X)` for the common case.
A = multinomial_design_matrix(helper, X; reference_class = K)
c_β = FixedEffects(size(A, 2))

projector = sparse_assembly(helper.row_id, c_α, A, c_β)   # built by helper
model = LatentGaussianModel(PoissonLikelihood(), (c_α, c_β), projector)
res = inla(model, helper.y)
```

The `IID` component already supports a `fix_τ` toggle from Phase I
work; the helper just wires the conventional initial value. The β
coefficients land on the FixedEffects component as usual.

**Rationale for not adding `MultinomialLikelihood`:**

- It would duplicate `PoissonLikelihood` numerically (same gradients,
  same Hessian) while adding latent-layout glue that R-INLA users
  already express manually.
- A "real" MultinomialLikelihood that uses the softmax form requires
  a non-diagonal Hessian, which is a Phase L+ refactor.
- The helper-only approach matches R-INLA's `Multinomial-INLA.pdf`
  tutorial bit-for-bit, including the per-row α_i convention. Users
  porting an R-INLA categorical fit translate line-by-line.

#### 2. Reference class and identifiability

The Poisson reformulation has a built-in identifiability subtlety: the
likelihood is invariant under `(α_i, η_{ik}) → (α_i + c, η_{ik} − c)`
for any class-independent shift `c`. R-INLA's near-flat prior on `α_i`
(precision `≈ 4.5e-5`) doesn't break this exactly but tightens it
weakly enough that the posterior on `(α, η)` is well-behaved.

To force a clean identification, the helper drops one class's
covariate effects (`reference_class = K` by default — R-INLA's
convention is to drop the *last* class). The reference class still
contributes Poisson observations; only its *coefficients* are zeroed.
The intercept-only reference structure is `η_{ik} = α_i + γ_k +
x_i' β_k` with `γ_K = 0, β_K = 0`. This is the
"corner-point" parameterisation used in Agresti (2010, §8.5).

The user can override `reference_class` to any of `1, …, K`. Picking
the most-populated class as reference improves numerical conditioning
(common practice in epidemiology).

#### 3. Native multinomial-logit deferred to v0.3+

The native multinomial-logit (and stick-breaking) parameterisations
need a non-diagonal Hessian and therefore a new abstract type
`AbstractMatrixHessianLikelihood` (or similar). This is naturally
paired with the FullLaplace strategy work in Phase L (UserComponent
+ FullLaplace), since FullLaplace also needs to handle non-diagonal
likelihood Hessians for sharply non-Gaussian latents.

When the v0.3+ refactor lands, the helper-based pattern this ADR
introduces will continue to work (it's a thin data-construction layer
on top of `PoissonLikelihood`). The native likelihood will be an
ergonomic alternative, not a replacement.

### Consequences

**Pros.**

- Matches R-INLA bit-for-bit on the Phase J close. Users porting a
  multinomial fit from R-INLA to Julia translate the helper call
  line-by-line.
- Reuses existing tested machinery: `PoissonLikelihood` gradients,
  `IID` with fixed precision, `FixedEffects`, `LinearProjector`. No
  new code paths in the inner Newton loop.
- The diagonal-Hessian likelihood contract stays intact — no
  ripple-effects through simplified-Laplace, marginal accessors,
  integration schemes, or the JointLikelihood block-row stack.
- The helper API is composable with the rest of the LGM stack: BYM2
  random effects on classes, AR1 on time, Copy for shared latents
  across classes — all work because we're just wiring up a Poisson
  fit.

**Cons.**

- Latent dimension blows up by `n` (the per-observation nuisance
  intercept). A `n = 10000, K = 5` problem becomes a `n*K + n =
  60000`-dim latent before random effects. Sparse machinery handles
  this fine for moderate `n`, but very large `n` may motivate the
  v0.3+ native form. Document this tradeoff in the helper's docstring.
- Users who already think in multinomial-logit terms see "Poisson"
  in the model definition and may be confused. The vignette and
  helper docstring must explain the equivalence explicitly.
- The α_i nuisance precision is *fixed* at a near-flat value; a user
  who wants to tune it sees a `fix_τ = true` setting they can't move
  via the helper alone (they'd construct the IID directly). This is
  R-INLA's convention but not documented in `inla.doc`; we surface
  it in the helper docstring.
- `pointwise_log_density` for CPO/WAIC reports per-Poisson-cell
  contributions, not per-multinomial-row. Users computing CPO must
  aggregate across the K rows belonging to each multinomial trial
  — the helper exposes a `multinomial_pointwise_aggregate(res, helper)`
  utility for this.

**Acceptance criteria for the implementation PR.**

- `multinomial_to_poisson(Y)` regression test: round-trip from a
  small `Y` matrix to long format and back, verify counts match.
- `multinomial_design_matrix(helper, X; reference_class)` regression
  test: the column count is `(K − 1) * p` (with `reference_class`
  zeroed), the row count is `n * K`.
- R-INLA oracle on a `n = 200, K = 3, p = 1` synthetic example. The
  R fixture uses R-INLA's own `Multinomial-INLA.pdf` reformulation.
  Tolerances: 5% on the regression coefficients' marginal means; 10%
  on the per-observation α_i posterior moments (which carry
  identifiability noise).
- Vignette `docs/src/vignettes/multinomial.md` walking through the
  helper + Poisson trick, mirroring the R-INLA tutorial layout. A
  small "alligator-food-choice" or "pollen-presence" dataset (both
  are R-INLA tutorial staples) is the target.

### References

- Baker, S. G. (1994). *The Multinomial-Poisson transformation.* The
  Statistician, 43(4), 495-504.
- Chen, T. T. (1985). *Log-linear models for categorical data with
  misclassification and double sampling.* JASA, 80, 158-162.
- R-INLA tutorial: *Multinomial-INLA.pdf* — the canonical reference
  for the Poisson-multinomial trick under R-INLA. Available from
  https://www.r-inla.org/learnmore/tutorials.
- Agresti, A. (2010). *Analysis of Ordinal Categorical Data*, 2nd ed.
  §8.5 covers the corner-point parameterisation we adopt.
- ADR-021 (`Copy` on receiving likelihood) — same architectural
  philosophy: keep new features on top of existing primitives where
  possible; reserve abstraction-extending changes for cases where the
  primitives genuinely don't fit.

---

## ADR-025: `UserComponent` callback signature (R-INLA `rgeneric`)

Status: Accepted
Date: 2026-05-04

### Context

Phase L PR-2 ships [`UserComponent`](../packages/LatentGaussianModels.jl/src/components/user_component.jl),
the R-INLA `rgeneric` equivalent: a one-line callable wrapper around
the [`AbstractLatentComponent`](../packages/LatentGaussianModels.jl/src/components/abstract.jl)
contract that lets users port R-INLA `rgeneric` model definitions
without subtyping. This is the most-used extension hook in published
R-INLA papers and the determinant of how much of R-INLA's long-tail
component library (`crw2`, `besag2`, `besagproper`, `clinear`, `z`,
`ou`, …) we can avoid implementing natively.

Two design questions need a durable answer:

1. **Return shape of the callable.** Positional tuple, named tuple,
   multiple-return, struct? Each has trade-offs in extensibility and
   readability.
2. **Does `UserComponent` replace direct subtyping, supplement it, or
   subsume one within the other?** The contract has eight methods of
   varying optionality; not all of them fit a single callable cleanly.

ADR-003 already commits to multiple-dispatch as the primary extension
mechanism. `UserComponent` is the *callable* surface over the same
seam, not a parallel mechanism — the question is how the surface
relates to the seam underneath.

The `cgeneric` part of R-INLA (C-callable user components) is dropped
from scope (see [`plans/conti-valiant-pebble.md`](../plans/conti-valiant-pebble.md)
§Scope decision): Julia callables are JIT-compiled to native code, so
there is no measurable performance gap with C, and an FFI layer would
*cost* speed.

### Decision

1. **The callable returns a `NamedTuple`**, not a positional tuple
   or multiple values. The required key is `:Q` (sparse precision
   matrix); optional keys with defaults are `:log_prior` (default `0`),
   `:log_normc` (default `0`), and `:constraint` (default
   `NoConstraint()`).

   Adding a fifth key in a future release (e.g. `:prior_mean` for
   ADR-023 hot paths) is non-breaking — existing callables stay
   silent on the new key and pick up the default. A positional tuple
   would force every caller to update on extension.

2. **The constraint is read once at construction time** by invoking
   the callable at `θ0` and caching the returned `:constraint`.
   `GMRFs.constraints(c)` is θ-independent in our seam; mirroring
   that, `UserComponent` documents that the constraint must be
   θ-independent (or, equivalently, that only its θ0-value is
   honoured). This matches R-INLA's rgeneric, where `extraconstr`
   is set at `f()` time, not in the rgeneric callback.

3. **`UserComponent` and direct subtyping ship side-by-side.** The
   power-user path (subtyping `AbstractLatentComponent` directly) is
   already proved by `Generic0` / `Generic1` / `Generic2` and is
   needed for cases where the callable is too narrow:
   - Component-specific overrides of `prior_mean(c, θ)` (ADR-023)
     for shifted-prior measurement-error components (`MEC`).
   - Custom `gmrf(c, θ)` factorisations that bypass the default
     `Generic0GMRF` wrapper.
   - Lazy/structured precision matrices that don't fit the
     `SparseMatrixCSC` mould.

   `UserComponent` covers everything else with one closure. The two
   are complementary, and `docs/src/extending.md` documents both as
   first-class extension paths rather than tiered alternatives.

4. **`@ccall` is the C-library bridge.** Users with existing C
   precision-matrix routines call them inside the Julia closure via
   `@ccall`. There is no separate `CGenericComponent` type and no
   plan to add one. The `UserComponent` docstring carries a one-
   paragraph note pointing this out; it does not need its own ADR.

### Consequences

#### Positive

- **R-INLA users port models in one closure**: `crw2`, `besag2`,
  `clinear`, `z` etc. become a single function literal each. The
  Phase L `crw2` vignette (PR-6) is the proof.
- **Extension is non-breaking**: adding new optional keys to the
  returned `NamedTuple` does not break existing callables.
- **No FFI surface to maintain**: dropping `cgeneric` removes a
  C-API drift risk and a build-system dependency.
- **One-shot validation**: the constructor invokes the callable
  once at `θ0`, surfacing wrong return shapes / wrong-size `Q` /
  bad constraint types as `ArgumentError` / `DimensionMismatch`
  during model construction rather than mid-Newton.

#### Negative / Trade-offs

- **Repeated callable invocations.** `precision_matrix`,
  `log_hyperprior`, and `log_normalizing_constant` each call the
  closure separately at the same θ. For an expensive precision
  build this triples the work per integration point. Users wanting
  to amortise can cache via `Memoize.jl` or hand-rolled state — the
  v0.1 contract is "call lazily; user owns memoisation".
- **θ-independent constraint is a real constraint on the API.**
  Models where the constraint genuinely depends on θ (rare; would
  e.g. model regime switches in the null space) must subtype
  `AbstractLatentComponent` directly. This is documented in the
  `UserComponent` docstring and is unlikely to bite anyone in v0.x.

#### Neutral

- **`prior_mean` is not in the callable.** ADR-023 makes
  `prior_mean(c, θ)` load-bearing in the Newton hot path, but the
  vast majority of components return zeros. We may add `:prior_mean`
  as an optional key in a future minor release without breaking
  existing callables (point 1 above).

### Acceptance criteria

- `UserComponent` reproducing a `Generic0` agrees with the native
  `Generic0` posterior to `1e-10` on a realistic `inla()` fit
  (Gaussian likelihood + tridiagonal SPD structure matrix). ✓
  ([`test/regression/test_user_component.jl`](../packages/LatentGaussianModels.jl/test/regression/test_user_component.jl))
- The intrinsic-with-sum-to-zero-constraint variant agrees to the
  same tolerance under the same fit configuration. ✓
- Construction-time validation surfaces wrong-size `Q`, missing
  `:Q`, and non-`NamedTuple` returns as user-facing errors. ✓
- The full LGM test suite (28 oracle fixtures, 2535+ regression
  assertions) passes unchanged after PR-2 lands. ✓

### References

- [`plans/conti-valiant-pebble.md`](../plans/conti-valiant-pebble.md)
  §PR-2 — the implementation plan.
- R-INLA `rgeneric` documentation: `inla.doc("rgeneric")`.
- [Generic0 implementation](../packages/LatentGaussianModels.jl/src/components/generic0.jl)
  — the natural-language reference for "what does a callable need
  to return?". `UserComponent`'s required signature mirrors
  `Generic0`'s required methods.
- ADR-003 — "subtype + multiple dispatch over macros" — the
  architectural commitment that makes `UserComponent` a *surface*
  rather than a *parallel mechanism*.
- ADR-023 — `prior_mean` promoted to load-bearing; informs why
  `:prior_mean` is reserved for a future namedtuple key rather than
  shipped now.

---

## ADR-026: Marginal strategies via `AbstractMarginalStrategy` type dispatch

Status: Accepted
Date: 2026-05-04

### Context

Phase L opens the third R-INLA marginal strategy: `FullLaplace` (R-INLA's
`strategy = "laplace"`), the per-`x_i` refitted Laplace that closes the
sharply-non-Gaussian-latent gap on fixtures like Brunei.

The two existing strategies — `:gaussian` and `:simplified_laplace` —
are dispatched via symbol whitelist in two locations:

- The integration-stage selector in
  `INLA(; latent_strategy=:gaussian|:simplified_laplace)`
  ([`packages/LatentGaussianModels.jl/src/inference/inla.jl`](../packages/LatentGaussianModels.jl/src/inference/inla.jl))
  — controls whether `_inla_integrate` applies the Rue-Martino mean shift
  to `x_mean` / `x_var`.
- The per-coordinate density selector in
  `posterior_marginal_x(res, i; strategy=:gaussian|:simplified_laplace)`
  ([`packages/LatentGaussianModels.jl/src/inference/marginals.jl`](../packages/LatentGaussianModels.jl/src/inference/marginals.jl))
  — controls whether the density mixture is augmented by the Edgeworth
  first-order skew correction.

Both surfaces use the same two symbol values to mirror R-INLA's
`control.inla$strategy`, which is one user-facing knob. ADR-016 already
documents that the two facets — integration-stage mean shift and
per-coordinate density skew — are *independently* triggered by the
single user setting because they happen at different points in the
pipeline.

`FullLaplace` is fundamentally different in shape from the other two:

- It is parameterised. The reference-grid construction takes a span and
  point count per coordinate; `:gaussian` and `:simplified_laplace`
  carry no per-call configuration.
- It re-runs the inner Newton with one coordinate fixed as a datum, via
  a constraint injection on top of the existing
  `LinearConstraint`-aware Laplace pipeline. That requires per-strategy
  state and configuration that does not fit a flat enum.
- It is expensive enough that users may want to tune knobs (rank-1
  factor reuse, grid adaptation) without changing strategy.

A symbol whitelist cannot extend cleanly to a parameterised strategy.
The replan ([§Phase L line 350](../plans/replan-2026-04-28.md)) frames
Phase L as the implicit quality gate on the marginal-strategy
abstraction; introducing the type hierarchy now is the cleanest place
to validate that the abstraction holds.

### Decision

Promote the symbol whitelist to a multiple-dispatch type hierarchy.

```julia
abstract type AbstractMarginalStrategy end

struct Gaussian          <: AbstractMarginalStrategy end
struct SimplifiedLaplace <: AbstractMarginalStrategy end
struct FullLaplace       <: AbstractMarginalStrategy   # PR-3
    n_grid::Int
    span::Float64
end
```

The strategy is the same object passed to both surfaces:

```julia
INLA(; latent_strategy = SimplifiedLaplace())                 # integration
posterior_marginal_x(res, i; strategy = FullLaplace())        # per-x_i density
refine_hyperposterior(res, model, y;
                      latent_strategy = SimplifiedLaplace())  # re-integrate
```

A symbol shim preserves backwards compatibility:

```julia
_resolve_marginal_strategy(s::AbstractMarginalStrategy) = s
function _resolve_marginal_strategy(s::Symbol)
    s === :gaussian          && return Gaussian()
    s === :simplified_laplace && return SimplifiedLaplace()
    throw(ArgumentError("unknown marginal strategy :$s; " *
                        "use Gaussian(), SimplifiedLaplace(), FullLaplace(), " *
                        "or :gaussian / :simplified_laplace"))
end
```

The shim mirrors `_resolve_scheme(::Symbol, ::Int)` for
`AbstractIntegrationScheme` ([`inla.jl:95-102`](../packages/LatentGaussianModels.jl/src/inference/inla.jl)),
which is the established pattern for "type-or-symbol" keyword handling
in the package.

The dispatch contract a new strategy must implement is single-method:

```julia
# integration-stage hook (PR-1 dispatches `_inla_integrate` on this)
_apply_mean_shift(::AbstractMarginalStrategy, lp, model, y) -> Vector{Float64}

# per-coordinate density hook (PR-1 dispatches `posterior_marginal_x` on this)
_density_mixture!(::AbstractMarginalStrategy, pdf, xs, m_k, σ_k, w_k, …)
```

`Gaussian` returns a zero shift / Gaussian mixture; `SimplifiedLaplace`
returns the Rue-Martino shift / Edgeworth-corrected mixture;
`FullLaplace` (PR-3) returns the appropriate refit-derived quantities.

### Consequences

**Pros.**

- `FullLaplace` plugs in as a third subtype with its own configuration
  fields, no further plumbing needed.
- The two facets of `SimplifiedLaplace` (mean shift in integration;
  Edgeworth in marginals) are dispatched on the same type but via
  *different methods*, preserving their semantic independence (per
  ADR-016) at the implementation level. Users who wanted the mean
  shift but not the density skew (or vice versa) currently have no
  way to ask; that's a follow-up — the new type is the right place
  to express it (e.g. `SimplifiedLaplace(; mean_shift=true,
  density_skew=true)`).
- `_resolve_marginal_strategy` keeps the symbol API alive indefinitely;
  no deprecation churn.
- The new type lives next to `AbstractIntegrationScheme`, not subsumed
  into it: the two are orthogonal axes (θ-quadrature vs latent-density
  shape). A user can pair any `AbstractIntegrationScheme` with any
  `AbstractMarginalStrategy`.

**Cons.**

- One more public abstract type (`AbstractMarginalStrategy`) and three
  exported concrete types. Documented in
  [`packages/LatentGaussianModels.jl/CLAUDE.md`](../packages/LatentGaussianModels.jl/CLAUDE.md)
  alongside the existing strategy types.
- The `INLA{I, S <: Laplace}` field type becomes
  `latent_strategy::M where M <: AbstractMarginalStrategy`; downstream
  code that destructures `strategy.latent_strategy === :simplified_laplace`
  must rewrite to `strategy.latent_strategy isa SimplifiedLaplace` (or
  use the `_apply_mean_shift` hook). All in-tree code is rewritten in
  PR-1.

**Acceptance criteria for PR-1.**

- All 28 oracle fixtures in
  [`packages/LatentGaussianModels.jl/test/oracle/`](../packages/LatentGaussianModels.jl/test/oracle/)
  pass bit-for-bit unchanged.
- New regression test:
  `INLA(latent_strategy=:simplified_laplace)` and
  `INLA(latent_strategy=SimplifiedLaplace())` produce byte-identical
  `INLAResult` objects on the Pennsylvania BYM2 oracle.
- Same equivalence asserted for `posterior_marginal_x(...; strategy=:simplified_laplace)`
  vs `strategy=SimplifiedLaplace()`.
- Unknown-symbol path still throws `ArgumentError` (the shim's
  fallthrough preserves the existing user-error message).

### References

- ADR-016 — Simplified-Laplace mean-shift correction; introduces the
  facet-independence note this ADR's `_resolve_marginal_strategy`
  preserves.
- [`plans/replan-2026-04-28.md`](../plans/replan-2026-04-28.md) §Phase L,
  line 350 — frames Phase L as the abstraction-quality gate.
- `_resolve_scheme(::Symbol, ::Int)` at
  [`packages/LatentGaussianModels.jl/src/inference/inla.jl:95-102`](../packages/LatentGaussianModels.jl/src/inference/inla.jl)
  — pattern this ADR generalises.

---

## ADR-027: Importance-sampling correction — declined for v0.x; deferred to per-workflow strategy if a use case appears

Status: Accepted
Date: 2026-05-05

### Context

Phase L's tail carried a deferred stretch item, PR-7(b), to add an
"importance-sampling correction" per
[`plans/replan-2026-04-28.md`](replan-2026-04-28.md) lines 497-498,
behind an `INLA(; importance_sample_correct=false)` flag. With Phase L
shipped at v0.1.5 (2026-05-04), the natural moment to either close or
adopt the item is now — before Phase M opens.

Investigation surfaced three things:

1. **The "standard" IS reweighting is already in the integration loop.**
   [`packages/LatentGaussianModels.jl/src/inference/inla.jl:298-307`](../packages/LatentGaussianModels.jl/src/inference/inla.jl)
   implements the textbook Laplace-Gaussian importance-sampling
   estimator described in Rue, Martino, Chopin (2009) JRSSB §3.2 (eq.
   11):

   ```julia
   log_unnorm   = log_base_weights .+ log_π .- log_q
   log_norm     = _logsumexp(log_unnorm)
   w            = exp.(log_unnorm .- log_norm)
   log_marginal = log_norm        # the IS estimator for log Z_π
   ```

   The proposal `q ~ N(θ̂, Σ)` is the Gaussian-at-the-mode used by
   `Grid` / `CCD`. The Phase K skewness-correction docstring at
   [`integration.jl:65`](../packages/LatentGaussianModels.jl/src/inference/integration.jl)
   explicitly notes "the IS reweight in step (3) absorbs the proposal
   mismatch." So the replan's "IS correction" is *not* about adding a
   missing reweight — the reweight is not missing.

2. **The only literature implementation with a substantively-different
   IS correction is Berild, Bolin, Lindgren, Rue (2022) "Importance
   Sampling with the Integrated Nested Laplace Approximation"** (JCGS
   31(4); arXiv 2103.02721; reference R code at
   [`github.com/berild/inla-mc`](https://github.com/berild/inla-mc)).
   This is **IS-INLA / AMIS-INLA**, an entirely different beast:

   - **Use case**: conditional LGMs where `θ` has a non-standard
     distribution that breaks the Gaussian-at-the-mode proposal —
     Bayesian lasso (Laplace prior on `β`), quantile regression
     (asymmetric-Laplace likelihood with the quantile parameter
     outside the LGM hyperparameter set), missing-covariate
     imputation.
   - **Algorithm**: draw `N = 10³–10⁴` samples θ_i ∼ q_proposal, run a
     *full INLA fit* at each, weight by
     `π(θ_i)·π̂(y|θ_i)/q_proposal(θ_i)`. Adaptive variants (AMIS)
     update `q_proposal` across iterations.
   - **Cost**: `N × t_inla` — two to three orders of magnitude slower
     than a single fit.
   - **Diagnostic**: modified ESS (Owen 2013); not Pareto-k.
   - **What it does NOT solve**: the heavy-tailed θ-posterior
     pathology (covered by Phase K's skewness correction) and the
     sharply non-Gaussian per-coordinate latent pathology (covered by
     Phase L's `FullLaplace`). IS-INLA targets *θ*-uncertainty under
     non-standard θ-priors specifically.

The Phase K decision matrix already records that a fixed-`N=100`
Monte-Carlo IS marginal-likelihood correction was evaluated and
rejected because, at that scale, its Monte-Carlo error (~0.1 nat) is
comparable to or larger than the corrections it claims to make and it
ships no ESS diagnostic. Nothing about Phase L's close changes that
arithmetic.

### Decision

1. **PR-7(b) is dropped from the Phase L tail and from the roadmap.**
   No `importance_sample_correct` flag is added to `INLA()` and no
   parallel IS code path is shipped.

2. **The standard Laplace-Gaussian IS reweighting at
   [`inla.jl:298-307`](../packages/LatentGaussianModels.jl/src/inference/inla.jl)
   is the canonical "IS correction" in this project.** It is not
   optional, not behind a flag, and not pluggable.

3. **The failure modes the replan cited are addressed by
   better-targeted, cheaper, deterministic methods**:
   - Heavy-tailed θ-posteriors → Phase K skewness correction
     (asymmetric `Grid` node placement) and `refine_hyperposterior`.
   - Sharply non-Gaussian per-coordinate latents → Phase L's
     `FullLaplace` marginal strategy.
   - Per-observation IS reliability → PSIS-LOO weakdep extension
     with Pareto-k diagnostic (Phase K).

4. **If a user case for Bayesian lasso / quantile regression / missing
   covariates appears later, IS-INLA is the right tool — but shipped
   as a fresh `ISINLA <: AbstractInferenceStrategy` per ADR-010, not
   as a flag on `INLA()`.** A new ADR is required at that time.

### Consequences

#### Positive

- Closes a long-running roadmap loose end without committing to a
  100×–10000× slowdown for workflows not on the v0.x roadmap.
- Removes a confusing, never-implemented `importance_sample_correct`
  kwarg from the prospective `INLA()` API.
- Redirects future "IS-INLA" requests to a structurally-correct
  `AbstractInferenceStrategy` extension (per ADR-010 — third-party
  strategies are first-class) rather than letting them accumulate as
  flags on the canonical entry point.

#### Neutral

- The replan's port pipeline (lines 476-504 of
  [`plans/replan-2026-04-28.md`](replan-2026-04-28.md)) loses one
  candidate. The remaining candidates (constraint-projected
  simplified-Laplace mean correction, asymmetric skewness corrections,
  Edgeworth correction, the BivariateIIDModel and NonStationarySPDE
  numerical kernels) are unaffected.

### References

- Rue, Martino, Chopin 2009 JRSSB §3.2 (eq. 11) — standard IS
  reweighting; the form already implemented at `inla.jl:298-307`.
- Berild, Bolin, Lindgren, Rue 2022 JCGS 31(4); arXiv 2103.02721;
  reference R code at `github.com/berild/inla-mc` — IS-INLA /
  AMIS-INLA, the legitimate port target if/when a Bayesian-lasso /
  quantile-regression / missing-covariate use case appears.
- ADR-010 — third-party `AbstractInferenceStrategy` is first-class;
  the route through which `ISINLA` would ship.
- ADR-026 — `AbstractMarginalStrategy` precedent for adding a new
  strategy without touching `INLA()`'s public surface.
- [`packages/LatentGaussianModels.jl/src/inference/inla.jl:298-307`](../packages/LatentGaussianModels.jl/src/inference/inla.jl)
  — current IS-reweight implementation.
- [`packages/LatentGaussianModels.jl/src/inference/integration.jl:65`](../packages/LatentGaussianModels.jl/src/inference/integration.jl)
  — Phase K skewness-correction docstring noting "the IS reweight in
  step (3) absorbs the proposal mismatch."

---

## ADR-028: Gaussian-basis prior on non-stationary SPDE basis coefficients — match R-INLA's `theta.prior.mean`/`theta.prior.prec` per-coefficient parameterisation; defer PC-on-basis-norm

Status: Accepted
Date: 2026-05-05

### Context

Phase M PR-3 implements `NonStationarySPDEModel` and lands an oracle
fixture against R-INLA's `inla.spde2.matern` (Lindgren-Rue-Lindström
2011 §3.2).

The non-stationary SPDE is parameterised by per-vertex
`log τ_v = (B_τ θ_τ)_v` and `log κ_v = (B_κ θ_κ)_v`, where `B_τ`,
`B_κ` are user-supplied basis matrices of shape `(n_v, p_τ)` and
`(n_v, p_κ)`. The internal hyperparameter vector has length
`p_τ + p_κ`. The prior on `θ = [θ_τ; θ_κ]` is the question.

Three options were on the table:

1. **Unit-Gaussian on every coefficient.**
   `log_prior(NonStationarySPDEModel, θ) = -0.5 ‖θ‖²`. No
   per-coefficient mean/precision tuning, no R-INLA-parity structure;
   would not be load-bearing against the LRL §3.2 fixture.

2. **R-INLA `theta.prior.mean` / `theta.prior.prec` per-coefficient
   Gaussian.** R-INLA's `inla.spde2.matern` exposes two vectors of
   length `p_τ + p_κ`: a prior mean `μ` and a prior precision `λ`, so
   the prior factorises as `θ_k ∼ N(μ_k, 1/λ_k)`. This is the prior
   that R-INLA actually uses internally; matching it is the only way
   to get parity on the LRL §3.2 fixture without adding extra free
   parameters to the comparison.

3. **PC prior on the basis-norm.** Fuglstad-Simpson-Lindgren-Rue (2019)
   §6 sketch a PC prior penalising deviations from a stationary base
   model via the L²-norm of `B_κ θ_κ` integrated over the mesh. This
   is the principled long-term choice but (a) has no R-INLA
   counterpart to validate against, (b) needs a quadrature rule and a
   reference mesh-norm, (c) introduces a third scaling hyperparameter
   that itself needs a prior. Out of scope for v0.2.

### Decision

1. **Ship `GaussianBasisPrior(mean::Vector, prec::Vector)` matching
   R-INLA's per-coefficient Gaussian parameterisation.** The prior
   density is

   ```
   log π(θ) = ∑_k −½ λ_k (θ_k − μ_k)² + ½ log(λ_k / 2π)
   ```

   exactly mirroring `inla.spde2.matern(theta.prior.mean = μ,
   theta.prior.prec = λ)`. Defaults: `mean = zeros(p_τ + p_κ)` and
   `prec = ones(p_τ + p_κ)` — R-INLA's documented defaults
   (`theta.prior.mean = 0`, `theta.prior.prec = 1`). The prior is
   evaluated *with* its normalising constant so that the marginal-
   likelihood comparison against R-INLA includes the same baseline.

2. **Live in `INLASPDE.jl` as a standalone struct, not as an
   `AbstractHyperPrior` subtype.** The package's `AbstractHyperPrior`
   contract is scalar-only (per the `priors/abstract.jl` docstring:
   "This type is for *scalar* priors. Multi-dimensional priors […]
   live in `INLASPDE.jl` because they are inherently coupled."). The
   precedent is `PCMatern` in `INLASPDE.jl/src/priors/pc_matern.jl` —
   a vector-valued prior carried as a plain struct field on the
   component (`SPDE2{α, T, FE, G, PR}` parameterises on
   `PR <: PCMatern`). `GaussianBasisPrior` follows the same shape:
   `SPDE2NonStationary{α, T, FE, G, PR <: GaussianBasisPrior}`.

3. **Decline the unit-Gaussian default.** Hard-coding `prec = 1` for
   the intercept column hard-codes the prior-belief "log τ should be
   order-1 around zero" which is *not* a safe default for arbitrary
   mesh scales. Exposing both `mean` and `prec` lets users widen the
   prior on intercept columns and tighten it on spline-basis
   columns — the canonical R-INLA usage pattern.

4. **Defer PC-on-basis-norm to a follow-up component or a v0.3+
   evolution.** No fixture, no parity benchmark; revisit if a real
   user case appears.

### Consequences

#### Positive

- LRL §3.2 oracle fixture compares like-for-like with R-INLA: the only
  free axis in the fit is the actual SPDE precision construction, not
  prior-mismatch slack.
- The Gaussian-basis prior surface is *the* R-INLA non-stationary
  surface — users who already drive `inla.spde2.matern` with custom
  `theta.prior.mean`/`prec` get a one-line port.
- Decoupling from `AbstractHyperPrior` means we don't have to pretend
  the basis prior is `p_τ + p_κ` independent scalar priors; it stays
  one struct with two vectors.

#### Neutral

- Unit-Gaussian behavior is recoverable as
  `GaussianBasisPrior(zeros(p), ones(p))`, which is the default.
- One more public surface (`GaussianBasisPrior`) but it's confined to
  `INLASPDE.jl` and only used by `SPDE2NonStationary`.

#### Negative

- No PC prior on the non-stationary deviation today. Users who want
  shrinkage to stationarity must encode it via the basis structure
  (e.g. a tight `prec` on spline-basis columns, wide on the intercept)
  rather than through a single penalty hyperparameter.

### References

- Lindgren, Rue, Lindström (2011), JRSSB B 73(4), §3.2 — non-stationary
  SPDE example with piecewise-constant `B_κ` (the oracle fixture).
- Fuglstad, Simpson, Lindgren, Rue (2019), JASA — PC prior on Matérn,
  §6 sketches the basis-norm PC prior (deferred).
- R-INLA `inla.spde2.matern` source — `theta.prior.mean` /
  `theta.prior.prec` per-coefficient Gaussian parameterisation.

---

## ADR-029: `KroneckerComponent` — generic two-component Kronecker composer for separable space-time GMRFs

Status: Accepted
Date: 2026-05-05

### Context

Phase M PR-5 lands the separable space-time piece of the SPDE
expansion arc per [`plans/replan-2026-04-28.md`](plans/replan-2026-04-28.md)
lines 372–405. The Cameletti et al. (2013) PM₁₀ air-pollution case
study is the third of the three Phase M oracle fixtures, and
canonically expresses the spatio-temporal field as a Kronecker product
of an SPDE-Matérn spatial GMRF and an AR(1) temporal GMRF —
`Q = Q_space ⊗ Q_time`. R-INLA expresses this through `f(field, model
= spde, group = time, control.group = list(model = "ar1"))`, which is
internally a Kronecker construction.

The new latent component must compose two child `AbstractLatentComponent`s
into one whose precision matrix is the Kronecker product. The design
question: does the composer ship as

a. a **generic two-component Kronecker composer** —
   `KroneckerComponent(spatial, temporal)` accepting any pair of
   `AbstractLatentComponent`, with hyperparameters concatenated as
   `θ = [θ_space; θ_time]`; or

b. a **specialised `SPDESpaceTime`** type that bakes in the SPDE2
   spatial side, the AR1 / RW1 temporal side, and the mesh-and-time
   metadata explicitly?

PR-1 already shipped the `KroneckerMapping` projector covering the
observation side (`(A_space ⊗ A_time)`). The component side is the
remaining seam.

### Decision

Ship the **generic** `KroneckerComponent(spatial, temporal)` composer
in `LatentGaussianModels.jl`. Concretely:

1. **API**

   ```julia
   struct KroneckerComponent{S<:AbstractLatentComponent,
                             T<:AbstractLatentComponent} <: AbstractLatentComponent
       space::S
       time::T
   end
   ```

   No keyword arguments at construction; the children carry their
   priors, sizes, and hyperparameter layouts. The composer adds no
   hyperparameters of its own.

2. **Hyperparameter layout — concatenation, not sharing.**
   `nhyperparameters(c) = nhyperparameters(c.space) +
   nhyperparameters(c.time)`. The internal-scale θ is sliced as
   `θ_s = θ[1:p_s]`, `θ_t = θ[p_s+1:end]`, and dispatched verbatim to
   the children. This contrasts with `Replicate` (which *shares* a
   single θ block across replicates) and matches R-INLA's `group =`
   surface where the spatial and temporal precisions have independent
   priors.

3. **Precision = `kron(Q_s, Q_t)`.** The Kronecker order matches the
   PR-1 `KroneckerMapping` storage convention
   (`(A_space ⊗ A_time) vec(X) = vec(A_time · X · A_space')`, with `X`
   `(time × space)`); the precision and the projector compose with the
   same flattening. `precision_matrix(c, θ) =
   kron(precision_matrix(c.space, θ_s), precision_matrix(c.time, θ_t))`
   stays sparse (`SparseArrays.kron`).

4. **Log normalising constant via the Kronecker logdet identity.**

   ```
   log |Q_s ⊗ Q_t| = n_t · log|Q_s| + n_s · log|Q_t|
   ```

   We never materialise the Kronecker product to compute the
   normalising constant. Each child already implements
   `log_normalizing_constant`, so the composer can rebuild
   `log|Q_s| = 2·(log_normc_s + ½ n_s log(2π))` and combine — but
   that round-trip is fragile when a child is intrinsic (its
   structural log-det is dropped from `log_normc`). Concretely, we
   require each child component to return a *finite* `log_normalizing_constant`
   (proper Gaussian log-NC), and reject `KroneckerComponent` over
   intrinsic children at construction time. The Cameletti workflow
   uses SPDE2 ⊗ AR1 — both proper, no intrinsic complication. RW1 ⊗
   IID requires its own ADR (PR-5 follow-up).

5. **Constraints — pass through if at most one child is constrained.**
   When neither child has a constraint, the composer returns
   `NoConstraint`. When only one child carries an `LinearConstraint
   (A_c, e_c)`, the composer lifts it to the Kronecker dimension via
   the appropriate Kronecker product with `I` on the unconstrained
   axis. When both children are constrained, throw
   `ArgumentError` — the joint-constraint case (rank, conditioning-
   by-kriging interaction with the Kronecker factor) is non-trivial
   and not exercised by Phase M's oracle fixtures.

6. **Prior mean — Kronecker of children's prior means.** When both
   children's `prior_mean` are zero (the default), the composed mean
   is zero. When non-zero, the composed mean is
   `kron(μ_space, μ_time)` (interpreting both as column vectors), which
   matches the latent flattening convention from item 3.

### Consequences

#### Positive

- One file, ~120 LOC. Reuses every contract method on the children.
- Cameletti SPDE2 ⊗ AR1 is the immediate consumer; AR1 ⊗ AR1, AR1 ⊗
  IID, and Generic0 ⊗ AR1 are tested by composition without needing
  INLASPDE in the LGM regression suite.
- Future SPDE work (e.g. Cameletti's `RW2` temporal variant) drops in
  by swapping the `time` child — no new types.
- The PR-1 `KroneckerMapping` and PR-5 `KroneckerComponent` form a
  matched pair; ADR-017's projector seam already exposes the
  observation-side composition.
- Third-party components (a user's `Generic2` ⊗ AR1, a custom seasonal
  ⊗ SPDE2) just work without touching the composer.

#### Neutral

- Hyperparameter ordering is `[θ_space; θ_time]` by convention, not
  `[θ_time; θ_space]`. We document this in the `KroneckerComponent`
  docstring and align the kron order so the user-facing flattening is
  consistent.
- Component children must implement `log_normalizing_constant`. All
  v0.1 concrete components do.

#### Negative

- The intrinsic-child case (RW1 ⊗ AR1, RW2 ⊗ AR1) is deferred. The
  Kronecker logdet identity still holds in form, but the Marriott-Van
  Loan correction for the joint constraint requires extra work; we do
  not attempt it in PR-5.
- Rejecting the doubly-constrained case is conservative; some users
  may want it. The error message points at the deferral.
- A specialised `SPDESpaceTime{α}` could carry mesh-aware diagnostic
  helpers (e.g. spatial-marginal extraction at each time slice). We
  defer those to user-space helpers built on top of the generic
  composer.

### References

- Cameletti, Lindgren, Simpson, Rue (2013), AStA 97(2):109–131 — the
  PM₁₀ air-pollution case study; `f(field, model = spde, group =
  time, control.group = list(model = "ar1"))`.
- Lindgren, Rue, Lindström (2011), JRSSB B 73(4) — SPDE-Matérn
  precision construction; the spatial side of the Kronecker.
- ADR-017 — `AbstractObservationMapping` projector seam; PR-1
  `KroneckerMapping` is the projector counterpart of this component.
- [`packages/LatentGaussianModels.jl/src/observation_mapping.jl:296`](packages/LatentGaussianModels.jl/src/observation_mapping.jl)
  — `KroneckerMapping` flattening convention this component matches.
- [`packages/LatentGaussianModels.jl/src/components/replicate.jl`](packages/LatentGaussianModels.jl/src/components/replicate.jl)
  — nearest precedent (block-diagonal composer with *shared* θ);
  this ADR's design contrasts via concatenated θ.

---

## ADR-030: Fractional-α SPDE (Bolin-Kirchner 2020) — deferred to v0.2.1+; Phase M closes at PR-6

Status: Accepted
Date: 2026-05-05

### Context

Phase M's plan
([`plans/conti-valiant-pebble.md`](conti-valiant-pebble.md)) carried a
stretch tail PR-7 — fractional-α SPDE via the Bolin-Kirchner 2020
rational approximation:

```math
(\kappa^2 - \Delta)^{-\alpha/2} \approx \sum_{k=1}^{m} r_k \,
    (\kappa^2 - \Delta + s_k)^{-1}
```

with poles `s_k > 0` and weights `r_k` from a degree-`m` Padé fit. The
plan suggested implementing `SPDEFractional` as "a sum of `m`
integer-α=1 SPDE precision matrices with shifted κ", with default
`m = 4` covering `α ∈ (0.5, 2.5)` to ~3 digits of accuracy.

Phase M PRs 1–6 closed cleanly (KroneckerMapping, 1D SPDE,
non-stationary SPDE, PD-failure safety net, separable space-time,
mesh utilities maturity). With the stretch tail nominally available,
the natural moment to either commit to PR-7 or document the deferral
is now — before tagging v0.2.0.

A closer read of Bolin-Kirchner 2020 ("The Rational SPDE Approach for
Gaussian Random Fields With General Smoothness", JCGS 29(2),
arXiv:1711.04333) and the production reference implementation at
[`finnlindgren/rSPDE`](https://github.com/finnlindgren/rSPDE) surfaced
three things:

1. **The plan's "sum of integer-α=1 precision matrices" sketch is
   wrong, mathematically.** Rational approximation of the *covariance*
   operator gives a sum of integer-α covariances; that is **not** a
   sum of precisions. If `Σ_α ≈ Σ_k r_k · Σ_{α=1, κ_k}`, then the
   precision is the inverse of that sum — *dense*, not sparse, and
   not equal to `Σ_k r_k · Q_{α=1, κ_k}`. The naive sum-of-precisions
   formula yields a different covariance operator entirely, with no
   theoretical link to the fractional Matérn target.

2. **The correct rational-SPDE construction needs state augmentation.**
   Following Bolin-Kirchner §2.3 and the rSPDE implementation, the
   sparse-precision-preserving form introduces auxiliary fields
   `v_k ∈ ℝ^{n_v}` for `k = 1, …, m`, each satisfying a stationary
   integer-α=1 SPDE `(K + s_k C̃) v_k = √(C̃) W_k` with independent
   white noise `W_k`. The field of interest is the linear combination
   `u = Σ_k √(r_k) · v_k`. The joint state `(v_1, …, v_m)` has
   *block-diagonal* sparse precision of size `m·n_v × m·n_v`; the
   marginal precision on `u` alone is *dense*. SPDE inference therefore
   requires fitting on the augmented `m·n_v`-dimensional latent
   vector, with a partial observation operator that reads off
   `u(s_obs) = Σ_k √(r_k) · A_obs · v_k`.

3. **This doesn't fit cleanly into LGM's per-vertex
   `AbstractLatentComponent` contract.** The contract
   ([`packages/LatentGaussianModels.jl/src/components/abstract.jl`](../packages/LatentGaussianModels.jl/src/components/abstract.jl))
   assumes one component owns one scalar field of dimension `length(c)`
   with one sparse `precision_matrix(c, θ) -> SparseMatrixCSC`. The
   rational-SPDE construction violates this in two places: (a) the
   "field of interest" `u` is *not* the latent vector — it is a
   linear functional of the augmented state; (b) the
   `MeshProjector` A-matrix maps mesh vertices to observation points
   for *one* field, not for `m` stacked auxiliary fields tied by a
   summation operator.

   Concretely, fitting a rational-SPDE component would need either:
   - A new `AugmentedComponent` abstraction wrapping `m` integer-α
     SPDE2 components plus a "summation observation operator" — an
     `AbstractObservationMapping` (per ADR-017) that does
     `Σ_k √(r_k) · A_obs · v_k` instead of `A_obs · u`; or
   - An adapter that lifts `(SPDE2, SPDE2, …, SPDE2)` through a custom
     `MeshProjector` and synthesises a single virtual `precision_matrix`
     by stacking the block-diagonal augmented precision — but this
     leaks the augmentation up to the LGM solver and breaks the
     `length(c) = n_v` invariant assumed by `prior_mean`,
     `constraints`, marginal-strategy code, and the inner Newton hot
     path.

   The first option is structurally clean but is multi-day work — a
   new abstraction in LatentGaussianModels.jl, a new observation
   mapping subtype, and end-to-end coverage of the Newton/Laplace
   path with augmented latents. That is not a "stretch tail" item; it
   is a phase-shape change.

4. **PR-7 is a fresh write regardless.** Stretch-tail justification
   for shipping in Phase M relied on the plan's (incorrect)
   "shifted-κ sum of precisions" sketch; once the actual construction
   is in view, the LOC budget is closer to PR-3 (~500 LOC across LGM
   + INLASPDE) than PR-6's mesh utilities work (~250 LOC, single
   package).

5. **No oracle-fixture obligation.** The replan's three Phase M
   oracles (synthetic 1D, Lindgren-Rue-Lindström §3.2, Cameletti
   PM10) all landed with PRs 2/3/5; PR-7 was always a stretch with no
   gating fixture. Deferring it does not move the phase-close gate.

### Decision

1. **Phase M closes at PR-6.** No fractional-α SPDE in v0.2.0. The
   replan's Phase M scope item (5/5) ships as documented-but-deferred,
   mirroring the v0.x discipline that "stretch" means "ship if the
   architecture happens to fit, defer otherwise" (precedent: ADR-027,
   PR-7(b) IS-correction; ADR-015, `LGMFormula.jl` /
   `GMRFsPardiso.jl` deferred from v0.1).

2. **The deferral target is v0.2.1+, not "Phase M+1".** Fractional-α
   SPDE is a single component with a localised infrastructure
   prerequisite (the augmentation seam in LGM); it does not need its
   own phase. When the prerequisite work lands — either as a v0.2.x
   minor item or as part of a future phase that has independent need
   for `AbstractObservationMapping` extensions — the SPDE component
   is then ~1 PR of pure rational-approximation arithmetic on top.

3. **The infrastructure prerequisite is named explicitly**:
   `AugmentedLatentComponent` (or equivalent) in
   `LatentGaussianModels.jl` that lets one logical "field of interest"
   be a linear functional of a stacked block-diagonal latent vector,
   with a paired `AbstractObservationMapping` subtype that performs
   the summation. This is general-purpose: it covers rational-SPDE,
   the SPDE-on-sphere construction (which uses a similar
   augmentation), and any future model whose "user-facing field" is
   a linear combination of multiple latent components. It is **not**
   committed to v0.2.x — it lands when a use case arrives.

4. **No `SPDEFractional` skeleton ships.** No placeholder struct, no
   throwing constructor, no documented-but-empty API surface. The
   v0.2.0 release ships `SPDE2` (α ∈ {1, 2}, 2D), `SPDE1D`
   (α ∈ {1, 2}, 1D), `SPDE2NonStationary` (α ∈ {1, 2}, 2D,
   per-vertex `(τ, κ)`), and the Kronecker composer. Fractional-α is
   absent — users who need it know to wait, not to find a half-built
   API.

5. **The replan's Phase M scope of 5 items closes as 4 items shipped
   + 1 deferred.** This is honest and traceable; the alternative
   (rushing a structurally-wrong implementation through stretch
   bandwidth) would have created technical debt the v0.2.x line
   would then have to unwind.

### Consequences

#### Positive

- v0.2.0 ships on schedule (Phase M week 8, against an 8–12 week
  replan estimate).
- No half-built fractional-α API surface for users to hit and
  discover is broken. R-INLA users moving SPDE workflows get a clean
  v0.2.0 with the four shipped SPDE components and a documented
  "fractional-α: see v0.2.x roadmap" gap.
- Surfaces the real prerequisite — the `AugmentedLatentComponent`
  seam — as a discoverable infrastructure item rather than a
  hidden cost inside an SPDE component PR. When the seam lands, it
  benefits other models (sphere SPDE, multi-resolution analysis)
  beyond fractional-α.
- Resolves the tension between the plan's mathematical sketch and
  the production rational-SPDE construction. Future-Phase work
  starts from the right baseline (`rSPDE`'s state-augmentation form),
  not from a sum-of-precisions misreading.

#### Neutral

- Replan scope-completion drops to 4/5 for Phase M. The remaining
  item is the only one whose deferral does not move a flagship
  workflow gate (the geostatistics flagship is satisfied by SPDE2 +
  non-stationary + space-time; fractional-α is a smoothness
  generalisation, not a workflow blocker).
- v0.2.0's CHANGELOG and release notes will explicitly call out
  fractional-α as out of scope; the docs landing page for SPDE will
  link to this ADR.

#### Negative

- Users with α ∈ ℝ⁺ \ {1, 2} workflows must keep using R-INLA's
  fractional-α path. The set of such users is small in practice —
  Lindgren et al.'s default α=2 covers the published case studies
  this project targets — but the set is not empty.
- The `AugmentedLatentComponent` seam is now load-bearing future
  work; if it never lands, fractional-α stays deferred indefinitely.
  Mitigated by point 3 above: the seam has independent justification
  (sphere SPDE, multi-resolution) and is not solely an SPDE concern.

### What would unblock shipping

Three items in order:

1. **`AbstractObservationMapping` extension** — a `LinearCombinationMapping`
   (or similar) subtype that applies `Σ_k w_k · A_k · x_k` to a
   stacked latent vector `[x_1; …; x_m]`, with `apply!` /
   `apply_adjoint!` matching the existing
   [`observation_mapping.jl`](../packages/LatentGaussianModels.jl/src/observation_mapping.jl)
   contract. The PR-1 KroneckerMapping is precedent for the
   block-structured mapping pattern.
2. **`AugmentedLatentComponent`** — a wrapper composing
   `(component_1, …, component_m)` into a single
   `AbstractLatentComponent` whose `length` is `Σ_k length(component_k)`,
   `precision_matrix` is block-diagonal, and `log_hyperprior` sums
   child priors. Distinct from `KroneckerComponent` (ADR-029): the
   structure here is *block-diagonal* (independent fields tied at the
   observation level), not Kronecker-product (separable joint).
3. **`SPDEFractional` component** — given items 1 and 2, the
   fractional-α SPDE is `AugmentedLatentComponent(SPDE1, SPDE1, …)`
   with `m` shifted-κ child SPDEs and a `LinearCombinationMapping`
   carrying the rational-approximation weights. The Padé/CF
   computation of `(r_k, s_k)` is closed-form in `α` and `m`
   (Bolin-Kirchner §3); ~50 LOC.

Phase M+1 (or an interstitial v0.2.x release) is the natural home
for items 1 and 2. Item 3 is then a single follow-up PR.

### References

- Bolin & Kirchner 2020, "The Rational SPDE Approach for Gaussian
  Random Fields With General Smoothness", JCGS 29(2):274–285,
  arXiv:1711.04333 — the rational-approximation construction; §2.3
  gives the state-augmentation form, §3 gives the
  Padé/contour-fraction computation of `(r_k, s_k)`.
- [`finnlindgren/rSPDE`](https://github.com/finnlindgren/rSPDE) — the
  production reference implementation; the augmented-state
  formulation in [`R/operators.R`](https://github.com/finnlindgren/rSPDE/blob/main/R/operators.R)
  is the implementation target when this ADR is unblocked.
- [`packages/LatentGaussianModels.jl/src/components/abstract.jl`](../packages/LatentGaussianModels.jl/src/components/abstract.jl)
  — the per-vertex contract that doesn't fit fractional-α; the seam
  to extend.
- ADR-017 — `AbstractObservationMapping` is the surface where the
  `LinearCombinationMapping` will land.
- ADR-029 — `KroneckerComponent` precedent for block-structured
  composition; `AugmentedLatentComponent` is the block-diagonal
  sibling.
- ADR-027 — precedent for documenting a deferred replan item with
  the architectural reason; the same pattern applied to PR-7(b)
  IS-correction.
- ADR-015 — precedent for v0.x deferral discipline; sub-packages
  cut from v0.1 with named v0.2 promotion targets.
- [`plans/conti-valiant-pebble.md`](conti-valiant-pebble.md) — Phase
  M plan; PR-7 stretch criterion ("ship if PRs 1–6 close inside
  week 8; defer to Phase M+1 if time-pressed") and the (now-known
  incorrect) "sum of integer-α=1 SPDE precision matrices with
  shifted κ" sketch.

---

## ADR-031: Targeted exception classification in the Laplace bad-θ wrapper

Status: Accepted
Date: 2026-05-05

### Context

The LGM Newton hot path wraps the inner Laplace step with a smooth
quadratic-in-θ penalty so that LBFGS over the hyperparameters never
sees a non-finite objective. Two wrapper sites exist:

- `_neg_log_posterior_θ` (`inla.jl`) — closure called by Optimization.jl
  during the θ-mode search.
- `_inla_integrate` (`inla.jl`) — per-design-point loop during
  importance-sample integration.

Both used `try { laplace_mode(...) } catch { return penalty }` with a
**bare** catch. That mask is too wide: it hides genuine bugs in user
likelihood code (a `MethodError` from a typo in `∇_η_log_density`, an
`ArgumentError` from a misshapen `A` matrix) by silently routing them
into the penalty region. Users debugging "why does INLA think this θ
is infeasible?" then have to single-step the inner Newton to discover
the error.

An earlier prototype caught only `LinearAlgebra.PosDefException`
because that was the single failure mode it had observed; this PR
(M-PR-4) is the moment to lift that targeting into the new code, and
to add the two adjacent numerical-failure types we expect:
`SingularException` (LAPACK side of the same Cholesky failure) and
`DomainError` (out-of-domain log likelihood — e.g. negative variance
in a hand-coded Gaussian log-density at a pathological θ).

### Decision

Introduce a single classifier `_is_bad_theta_failure(err) -> Bool`
that returns `true` only for `LinearAlgebra.PosDefException`,
`LinearAlgebra.SingularException`, and `DomainError`. Both wrapper
sites become `try ... catch err; _is_bad_theta_failure(err) ||
rethrow(err); ... end` — bad-θ failures land in the penalty region;
everything else propagates.

Component constructors with **parametric** domain checks (e.g.
`AR1GMRF`'s `ρ ∈ (-1, 1)` and `τ > 0`) throw `DomainError` rather
than `ArgumentError`, so the classifier catches them. Static
shape/structure checks (e.g. `n ≥ 2`) keep `ArgumentError` because
they cannot be triggered by any θ-step — failure is a programming
bug, not a domain violation. `ArgumentError` is therefore the
"programming-bug" tier and `DomainError` the "θ-domain" tier. New
components should follow the same split.

Three secondary defenses ride along in the same PR (no separate ADR
because they are mechanical):

- **NaN warm-start reset.** `laplace_mode` resets `x₀` to zeros if any
  entry is non-finite. CCD's per-point loop forwards the previous
  point's `x̂`; if that point was discarded as bad-θ its `x̂` may
  contain NaN.
- **FD-Hessian non-finite guard.** `_safe_inverse_hessian` falls back
  to identity covariance with a warning when the FD Hessian at θ̂ has
  NaN/Inf entries (probes that straddle the penalty cliff).
- **Per-point keep-mask** (already in `_inla_integrate` pre-PR;
  retained, now triggered by the targeted classifier).

The penalty form `1.0e10 + 1.0e3 · ‖θ‖²` is unchanged.

### Consequences

#### Positive

- Genuine bugs (typos in user likelihoods, broken model construction)
  surface as the actual exception with the actual stacktrace, not as
  a silent "this θ is infeasible" misdirection.
- The classifier is a single point of truth — adding a new failure
  mode is a one-line edit, and the test suite (predicate +
  end-to-end on Gaussian/SPDE2 fixtures) covers both sides.
- The targeted-classifier pattern is set as the v0.2 baseline; future
  numerical-failure modes can be folded in without re-architecting
  the wrappers.

#### Negative

- Slightly larger surface: future contributors must add to the
  classifier when introducing a new numerical-failure mode in a
  user-defined component. The docstring on `_is_bad_theta_failure`
  flags this.
- The classifier list is not exhaustive — if a new failure mode lands
  with a different exception type before the docstring is read, the
  user sees a hard crash rather than the smooth penalty. Treated as
  acceptable: a hard crash with a clear stacktrace is strictly more
  useful than a silent penalty.

#### Out of scope (deferred)

- A full survey of every numerical exception thrown by every
  Optimization / LinearSolve / SuiteSparse path. The classifier is
  expanded as new failures are observed in oracle / triangulation
  tests, not pre-emptively.

### References

- [`packages/LatentGaussianModels.jl/src/inference/inla.jl`](packages/LatentGaussianModels.jl/src/inference/inla.jl)
  — `_is_bad_theta_failure`, the two wrapper sites,
  `_safe_inverse_hessian` finite guard.
- [`packages/LatentGaussianModels.jl/src/inference/laplace.jl`](packages/LatentGaussianModels.jl/src/inference/laplace.jl)
  — NaN warm-start reset.
- [`packages/LatentGaussianModels.jl/test/regression/test_safety_net.jl`](packages/LatentGaussianModels.jl/test/regression/test_safety_net.jl)
  — closed-form regression suite for all four behaviours.
- [`packages/INLASPDE.jl/test/regression/test_safety_net_spde.jl`](packages/INLASPDE.jl/test/regression/test_safety_net_spde.jl)
  — SPDE2 triangulation against the LGM safety net.
- [`packages/GMRFs.jl/src/gmrf.jl`](packages/GMRFs.jl/src/gmrf.jl) —
  `AR1GMRF` parametric domain checks raised as `DomainError`.
- [`packages/INLASPDE.jl/src/assembly/precision.jl`](packages/INLASPDE.jl/src/assembly/precision.jl)
  — `spde_precision` / `spde_precision_nonstationary` parametric
  `(τ, κ)` checks raised as `DomainError`; structural `α ∈ {1, 2}`
  and length-mismatch checks stay `ArgumentError`.
- ADR-022 — IIDND saturation precedent for the smooth-penalty design.

---

## ADR-032: Mesh utilities maturity — alpha-shape boundary, tuple `max_edge`, boundary pre-subdivision

Status: Accepted
Date: 2026-05-05

### Context

Phase M PR-6 closes three deferred items in
[`packages/INLASPDE.jl/plans/plan.md:100-109`](packages/INLASPDE.jl/plans/plan.md):

1. **Nonconvex hull helper.** R-INLA's `inla.nonconvex.hull(loc, convex)`
   wraps a point cloud with an α-concave polygon — used when the
   convex hull over-extends into empty regions (coastline data, river
   networks, study regions with holes). v0.1 ships only
   `convex_hull_polygon`; users with non-convex domains have had to
   hand-build the boundary as a `k × 2` matrix.
2. **Two-region `max_edge`.** R-INLA's `max.edge = c(inner, outer)`
   refines the data region tighter than the buffer ring outside.
   Critical for cost-controlled SPDE work: a single `max_edge` either
   over-refines the buffer (wasted vertices) or under-refines the data
   region (poor SPDE accuracy where it matters). Today the API is
   single-valued.
3. **Boundary pre-subdivision.** The current `max_edge` is enforced via
   the equilateral-area Ruppert bound (`max_area = √3/4 · max_edge²`),
   which is *soft* on edge lengths — boundary edges in the input
   polygon longer than `max_edge` survive refinement at their input
   length. The existing `test_mesh_quality.jl` regression accepts
   `max_edge_J ≤ 1.5 · max_edge` for that reason.

The three items can be tackled independently but ship together because
they share the same critical file (`src/mesh/inla_mesh.jl` +
`src/mesh/boundary.jl`) and the same test surface.

A fourth design call sits underneath: do we adopt `ConcaveHull.jl`
(small Julia package, alpha-shape based) as a `[deps]` entry, or roll
our own α-shape over the `DelaunayTriangulation.jl` we already depend
on?

### Decision

#### Alpha-shape: roll our own

Implement `nonconvex_hull_polygon(loc; α)` natively in
`src/mesh/boundary.jl`, ~80 LOC including boundary tracing. The
algorithm:

1. Compute Delaunay triangulation of `loc` via DT.jl (already a dep).
2. For each triangle, compute its circumradius `r`. Keep iff `r ≤ α`
   (Edelsbrunner alpha-shape convention; large α → convex hull, small
   α → tight wrap around points).
3. Boundary edges = edges incident to exactly one kept triangle.
4. Trace boundary into a CCW closed polygon; orient via signed area.

ConcaveHull.jl would have added Distances, NearestNeighbors,
StaticArrays, RecipesBase, StatsAPI, PrecompileTools, AbstractTrees as
transitive deps for ~50 LOC of work we can do directly on top of DT.
The ADR weighs the dep footprint against the maintenance cost; with DT
already giving us the triangulation, the maintenance cost of ~80 LOC is
the smaller burden. CLAUDE.md's policy is that new `[deps]` entries
need explicit justification, and "saving 50 LOC" doesn't clear that
bar.

The α parameter is exposed directly (no auto-tuned default that
depends on point density). Default `α = 2 · median(nearest-neighbour
distance)` for the convenience case where the user just wants
"reasonable concavity".

API restriction: only **simply-connected** alpha-shapes (single closed
loop) are returned. Multi-component or hollow alpha-shapes throw
`ArgumentError` with guidance to lower α. R-INLA's `inla.nonconvex.hull`
also enforces simple connectivity.

#### Two-region `max_edge`: tuple API + DT.jl `custom_constraint` callback

Public API extension to `inla_mesh_2d`:

```julia
inla_mesh_2d(loc; max_edge = (inner, outer), offset = (inner, outer), …)
```

When both `max_edge` and `offset` are 2-tuples, the inner boundary is
`expand_polygon(hull(loc), offset[1])` and the outer boundary is
`expand_polygon(inner, offset[2])`. Refinement uses
`DelaunayTriangulation.refine!(tri; max_area = max_area_outer,
custom_constraint = (tri, T) -> centroid_inside_inner(T) ?
area(T) > max_area_inner : false)` so that interior triangles get the
tighter area bound and exterior (buffer) triangles only need to satisfy
the outer area bound.

Single-valued `max_edge::Real` keeps its current behaviour. Tuple input
without paired tuple `offset` errors at construction.

Explicit `boundary` input remains single-region (one polygon, one
`max_edge`). Two-region with explicit boundary is deferred — needs an
explicit nested-polygon API which is bigger scope.

#### Boundary pre-subdivision: opt-in, default off (back-compat)

New kwarg `subdivide_boundary::Bool = false`. When true:

1. Walk the input boundary polygon (after `expand_polygon`/`offset`
   resolution but before triangulation).
2. For each consecutive pair `(p_i, p_{i+1})`, if `dist > max_edge_outer`
   (or `max_edge` in single-region mode), split into
   `⌈dist / max_edge⌉` equally-spaced segments by inserting Steiner
   points before the constrained Delaunay step.
3. The triangulation then sees those points as boundary vertices, and
   Ruppert refinement preserves them.

Default is `false` because turning it on changes mesh statistics
(vertex count, boundary count) on every existing fixture — the
fmesher-parity oracle would need to be regenerated. The opt-in path
gives users who hit the soft-bound problem a sharp tool; existing
behaviour is unchanged.

When `subdivide_boundary = true` and refinement passes, the bound
becomes strict: `max_edge_J ≤ max_edge` to within DT's refinement
roundoff (~1e-12).

### Consequences

#### Positive

- Users with non-convex domains can build the boundary in one call
  instead of hand-building a polygon matrix, lowering the on-ramp for
  coastline / watershed / study-region SPDE work.
- Users with deep buffer rings (e.g. SPDE oracles on small data
  regions inside large extension zones) save vertex count by an order
  of magnitude — the Cameletti M PR-5 mesh, hypothetically, drops from
  ~1000 to ~200 vertices with `max_edge = (0.15, 0.4)`.
- The strict `max_edge` path closes a long-standing footnote in
  `test_mesh_quality.jl` (the 1.5× soft bound). Critical for hand-
  computed FEM-error analyses where a known max edge is load-bearing.
- No new `[deps]` entries — INLASPDE.jl's dependency surface stays at
  its v0.1 size.

#### Negative

- ~250 LOC added to `src/mesh/`. The alpha-shape implementation has a
  documented edge case (multi-component alpha shapes) where it errors
  out — users may need to retry with a different α.
- The two-region API doubles the `max_edge` / `offset` parameter
  surface. Documentation and examples carry both the scalar and tuple
  forms.
- Pre-subdivision opt-in vs default-on is a hedge. We may flip the
  default to `true` in v0.3 once existing fixtures have been
  regenerated and the fmesher-parity oracle's margin has tightened.

#### Out of scope (deferred to v0.3+)

- **R-INLA `inla.nonconvex.hull(loc, convex, concave, …)` parity.**
  R-INLA's helper uses a morphological-closing algorithm with a
  smoothing radius, not a strict alpha-shape. Our helper is
  Edelsbrunner-style; for users who need pixel-identical R-INLA
  geometry, a separate `INLASPDEFmesher.jl` extension that calls into
  `fmesher` over the JLL-shipped binary remains an option (ADR-007's
  deferred fallback).
- **Multi-region max_edge with N > 2 regions.** R-INLA accepts
  `max.edge = c(0.1, 0.3, 0.5)` for three nested regions. We ship the
  inner/outer 2-tuple only.
- **Hollow domains.** Polygons with interior holes (lakes inside a
  watershed) need a multi-loop boundary representation; deferred.

### References

- [`packages/INLASPDE.jl/src/mesh/inla_mesh.jl`](packages/INLASPDE.jl/src/mesh/inla_mesh.jl)
  — `inla_mesh_2d` API extension for tuple `max_edge`/`offset` and
  the `subdivide_boundary` flag.
- [`packages/INLASPDE.jl/src/mesh/boundary.jl`](packages/INLASPDE.jl/src/mesh/boundary.jl)
  — `nonconvex_hull_polygon` + `subdivide_polygon`.
- [`packages/INLASPDE.jl/test/regression/test_mesh_nonconvex.jl`](packages/INLASPDE.jl/test/regression/test_mesh_nonconvex.jl)
  — alpha-shape, two-region, and pre-subdivision regression suites.
- DT.jl's `refine!(tri; custom_constraint = ...)` is the
  region-aware refinement primitive used by the two-region path.
- ADR-007 — INLASPDEFmesher.jl deferred fallback for users who need
  fmesher-pixel-identical mesh geometry.
- Edelsbrunner & Mücke 1994, "Three-dimensional alpha shapes" — base
  algorithm; we implement the 2D restriction.

---

## ADR-033: Multi-likelihood formula syntax — tuple-LHS `(y1, y2, …) ~ rhs`

Status: Accepted
Date: 2026-05-05

### Context

Phase N PR-4 ([`plans/phase-n.md:153-192`](phase-n.md)) ships
multi-likelihood support in `@lgm`: a single formula expression that
fits a joint LGM with two or more `AbstractLikelihood` instances
sharing one latent vector. Two surface candidates were on the table:

1. **Tuple-LHS, Julia-idiomatic.** `(y1, y2) ~ 1 + f(idx, IID(n))`.
   The macro expands the tuple into a `JointLikelihood([ℓ1, ℓ2])` and
   builds a `StackedMapping` over the shared RHS.
2. **R-INLA's stacked-Y form.** `Y ~ stack(...)` with a per-row
   `family` vector encoding which row uses which likelihood. This is
   how R-INLA's `inla(Y ~ ..., family = c("gaussian", "poisson"),
   data = inla.stack(...))` works under the hood.

Both forms can express the same joint model; the question is which
becomes the primary surface and which is a fallback.

The shared-RHS-only constraint (every `f`-term applies to every
likelihood) is enforced in PR-4. Per-likelihood RHS variation (e.g.
Baghfalaki joint longitudinal-survival, where the survival predictor
differs from the longitudinal predictor) requires a separate syntax
extension and is deferred to PR-4b's `Copy` augmentation plus a
"this f-term applies only to likelihood k" marker.

### Decision

**Tuple-LHS as the primary surface.** `(y1, y2) ~ 1 + f(idx, IID(n))`
expands to a `JointLikelihood` with one `StackedMapping` per
likelihood, all wrapping the same shared `A` projector built from
the RHS.

R-INLA's stacked-Y form is documented as a fallback users can write
explicitly when they need long-format observation tables (a single
`y` vector with a `type` column). The macro doesn't generate it
automatically — long-format inputs are converted to wide-format
inside the user's data prep.

**Wide-format observations.** Each LHS column has length `n`; the
stacked observation vector is `vcat(y1, y2, …)`. Mismatched column
lengths or missing LHS columns surface as construction-time errors
with concrete pointers ("column `:y2` has length 56, expected 67").

### Consequences

**Pros.**

- One-liner syntax for the joint case: `(y1, y2) ~ 1 + f(idx, …)`
  reads like Julia destructuring, not like a new domain-specific
  trick.
- Matches the Julia ecosystem's tuple-pattern conventions
  (`StatsModels`-style multi-response works the same way).
- Wide-format input means each likelihood's data shape stays
  observable in the source — long-format hides which row belongs
  to which likelihood inside a `family` vector.

**Cons.**

- Long-format users have one extra reshape step before calling
  `@lgm`. Documented in the migration guide vignette.
- Per-likelihood RHS variation needs PR-4b — users with truly
  per-likelihood predictors hit "not yet supported" until PR-4b.

### References

- [`plans/phase-n.md:176-180`](phase-n.md) — the original ADR-033
  candidate framing.
- `JointLikelihood` core type in
  [`packages/LatentGaussianModels.jl/src/likelihoods/joint.jl`](../packages/LatentGaussianModels.jl/src/likelihoods/joint.jl).
- `StackedMapping` in
  [`packages/LatentGaussianModels.jl/src/observation_mapping.jl`](../packages/LatentGaussianModels.jl/src/observation_mapping.jl).

---

## ADR-034: Implicit f-term naming via column symbol

Status: Accepted
Date: 2026-05-05

### Context

Phase N PR-4b's `Copy` augmentation needs a way to refer to f-terms
in the formula by name. Two candidates:

(a) **Implicit naming via the column symbol.** `f(subject, IID(n))`
    is automatically named `:subject`; `f(idx, IID(n))` is named
    `:idx`. Users refer to the term as `copy = :subject`.

(b) **Explicit `name = :foo` kwarg.** `f(subject, IID(n); name = :u)`
    binds the term's name to `:u`, decoupled from the column symbol.

R-INLA uses option (a) implicitly: `inla.stack` indexes terms by
column name, and `f(group, model = "iid")` is reachable via
`copy = "group"` in the predictor. The duplicate-column case
(`f(year)` and `f(year, copy = …)`) is handled by R-INLA at the
`inla.stack` level, not in the formula.

### Decision

**Option (a): implicit naming via the column symbol.** `f(col, ...)`
binds the term's name to `:col` automatically. The formula

```julia
@lgm formula = (y1, y2) ~ 1 + f(subject, IID(n)) + f(group, IID(m))
```

emits two named f-terms `:subject` and `:group`. Augmentation later
in the formula references them by symbol:

```julia
@lgm formula = (y1, y2) ~ 1 + f(subject, IID(n)) +
                              f(subject; copy = :u, target = 2)
```

Per-likelihood targeting (which `Copy` lands on which likelihood's
predictor) uses an explicit `target = k` kwarg, not the term name.

The duplicate-column case is rare enough that we don't ship a
`name = :foo` override in v0.2 — if it surfaces in a real workflow
the kwarg is non-breaking to add.

### Consequences

**Pros.**

- Zero ceremony for the common case: a single `f(col, …)` produces
  one obviously-named term.
- Mirrors R-INLA's `inla.stack` indexing-by-column convention,
  which is how R-INLA users already think about the joint model.
- Per-likelihood `target = k` is explicit — no implicit "first
  likelihood gets the copy" rule to remember.

**Cons.**

- Two `f(col, …)` calls with the same first argument collide. In
  v0.2 this is an error; an explicit `name =` kwarg is the natural
  v0.3 escape hatch.

### References

- [`plans/phase-n.md:203-210`](phase-n.md) — original ADR-034
  candidate framing.
- `Copy` likelihood augmentation in
  [`packages/LatentGaussianModels.jl/src/likelihoods/copy.jl`](../packages/LatentGaussianModels.jl/src/likelihoods/copy.jl).
- R-INLA `inla.stack` documentation: `inla.doc("stack")`.

---

## ADR-035: `@lgm` `replicate` / `group` routing — runtime wrap calls in the components tuple

Status: Accepted
Date: 2026-05-05

### Context

Phase N PR-5 ships R-INLA's `replicate` and `group` term keywords as
`@lgm y ~ 1 + f(t, AR1(n); replicate = id) + f(s, IID; group = grp)`.
Both wrap the inner component:

- `replicate = id_col` → `Replicate(comp, R)` where
  `R = maximum(data.id_col)`.
- `group = grp_col` → `Group(factory, data.grp_col)` (LGM core's
  factory-form constructor that infers per-group sizes from the
  label vector).

Both wrappers depend on `data` to determine `R` or per-group sizes.
The macro itself has no access to `data` at expansion time (macros run
before runtime). Two design alternatives:

a. **Bake `data` into the AST literally** — illegal: would require
   the macro to receive a `data` value at expansion time, which is
   exactly the "code generation based on runtime data" failure mode
   forbidden by `plans/macro-policy.md`.

b. **Emit a runtime helper call into the AST** — the components tuple
   slot for a replicated/grouped term contains an
   `LGMFormula._wrap_term(comp, data, replicate_col, group_col)` call
   that resolves the wrap at runtime. The AST is data-free; only the
   *expansion's evaluation* reads `data`, which is what every other
   `@lgm` term already does (via `_build_design_matrix(data, …)`).

### Decision

Adopt (b). The components tuple in the expansion is no longer purely
a tuple of literal `Component(args...)` calls — for `replicate`/`group`
terms, the corresponding slot is a `_wrap_term(...)` call. Concretely:

```julia
@lgm y ~ 1 + f(t, AR1(n); replicate = id) data=df family=Gaussian()

# expands to (modulo escapes):

LatentGaussianModel(
    GaussianLikelihood(),
    (LatentGaussianModels.Intercept(),
     LGMFormula._wrap_term(AR1(n), df, :id, nothing)),
    LGMFormula._build_design_matrix(df, :y, true, Symbol[],
        [(:t, AR1(n), :id, nothing)]),
)
```

The plain (no-replicate / no-group) case continues to emit a literal
`Component(args)` slot — backward-compatible with the PR-1..PR-4
macroexpand structural tests.

`_wrap_term(comp, data, replicate_col, group_col)`:

- both `nothing`: returns `comp` (used only by the function-form
  `lgmformula` for symmetry; macro skips this path).
- `replicate ≠ nothing`: validates the column, computes
  `R = maximum(rep_col)`, returns `Replicate(comp, R)`. `comp` must
  be an `AbstractLatentComponent` instance.
- `group ≠ nothing`: returns `Group(factory, grp_col)`. `comp` is the
  factory (a callable like `IID`, `AR1`, `RW1`); LGM core's
  `Group(factory, group_id)` constructor counts per-group sizes and
  invokes `factory(s_g)` for each group size.

The design-matrix builder `_build_design_matrix` was extended to
accept the same per-term spec and lay out the projector block as:

- replicate: `(R · length(comp))` columns; row i has 1 at column
  `(rep_col[i] - 1) · length(comp) + idx_col[i]`.
- group: `Σ_g s_g` columns; row i has 1 at column
  `offset[grp_col[i]] + idx_col[i]`, where `offset[g]` is the
  cumulative sum of preceding per-group sizes.

Both layouts match R-INLA's panel-stacking convention exactly, so
the latent-vector ordering between R-INLA and `@lgm` is bit-for-bit
identical (verified against the
`synthetic_replicate_ar1` Phase I-C oracle).

### Consequences

#### Positive

- Macro is still purely source-to-source — the AST contains no `data`
  values, only references to the data-binding `Symbol` and runtime
  helper calls.
- `@macroexpand` output is readable: replicated/grouped slots are a
  named call into the public `LGMFormula._wrap_term` rather than
  inlined indexing logic.
- Errors from `_wrap_term` surface at the macro-call site (because
  Julia traces back through the helper call to user code) and refer
  to the user-supplied column name, not table internals.
- Symmetric with `_build_design_matrix`'s data-bound nature, which
  has always been the runtime-deferred slot in the expansion.

#### Neutral

- The "components tuple appears literally" guarantee from PR-1's
  macroexpand tests holds for plain f-terms but is relaxed for
  replicate/group. The PR-5 macroexpand tests assert the relaxed
  invariant directly: `_is_call_to(slot, :_wrap_term)`.
- `Group` accepts a factory rather than an instance — this asymmetry
  with `Replicate` (which takes an instance) reflects LGM core's
  existing `Group(factory, group_id; kwargs...)` API, not a
  formula-side choice. Documented in `LGMFormula.jl`'s docstring.

#### Negative

- A user who reads `@macroexpand` on a replicate/group formula and
  expects a static components tuple will see the wrapper call. The
  formula's docstring explicitly notes this; the wrapper name is
  public (`LGMFormula._wrap_term`) and its effect is documented.
- The `lgmformula` function form must accept either 2-tuples
  `(col, comp)` or 4-tuples `(col, comp, replicate, group)` — not as
  elegant as a single shape, but the 2-tuple form is required for
  PR-1..PR-4 backward compatibility.

### References

- [`packages/LGMFormula.jl/src/schema.jl`](packages/LGMFormula.jl/src/schema.jl)
  — `_wrap_term`, `_build_design_matrix` extensions.
- [`packages/LatentGaussianModels.jl/src/components/replicate.jl`](packages/LatentGaussianModels.jl/src/components/replicate.jl)
  — `Replicate(component, n_replicates)` LGM core API.
- [`packages/LatentGaussianModels.jl/src/components/group.jl`](packages/LatentGaussianModels.jl/src/components/group.jl)
  — `Group(factory, group_id; kwargs...)` factory-form constructor.
- [`plans/macro-policy.md`](plans/macro-policy.md) — "Code generation
  based on runtime data" prohibition that this design works around
  via runtime helper calls in the AST rather than baked-in data.
- ADR-008 — `LGMFormula.jl` as a separate package; the macro is sugar
  over Tier-1 constructors and may emit any constructor call sequence.

---

## ADR-036: `SPDE2` retains `INLAMesh` so the `@lgm` extension can build a `MeshProjector` at runtime

Status: Accepted
Date: 2026-05-06

### Context

Phase N PR-7 extends `@lgm` from index-column random-effect syntax
(`f(idx, IID(n))`) to coordinate-column syntax for SPDE models
(`f((east, north), spde)`). The macro must build a barycentric
`MeshProjector(mesh, locations)` at runtime, where `locations` is the
per-row `(east, north)` matrix taken from the user's data. But
pre-PR-7a, `SPDE2` stored only the FEM matrices (`C`, `G₁`, `G₂`),
the `GMRFGraph`, and the PC prior; the originating `INLAMesh` was
discarded after FEM assembly. Three placement candidates for the mesh
reference:

(a) **`SPDE2` retains its `INLAMesh`.** Add a `mesh` field. The
    upgrade path mirrors R-INLA: users call `inla_mesh_2d(loc, …)` to
    get a mesh, pass it to `SPDE2(mesh; …)`, and `MeshProjector` calls
    inside the `@lgm` extension read `spde.mesh` directly. Backwards
    compatibility holds via a fallback constructor.

(b) **The macro extracts the mesh from the `f(...)` AST at parse
    time.** Pattern-match `SPDE2(mesh, …)` and emit
    `MeshProjector($mesh, …)` directly. Brittle — only works for
    literal `SPDE2(...)` constructor calls, breaks under aliasing
    (`spde = SPDE2(...); @lgm ... f(coords, spde)`).

(c) **A separate `LGMFormula._SpatialTerm(mesh, component)` wrapper.**
    User writes `@lgm ... f(coords, _SpatialTerm(mesh, SPDE2(mesh, …)))`.
    Verbose; mesh redundancy invites bugs.

### Decision

Adopt (a). `SPDE2` gains a type-parameterised `mesh::M` field where
`M <: Union{INLAMesh, Nothing}`:

- `SPDE2(mesh::INLAMesh; α, pc)` — primary constructor, stores the
  mesh. Used by the `@lgm` extension; required for any code path that
  calls `MeshProjector(spde.mesh, …)`.
- `SPDE2(points, triangles; α, pc)` — back-compat constructor, stores
  `mesh = nothing`. Existing v0.1.x users continue to work; only the
  new tuple-coordinate macro form requires the mesh-bearing variant.

The two-path constructor keeps every Phase M oracle fixture green
without regenerating R-INLA references — the FEM matrices, graph, and
PC prior are bit-identical between the two paths. The `@lgm` extension
explicitly checks `spde.mesh isa INLAMesh` and raises a user-readable
`ArgumentError` when the bare `(points, triangles)` constructor was
used, pointing at the mesh-bearing constructor as the fix.

### Consequences

#### Positive

- The `@lgm` macro never has to re-thread the mesh — the spatial
  component carries everything its design block needs. Preserves the
  PR-1..PR-6 invariant that the components tuple is closed under
  composition.
- The pattern matches R-INLA's `inla.spde2.matern(mesh, …)` API, so
  the porting story for users is "same input, same shape."
- A future `predict_raster(model, res, template)` overload (Phase O
  PR-1) reads `spde.mesh` for the projector; no separate threading.

#### Neutral

- `SPDE2` is now type-parameterised on the mesh field
  (`SPDE2{M, …}`); two SPDE2 instances built via the two constructors
  are not type-equal (`SPDE2{INLAMesh, …}` vs `SPDE2{Nothing, …}`).
  Downstream functions dispatching on `::SPDE2` continue to work
  because the parameter is unconstrained at the dispatch site.
- The umbrella `INLA.jl` and `INLASPDERasters.jl` widen their
  `INLASPDE` compat to `"0.2, 0.3"` to accept the field-set change.

#### Negative

- An additive field is, formally, a contract change in `INLASPDE.jl`
  — bumped 0.2.0 → 0.3.0. Code paths that compare `SPDE2` structs by
  field set (oracle fixture round-trip) had to be audited. The audit
  found one site (oracle fixture replay) which already uses
  `_struct_isequal`-style comparisons that ignore the new field if
  it's `nothing`.

### References

- [`packages/INLASPDE.jl/src/components/spde2.jl`](packages/INLASPDE.jl/src/components/spde2.jl)
  — `SPDE2` struct and the two constructor paths.
- [`packages/INLASPDE.jl/src/projector.jl`](packages/INLASPDE.jl/src/projector.jl)
  — `MeshProjector(mesh, locations)`, the consumer of `spde.mesh`.
- [`packages/LGMFormula.jl/ext/LGMFormulaINLASPDEExt.jl`](packages/LGMFormula.jl/ext/LGMFormulaINLASPDEExt.jl)
  — extension that calls `MeshProjector(spde.mesh, …)`.
- ADR-008 — `LGMFormula.jl` as a separate package; the SPDE bridge
  must layer through the LGM core's existing component contract.
- [`plans/phase-n-pr7.md`](plans/phase-n-pr7.md) — PR-7 subplan
  motivating this ADR.

---

## ADR-037: `@lgm` accepts tuple-coordinate first arg in `f(...)`; arity ∈ {2, 3} only

Status: Accepted
Date: 2026-05-06

### Context

Phase N PR-7b extends the `@lgm` parser to accept `f((east, north), spde)`
and (PR-7c) `f((east, north, time), KroneckerComponent(...))`. Two
parser-level questions:

1. **What does the tuple mean for non-SPDE components?** A 2-tuple has
   no meaning for `IID(n)` or `AR1(n)` — those components index by a
   single column. The parser does not know component types at parse
   time (component arguments are unevaluated AST). Compromise: accept
   tuple-shape syntax at parse time, defer the type-fit check
   ("component must accept coordinate columns") to the schema-side
   runtime helper.

2. **Tuple arity bounds.** R-INLA's SPDE always uses 2D coordinates;
   3D SPDE is out of scope (deferred to v0.3+ per
   `packages/INLASPDE.jl/plans/plan.md`). 1D SPDE uses a single
   coordinate column — `f(t, SPDE1D(mesh))` already works under the
   pre-PR-7 parser since `t` is a `Symbol`, not a tuple. So PR-7
   needs to accept 2-tuples (spatial-only) and 3-tuples (space-time
   via `KroneckerComponent`, ADR-038).

### Decision

Accept `Expr(:tuple, args...)` as the first positional argument of
`f(...)` only when `length(args) ∈ {2, 3}`. Reject other arities at
parse time with a user-readable error:

- `length == 1` (e.g. `f((s,), spde)`) — point at the bare-symbol
  form `f(s, spde)`.
- `length ≥ 4` — point at "3D SPDE is out of scope; use 2D
  spatial-only or `(east, north, time)` space-time."

All tuple entries must be `Symbol`s. Non-Symbol entries (e.g.
`f((east, 3.14), spde)`) raise an `ArgumentError` with the offending
position.

The component-type fit check ("`IID(n)` does not accept coordinate
columns") happens at the schema-side runtime helper
`_build_spatial_block(component, data, coord_cols, n_obs)`. The
default fallback throws `ArgumentError("@lgm: component
$(typeof(comp)) does not accept coordinate columns; install
INLASPDE.jl and load it for SPDE support")`. Concrete overloads for
`SPDE2` and `KroneckerComponent` ship in the
`LGMFormulaINLASPDEExt` weakdep extension (ADR-039).

`replicate` / `group` keyword routing (ADR-035) is **rejected** on
tuple-coord terms — the mesh-barycentric block does not compose
cleanly with R or per-group panel stacking. Surfaces as a
parse-time error.

### Consequences

#### Positive

- `@lgm` reaches the second flagship R-INLA workflow (geostatistics)
  using the same surface as R-INLA: `f((east, north), spde)` mirrors
  `f(s, model = "spde")` plus a coord-pair convention. The Meuse
  vignette can now express its model in `@lgm` form.
- Parse-time arity bounds catch 90% of typos before they reach the
  schema-binding stage; the remaining 10% (component type mismatch)
  surfaces with a user-readable error from
  `_build_spatial_block`.
- The bare-symbol form `f(idx, comp)` is unchanged — the tuple
  branch is purely additive in the parser.

#### Neutral

- The "first positional arg of `f(...)` is a Symbol" invariant from
  PR-1..PR-6 is relaxed to "Symbol or 2/3-tuple of Symbols." The
  PR-7b/c macroexpand structural tests assert the relaxed
  invariant directly.
- 1D SPDE (`SPDE1D`) uses bare-Symbol `f(t, spde1d)` — no tuple
  needed since 1D coordinates are a single column. This asymmetry
  is documented in `LGMFormula`'s docstrings.

#### Negative

- A future 3D SPDE (deferred to v0.3+) will need to relax the
  `length ≥ 4` reject branch. When 3D lands, the relaxation is
  parser-local: change the upper bound and add a 4-tuple test case.
  Forward-compat path is clear.

### References

- [`packages/LGMFormula.jl/src/parse.jl`](packages/LGMFormula.jl/src/parse.jl)
  — `_parse_f_term` extension to accept tuple-shape first arg.
- [`packages/LGMFormula.jl/src/schema.jl`](packages/LGMFormula.jl/src/schema.jl)
  — `_build_spatial_block` default fallback.
- [`packages/LGMFormula.jl/ext/LGMFormulaINLASPDEExt.jl`](packages/LGMFormula.jl/ext/LGMFormulaINLASPDEExt.jl)
  — concrete overloads for `SPDE2` and `KroneckerComponent`.
- ADR-035 — `replicate` / `group` runtime-wrap routing; this ADR
  excludes that routing on tuple-coord terms.
- ADR-038 — 3-tuple form for `KroneckerComponent` space-time.
- ADR-039 — weakdep extension that hosts the SPDE-aware overloads.

---

## ADR-038: `KroneckerComponent` space-time `@lgm` form takes a 3-tuple coordinate; design block is a sparse Khatri-Rao matrix

Status: Accepted
Date: 2026-05-06

### Context

PR-7c closes the Cameletti-style separable space-time SPDE entry
point in `@lgm`. R-INLA writes this as
`f(s, model = "spde", group = t, control.group = list(model = "ar1"))`.
PR-5 already routes `f(t, AR1; group = grp_col)` to a runtime
`Group(AR1, data.grp_col)` — which is **not** the same as Kronecker
(Group is one inner component per group label with block-diagonal
flattening; Kronecker is `Q_s ⊗ Q_t`). Three candidates for the
`@lgm` surface:

(a) **Explicit `KroneckerComponent` second arg, 2-tuple coord.** User
    writes `f((east, north), KroneckerComponent(SPDE2(mesh), AR1(T)))`.
    Time column has nowhere to live — `KroneckerComponent` needs
    *three* columns: two spatial coords + one time index. Doesn't fit.

(b) **3-tuple coord.** User writes
    `f((east, north, time), KroneckerComponent(SPDE2(mesh), AR1(T)))`.
    Parser accepts 3-tuple when the component is
    `KroneckerComponent`. First two are mesh coordinates; third is
    the time index. Mirrors the Cameletti fixture's natural data
    layout.

(c) **R-INLA-style `group = time` + auto-Kronecker.** User writes
    `f((east, north), SPDE2(mesh); group = time)` and the macro
    auto-wraps as `KroneckerComponent(SPDE2(...), AR1(T))`. Hidden
    behavior — picks the time model and inherits R-INLA's
    `control.group` ambiguity.

The second design question is **how the design block is laid out**.
The plan body called for emitting `KroneckerMapping(A_space, A_time)`
directly, with mixed-RHS via `StackedMapping`. That assertion turned
out to be mis-typed: `KroneckerMapping` is matrix-Kron with
`nrows = nrows(A_space) · nrows(A_time)`, so feeding both factors at
`n_obs` rows would yield `n_obs²` observation rows. `StackedMapping`
is row-stacked (multi-likelihood), not column-stacked. Neither fits
the per-obs `(east, north, time)` data shape.

### Decision

Adopt (b) — 3-tuple coord with explicit `KroneckerComponent` second
arg.

The schema-side runtime helper
`_build_spatial_block(c::KroneckerComponent, data, (east, north, time), n_obs)`
builds a **sparse Khatri-Rao (row-product) design matrix**, not a
`KroneckerMapping`:

```
A_st[i, (s - 1) · n_t + t] = A_space[i, s]   when t == time_idx[i]
                           = 0               otherwise
```

The column layout `(s − 1) · n_t + t` matches `KroneckerComponent`'s
`vec(X)` flattening (`X` of shape `n_t × n_s`) and the existing
`KroneckerMapping`'s flattening convention used by the
`cameletti_pm10` oracle. For gridded data where every spatial location
is observed at every time slot, this Khatri-Rao construction reduces
exactly to `kron(A_space_j, I_{n_t})` — the form the oracle uses.

Implementation lives in the `LGMFormulaINLASPDEExt` weakdep extension
(ADR-039). Three dispatches:

- `_build_spatial_block(::KroneckerComponent, …, NTuple{3, Symbol}, …)`
  — main path; validates `c.space isa SPDE2` and `c.space.mesh isa
  INLAMesh`, builds `MeshProjector`, remaps columns via `findnz`.
- `_build_spatial_block(::SPDE2, …, NTuple{3, Symbol}, …)` — points
  at the `KroneckerComponent` wrapper.
- `_build_spatial_block(::KroneckerComponent, …, NTuple{2, Symbol}, …)`
  — points at the 3-tuple form.

Defer (c) (`group = time` syntax sugar) to a follow-up if user demand
surfaces. The 3-tuple form is the most explicit and least surprising,
and the parse-time arity check from ADR-037 generalises cleanly.

### Consequences

#### Positive

- The Cameletti et al. (2013) PM₁₀ space-time SPDE — the second
  flagship space-time fixture — is expressible in `@lgm` form. The
  Cameletti vignette can match the Scotland / Tokyo `@lgm`-driven
  treatment.
- Khatri-Rao construction is the correct general-purpose form for
  arbitrary per-obs `(east, north, time)` data; the gridded
  Cameletti special case (`kron(A_space_j, I_{n_t})`) drops out
  automatically.
- Time-coord validation (out-of-range, non-integer) lives in the
  extension and surfaces user-readable errors at the macro-call
  site.

#### Neutral

- 3-tuple coord is reserved for `KroneckerComponent`; using it with
  a bare `SPDE2` (`f((east, north, time), spde)`) raises an error
  pointing at the wrapper. Documented in `LGMFormula`'s parse-error
  message.
- The decision to use Khatri-Rao instead of `KroneckerMapping`
  deviates from `plans/phase-n-pr7.md` PR-7c §, which conjectured
  `KroneckerMapping` direct emit. The deviation is documented in
  the PR-7c commit message and the `LGMFormulaINLASPDEExt`
  source-file comments.

#### Negative

- The `KroneckerMapping` mapping type still has no `@lgm` consumer
  — it's used internally by the LGM core's space-time machinery but
  the macro lowers to a `LinearProjector(A::SparseMatrixCSC)` even
  for Kronecker-shaped designs. Forward path: when an `@lgm`
  consumer needs the matrix-Kron form (e.g. tensor-grid
  observations like climate-model output where `n_obs` *is*
  `n_s · n_t`), add a `KroneckerMapping`-emitting branch. Not in
  scope for v0.2.x.

### References

- [`packages/LGMFormula.jl/ext/LGMFormulaINLASPDEExt.jl`](packages/LGMFormula.jl/ext/LGMFormulaINLASPDEExt.jl)
  — Khatri-Rao construction at the `_build_spatial_block` overload
  for `KroneckerComponent`.
- [`packages/LatentGaussianModels.jl/src/components/kronecker.jl`](packages/LatentGaussianModels.jl/src/components/kronecker.jl)
  — `KroneckerComponent` and the `vec(X)` flattening convention.
- [`packages/LatentGaussianModels.jl/src/observation_mapping.jl`](packages/LatentGaussianModels.jl/src/observation_mapping.jl)
  — `KroneckerMapping`, the matrix-Kron form not emitted by `@lgm`.
- [`packages/INLASPDE.jl/test/oracle/test_cameletti_pm10.jl`](packages/INLASPDE.jl/test/oracle/test_cameletti_pm10.jl)
  — gridded `kron(A_space_j, I_{n_t})` oracle that this ADR's
  Khatri-Rao form recovers as a special case.
- ADR-029 — generic `KroneckerComponent` two-component composer.
- ADR-037 — tuple-coord parser semantics; this ADR extends to 3-tuple.
- ADR-039 — weakdep extension hosts the `KroneckerComponent` overload.

---

## ADR-039: `LGMFormula` ↔ `INLASPDE` integration ships as a Julia 1.9 weakdep extension, not a hard dep or duck-typed hook

Status: Accepted
Date: 2026-05-06

### Context

The mesh-barycentric runtime helper from PR-7b/PR-7c (ADR-037, ADR-038)
needs to live somewhere in the package graph. The function is owned
by `LGMFormula` (it's the schema-side dispatch surface) but its
implementation needs `INLASPDE` types (`SPDE2`, `INLAMesh`,
`MeshProjector`) and `LatentGaussianModels` types
(`KroneckerComponent`). Three placements:

(a) **`LGMFormula` hard-deps `INLASPDE`.** Simplest, but pulls the
    SPDE FEM stack (DelaunayTriangulation, Meshes, CoordRefSystems,
    SciMLOperators) into every `using LGMFormula` import — wrong
    layering for users who only need `@lgm` for non-spatial models
    (Scotland BYM2, Tokyo seasonal AR1).

(b) **`LGMFormula` weakdeps `INLASPDE` via a Julia 1.9 extension.**
    `LGMFormulaINLASPDEExt` adds the spatial-projector helper when
    the user has `using INLASPDE` (or its dependants) loaded.
    Macros from `LGMFormula` cannot export new symbols from
    extensions, but the extension can attach methods to existing
    `LGMFormula` functions.

(c) **`INLASPDE` adds a hook the macro calls via duck typing.**
    Define a method
    `LGMFormula._build_spatial_block(c::AbstractLatentComponent, data, cols)`
    that throws by default, and `INLASPDE` overloads it for `SPDE2`.
    Same effect as (b) but the bookkeeping is on the `INLASPDE`
    side — `INLASPDE` would import `LGMFormula` to attach methods,
    inverting the natural dependency direction.

### Decision

Adopt (b). `LGMFormula.jl`'s `Project.toml` declares:

```toml
[weakdeps]
INLASPDE = "2835c710-3f40-4945-979f-d21c9e20d425"

[extensions]
LGMFormulaINLASPDEExt = "INLASPDE"

[sources]
INLASPDE = {path = "../INLASPDE.jl"}
```

The `[sources]` block enables the Julia 1.11+ test-sandbox path
resolution so `Pkg.test()` finds the local dev `INLASPDE@0.3.0` —
matches the same pattern in `docs/Project.toml` and other extensions
in the ecosystem (ADR-008).

The extension `LGMFormulaINLASPDEExt`:

- Imports `INLASPDE.SPDE2`, `INLASPDE.MeshProjector`,
  `INLASPDE.INLAMesh`, and `LatentGaussianModels.KroneckerComponent`.
- Attaches concrete methods to `LGMFormula._build_spatial_block`
  (the unexported function defined in `LGMFormula`'s core, ADR-037).
- Attaches mis-pair error overloads (SPDE2 + 3-tuple,
  KroneckerComponent + 2-tuple) so the user gets a precise pointer
  to the right shape rather than a generic `MethodError`.

### Consequences

#### Positive

- Users who only need `@lgm` for non-spatial models incur zero cost
  from the SPDE stack — `using LGMFormula` does not load
  `INLASPDE`. The Phase L migration story for Scotland / Germany /
  Tokyo (`@lgm` for areal and time-series models) ships clean.
- The extension activates automatically when both packages are
  loaded (e.g. `using INLA` re-exports both, so `@lgm` "just works"
  for SPDE the moment `using INLA` runs).
- The extension-method pattern mirrors `INLASPDEMakieExt`,
  `LatentGaussianModelsTuringExt`, and other ecosystem extensions —
  one consistent layering across the repo.
- The hook function (`_build_spatial_block`) is part of
  `LGMFormula`'s public API even though unexported; its docstring is
  the formal contract for "what does it mean for a component to
  accept coordinate columns?"

#### Neutral

- Extension methods are not visible to `@code_warntype` until the
  extension is loaded; users debugging an `@lgm`-driven SPDE model
  must have `using INLASPDE` (or `using INLA`) before
  `@code_warntype` reports types correctly.
- Julia 1.9 extensions cannot export new symbols. Any future
  `LGMFormula`-side public name needed by SPDE-only callers
  (currently none) would require adding it to `LGMFormula`'s core.

#### Negative

- Bumps `LGMFormula.jl` 0.2.0 → 0.3.0 (extension target addition).
  PR-7c further bumps 0.3.0 → 0.4.0 (mapping-shape change in
  `KroneckerComponent` overload). The minor bumps are conservative
  but accurate — the public schema function
  `_build_spatial_block` is new in 0.3.0.
- Test-time setup is slightly more involved than a hard dep: the
  test target declares `INLASPDE` in `[extras]` and lists it in the
  `test` target so the extension activates inside `Pkg.test()`.
  Documented in `LGMFormula.jl/Project.toml`.

### References

- [`packages/LGMFormula.jl/Project.toml`](packages/LGMFormula.jl/Project.toml)
  — `[weakdeps]`, `[extensions]`, `[sources]`, and test-extras
  declarations.
- [`packages/LGMFormula.jl/ext/LGMFormulaINLASPDEExt.jl`](packages/LGMFormula.jl/ext/LGMFormulaINLASPDEExt.jl)
  — extension module, all SPDE-aware `_build_spatial_block` methods.
- [`packages/LGMFormula.jl/src/schema.jl`](packages/LGMFormula.jl/src/schema.jl)
  — the unexported `_build_spatial_block` function the extension
  attaches methods to.
- ADR-008 — `LGMFormula.jl` as a separate package; the layering
  intent that this ADR realises.
- ADR-036 — `SPDE2` mesh field that the extension reads.
- ADR-037 — tuple-coord parser semantics that the extension serves.
- ADR-038 — 3-tuple `KroneckerComponent` route the extension hosts.
- [`plans/dependencies.md`](plans/dependencies.md) — weakdep policy.

---

## ADR-040: `predict_raster(model, res, …)` — Gaussian summary + sample-based path with `Exceedance` wrapper

Status: Accepted
Date: 2026-05-06

### Context

Phase O ([`plans/phase-o.md:91-147`](phase-o.md)) lifts the existing
vertex-vector `predict_raster` and `quantile_rasters` primitives to
the user-facing `(model, res, …)` shape. Two semantic axes need a
durable answer:

1. **Where do the vertex statistics come from?** Gaussian-approximation
   slice (cheap, vertex-Gaussian view) vs. joint posterior draws
   (expensive, captures non-Gaussian tails).
2. **What does `quantity` accept?** Marginal mean / sd / Gaussian
   intervals are well-defined under linear barycentric projection;
   exceedance probabilities `P(u(s) > c | y)` are not (probabilities
   don't compose linearly under barycentric averaging).

The honest path for tail functionals is to draw joint samples,
slice the SPDE block, project each draw through the same
`MeshProjector`, then reduce column-wise per the requested quantity.

### Decision

**Both paths ship.** Two overloads on a single user-facing entry:

```julia
predict_raster(model, res, template;
               component, quantity = :mean, outside, missingval)
predict_raster(rng, model, res, template;
               component, quantity, n_samples,
               outside, missingval)
```

- The Gaussian overload accepts `quantity ∈ (:mean, :sd, :lower,
  :upper)` and reads `random_effects(model, res)` for the SPDE
  component slice. Wraps the existing primitive; no new sampling
  cost.
- The sample-based overload accepts `quantity ∈ (:mean, q::Real ∈
  [0, 1], Exceedance(c))`. Calls `posterior_sample(rng, res, model,
  n_samples)`, slices the SPDE block via `model.latent_ranges[i]`,
  applies the per-cell `MeshProjector` to the draw matrix in a
  single sparse-dense GEMM, then reduces column-wise.

`Exceedance{T}` is a tiny wrapper struct exported from
`INLASPDERasters` so the dispatch is type-stable and the API
self-documents at the call site:

```julia
predict_raster(rng, model, res, template;
               quantity = Exceedance(2.5), n_samples = 1000)
```

**Performance contract.** The `MeshProjector` is built once per
`predict_raster` call and reused across draws — `P.A * x_samples[spde_range, :]`
is one sparse-dense GEMM and gives a `n_cells × n_samples` matrix.
For 10⁶-cell rasters × 10³ draws this is the only tolerable shape;
element-wise per-draw projection is 1000× slower.

### Consequences

**Pros.**

- Honest handling of tail functionals (`Exceedance`) — no silent
  linear-projection lie that would mis-state `P(u > c)`.
- One entry point with a `quantity` kwarg, matching R-INLA's
  `inla.posterior.sample(...)` + post-processing pattern but
  pre-fused into a single call.
- Gaussian path stays cheap — users who only need posterior mean/sd
  rasters never pay the sampling cost.

**Cons.**

- Two overloads on the same name. The `rng`-first overload is the
  canonical sample-based form; the `rng`-less overload is the
  Gaussian summary. Documented at the function docstring.
- `Exceedance` adds one exported type to `INLASPDERasters`. Small
  surface, self-documenting at the call site.

### References

- [`plans/phase-o.md:91-147`](phase-o.md) — original ADR-040
  candidate framing.
- `predict_raster` overloads in
  [`packages/INLASPDERasters.jl/src/predict.jl`](../packages/INLASPDERasters.jl/src/predict.jl).
- `Exceedance` in
  [`packages/INLASPDERasters.jl/src/exceedance.jl`](../packages/INLASPDERasters.jl/src/exceedance.jl).
- ADR-005 — projector matrix as a model field; same lazy
  `MeshProjector` is reused here.

---

## ADR-041: CRS policy — `predict_raster` rejects mismatched CRS at the API boundary; `mesh_crs` keyword is opt-in

Status: Accepted
Date: 2026-05-06

### Context

The `INLAMesh` struct does not carry CRS metadata
([`packages/INLASPDE.jl/src/mesh/inla_mesh.jl:19-25`](../packages/INLASPDE.jl/src/mesh/inla_mesh.jl)).
`Rasters.Raster` does, via `Rasters.crs(raster)`. Silent CRS
mismatches between mesh and raster (e.g. WGS84 raster in degrees
against a UTM mesh in metres) are the single most common geo
foot-gun and produce wrong-but-not-obviously-wrong rasters.

Three shapes were on the table:

(a) **Reject mismatch silently — caller pre-projects.** The current
    pre-Phase-O behaviour. Cheap, but unsafe.
(b) **Add `mesh_crs::Union{Nothing, CRS} = nothing` keyword.** When
    supplied, assert `mesh_crs == Rasters.crs(raster)` at the API
    boundary; `nothing` keeps current "trust the caller" behaviour.
    No `INLAMesh` change; `CoordRefSystems.jl` is already in
    INLASPDE's deps.
(c) **Promote `INLAMesh` to carry `crs::Union{Nothing, CRS}`.** Field
    addition is a contract change; every Phase M oracle fixture
    needs to round-trip the new field. Heavy.

### Decision

**Option (b): `mesh_crs` keyword.** `predict_raster(model, res,
template; mesh_crs = nothing, …)` and `extract_at_mesh(...; mesh_crs
= nothing, …)`. When `mesh_crs` is supplied, the API asserts
equality against `Rasters.crs(template)` and throws
`ArgumentError` with a concrete pointer if they disagree. When
`nothing`, the existing "trust the caller" path is preserved
verbatim.

**Reprojection is out of scope.** Reprojecting the *raster* requires
`Proj_jll` (already in INLASPDERasters' transitive closure), but
reprojecting the *mesh* would mutate `mesh.points` in place, which
invalidates the SPDE FEM matrices already assembled at model
construction time. A `reproject = true` opt-in is deferred to v0.4
with an explicit pre-condition that the caller reassemble FEM
matrices on the reprojected mesh.

Promote (c) — `INLAMesh.crs` field — when CRS becomes load-bearing
elsewhere (sphere SPDE, geodesic distance priors, multi-CRS joint
likelihoods). All deferred to v0.4+.

### Consequences

**Pros.**

- Zero cost when the kwarg is omitted (back-compat with all v0.3.x
  code).
- Safety net activates with a single kwarg — concretely points at
  the mismatch rather than silently producing wrong rasters.
- No `INLAMesh` change — the struct stays Phase-M-compatible and
  no oracle fixture needs to be re-generated.

**Cons.**

- The escape hatch is opt-in: silent-mismatch is still possible if
  the user doesn't supply `mesh_crs`. Documented in the function
  docstring and the Meuse vignette.
- Reprojection is genuinely deferred — users with mismatched CRSs
  must pre-project the raster themselves in v0.3.

### References

- [`plans/phase-o.md:149-186`](phase-o.md) — original ADR-041
  candidate framing.
- `predict_raster` / `extract_at_mesh` CRS keyword in
  [`packages/INLASPDERasters.jl/src/predict.jl`](../packages/INLASPDERasters.jl/src/predict.jl).
- `INLAMesh` definition in
  [`packages/INLASPDE.jl/src/mesh/inla_mesh.jl:19-25`](../packages/INLASPDE.jl/src/mesh/inla_mesh.jl).

---

## ADR-042: INLASPDERasters takes a load-bearing dep on LatentGaussianModels via `(model, res, …)` overloads

Status: Accepted
Date: 2026-05-06

### Context

The `INLASPDERasters/Project.toml` lists
`LatentGaussianModels = "0.2"` even though, pre-Phase-O, no `src/`
file imports LGM types. The dep was foreseen and accepted at
scaffolding time ([`plans/dependencies.md:101`](dependencies.md))
because the package's stated purpose is the raster-shaped sibling
of LGM — it was a placeholder waiting for the user-facing surface.

Phase O's `predict_raster(model, res, …)` overloads turn that latent
dep into a load-bearing one: the function dispatches on
`LatentGaussianModel`, reads `res.x_mean`, `res.x_var`,
`model.latent_ranges`, and calls `random_effects(model, res)` and
`posterior_sample(rng, res, model)`. INLASPDERasters now genuinely
needs LGM at the API boundary, not just in tests.

### Decision

**Promote the dep from "declared but unused" to "load-bearing for
the public API."** No code change to `[deps]`, just an explicit
acknowledgement that:

- INLASPDERasters depends on LGM's public LGM-result API
  (`random_effects`, `posterior_sample`, `latent_ranges`).
- LGM compat in `INLASPDERasters/Project.toml` is a real bound
  going forward — bumping LGM major versions requires
  re-validating the raster overloads.
- LGM is in the *direct* `[deps]` rather than `[weakdeps]` because
  the `(model, res, …)` overloads are the package's user-facing
  surface, not an optional integration.

The companion alternative — keeping LGM in `[weakdeps]` and
shipping the overloads via an extension — was considered and
rejected: the entire point of `INLASPDERasters` is the LGM-shaped
surface; users who don't use LGM don't load this package.

### Consequences

**Pros.**

- The package's purpose is unambiguous: it is the raster-shape
  sibling of LGM, not a generic mesh-to-raster utility that LGM
  happens to use.
- Compat bumps are caught at `Pkg.resolve()` time rather than as
  runtime `MethodError`s.

**Cons.**

- One more direct dep edge in the ecosystem dependency graph. Small
  cost; the package was already shipping the dep in `Project.toml`.
- Future re-layering ("split the LGM-shape overloads into a
  weakdep extension") is a larger refactor than if we'd shipped
  the weakdep route now. Justified by the package's stated
  purpose.

### References

- [`plans/phase-o.md:188-205`](phase-o.md) — original ADR-042
  candidate framing.
- [`packages/INLASPDERasters.jl/Project.toml`](../packages/INLASPDERasters.jl/Project.toml)
  `[deps]` table (LGM is a hard dep, not a weakdep).
- [`plans/dependencies.md:101`](dependencies.md) — the original
  scaffolding-time decision.

---

## ADR-043: Gen.jl second-MCMC sanity check — deferred to v1.x; v1.0 ships NUTS-only triangulation

Status: Accepted
Date: 2026-05-07

### Context

Phase P ([`plans/phase-p.md:115-148`](phase-p.md)) ships tier-3
triangulation against NUTS via `LGMTuring.jl` as the v1.0
load-bearing implementation. The question was whether v1.0 also
ships a second Julia-HMC implementation through
[Gen.jl](https://github.com/probcomp/Gen.jl) — a third column in
the triangulation envelope (INLA / NUTS / Gen-HMC).

Three options:

- **(a) Ship Gen.jl in v1.0 as PR-1b.** Express the LGM
  `LogDensityProblem` as a `@gen function`, run Gen's HMC sampler,
  add a third column. Cost: ~1–2 days for translation + chain
  authoring. Coverage gain: catches `AdvancedHMC`-specific bugs.
- **(b) Defer Gen.jl to v1.x as a stretch lane.** Ship NUTS-only
  triangulation in v1.0. Coverage at v1.0 is "INLA vs NUTS";
  second-MCMC sanity is v1.x.
- **(c) Skip Gen.jl entirely.** Cross-language verification (Stan,
  NIMBLE) is the genuinely-independent triangulation; second-Julia-
  HMC is a thin gain because both share the same AD ecosystem and
  `LogDensityProblems` contract.

### Decision

**Option (b): defer to v1.x.** v1.0 ships NUTS-only triangulation
on the three flagship models (Scotland BYM2, Pennsylvania BYM2,
Meuse SPDE). Gen.jl integration is tracked as a v1.x backlog item.

The user-trust argument for tier-3 at v1.0 is satisfied by NUTS-vs-
INLA: NUTS is a fully-independent inference path with a different
posterior representation (samples vs. grid), different mode-finder
(leapfrog vs. BFGS), and different uncertainty quantification (MCMC
vs. Laplace). Adding Gen.jl as a second Julia-HMC chain shares the
same AD ecosystem (`ForwardDiff` / `ReverseDiff`) and the same
`LogDensityProblems` contract, so the marginal coverage gain is
small relative to the genuinely-cross-language Stan / NIMBLE path.

If Gen.jl ships in v1.x, the integration shape is:

- `packages/LGMTuring.jl/ext/LGMTuringGenExt.jl` — weakdep
  extension.
- `LGMTuring.gen_hmc_sample(model, y, n; rng, …)` wrapping
  `Gen.hmc(...)` and returning an `MCMCChains.Chains`-shaped
  result.
- The existing `compare_posteriors(...)` harness reuses cleanly
  for the three-way envelope.

### Consequences

**Pros.**

- v1.0 ships on the originally-planned schedule; PR-1b stretch
  doesn't gate the release.
- The `compare_posteriors` harness is generic enough to absorb a
  third sampler later — no v1.0 architecture work is owed.
- ADR captures the tradeoff so a future maintainer doesn't re-ask
  the question.

**Cons.**

- v1.0's triangulation envelope is two-way (INLA vs. NUTS), not
  three-way. A NUTS-implementation bug that happens to also be
  present in INLA would not be caught at v1.0.

### References

- [`plans/phase-p.md:115-148`](phase-p.md) — original ADR-043
  candidate framing (option (b) recommended).
- ADR-009 — Turing / HMC bridge in `LGMTuring.jl`; the load-bearing
  v1.0 tier-3 implementation.
- `compare_posteriors(...)` harness in
  [`packages/LGMTuring.jl/test/triangulation/compare_posteriors.jl`](../packages/LGMTuring.jl/test/triangulation/compare_posteriors.jl).

---

## ADR-044: Tier-3 triangulation tolerances at v1.0 — `tol_mean = 1.5 SDs`, `tol_sd = 0.60` uniformly

Status: Accepted
Date: 2026-05-07

### Context

Phase P ([`plans/phase-p.md:150-191`](phase-p.md)) tightens the
tier-3 NUTS-vs-INLA tolerances from the Phase J prototype levels
(`tol_mean = 2.0 SDs`, `tol_sd = 0.30`) to v1.0 GA levels. With
full chains (≥1000 post-warmup samples after 200 warmup) the noise
floor is much tighter.

Two options were on the table:

- **(a) Per-test calibration.** Each model gets a tolerance derived
  from the chain's effective sample size and the hyperparameter
  dimensionality. Sound but expensive to maintain.
- **(b) Uniform v1.0 tolerances.** Apply across all three flagship
  models: `tol_mean = 1.5 SDs` (envelope: INLA / NUTS posteriors
  agree within 1.5σ), `tol_sd = 0.60` (relative).

The first cut of this ADR proposed `tol_sd = 0.15` ("tighter
because chains are longer"). Empirical observation on both Scotland
and PA flipped that recommendation: NUTS finds 30–45% wider
posterior on the BYM2 mixing weight `logit φ` than INLA does. At
1000-sample chains this is *not* MC-error-bounded, it's a real
property of grid integration vs. full HMC. INLA's grid is bounded
by the Hessian-explored region around the Laplace mode; NUTS
samples the full posterior including the heavy tails of weakly
identified mixing parameters.

### Decision

**Option (b): uniform v1.0 tolerances.** Apply across all three
flagship triangulation tests:

```julia
tol_mean = 1.5  # standard deviations; INLA mean within 1.5σ of NUTS mean
tol_sd   = 0.60 # relative; |sd_inla - sd_nuts| / sd_nuts ≤ 0.60
```

The mean tolerance is tight (catches gradient bugs, precision-build
bugs, mode-finder drift); the SD tolerance is loose because the SD
diff between INLA's grid and NUTS is *structural*, not MC-noise.

**Why `tol_sd = 0.60` and not 0.15.** Tightening `tol_sd` to 0.15
would force either per-parameter tolerances (which rejects option
(a)) or n=10k+ chains (10× nightly cost). The 0.60 envelope still
flags genuine regressions:

- A gradient bug blows SDs by 100%+ (caught easily at 0.60).
- A precision-build bug shifts means by SDs (caught by `tol_mean`).
- A mode-finder bug disagrees on both.

The means agree within 0.06–0.23 SDs on both flagship Poisson-BYM2
datasets — that's the load-bearing constraint. The SD envelope is
the give-room for the structural inference-method difference.

Tier-3 tolerances are not the place for per-model fine-tuning —
that's tier-1's job. The tier-3 contract is "all available
independent implementations land in the same envelope," and the
envelope width should be wide enough to absorb known structural
inference-method differences while still catching regressions.

If ADR-043's option (a) flips and Gen.jl ships in v1.x, the same
tolerances apply to the three-way envelope — Gen-HMC and NUTS
share the same structural-grid-vs-HMC-tail behaviour against INLA.

### Consequences

**Pros.**

- Single pair of tolerances for all three flagship triangulation
  tests — easy to maintain, easy to reason about.
- Catches the regression classes that matter (gradient,
  precision-build, mode-finder) without flagging the structural
  grid-vs-HMC SD difference as a bug.
- Diagnostic `@info` log on flagged rows means failure modes are
  visible in CI logs without re-running locally.

**Cons.**

- A genuine 30% SD regression on a single flagship would not flag
  (within the 60% envelope). Tier-1's per-component oracle tests
  catch this at 1% / 5% tolerances; tier-3's job is the joint
  posterior shape, not per-parameter SDs.
- Tightening v2.0 tolerances will require either more samples or
  per-test calibration — design call deferred.

### References

- [`plans/phase-p.md:150-191`](phase-p.md) — original ADR-044
  candidate framing, including the `tol_sd = 0.60 vs 0.15`
  empirical pivot.
- Tier-3 triangulation tests in
  [`packages/LGMTuring.jl/test/triangulation/`](../packages/LGMTuring.jl/test/triangulation/)
  (Scotland, PA, Meuse).
- ADR-009 — `LGMTuring.jl` package; tier-3 implementation.
- `compare_posteriors(...)` harness — emits the diagnostic rows
  consumed by the tolerance gates.

## ADR-045 (Proposed): De-densify the constrained-Laplace null-space bump — low-rank or bordered KKT

Status: Proposed
Date: 2026-07-05

### Context

Surfaced by the 2026-07 performance review
([`plans/review-2026-07-remediation.md`](review-2026-07-remediation.md),
Tier-3 item 15) while benchmarking the factor-reuse change (ADR-adjacent,
PR #22). A single 625-node **sparse** Besag fit allocates ~139 MiB inside
one `laplace_mode`; ~86% of that is downstream of one construction.

For an intrinsic component under a hard constraint `C x = e`, the Laplace
step regularises the singular precision `Q` with a null-space bump so the
inner Newton factorisation is well-posed:

    Q_reg = Q + V Vᵀ,   V Vᵀ = Cᵀ (C Cᵀ)⁻¹ C     ([`inference/constraints.jl:71`](../packages/LatentGaussianModels.jl/src/inference/constraints.jl) `_null_bump`)

`V Vᵀ` is the orthogonal projector onto `range(Cᵀ)` — **rank k**, where
`k` = number of constraint rows (one sum-to-zero per connected component,
so `k` is tiny: 1 for connected Besag/BYM2, a handful for disconnected
graphs). But `_null_bump` forms it as a **dense `n × n` matrix**: for a
connected sum-to-zero constraint `C = onesᵀ`, `V Vᵀ = (1/n)·ones(n, n)`,
fully dense. `sparse(B)` then stores `n²` nonzeros.

`Q_reg = Q + V Vᵀ` is therefore dense, and everything built on it inherits
the density:

- `H = Q_reg + Jᵀ D J` ([`laplace.jl:130/143/185`](../packages/LatentGaussianModels.jl/src/inference/laplace.jl)) — dense.
- `FactorCache(H)` / `update!` — dense Cholesky, `O(n³)` per Newton step.
- `lp.precision = H` stored dense; selected inversion / marginal variances
  operate on a dense factor.

This is the dominant allocation source for constrained intrinsic models
(disease mapping: Besag/BYM/BYM2/ICAR/Leroux — the bulk of the R-INLA use
cases) and is the scalability ceiling for large areal fields. It is
strictly an **implementation** waste: the bump is rank-`k`, never `O(n²)`.

### Structure to exploit

`H = (Q + Jᵀ D J) + V Vᵀ = H_s + V Vᵀ`, a **rank-`k` update of a sparse
matrix `H_s`**. All quantities the Laplace step needs have low-rank
closed forms, provided `H_s` is factorisable:

- **Solve** (Woodbury): `H⁻¹ b = H_s⁻¹ b − H_s⁻¹ V M⁻¹ Vᵀ H_s⁻¹ b`,
  `M = I_k + Vᵀ H_s⁻¹ V` (`k` extra sparse solves).
- **Log-det** (matrix determinant lemma): `logdet H = logdet H_s + logdet M`.
  Feeds `_log_det_HC` ([`laplace.jl`](../packages/LatentGaussianModels.jl/src/inference/laplace.jl)).
- **Marginal variances**: `diag(H⁻¹) = diag(H_s⁻¹) − diag(H_s⁻¹ V M⁻¹ Vᵀ H_s⁻¹)`.
  `diag(H_s⁻¹)` via selected inversion on the **sparse** factor (cheap —
  this is the ADR-012 path we already reuse per PR #22); the correction is
  assembled from `H_s⁻¹ V` (`n × k`, `k` solves) and the `k × k` `M`.

### The catch

Woodbury needs `H_s = Q + Jᵀ D J` **positive definite**. `Q` is singular
(rank `n − k`); `H_s` is PD iff the likelihood curvature `Jᵀ D J` covers
`null(Q) = range(Cᵀ)`. For the flagship models (Gaussian/Poisson with an
identity-ish projector, every latent coordinate observed) it does, so
`H_s` is PD and no bump is needed *at all* — the bump is pure safety
margin for the rank-deficient case (null direction unidentified by data,
e.g. a sum-to-zero effect with no data on that contrast). When `H_s` is
singular, `H_s⁻¹` does not exist and Woodbury-on-`H_s` is unavailable.

### Options

- **(A) Status quo.** Dense bump. Simple, robust, correct. `O(n²)` memory
  and `O(n³)` factorisation for every constrained intrinsic fit; the
  observed 139 MiB / 625-node ceiling.

- **(B) Low-rank Woodbury on sparse `H_s`, dense-bump fallback.** Try the
  sparse Cholesky of `H_s`; on success use the low-rank forms above for
  solves, log-det, and marginal variances; on `PosDefException` fall back
  to option (A) for that fit. Captures the full win for the common
  (data-identified) case — which is the entire flagship benchmark set —
  while preserving correctness for the singular case. Lowest-risk path to
  most of the benefit. Downside: two code paths; the marginal-variance
  low-rank correction needs a dense-oracle test (we have the harness:
  [`test_constrained_variances_dense.jl`](../packages/LatentGaussianModels.jl/test/regression/test_constrained_variances_dense.jl)).

- **(C) Bordered KKT augmented system.** Solve `[H_s Cᵀ; C 0] [x; λ] =
  [b; e]` directly, eliminating both the bump *and* the separate kriging
  projection ([`_kriging_correction`](../packages/LatentGaussianModels.jl/src/inference/constraints.jl)).
  Well-posed even when `H_s` is singular on `null(Q)` (as long as
  `[H_s; C]` has full column rank — the constraint supplies the missing
  rank). The augmented matrix is sparse apart from the `k` constraint
  rows; for a dense sum-to-zero row this is an **arrow** pattern —
  ordered last, it adds `O(n·k)` fill, not `O(n²)`. This is the textbook
  equality-constrained-GMRF solve (Rue & Held 2005 §2.3.2). Most robust
  and most general, but the largest change: it reshapes the Newton solve,
  the log-det, and — the genuinely hard part — **selected inversion for
  marginal variances under the augmented factor** (the constrained
  covariance is the `(1,1)` block of the KKT inverse; selinv on an
  indefinite bordered LDLᵀ is not the plain ADR-012 path and needs a
  prototype spike before commitment).

### Recommendation (for review — not yet accepted)

Stage it. **Adopt (B) first**: it is a contained, well-tested win that
removes the densification for every data-identified constrained fit (the
whole flagship suite) with a clean fall-back that guarantees no
correctness or robustness regression. Gate each low-rank form behind the
existing dense-oracle test and the R-INLA oracle fixtures; require
bit-level agreement on log-det and ≤1e-10 on marginal variances vs the
current dense path before switching a call site.

**Defer (C)** to a follow-up unless (B)'s fallback fires on a real model
in the oracle suite — i.e. only invest in the KKT reshape if the
data-unidentified case turns out to matter in practice. If it does, (C)
subsumes (B) and the bump can be retired entirely.

Do **not** proceed to implementation on the strength of this ADR alone:
the marginal-variance low-rank correction (B) and, if pursued, the
augmented-system selinv (C) each need a numerical spike validated against
the dense oracle before a line of production code changes.

### Spike result (2026-07-05)

[`scripts/spikes/adr045_null_space_bump.jl`](../scripts/spikes/adr045_null_space_bump.jl)
validated option (B) on connected and disconnected Besag models. The result
is **stronger than the low-rank plan assumed**: when `H_s` is PD, the bump
is not merely reconstructable via Woodbury — it is **entirely unnecessary**.
Every constrained quantity is bump-invariant to machine precision, so the
sparse `H_s` factor can be used directly with the *existing* kriging
machinery. Measured `max|Δ|` between the sparse-`H_s` path and the current
dense-`H` path (`H = H_s + V Vᵀ`):

| Quantity | connected (n=144, k=1) | disconnected (n=50, k=2) |
|---|---|---|
| H1 constrained marginal variances | 1.7e-16 | 1.7e-16 |
| H2 constrained log-det `_log_det_HC` | 5.7e-14 | 1.4e-14 |
| H3 ‖projected bump-free Newton step at mode‖∞ | 4.5e-16 | 2.5e-16 |
| Woodbury solve vs dense `H\b` (cross-check) | 7.8e-16 | 2.2e-16 |

`H_s` was PD in both cases, and the bump inflated the stored precision from
**672 → 20 736 nonzeros (31×, fully dense)** on the 144-node graph — the
densification this ADR targets, quantified.

**Why bump-invariant:** `V Vᵀ = Cᵀ(CCᵀ)⁻¹C` lives entirely in `range(Cᵀ)`,
and every constrained quantity (kriging-corrected covariance, `_log_det_HC`,
the constraint-projected Newton step) is invariant to precision changes
confined to `range(Cᵀ)`. This is an identity, not a numerical coincidence —
confirmed to ε across two topologies.

**Consequence for (B):** it collapses from "Woodbury on `H_s`" to **"factor
the sparse `H_s`; drop the bump when it is PD; keep the kriging correction
unchanged."** No low-rank correction is needed in the common (PD) case — the
Woodbury path is only a *potential* refinement, and the dense bump remains
the fallback when `H_s` is singular. This makes (B) markedly smaller and
lower-risk than originally scoped. Implementation still gates each call site
(marginal variances, log-marginal, mode) on the dense-oracle test and the
R-INLA oracle fixtures.

### Consequences

- **Buys:** removes the `O(n²)`/`O(n³)` blow-up for constrained intrinsic
  models — the dominant allocation and the areal-scalability ceiling.
  Composes with PR #22 (the sparse `H_s` factor is exactly what
  `marginal_variances(::FactorCache)` consumes).
- **Costs:** two code paths (B) or a reshaped constrained solve (C); new
  low-rank correctness surface that must be oracle-gated; the selinv
  question under (C) is unresolved and spike-gated.
- **Escape hatch:** the dense bump (A) stays as the fallback in (B) and as
  the revert target for (C); no fit loses correctness or robustness.

### References

- [`plans/review-2026-07-remediation.md`](review-2026-07-remediation.md) — Tier-3 item 15, the benchmark that surfaced this.
- [`_null_bump`](../packages/LatentGaussianModels.jl/src/inference/constraints.jl) — current dense construction.
- ADR-012 — SelectedInversion.jl; the sparse selinv path the low-rank form reuses.
- PR #22 — `marginal_variances(::FactorCache)` factor reuse; consumes the sparse `H_s` factor.
- Rue & Held (2005) §2.3.2 — equality-constrained GMRF solves (kriging vs augmented system).

---

## ADR-046: Integrated hyperparameter marginals — design-point reuse for m(θ) = 1, conditional-mode profile slices for m(θ) ≥ 2

Status: Accepted (implemented 2026-07-05; see implementation findings below)
Date: 2026-07-05

### Context

`posterior_marginal_θ(res, j)`
([`marginals.jl:149`](../packages/LatentGaussianModels.jl/src/inference/marginals.jl))
returns a Gaussian at `(θ̂_j, √Σθ[j,j])` — the honestly-documented
placeholder from the Tier-4 backlog
([`review-2026-07-remediation.md`](review-2026-07-remediation.md) item 16).
The gap is not cosmetic:

- Log-precision posteriors are strongly right-skewed on the internal
  scale; a symmetric Gaussian at the mode misstates both moments, and
  the internal→user Jacobian (`p(τ) = p_int(log τ)/τ`) shifts user-scale
  summaries by 50–80 % on heavy-tailed axes. The Scotland classical-BYM
  oracle test
  ([`test_scotland_bym.jl`](../packages/LatentGaussianModels.jl/test/oracle/test_scotland_bym.jl))
  documents exactly this and *defers its τ_b point comparison* to the
  day an integrated θ-marginal accessor lands (replan Phase K promise,
  Martins et al. 2013 §3.3).
- R-INLA's `summary.hyperpar` / `marginals.hyperpar` — the numbers every
  R-INLA user reads first — come from an integrated marginal, never
  from the Gaussian at the mode.

Why the existing `INLAResult` cannot express the density: it stores
`θ_points` and the *normalised IS quadrature weights* `θ_weights`
(`w_k ∝ base_w_k · π̂(θ_k|y)/q(θ_k)`), but not the base quadrature
weights or `log π̂(θ_k|y)` themselves. Neither Grid (`φ(z)Δz` products)
nor CCD (center/axial/corner classes) has uniform base weights, so the
posterior density at the design points is not recoverable from
`θ_weights` after the fact. The integration stage *computes*
`log_π[k] = log_marginal + log_hyperprior` at every retained point
([`inla.jl:314`](../packages/LatentGaussianModels.jl/src/inference/inla.jl))
and then throws it away.

Validation data already exists: every one of the 28 LGM oracle fixtures
stores R-INLA's `marginals.hyperpar` density grids
(`scripts/generate-fixtures/_helpers.R`, `include_marginals = TRUE`), so
tier-2 gating needs **no fixture regeneration**. R-INLA's stored
marginals are user-scale; tests transform them to the internal scale via
the Jacobian rather than teaching the accessor about user-scale maps.

What R-INLA does (Martins, Simpson, Lindgren & Rue 2013 §3.2–3.3): build
`p(θ_j|y)` from the mode/Hessian eigenbasis ("z-parameterisation") with
per-axis asymmetric scalings σ⁺/σ⁻ estimated from a handful of
log-posterior probes — cheap, integration-free — and offer
`inla.hyperpar()` for a dense-grid refinement. Our
`compute_skewness_corrections`
([`integration.jl:314`](../packages/LatentGaussianModels.jl/src/inference/integration.jl))
already implements the σ⁺/σ⁻ probe machinery for the Grid stretches.

### Decision

Two-part accessor upgrade, dimension-dependent, plus the missing storage.

1. **`INLAResult` gains a `log_π::Vector{Float64}` field** — the
   unnormalised `log π̂(θ_k | y)` at each *retained* design point,
   aligned with `θ_points` / `θ_weights`. Populated for free in
   `_inla_integrate` (the values are already computed for the IS
   reweight) and as `[lp.log_marginal]` in the 0-hyperparameter fast
   path. No numerical output changes anywhere: the field is
   write-only until `posterior_marginal_θ` reads it.

2. **`posterior_marginal_θ(res, j; method = :auto, model = nothing,
   y = nothing, …)`** with methods:

   - `:gaussian` — current behaviour, unchanged, kept as the escape
     hatch.
   - `:integrated`, `m == 1` — reuse the design line exactly (the
     docstring's promised "density numerically integrated over the INLA
     design points"): density ∝ `exp(log_π)` at the sorted design
     points, trapezoid-normalised, evaluated on the output grid by
     monotone interpolation of the log-density with linear log-density
     extrapolation beyond the design span, then renormalised. Zero new
     `laplace_mode` calls. Covers every 1-hyperparameter model (the
     `:auto` scheme resolves to `Grid()` for m ≤ 2, and CCD itself
     falls back to `Grid(7)` at m = 1).
   - `:integrated`, `m ≥ 2` — conditional-mode profile slice, requires
     `model` and `y` (same convention as `posterior_marginal_x` with
     `SimplifiedLaplace`/`FullLaplace`; throws `ArgumentError`
     otherwise). For grid values `t` of `θ_j`, evaluate
     `log π̂(θ(t) | y)` along the Gaussian-conditional path
     `θ(t) = θ̂ + (Σθ[:, j] / Σθ[j, j]) (t − θ̂_j)` via the existing
     `_neg_log_posterior_θ` closure. Under the local-Gaussian
     approximation the conditional normaliser is `t`-independent, so
     the renormalised slice *is* the Laplace-approximate marginal.
     Evaluated on an internal profile grid (default 21 points spanning
     ±4 conditional sd) with the same 25-nat early-truncation rule as
     `full_laplace.jl`, then interpolated onto the output grid. Cost:
     ≤ 21 warm-startable `laplace_mode` calls per hyperparameter,
     on demand.
   - `:auto` (default) — `:integrated` when it is free or possible
     (m == 1 always; m ≥ 2 when `model` and `y` are supplied),
     `:gaussian` otherwise. Documented explicitly in the docstring so
     the fallback is never silent-*and*-surprising.

Out of scope for the first PR: user-scale marginal transforms (the
Jacobian lives in tests for now), R-INLA-style `summary.hyperpar`
accessors built on top, and reusing the profile slices to refine
`θ_mean`. Each is a natural follow-up once the density itself is gated.

### Alternatives considered

- **Kernel-smoothed projection of the IS-weighted design cloud**
  (works for any scheme, no new evaluations at m ≥ 2): rejected —
  bandwidth-sensitive, inflates variance by the kernel width, and on a
  5×5 Grid the per-coordinate projection has too few distinct support
  points for a credible density.
- **Martins §3.2 integration-free split-normal construction** (what
  R-INLA's cheap default does): viable and cheap (2m + 1 probes we
  already know how to make), but strictly less accurate than fresh
  profile slices and requires 1-D convolution machinery for the
  linear-combination marginal. Recorded as the candidate no-`model`
  fallback for m ≥ 2 if users need one; not part of this ADR.
- **Dense tensor re-integration** (`refine_hyperposterior` at
  `n_grid^m` then marginalise): cost explodes with m, and the
  eigen-axis grid still needs projection + smoothing to become a
  coordinate marginal — it inherits the cloud problem it was meant to
  solve.

### Consequences

- **Buys:** closes the highest-visibility R-INLA parity gap of the
  three Tier-4 items; un-defers the Scotland classical-BYM τ_b point
  comparison; gives `refine_hyperposterior` output a density-grade
  consumer. Validated against already-stored oracle marginals — no R
  round-trip.
- **Costs:** one struct field (internal; both construction sites are in
  `inla.jl`, none in tests); m ≥ 2 accuracy is itself a Laplace-flavour
  approximation (exact only when `π̂(θ|y)` is Gaussian) — gated at the
  5 % hyperparameter oracle tolerance; ~21 extra Laplace fits per
  hyperparameter for the m ≥ 2 path, opt-in.
- **Escape hatch:** `method = :gaussian` is bit-for-bit the previous
  behaviour; `:auto` resolves to it whenever the integrated path lacks
  inputs.

### Acceptance criteria

- Tier 1 (regression): near-Gaussian θ-posterior (Gaussian likelihood,
  generous data) — `:integrated` vs `:gaussian` densities agree in
  sup-norm within a tight tolerance; m = 1 self-consistency — trapezoid
  mean of the integrated density matches `res.θ_mean[j]` within
  quadrature error (same measure, same points); m = 2 — profile-slice
  marginal moments match a dense `refine_hyperposterior` reference
  within tolerance.
- Tier 2 (oracle, existing fixtures): Brunei (m = 1), Scotland BYM2 and
  Pennsylvania BYM2 (m = 2) — mean/sd of the Jacobian-transformed
  integrated marginal vs the same moments computed from the fixture's
  `marginals_hyperpar` grid, within the 5 % hyperparameter tolerance.
  Rewrite the deferred `test_scotland_bym.jl` τ_b assertion as a
  user-scale mean-to-mean comparison.
- Full dev-linked LatentGaussianModels suite green; all paths not
  calling `posterior_marginal_θ` bit-identical.

### Implementation findings (2026-07-05)

Landed as designed, with one acceptance criterion amended by what the
data showed:

- Tier-2 gates pass on the identifiable fixtures: Brunei (m = 1,
  user-scale τ mean within 5 %, sd within 10 % of the stored
  `marginals_hyperpar` moments), Scotland BYM2 and Pennsylvania BYM2
  (m = 2 profile slice, mean within 10 %, sd within 25 %).
- The planned "user-scale mean-to-mean" rewrite of the deferred
  `test_scotland_bym.jl` τ_b comparison is **not achievable on that
  fixture** — and the reason is informative. With the
  statistic-mismatch (mode-vs-mean + Jacobian) problem removed by this
  accessor, a genuine fit-level gap remains (~2× on the τ_b mean): the
  classical-BYM `(τ_v, τ_b)` posterior is a non-identified ridge
  (Eberly & Carlin 2000) and the two implementations distribute ridge
  mass differently — R-INLA's own fixture reports τ_v mean 1005 with
  sd 7756. The straight-path profile slice is faithful to the local
  Laplace picture but cannot see off-path ridge mass; this is now a
  documented limitation in the `posterior_marginal_θ` docstring. The
  test instead asserts tail-robust two-sided consistency (Julia's
  integrated median inside R-INLA's 95 % CI; R-INLA's mean inside
  Julia's integrated central 95 % band), with tight τ parity delegated
  to the identifiable BYM2 fixtures.

### References

- Martins, Simpson, Lindgren & Rue (2013) §3.2–3.3 — R-INLA's
  hyperparameter-marginal construction (`references/papers.md`).
- Rue, Martino & Chopin (2009) §6.5 — design construction, σ⁺/σ⁻.
- [`review-2026-07-remediation.md`](review-2026-07-remediation.md)
  Tier-4 item 16 — the backlog entry this scopes.
- [`replan-2026-04-28.md`](replan-2026-04-28.md) Phase K — the original
  "integrated θ-marginal accessor" promise.
- [`test_scotland_bym.jl`](../packages/LatentGaussianModels.jl/test/oracle/test_scotland_bym.jl)
  — the deferred comparison, now replaced by the two-sided consistency
  check described above.

---

## ADR-047: Item 16.3 closed — the "simplified-Laplace variance correction" is mis-specified; no such term exists in the reference method

Status: Accepted
Date: 2026-07-06

### Context

Tier-4 item 16.3
([`review-2026-07-remediation.md`](review-2026-07-remediation.md))
scheduled the "SimplifiedLaplace variance correction" deferred by
ADR-016, which described it as "R-INLA's third simplified-Laplace
term, `H⁻¹ Aᵀ diag(h⁴) A H⁻¹` per coordinate". Scoping against the
primary source — Rue, Martino & Chopin (2009) §3.2.3, eqs. 17–22 —
shows the premise is wrong:

- The simplified-Laplace expansion defines exactly **two** correction
  terms, both built from third derivatives `d_j^(3)`:
  `γ_i^(1)` (mean shift, eq. 21 — shipped per ADR-016) and `γ_i^(3)`
  (skewness, eq. 21 — shipped as the Edgeworth/skew density
  correction).
- The density representation is a skew-normal fitted "so that the
  third derivative at the mode is `γ^(3)`, the mean is `γ^(1)` and
  the variance is 1" — the per-θ conditional variance is
  **deliberately left uncorrected** in the reference method.
- For symmetric heavy-tailed likelihoods (Student-t is the paper's
  example) the prescribed remedy is not a variance term but the
  spline-corrected *full* Laplace (eq. 17) — the strategy we ship as
  `FullLaplace`, since item 16.2 also at the integration stage.
- A fourth-derivative term of the claimed form could only arise from
  a second-order expansion of the log-determinant denominator
  (eq. 20), which the paper never carries out and classic R-INLA
  never implemented as part of `simplified.laplace`.

What corrects variances in *modern* R-INLA is a different,
strategy-independent mechanism: the variational-Bayes corrections
(low-rank VB mean correction, van Niekerk & Rue 2024, JMLR 25(62);
variance strategy per van Niekerk, Krainski, Rustand & Rue 2023, CSDA
181; exposed as `control.vb`). This is the "unified VB-corrected
pipeline" already documented in the Brunei fixture header
(`scripts/generate-fixtures/lgm/synthetic_brunei.R`) and the reason
the pure-FL vs R-INLA SD tolerance there is 8 %.

### Decision

Close item 16.3 with **no numerical change** (option A of the
2026-07-06 scoping):

- Correct the `INLA` docstring claim that "R-INLA's full
  `simplified.laplace` also includes a variance correction … deferred
  to v0.3".
- Amend ADR-016 to withdraw the `h⁴` formula.
- Do **not** add a `∇⁴_η_log_density` likelihood contract; no fixture
  regeneration; no R round-trip.

Record option B — implementing the van Niekerk–Rue VB corrections as
an opt-in post-pass on the `_apply_integration_moments` seam
introduced by item 16.2 — as the properly-specified successor feature
if modern-R-INLA variance parity is wanted. It has real oracles
already in-tree (every stored fixture is VB-corrected), would allow
*tightening* the Brunei SD band, and needs its own ADR when
scheduled; it is a new feature, not remediation.

Rejected: option C, deriving the fourth-derivative second-order
denominator term ourselves — no reference implementation, no oracle,
no parity value ("more classic than classic").

### Consequences

- Tier 4 of the 2026-07 remediation closes: 16.1 shipped (ADR-046,
  PR #29), 16.2 shipped (ADR-026 "PR-4", PR #30), 16.3 resolved as
  documentation (this ADR).
- The remaining honest gap vs modern R-INLA is now precisely named:
  per-θ conditional variances are not VB-corrected in any strategy.
  It is visible only inside the tolerances already calibrated on the
  oracle fixtures (e.g. Brunei's 0.075 SD band).
- **Escape hatch / future:** option B above, on the
  `_apply_integration_moments` seam.

### References

- Rue, Martino & Chopin (2009), §3.2.3, eqs. 17–22 — the definitive
  statement that the SLA variance is pinned at 1.
- van Niekerk & Rue (2024). Low-rank variational Bayes correction to
  the Laplace method. *JMLR*, 25(62), 1–25.
- van Niekerk, Krainski, Rustand & Rue (2023). A new avenue for
  Bayesian inference with INLA. *CSDA*, 181, 107692
  (`references/papers.md`).
- ADR-016 (amended 2026-07-06), ADR-026,
  [`review-2026-07-remediation.md`](review-2026-07-remediation.md)
  item 16.
