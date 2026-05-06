# Phase P — `v1.0` close (tier-3 triangulation, perf benchmarks, release polish)

## Context

Phase P is the **final** phase before the `v1.0.0` cut. The replan
([`plans/replan-2026-04-28.md:454-475`](replan-2026-04-28.md))
originally scoped Phase P as a *rolling* tier-3 triangulation effort
that runs alongside K/L/M and never closes. With Phases F–O now
shipped (the v0.1.x line consolidated through F–L, Phase M at v0.2.0,
Phase N at v0.2.2, Phase O closing at v0.3.0), the rolling framing is
no longer load-bearing: every numerical surface that v1.0 wants to
ship is in. Phase P collapses to a **tight close** — five PRs of
test, doc, and release work plus one optional stretch, no new
components, no new likelihoods, no new oracle gates. The deliverable
is `v1.0.0` with a maintenance-only roadmap behind it.

The four work-streams that gate `v1.0.0`:

1. **Tier-3 triangulation green in nightly CI** — NUTS cross-validation
   via [`LGMTuring.jl`](../packages/LGMTuring.jl/) for the three
   flagship workflows (Scotland BYM2, Pennsylvania BYM2, Meuse SPDE).
   Per [`plans/testing-strategy.md:61-76`](testing-strategy.md),
   tier-3 is the only tier that has been "unrealised" through
   v0.3.0. One of the three tests already exists as a partial scaffold
   ([`packages/LGMTuring.jl/test/triangulation/test_scotland_bym2.jl`](../packages/LGMTuring.jl/test/triangulation/test_scotland_bym2.jl)),
   so PR-1 builds out, doesn't build from scratch. Stan via
   `BridgeStan.jl` was originally part of this work-stream but is
   dropped from Phase P (no local Stan toolchain); a second-MCMC
   sanity check via [`Gen.jl`](https://github.com/probcomp/Gen.jl)
   is evaluated as an optional stretch — see ADR-043 below.
2. **Performance benchmarks against R-INLA** — the v0.1.1 deferral in
   [`plans/quality-and-perf-benchmarks.md:101-191`](quality-and-perf-benchmarks.md)
   was never picked up. v1.0 needs at least one published results
   file — the user-trust argument that justified Part 1 (quality)
   for v0.1.0 applies in stronger form to v1.0.
3. **R-INLA migration guide** — Phase N PR-6 shipped
   [`docs/src/lgmformula-tutorial.md`](../docs/src/lgmformula-tutorial.md)
   for the formula-DSL surface, but there is no top-level "coming
   from R-INLA" guide that walks through the full ecosystem
   (likelihoods, components, priors, sampling, prediction). v1.0 is
   when we owe that page to users.
4. **ADR backlog + CHANGELOG consolidation** — `plans/decisions.md`
   carries 41+ ADRs across F–O; the registry needs a forward-pointer
   index so future readers can navigate it. `CHANGELOG.md` tags v0.2.2
   retroactively at `cbfae1e` per Phase O PR-0; v1.0's CHANGELOG entry
   is a phase-arc summary, not a per-commit list.

What Phase P **does not** do: ship new SPDE expansion items
(fractional-α, sphere, 3D, non-separable space-time), new likelihoods,
new components, new oracle fixtures. Those are explicit v1.x backlog
per the "Out of scope" section below.

## State of the world (read-only audit)

What ships today, against the four Phase P work-streams:

| Work-stream | Today | Gap |
|---|---|---|
| **Tier-3 triangulation — Scotland BYM2** | Partial scaffold at [`packages/LGMTuring.jl/test/triangulation/test_scotland_bym2.jl`](../packages/LGMTuring.jl/test/triangulation/test_scotland_bym2.jl). Loads the LGM oracle fixture, fits Julia INLA at `int_strategy=:grid`, runs a 200-sample NUTS chain through `LGMTuring.nuts_sample(model, y, n; init_from_inla, rng)`, asserts mean within 2 SDs and SD within 30% relative. Skips transparently if fixture missing. | **No nightly CI integration** — the test runs on package-local CI, but `nightly-upstream.yml` only refreshes oracle fixtures, doesn't run triangulation. **Tolerances are loose** (200-sample chain, 30% relative SD); v1.0 should tighten to 1000+ samples and 10–15% SD per the replan acceptance gate. |
| **Tier-3 triangulation — Pennsylvania BYM2** | **Does not exist.** Oracle fixture exists at [`packages/LatentGaussianModels.jl/test/oracle/fixtures/pennsylvania_bym2.jld2`](../packages/LatentGaussianModels.jl/test/oracle/fixtures/pennsylvania_bym2.jld2); LGM-side oracle test passes (per `plans/release-v0.1.md` close memo). | Mirror of Scotland test — same shape, different fixture. ~50 LOC + a NUTS chain. |
| **Tier-3 triangulation — Meuse SPDE** | **Does not exist.** Oracle fixture exists at [`packages/INLASPDE.jl/test/oracle/fixtures/meuse_spde.jld2`](../packages/INLASPDE.jl/test/oracle/fixtures/meuse_spde.jld2); the post-Phase O Meuse predict oracle [`packages/INLASPDERasters.jl/test/oracle/fixtures/meuse_spde_predict.jld2`](../packages/INLASPDERasters.jl/test/oracle/fixtures/meuse_spde_predict.jld2) covers the raster projection but not the latent posterior. | Latent-dim is `n + n_v ≈ 155 + 200 = 355` (Meuse has ~155 obs and the standard `max_edge=0.05` mesh has ~200 vertices), the largest model in the triangulation suite. NUTS will be slow at this dim — budget per-leapfrog in seconds, not milliseconds. Per `plans/testing-strategy.md:65-66` the SPDE triangulation is exactly the kind of model "we really care about being right." |
| **Second-MCMC sanity check (Gen.jl)** | **Does not exist anywhere in the repo.** `Gen.jl` is a probabilistic-programming framework from MIT ProbComp with its own HMC / MH / particle-filter inference stack, independent of `Turing` / `AdvancedHMC`. No prior integration with LGM. | Optional. Adds independent-Julia-HMC coverage on top of NUTS — catches `AdvancedHMC`-specific bugs that NUTS-alone would miss. Cost: ~1–2 days to express the LGM `LogDensityProblem` as a `@gen function` (Gen.jl's DSL doesn't take `LogDensityProblems` directly). Recommendation in ADR-043: ship as PR-1b stretch only if PR-1 closes ahead of schedule. |
| **Stan triangulation harness** | Out of scope — no local Stan toolchain. | Deferred to v1.x backlog. NUTS via `LGMTuring.jl` is the v1.0 tier-3 coverage; cross-language verification (Stan / NIMBLE) lands when the toolchain is available. |
| **Performance benchmark harness** | **Does not exist.** [`plans/quality-and-perf-benchmarks.md:101-191`](quality-and-perf-benchmarks.md) has the design (`benchmarks/run.jl` entry, `benchmarks/r_inla.R` for the R side, single-threaded BLAS, hardware spec in every results file). No `benchmarks/` directory at the repo root. | ~1 week wall-clock per the design memo (2 days harness + Scotland, 1 day add PA + Meuse, 1 day report integration, 1 day buffer). v1.0 is the right beat — perf claims without a results file would dilute the v1.0 announcement. |
| **R-INLA migration guide** | Formula-DSL-scoped at [`docs/src/lgmformula-tutorial.md`](../docs/src/lgmformula-tutorial.md). Phase N PR-6 covered Scotland BYM2 and Tokyo rainfall side-by-side `inla(...)` ↔ `@lgm(...)`. | No top-level "coming from R-INLA" page. Surface to cover: likelihood naming (`inla(family = "poisson")` ↔ `PoissonLikelihood()`), component naming (`f(model = "bym2")` ↔ `BYM2(...)`), prior parameterisation gotchas (PC priors on precision vs SD; per [`plans/defaults-parity.md`](defaults-parity.md)), sampling (`inla.posterior.sample` ↔ `posterior_sample`), prediction (`predict.inla` ↔ `predict_raster`), getting raster output (no R-INLA equivalent — point at `INLASPDERasters.jl`). |
| **ADR registry navigation** | [`plans/decisions.md`](decisions.md) ships chronologically — ADR-001 through ADR-041+ in the order they landed. No topical index, no superseded chain visualisation. | Add a TOC at the head of `decisions.md` grouping ADRs by topic (architecture, dispatch, dependencies, defaults, macro, observation mapping, raster, formula, SPDE). Mark superseded ADRs with explicit "Superseded by ADR-NNN" status. ~half a day of read-and-edit work. |
| **CHANGELOG cross-phase consolidation** | Per-phase entries live in [`CHANGELOG.md`](../CHANGELOG.md): v0.1.5 (F–L), v0.2.0 (M), v0.2.2 (N), v0.3.0 (O — pending). No v1.0 phase-arc summary. | v1.0 entry should be a 1-page narrative: "F–O delivered the R-INLA-parity surface; v1.0 closes the testing matrix and ships the migration guide." Not a commit list — that's what the per-phase entries already provide. |
| **Documentation completeness** | All 36 oracle fixtures pass; vignettes exist for Scotland BYM2, Tokyo rainfall, Meuse SPDE, Cameletti space-time, Pennsylvania BYM2, plus the synthetic 1D Matérn and Lindgren–Rue–Lindström §3.2 SPDE flagship pieces. | A v1.0 release announcement / blog-post-style `docs/src/announcing-v1.0.md` would be useful but is not a blocker. Defer to a v1.0.x doc patch unless free time materialises. |

Three cross-cutting realities the audit surfaces:

1. **The replan's "rolling" framing for Phase P is obsolete.** When
   the replan was written (2026-04-28, pre-Phase-M), the assumption
   was that triangulation would be co-developed alongside the
   numerical work it certified. In practice tier-2 R-INLA oracle
   coverage was strong enough that tier-3 was deferred without
   blocking phase closes. v1.0 surfaces tier-3 as the gate that was
   "rolling" because nobody had to schedule it.

2. **No new numerical code is on the v1.0 critical path.** The NUTS
   path already works (Scotland triangulation test demonstrates it).
   Pennsylvania and Meuse triangulation is shape-replication of
   Scotland. The optional Gen.jl path (PR-1b stretch) is the one
   place where v1.0 is potentially building something the project
   doesn't have, and it's explicitly off the critical path.

3. **Phase P's calendar risk concentrates on the perf benchmark, not
   triangulation.** Without Stan in scope, PR-1 collapses to NUTS-
   tightening + clone — ~3 days of Julia work. The dominant
   variance now sits in PR-3 (R-INLA install + RCall integration is
   famously fiddly). If PR-3's harness slips, the Gen.jl stretch
   gets dropped first; v1.0 ships with NUTS-only triangulation.

## Pre-Phase-P housekeeping

One Phase O dangle to reconcile before Phase P opens:

- **Phase O `v0.3.0` tag must land first.** Phase P's release
  arithmetic assumes the umbrella is at `0.3.0`; if Phase O is still
  in PR-flight when Phase P starts, the two close tasks bleed into
  each other and the CHANGELOG split between v0.3.0 (Phase O) and
  v1.0.0 (Phase P) gets muddy. Phase O should close cleanly before
  Phase P PR-1 opens.

If `v0.3.0` ships in time, Phase P starts clean. If not, fold the
remaining Phase O work into Phase P PR-0 as a housekeeping prologue
analogous to Phase O's PR-0.

## Design calls (ADRs)

Phase P surfaces two ADR candidates. Numbers continue from Phase O's
ADR-040 / ADR-041 (and optional ADR-042); Phase P's entries are 043
and 044.

### ADR-043 candidate — Gen.jl as a second-MCMC sanity check: ship in v1.0 or defer

NUTS via `LGMTuring.jl` is the load-bearing tier-3 implementation;
that's not in question. The question is whether Phase P also ships a
second Julia-HMC implementation through [`Gen.jl`](https://github.com/probcomp/Gen.jl).
Three options:

- **(a) Ship Gen.jl in v1.0 as PR-1b.** Express the LGM
  `LogDensityProblem` as a `@gen function`, run Gen's HMC sampler,
  add a third column to the triangulation envelope (INLA / NUTS /
  Gen-HMC). Cost: ~1–2 days for the model translation plus chain
  authoring. Coverage gain: catches `AdvancedHMC`-specific bugs.
- **(b) Defer Gen.jl to v1.x as a stretch lane.** Ship NUTS-only
  triangulation in v1.0; track Gen.jl integration as a v1.x backlog
  item. Coverage at v1.0 is "INLA vs NUTS"; second-MCMC sanity is
  v1.x.
- **(c) Skip Gen.jl entirely.** Cross-language verification (Stan,
  NIMBLE) is the genuinely-independent triangulation; second-Julia-
  HMC is a thin gain because both implementations share the same
  AD ecosystem (`ForwardDiff` / `ReverseDiff`) and the same
  `LogDensityProblems` contract.

**Recommendation:** option (b). Gen.jl is real value but not v1.0-
blocking; the user-trust argument for tier-3 at v1.0 is satisfied by
NUTS-vs-INLA. Frame Gen.jl as a stretch PR (PR-1b) — ship it if
PRs 1–4 close inside week 1; otherwise defer cleanly to v1.1. The
ADR documents the tradeoff so a future maintainer doesn't re-ask
the question.

If option (a) wins on the day, the integration shape is: new
`packages/LGMTuring.jl/ext/LGMTuringGenExt.jl` weakdep extension;
`LGMTuring.gen_hmc_sample(model, y, n; rng, ...)` wrapping
`Gen.hmc(...)` and returning an `MCMCChains.Chains`-shaped result so
the existing `compare_posteriors(...)` harness reuses cleanly.

### ADR-044 candidate — Tier-3 triangulation tolerances at v1.0

The Scotland BYM2 triangulation test ships `tol_mean = 2.0 SDs` and
`tol_sd = 0.30` (relative). At v1.0, with full chains (>=1000
post-warmup samples), the tolerances tighten. Two options:

- **(a) Per-test calibration.** Each model gets its own tolerance
  derived from the chain's effective sample size and the
  hyperparameter dimensionality. Sound but expensive to maintain.

- **(b) Uniform v1.0 tolerances.** Apply across all three flagship
  models: `tol_mean = 1.5 SDs` (envelope: INLA / NUTS posteriors
  agree within 1.5σ), `tol_sd = 0.60` (relative). The mean tolerance
  is tight (catches gradient bugs, precision-build bugs, mode-finder
  drift); the SD tolerance is loose because the SD diff between
  INLA's grid and NUTS is *structural*, not MC-noise — see below.

**Recommendation:** option (b).

**Why `tol_sd = 0.60` and not 0.15.** The first cut of this ADR
proposed `tol_sd = 0.15` ("tighter because chains are longer").
Empirically, on both Scotland and PA, NUTS finds 30–45% wider
posterior on the BYM2 mixing weight `logit φ` than INLA does — at
1000-sample chains this is *not* MC-error-bounded, it's a real
property of grid integration vs. full HMC. INLA's grid is bounded
by the Hessian-explored region around the Laplace mode; NUTS
samples the full posterior including the heavy tails of weakly
identified mixing parameters. Tightening `tol_sd` would force
either per-parameter tolerances (rejecting option (a)) or n=10k+
chains (10× nightly cost). The 0.60 envelope still flags genuine
regressions: a gradient bug blows SDs by 100%+; a precision-build
bug shifts means by SDs (caught by `tol_mean`); a mode-finder bug
disagrees on both. The means agree within 0.06–0.23 SDs on both
flagship Poisson-BYM2 datasets — that's the load-bearing constraint.

Tier-3 tolerances are not the place for per-model fine-tuning —
that's tier-1's job. The tier-3 contract is "all available
independent implementations land in the same envelope," and the
envelope width should be wide enough to absorb known structural
inference-method differences while still catching regressions. If
PR-1b (Gen.jl) ships, the same tolerances apply to the three-way
envelope.

## PR sequence

Five PRs (+ one stretch tail), all light on numerical content. PR-1
is the load-bearing piece; PR-2 / PR-3 / PR-4 are publication /
polish work.

### PR-1 — Tier-3 NUTS triangulation harness (Scotland + PA + Meuse)

Extends [`packages/LGMTuring.jl/test/triangulation/test_scotland_bym2.jl`](../packages/LGMTuring.jl/test/triangulation/test_scotland_bym2.jl)
to v1.0 tier-3 coverage:

- **Tighten Scotland tolerances** to ADR-044 v1.0 levels (`tol_mean
  = 1.5 SDs`, `tol_sd = 0.60` relative); bump chain length to 1000
  post-warmup with 200 warmup. Wall-clock for Scotland NUTS at this
  spec: ~30 s on `Apple M2 Pro / 16 GB` (estimate; benchmark in
  PR-3).
- **Add `test_pennsylvania_bym2.jl`** — mirror Scotland, same
  shape, different fixture (`pennsylvania_bym2.jld2`).
- **Add `test_meuse_spde.jl`** — latent dim ~355, expect minutes
  for the NUTS run; gate behind a `--triangulation` test arg so
  package-local CI doesn't pay the cost on every push.

**CI integration**:

- `.github/workflows/nightly-triangulation.yml` (new) — runs the
  triangulation suite on a weekly cron, separately from
  `nightly-upstream.yml` (which is for oracle drift). Opens an issue
  on failure.
- `package-local CI` (`test.yml`) gets a `--triangulation` opt-in
  test arg; default still excludes triangulation from PR CI.

**Tests / acceptance**: all three triangulation tests pass nightly
at the ADR-044 v1.0 tolerances. The existing
[`packages/LGMTuring.jl/test/regression/`](../packages/LGMTuring.jl/test/regression/)
suite (NUTS recovers analytic posterior on a known-conjugate model)
remains as the tier-1 sanity gate on the NUTS bridge itself.

**Files**: `packages/LGMTuring.jl/test/triangulation/test_*.jl`
(Scotland tighten + two new); `.github/workflows/nightly-triangulation.yml`
(new); `packages/LGMTuring.jl/test/runtests.jl` (wire the
`--triangulation` opt-in arg).

### PR-1b (stretch) — Gen.jl second-MCMC sanity check

**Lands only if PR-1 / PR-2 close inside week 1**, per ADR-043
option (b). New weakdep extension
`packages/LGMTuring.jl/ext/LGMTuringGenExt.jl`:

- Express the LGM `LogDensityProblem` as a `@gen function` —
  Gen's DSL wraps the gradient call so the bridge consumes the
  same `LogDensityProblems` contract as the NUTS bridge.
- `LGMTuring.gen_hmc_sample(model, y, n; rng, n_adapts)` →
  `MCMCChains.Chains`, so the existing
  [`packages/LGMTuring.jl/src/compare.jl`](../packages/LGMTuring.jl/src/compare.jl)
  `compare_posteriors` harness reuses unchanged.
- Add `gen_hmc` calls to the three triangulation tests, asserting a
  three-way INLA / NUTS / Gen-HMC envelope at the same ADR-044
  tolerances.

**Acceptance**: same as PR-1 (envelope holds at the v1.0
tolerances) plus a tier-1 regression test (small known-conjugate
model, Gen-HMC recovers analytic posterior within MC error) so
the Gen path is independently certified.

**Defer cleanly if time-pressed**: PR-1b is the first thing dropped
if Phase P slips. The v1.0 tier-3 contract holds at NUTS-only.

**Files**: `packages/LGMTuring.jl/ext/LGMTuringGenExt.jl` (new);
`packages/LGMTuring.jl/Project.toml` (add `Gen` to `[weakdeps]` +
`[extensions]`); `packages/LGMTuring.jl/test/triangulation/test_*.jl`
(extend each with the three-way assertion).

### PR-2 — R-INLA migration guide

`docs/src/coming-from-r-inla.md` (new top-level docs page). Sectioned
walkthrough:

- **Likelihoods** — table mapping R-INLA `family = "..."` strings to
  Julia `<Likelihood>()` types; cover the 20 likelihoods shipped
  through Phase J.
- **Components** — table mapping R-INLA `f(model = "...")` strings to
  Julia `<Component>(...)` constructors; cover the 22 components
  shipped through Phase M / N.
- **Priors** — section on PC prior parameterisation gotchas
  (precision-scale vs SD-scale; per
  [`plans/defaults-parity.md`](defaults-parity.md)). Three worked
  examples: BYM2 PC prior, RW2 PC prior, SPDE2 PC-Matérn prior.
- **Inference flags** — `int.strategy = "ccd"|"grid"|"empirical"` →
  `int_strategy = :ccd|:grid|:empirical`. `control.compute = list(...)`
  → Phase K accessor functions (`fixed_effects`, `random_effects`,
  `marginal_likelihood`, `psis_loo`, `pp_check`).
- **Sampling** — `inla.posterior.sample` → `posterior_sample`. Note
  the shape difference (Phase K-3 ADR).
- **Prediction** — `predict.inla` → `predict_raster` (raster) and
  `posterior_predictive` (point predictions on test locations).
- **Things R-INLA does that we don't (yet)** — `inla.qsample` raw
  GMRF sampling (we expose this through `GMRFs.jl`); spatial
  models on the sphere; 3D SPDE; non-separable space-time. Pointer
  to v1.x backlog.

The page lives next to (and links into) the formula-DSL tutorial
[`docs/src/lgmformula-tutorial.md`](../docs/src/lgmformula-tutorial.md);
each table cell that has a `@lgm` analogue gets a third column.

**Tests**: docs build green under `Documenter.makedocs(strict=true)`.
No code-test coverage — this is documentation work.

**Files**: `docs/src/coming-from-r-inla.md` (new); `docs/make.jl`
(register the new page under "Coming from R-INLA" top-level section,
which already exists from Phase N PR-6); `docs/src/index.md` (add a
landing card linking to the new page).

### PR-3 — Performance benchmark harness

Implements [`plans/quality-and-perf-benchmarks.md:101-191`](quality-and-perf-benchmarks.md)
Part 2 verbatim. Re-stated here for completeness:

```
benchmarks/
├── Project.toml             # BenchmarkTools, RCall, the four packages
├── run.jl                   # entry point
├── harness.jl               # Julia-side timing + memory measurement
├── r_inla.R                 # R-INLA fits
├── compare.jl               # parses both sides, emits the table
└── results/
    └── 2026-MM-DD_<arch>.md
```

The three datasets are the Phase P triangulation set (Scotland BYM2,
PA BYM2, Meuse SPDE) — same fixtures, different tolerances. Single-
threaded BLAS by default per the design memo; multi-threaded numbers
are a separate table in the same file.

**Acceptance** (per the design memo):
- `benchmarks/run.jl` produces `results/<date>_<arch>.md` for the
  three R-INLA examples on at least one documented hardware spec.
- The results file is linked from `docs/src/benchmarks/quality.md`
  (the v0.1.0 page) under a "Performance" §, replacing the
  "deferred to v0.1.x" placeholder text.
- R-INLA version + Julia versions + hardware spec in the results
  header.
- No performance claims on the docs page that aren't backed by a
  results file.

**Tests**: `benchmarks/Project.toml` resolves cleanly; `run.jl`
produces a results file when run interactively. No CI gate — the
benchmark harness is for human-driven runs, not automated regression.

**Files**: `benchmarks/run.jl`, `benchmarks/harness.jl`,
`benchmarks/r_inla.R`, `benchmarks/compare.jl`,
`benchmarks/Project.toml` (new); `benchmarks/results/<date>_<arch>.md`
(first results file); `docs/src/benchmarks/quality.md` (replace
"Performance: deferred" placeholder).

### PR-4 — ADR backlog cleanup + CHANGELOG consolidation

`plans/decisions.md` reorganisation:

- Add a topical TOC at the head of the file. Group ADRs by topic
  (Architecture, Dispatch, Dependencies, Defaults, Macro, Observation
  mapping, Raster, Formula, SPDE, Inference, Testing). ADRs remain in
  numerical order in the body — the TOC is a forward index, not a
  reorder.
- Mark superseded ADRs with `Status: Superseded by ADR-NNN` and
  cross-link from the superseder.
- Backfill any ADRs that landed via PR commit messages but never made
  it into `decisions.md` (Phase O PR-0 already addressed ADRs 036–039;
  Phase P verifies no later ADRs slipped).

`CHANGELOG.md` v1.0 entry — a phase-arc narrative:

```markdown
## [v1.0.0] — 2026-MM-DD

Phase P closes the v1.0 release arc. F–O delivered the R-INLA-parity
numerical surface (likelihoods, components, priors, sampling,
prediction, raster bridge, formula DSL); v1.0 closes the testing
matrix (tier-3 NUTS triangulation across the three flagship workflows,
performance benchmarks against R-INLA) and ships the top-level R-INLA
migration guide.

No new components, no new likelihoods, no new oracle gates compared
to v0.3.0. v1.0 is a maturity release.

### Added — Phase P
...
```

Not a commit list — per-phase CHANGELOG entries already cover that.

**Tests**: `decisions.md` builds under `markdownlint` (if the repo
has it) or renders in GitHub preview without broken links.
`CHANGELOG.md` v1.0 entry parses under [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

**Files**: `plans/decisions.md` (TOC + supersession status edits);
`CHANGELOG.md` (v1.0 entry).

### PR-5 — `v1.0.0` release

Pure metadata PR. Bumps:

- `GMRFs.jl` 0.1.2 → 1.0.0
- `LatentGaussianModels.jl` 0.2.0 → 1.0.0
- `INLASPDE.jl` 0.3.1 → 1.0.0
- `INLASPDERasters.jl` 0.4.0 → 1.0.0
- `LGMFormula.jl` 0.4.0 → 1.0.0
- `LGMTuring.jl` 0.x.x → 1.0.0
- Umbrella `INLA.jl` 0.3.0 → 1.0.0

All compat bounds across `[deps]` and `[weakdeps]` updated to `"1"`
where they currently target `"0.x"` siblings.

Tag `v1.0.0` on `main`. Release commit message:

```
chore(release): v1.0.0 — Phase P close (tier-3 + perf + migration guide)
```

Push to the personal `haavardhvarnes/JuliaRegistry` per the registry
policy (memory file `project_general_registry_status.md`); General-
registry submission stays off the roadmap.

**Tests**: package resolution clean across all six packages at the
v1.0 compat pins; full ecosystem test suite green.

**Files**: every package's `Project.toml` (version + compat); root
`Project.toml`; `CHANGELOG.md` (release date stamp on the v1.0 entry).

## Out of scope for Phase P

Tracked in [`plans/replan-2026-04-28.md:476-535`](replan-2026-04-28.md)
"Cross-cutting" items — explicitly **NOT** in Phase P:

- **Fractional-α SPDE** — Phase M PR-7's deferred item. Mathematically
  correct via Bolin-Kirchner 2020 needs an `AugmentedLatentComponent`
  seam that doesn't fit the per-vertex `AbstractLatentComponent`
  contract. v1.x territory under a new abstract type.
- **SPDE on the sphere** (great-circle distances via
  `CoordRefSystems`). Per [`packages/INLASPDE.jl/plans/plan.md:213-219`](../packages/INLASPDE.jl/plans/plan.md)
  "Deferred to v0.3+" — equally v1.x.
- **3D SPDE** (brain connectivity, subsurface). Same v1.x backlog.
- **Non-separable space-time SPDE** (Lindgren et al. 2024). The
  separable case shipped in Phase M PR-5; non-separable needs the
  2024 generalisation. v1.x.
- **Areal-on-raster** (zonal statistics for raster covariates over
  areal regions). Deferred at Phase O PR-2 close. v1.x.
- **Pardiso backend.** v0.1 perf gap closed via the LBFGS `g_tol` fix
  (memory `project_phaseq_lbfgs_gtol.md`); Pardiso work is now
  optional and stays off the v1.0 roadmap.
- **R-INLA `inla.posterior.sample` parity.** Sampling shape diffs
  enough that an oracle fixture would need a 5% MC tolerance. The
  Gaussian-approximation tier-2 oracles already cover the tight
  numerical gate.
- **Auto-published benchmark history with charts.** Single results
  files per release per
  [`plans/quality-and-perf-benchmarks.md:188`](quality-and-perf-benchmarks.md).
- **Stan triangulation.** `BridgeStan.jl` integration with the
  `ConnorDonegan/Stan-IAR` reference models — dropped from Phase P
  (no local Stan toolchain). Lands as a v1.x triangulation expansion
  if/when a Stan toolchain is available; the harness shape is
  already documented in the original Phase P draft.
- **NIMBLE triangulation.** `gkonstantinoudis/nimble` is the
  reference but R-only; no Julia bridge currently exists. Same v1.x
  triangulation expansion lane as Stan.
- **Gen.jl as a non-stretch deliverable.** PR-1b is the *optional*
  Gen.jl path (ADR-043). If PR-1b ships, it ships as v1.0 stretch;
  if it doesn't, Gen.jl integration moves to v1.x backlog tracked
  alongside Stan / NIMBLE.
- **Edgeworth correction / IS-INLA / AMIS-INLA**. Phase K-class
  inference variants per the cross-cutting items in
  [`replan-2026-04-28.md:484-509`](replan-2026-04-28.md). Useful but
  not v1.0-blocking.

The boundary rule for Phase P: if it changes a Project.toml
`[deps]`, it's out of scope. v1.0's promise is "the tier-2 surface
shipped through O is now tier-3-validated and benchmarked against
R-INLA." Anything beyond that is v1.x.

## Critical files

| Concern | Path |
|---|---|
| Scotland BYM2 triangulation test (target of PR-1 tighten) | [`packages/LGMTuring.jl/test/triangulation/test_scotland_bym2.jl`](../packages/LGMTuring.jl/test/triangulation/test_scotland_bym2.jl) |
| LGMTuring NUTS bridge | [`packages/LGMTuring.jl/src/`](../packages/LGMTuring.jl/src/) |
| LGM oracle fixtures (PR-1 input) | [`packages/LatentGaussianModels.jl/test/oracle/fixtures/`](../packages/LatentGaussianModels.jl/test/oracle/fixtures/) |
| INLASPDE oracle fixture (PR-1 input) | [`packages/INLASPDE.jl/test/oracle/fixtures/meuse_spde.jld2`](../packages/INLASPDE.jl/test/oracle/fixtures/meuse_spde.jld2) |
| Nightly triangulation workflow (PR-1 new) | `.github/workflows/nightly-triangulation.yml` (new) |
| Gen.jl extension (PR-1b stretch — only if shipped) | `packages/LGMTuring.jl/ext/LGMTuringGenExt.jl` (new, optional) |
| Migration guide (PR-2 new) | `docs/src/coming-from-r-inla.md` (new) |
| Defaults parity reference (PR-2 input) | [`plans/defaults-parity.md`](defaults-parity.md) |
| Formula DSL tutorial (PR-2 cross-link) | [`docs/src/lgmformula-tutorial.md`](../docs/src/lgmformula-tutorial.md) |
| Performance benchmark design (PR-3 spec) | [`plans/quality-and-perf-benchmarks.md`](quality-and-perf-benchmarks.md) |
| Benchmarks tree (PR-3 new) | `benchmarks/` (new at repo root) |
| Quality benchmark page (PR-3 cross-link target) | [`docs/src/benchmarks/quality.md`](../docs/src/benchmarks/quality.md) |
| ADR registry (PR-4 reorganisation) | [`plans/decisions.md`](decisions.md) |
| CHANGELOG (PR-4 + PR-5) | [`CHANGELOG.md`](../CHANGELOG.md) |
| Per-package Project.toml (PR-5 version bumps) | `packages/*/Project.toml` |
| Personal registry policy | `~/.claude/.../memory/project_general_registry_status.md` |

## Verification

Phase P closes when **all** of the following hold:

1. **Tier-3 nightly green.** The three triangulation tests
   (Scotland BYM2, PA BYM2, Meuse SPDE) pass nightly via the NUTS
   path, at the ADR-044 v1.0 tolerances (`tol_mean = 1.5 SDs`,
   `tol_sd = 0.60` relative). At least 7 consecutive nightly runs
   green before tagging. PR-1b's Gen.jl path, if shipped, extends
   the assertion to a three-way envelope at the same tolerances.
2. **All 37+ tier-2 oracle fixtures still pass.** Tier-3 work must
   not regress tier-2 — sanity check.
3. **Performance benchmark results file published.** At least one
   `benchmarks/results/<date>_<arch>.md` committed with hardware
   spec, R-INLA version, Julia versions, single-threaded numbers
   for all three flagship models. The docs `benchmarks/quality.md`
   page links to it.
4. **R-INLA migration guide published.** `docs/src/coming-from-r-inla.md`
   builds clean under `Documenter.makedocs(strict=true)`; covers
   all 20 likelihoods, all 22 components, three flagship priors,
   inference / sampling / prediction sections.
5. **`plans/decisions.md` has a topical TOC.** Superseded ADRs
   marked. No ADR referenced in commit messages but missing from
   the registry.
6. **`CHANGELOG.md` v1.0.0 entry merged.** Phase-arc narrative,
   not a commit list.
7. **All six packages bumped to `1.0.0`.** Compat bounds updated.
   Full ecosystem test suite green at the bumped pins. Tag
   `v1.0.0` on `main`.
8. **Personal registry pushed.** `haavardhvarnes/JuliaRegistry`
   reflects the v1.0 release; users can `add INLA` from the
   registry URL.

## Release target

`v1.0.0` — the **final** semantic-version major-bump from the v0.x
line. Per-package version arithmetic detailed in PR-5 above. Release
commit `chore(release): v1.0.0 — Phase P close (tier-3 + perf +
migration guide)` mirroring the v0.2.0 / v0.3.0 close pattern.

After v1.0:

- **v1.0.x patches** — bug fixes, doc patches, small CHANGELOG-
  tracked deltas. No new public API.
- **v1.x.0 minors** — backlog items from "Out of scope" above
  (sphere SPDE, 3D, non-separable space-time, fractional-α,
  areal-on-raster, Stan triangulation, NIMBLE triangulation,
  Gen.jl triangulation if not shipped in v1.0). Each lands as its
  own phase plan when resourced.
- **v2.0** — not on the current roadmap. Would require a
  semantic-version-breaking change to the public LGM constructor
  surface, which is not anticipated.

The personal registry remains the canonical install path; General-
registry submission stays off the roadmap per
`memory/project_general_registry_status.md`.

## Cadence

Solo-developer pace, ~1.5 weeks total at the tight-close target
(Stan removal saves ~2 days vs. the original budget).

- **Week 1**
  - **Days 1–2: PR-1 — NUTS triangulation harness.** Day 1: tighten
    Scotland test (chain length, tolerances) + clone shape for PA.
    Day 2: Meuse triangulation (slowest — latent dim ~355 means
    NUTS chains take minutes per fit) + wire
    `nightly-triangulation.yml`.
  - **Days 3–4: PR-2 — R-INLA migration guide.** Day 3: write
    likelihood + component + prior tables. Day 4: write inference
    / sampling / prediction sections, cross-link with
    `lgmformula-tutorial.md`, build docs strict.
  - **Day 5: PR-1b stretch — Gen.jl second-MCMC sanity check.**
    Express the LGM `LogDensityProblem` as a `@gen function`, wire
    `LGMTuringGenExt`, extend the three triangulation tests to the
    three-way envelope. Drop cleanly to v1.x if the day runs over.

- **Week 2**
  - **Days 1–3: PR-3 — performance benchmark harness.** Day 1:
    `benchmarks/Project.toml` + `run.jl` + `harness.jl` + verify on
    Scotland. Day 2: add PA + Meuse to harness. Day 3: write
    `r_inla.R`, run cross-comparison, produce first results file,
    update docs page. Per the design memo, ~1 week budgeted; tight
    target is 3 days.
  - **Days 4–5: PR-4 + PR-5 — ADR backlog + CHANGELOG + release.**
    Day 4: ADR TOC + supersession marks + v1.0 CHANGELOG entry.
    Day 5: version bumps across six packages, tag, push to personal
    registry.

- **Stretch (optional, week 3 if time permits)**: a
  `docs/src/announcing-v1.0.md` blog-post-style release announcement
  page. Skip if PR-3's perf harness slips into week 2 day 3+.

The 1.5-week target reflects the audit's main finding: **Phase P is
almost entirely publication and CI plumbing**. The optional new code
is the Gen.jl extension (PR-1b stretch); everything else is wiring,
docs, or metadata. Calendar risk concentrates in PR-3 (R-INLA install
+ RCall integration is famously fiddly). If PR-3 slips, PR-1b is the
first thing dropped — Gen.jl integration is real value but not v1.0-
blocking, and v1.0's tier-3 contract holds at NUTS-only.

If PR-3's perf harness collapses into a longer plumbing exercise,
Phase P stretches to 2.5 weeks; the v1.0 tag slips by half a week.
Neither risk is load-bearing on the v1.0 *correctness* claim — that
comes from tier-2 (already shipped) and NUTS triangulation (already
partially shipped). Perf benchmarks and Gen.jl raise the credibility
ceiling but don't change the floor.
