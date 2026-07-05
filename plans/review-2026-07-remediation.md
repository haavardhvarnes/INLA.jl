# Review remediation plan — 2026-07

> Source: independent code review (2026-07-05) cross-checked against an
> external Codex review of an earlier tree. All findings re-verified against
> `main` at HEAD `85eba20`. This document is the actionable remediation
> backlog; ADR-worthy decisions arising from it get their own entry in
> `decisions.md`.

## Context

Two independent reviews (Codex, and a fresh in-repo pass) converged on the
same conclusion: the architecture and numerics are strong, but a cluster of
**release-polish defects** actively contradict the shipped `v1.0.0` claim,
and one **error-handling anti-pattern** is genuinely dangerous. Codex's
concrete claims all reproduced on the current tree; the independent pass
widened one finding (bare `catch` is in 3 sites, not 1), corrected the fix
for another (`int.strategy="eb"` maps to the existing `EmpiricalBayes`
strategy, not a missing `:empirical` symbol), and separated documented
roadmap items from real defects.

## Verified findings (with file:line)

- Version/doc drift: `README.md:16` says `v0.1.1`; packages + `CHANGELOG.md`
  say `1.0.0`. README also calls `INLASPDERasters` "scaffolding" (786 LOC of
  real code) and `benchmarks/` a "placeholder" (real harness).
- Missing extension sources: `LatentGaussianModels.jl/Project.toml:30-31`
  declares `LGMHCubatureExt` + `LGMIntegralsExt`; `ext/` has only Makie +
  PSIS. Breaks full `Pkg.precompile()`.
- Benchmark CI stale: `.github/workflows/benchmark.yml` uses Julia `1.10`
  (packages require 1.12) and targets `benchmarks/benchmarks.jl`; real
  entrypoint is `benchmarks/run.jl`.
- CI matrix gaps: `test.yml` omits `LGMFormula.jl` (v1.0.0, 4 src files) and
  `INLA.jl`; the justifying comment ("scaffolding without sources") is false
  for LGMFormula.
- Bare `catch` swallowing bugs — **3 sites**:
  `inference/empirical_bayes.jl:50` (worst — optimizer objective),
  `inference/full_laplace.jl:149`, `inference/log_density.jl:61`. Correct
  pattern already exists: `_is_bad_theta_failure` + `rethrow`
  (`inference/inla.jl:160`).
- Migration guide `docs/src/coming-from-r-inla.md:214` maps `int.strategy="eb"`
  → unsupported `int_strategy = :empirical`. Correct target is the existing
  `EmpiricalBayes()` strategy (`inference/empirical_bayes.jl`).
- Two benchmark dirs: `bench/` and `benchmarks/`; CI references neither.
- Allocation-heavy hot path (roadmap Phase 6): `laplace.jl:75` materializes
  `as_matrix` despite `apply!`/`apply_adjoint!` existing; `_symmetrize!`
  (`laplace.jl:266`) returns `(H+H')./2`; `marginal_variances` re-factorizes
  instead of reusing `LaplaceResult.factor`.
- Deeper correctness / test gaps: `_safe_inverse_hessian` eigenvalue floor
  leans on IS reweighting untested at θ boundaries; constraint-contract
  violations surface as opaque `PosDefException`; no Newton-convergence
  assertions or integration-scheme cross-checks in tests.

## Not defects (documented roadmap — do not treat as bugs)

- Gaussian `posterior_marginal_θ` (integrated design = future work).
- FullLaplace integration-stage means stay Gaussian until "PR-4".
- SimplifiedLaplace variance correction deferred.
  These are honestly labeled in docstrings/ADRs and parked in ROADMAP
  Phase 5/6.

## Remediation tiers

### Tier 0 — release-blocking (contradict the v1.0 claim)
1. Create `ext/LGMHCubatureExt.jl` + `ext/LGMIntegralsExt.jl`, or drop the
   `[extensions]`/`[weakdeps]` declarations if the features aren't wired.
2. Fix `benchmark.yml`: Julia `1.10`→`1`; target `benchmarks/run.jl`.
3. Reconcile README/docs version + kill "scaffolding"/"placeholder" language.
4. Add `LGMFormula.jl` + `INLA.jl` to CI matrix; delete false comment.

### Tier 1 — correctness & honesty (low effort)
5. Replace all 3 bare `catch` with the `_is_bad_theta_failure`/`rethrow`
   pattern. Prioritize `empirical_bayes.jl`.
6. Fix migration guide: `int.strategy="eb"` → `EmpiricalBayes()`.
7. ~~Consolidate `bench/` + `benchmarks/`~~ — **withdrawn.** On closer
   inspection these are distinct: `bench/` is the R-INLA *parity/correctness*
   reproducer (`oracle_compare.jl`), `benchmarks/` is the *performance*
   harness (`run.jl`). Not duplication. Only the stale README description of
   `benchmarks/` needed fixing (done in item 3).

### Tier 2 — test hardening  ✅ done (branch `test/review-2026-07-tier2`)
8. `test_newton_convergence.jl` — asserts every design-point Laplace in an
   `inla` fit converges within budget (Gaussian + Poisson), plus a robustness
   sweep fitting at extreme fixed θ (log-precision −6…+9).
9. `test_integration_scheme_consistency.jl` — Grid/CCD/GaussHermite agree on
   θ_mean, x_mean, and log-marginal for a 2-hyperparameter Gaussian+IID model
   (non-vacuous: 25/9/25 design points respectively).
10. `test_constrained_variances_dense.jl` — validates
   `_constrained_marginal_variances` against a full dense kriging oracle on a
   non-diagonal Besag precision, for a connected path graph (1 constraint) and
   a disconnected 2-component graph (one sum-to-zero per component).

### Tier 3 — performance (roadmap Phase 6)
11. ✅ Reuse `LaplaceResult.factor` for selected inversion (PR #22). 6.85×
    faster / ~34 MiB/call less on a 625-node constrained Besag.
12. ✅ Make `_symmetrize!` genuinely in-place. Zero-alloc two-pass average of
    the transpose pair, guarded by a structural-symmetry check with the
    `(H+H')/2` fallback; bit-identical (2 → 0 allocations per Newton step).
13. Fixed-sparsity-pattern numeric updates for `joint_precision`.
14. Hot path via `apply!`/`apply_adjoint!`; longer term a `LaplaceWorkspace`.
15. (new, from #22 benchmarking) The intrinsic null-space bump `V Vᵀ`
    densifies the regularised precision — biggest structural win. **Design
    written: ADR-045 (Proposed)** in `plans/decisions.md` — recommends
    staged low-rank Woodbury (option B) with a dense-bump fallback, KKT
    (option C) deferred. Implementation is spike-gated pending review.

### Tier 4 — feature completion (roadmap)
16. FullLaplace integration means; SimplifiedLaplace variance; integrated
    θ marginals.

## Execution status

- [x] Tier 0 (1–4) — done on branch `fix/review-2026-07-tier0-1`.
  1. `ext/LGMHCubatureExt.jl` + `ext/LGMIntegralsExt.jl` scaffold stubs.
  2. `benchmark.yml`: Julia `1`, targets `run.jl`, guards on `Rscript`.
  3. README version + scaffolding/placeholder language fixed.
  4. CI matrix adds `LGMFormula.jl` + `INLA.jl`; dev-link seed now includes
     in-repo weakdeps/extras (supersedes fragile relative `[sources]`).
- [x] Tier 1 (5–6) — done. All 3 bare catches narrowed; migration guide fixed.
- [~] Tier 1 (7) — withdrawn (see above).
- [x] Tier 2 (8–10) — done on branch `test/review-2026-07-tier2`. Three new
  regression files, wired into `runtests.jl`; all pass standalone (35 tests).
- Tiers 3–4: backlog.
