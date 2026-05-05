# Changelog

All notable changes to this repository are documented here. Format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versions follow [SemVer](https://semver.org/spec/v2.0.0.html).

## [v0.1.5] — 2026-05-04

Phase L is the marquee deliverable: `UserComponent` (R-INLA `rgeneric`
equivalent) lands as the structural extension hook for downstream users,
and `FullLaplace` (R-INLA's `strategy = "laplace"`) closes the
sharply-non-Gaussian-latent gap on the new Brunei oracle. Phases I-A tail
(IID3D / PCCor1 / Generic2), I-B (MEB / MEC), I-C (Replicate / Group),
J (six new likelihoods plus PCGevtail / BetaPrior), and K (posterior
tooling — sampling, predictive, PSIS-LOO, `refine_hyperposterior`,
`pp_check`) form the body. Patch release on `LatentGaussianModels.jl`
and the `INLA.jl` umbrella; `GMRFs.jl`, `INLASPDE.jl`, and
`INLASPDERasters.jl` are unchanged at v0.1.1. Oracle suite expands from
21 to 28 fixtures.

### Added

**Phase L — `UserComponent` and `FullLaplace`**

- **`UserComponent(callable; n, n_hyper)`** (PR-2, ADR-025) — the R-INLA
  `rgeneric` equivalent: a one-line callable wrapper around the
  `AbstractLatentComponent` contract that lets users port `rgeneric`
  models without subtyping. The callable returns a `NamedTuple` with
  required key `:Q` (sparse precision) and optional keys `:log_prior`,
  `:log_normc`, `:constraint`, all with defaults so future-extending
  the namedtuple is non-breaking. Constraint is read once at
  construction (θ-independent, mirroring R-INLA `rgeneric`'s
  `extraconstr`); precision and log-NC are re-evaluated at every θ.
  `cgeneric` is deliberately out of scope — users with C precision
  routines call them via `@ccall` inside the Julia closure.
  Implementation:
  [`src/components/user_component.jl`](packages/LatentGaussianModels.jl/src/components/user_component.jl).
- **`FullLaplace <: AbstractMarginalStrategy`** (PR-3 / PR-4) — R-INLA's
  `strategy = "laplace"`, the per-`x_i` refitted Laplace that closes
  the sharply-non-Gaussian-latent gap. Implemented via constraint
  injection: each grid point stacks `e_i^T x = a` onto the model-level
  constraint and reuses the existing kriging machinery in
  `inference/laplace.jl`. PR-4 ships warm-start Newton (each grid
  point bootstraps from the previous fit's mode, collapsing inner
  iterations from ~5–10 cold-start to 1–2) and adaptive truncation
  (sweep stops once running per-θ log-density drops 25 nats below
  the running max — a `< 1.4e-11` factor of the peak). Implementation:
  [`src/inference/full_laplace.jl`](packages/LatentGaussianModels.jl/src/inference/full_laplace.jl);
  bench harness at
  [`bench/full_laplace.jl`](packages/LatentGaussianModels.jl/bench/full_laplace.jl).
- **Brunei oracle** (PR-5,
  [`test/oracle/test_synthetic_brunei.jl`](packages/LatentGaussianModels.jl/test/oracle/test_synthetic_brunei.jl)) —
  Phase L acceptance gate: Poisson + Generic0 (RW1 + sum-to-zero) on a
  low-count, sharply non-Gaussian regime (n=24, E_i=1, most
  y_i ∈ {0, 1, 2, 3}), with a tight `loggamma(100, 100)` prior pinning
  τ ≈ 1. Asserts per-coordinate `FullLaplace` mean/sd against R-INLA's
  `strategy = "laplace"` reference fit, hyperposterior τ-mean within
  Phase F's 20% relative band, mlik finite within 5%, and FL ≠ SL on
  at least one coordinate (regression guard for the constraint-injection
  seam). Modern R-INLA's post-Laplace VB correction tightens SDs by a
  few percent — Julia's pure FL leaves this systematic offset
  uncorrected, so the SD tolerance is widened accordingly.
- **`crw2` `UserComponent` tutorial** (PR-6) at
  [`docs/src/vignettes/rgeneric-tutorial.md`](docs/src/vignettes/rgeneric-tutorial.md) —
  end-to-end implementation of R-INLA's `model = "crw2"` (continuous-time
  RW2 on irregularly-spaced knots) as a 30-line `UserComponent`. Closed-
  form precision `Q(τ; s) = τ · Dᵀ W D` derived from the standardised
  second-difference representation, verified against the built-in
  `RW2(n)` on a regular grid (max abs diff = machine zero), and pinned
  against R-INLA at Phase H tolerances on a synthetic
  `1 + s + f(idx, model = "crw2")` fit. Companion extension guide at
  [`docs/src/extending.md`](docs/src/extending.md) documents the two
  paths (`UserComponent` vs subtyping `AbstractLatentComponent`
  directly) as first-class alternatives.

**Phase K — diagnostics and posterior tools**

- **`posterior_sample(rng, res, model; n_samples)`** (PR-2) — joint
  draws of `(x, θ) ∼ π̂(· | y)` from the CCD/grid mixture. Pick `θ_k`
  with weight `res.θ_weights`, draw `x | θ_k ∼ N(mode_k, H_k⁻¹)` via
  the cached Cholesky (with hard-linear-constraint kriging correction
  when applicable). Complements the existing `posterior_samples_η`
  (which returns `η = A x` and is the right tool for WAIC / CPO / PIT)
  with the joint-`x` building block needed for random-effect contrast
  posteriors and Stan/NUTS triangulation.
- **`refine_hyperposterior(res, model, y; n_grid = 15, span = 4.0,
  skewness_correction = true)`** (PR-3) — R-INLA `inla.hyperpar`
  equivalent. Re-runs the integration stage against a denser `Grid`
  design without redoing the outer LBFGS / FD-Hessian pass; reuses
  `res.θ̂` and `res.Σθ` from the input fit. Closes the
  IntegratedNestedLaplace.jl-documented heavy-tail undersampling
  failure mode (Brunei pathology). Defaults differ from `inla()`'s
  (5×5 grid, no skewness correction) because the user calls
  `refine_hyperposterior` precisely when the default is too coarse.
  Refactors `fit(::INLA)` to factor the integration stage into a
  private `_inla_integrate` helper shared with `refine_hyperposterior`.
- **`posterior_predictive(rng, res, model; n_samples)`** (PR-4) —
  linear-predictor posterior predictive built on top of
  `posterior_sample`. Returns `(x, θ, η)` with
  `η[:, s] = mapping * draws.x[:, s]`; mapping accepts any
  `AbstractObservationMapping` (`LinearProjector`, `IdentityMapping`,
  `StackedMapping`) or a raw `AbstractMatrix` (auto-wrapped). The
  ADR-017 observation-mapping seam from Phase G is the load-bearing
  abstraction.
- **`posterior_predictive_y` + `sample_y` likelihood hook** (PR-6) —
  extends the predictive with response-scale draws `y_rep ∼ p(y | η, θ)`.
  `sample_y(rng, ℓ, η, θ)` is a new method on `AbstractLikelihood`;
  closed-form samplers ship for `GaussianLikelihood`,
  `PoissonLikelihood`, `BinomialLikelihood`,
  `NegativeBinomialLikelihood{LogLink}`, and `GammaLikelihood{LogLink}`.
  The default raises `ArgumentError`, so unsupported families
  (censored survival, zero-inflated) error explicitly.
  `pp_check(rng, res, model, y_obs; n_samples = 400)` is the thin
  convenience wrapper that returns `(y, y_rep::Matrix)` ready for
  Makie `density!` overlays.
- **PSIS-LOO** (PR-5) via the
  [`PSIS.jl`](https://github.com/arviz-devs/PSIS.jl) weakdep extension
  (`LGMPSISExt`). Pareto-smoothed importance sampling estimate of LOO
  elpd, supplementing the harmonic-mean `cpo` from Phase H. Surfaces
  per-observation `pareto_k` as a reliability diagnostic — values
  above 0.7 flag IS failure that the harmonic-mean `cpo` silently
  absorbs (Vehtari et al. 2017, 2024). `psis_loo` is exported as a
  no-method stub from
  [`src/inference/diagnostics.jl`](packages/LatentGaussianModels.jl/src/inference/diagnostics.jl);
  the extension fires on `using PSIS`.
- **Asymmetric per-axis skewness corrections for `Grid`** (PR-1) —
  optional `stdev_corr_pos` / `stdev_corr_neg` on `Grid` and the new
  `compute_skewness_corrections(log_post, θ̂, Σ)` helper implement
  Rue-Martino-Chopin (2009) §6.5: along each eigen-axis of the
  Gaussian approximation, the design point spacing is stretched on
  each side to match the actual log-posterior drop, so heavy-tailed
  log-precision posteriors get nodes placed where the mass lives.
  Stretches floored at 0.05, capped at 5.0. Quadrature weights
  remain standard-normal — IS reweight in the surrounding fit
  absorbs the proposal mismatch. `INLA(skewness_correction = true)`
  is the opt-in flag; default is `false` to match R-INLA's
  `control.inla$skew.corr.positive = 1.0`.

**Phase J — six new likelihoods, plus PCGevtail / BetaPrior**

- **`BetaLikelihood`** (PR-1) — R-INLA's mean-dispersion form
  `y | μ, φ ~ Beta(μφ, (1−μ)φ)` with `μ = expit(η)`, single
  hyperparameter `θ = log φ`. Closed-form ∇_η and ∇²_η via digamma /
  trigamma. Default hyperprior `GammaPrecision(1, 0.01)` matches
  R-INLA's `family = "beta"` `loggamma(1, 0.01)` bit-for-bit. Only
  `LogitLink` accepted; `LogLink`/`IdentityLink` raise `ArgumentError`.
  Implementation:
  [`src/likelihoods/beta.jl`](packages/LatentGaussianModels.jl/src/likelihoods/beta.jl).
- **`BetaBinomialLikelihood`** (PR-2) — R-INLA's mean-overdispersion
  form `y | n, μ, ρ ~ BetaBinomial(n, μs, (1−μ)s)` with
  `s = (1−ρ)/ρ`, `θ = logit(ρ)`. Closed-form derivatives via digamma /
  trigamma. Default hyperprior `GaussianPrior(0, √2)` matches R-INLA's
  `family = "betabinomial"` default. Implementation:
  [`src/likelihoods/betabinomial.jl`](packages/LatentGaussianModels.jl/src/likelihoods/betabinomial.jl).
- **`StudentTLikelihood`** (PR-3) — scaled Student-t observation model
  `y = η + ε/√τ`, `ε ~ Student-t(ν)`, two hyperparameters
  `θ = (log τ, log(ν − 2))`. The `ν > 2` floor matches R-INLA's
  `family = "T"` and ensures finite variance. Closed-form score and
  Hessian via the standard t-density expression. Implementation:
  [`src/likelihoods/studentt.jl`](packages/LatentGaussianModels.jl/src/likelihoods/studentt.jl).
- **`SkewNormalLikelihood`** (PR-4) — R-INLA `family = "sn"`
  parameterisation with `θ[1] = log τ` and `θ[2] = logit-skew` mapped
  to standardised skewness `γ = 0.988 · tanh(θ[2] / 2) ∈ (−0.988,
  0.988)`. Closed-form gradient and (diagonal) Hessian via the inverse
  Mills ratio `λ(t) = φ(t)/Φ(t)`, with `Distributions.logpdf` /
  `logcdf` for tail stability. Recovers Gaussian exactly at `γ = 0`.
  Implementation:
  [`src/likelihoods/skewnormal.jl`](packages/LatentGaussianModels.jl/src/likelihoods/skewnormal.jl).
- **`GEVLikelihood`** (PR-5) — R-INLA `family = "gev"`
  `F(y) = exp(− [1 + ξ √(τ s) (y − η)]^(−1/ξ))` with `θ[1] = log τ`,
  `θ[2] = ξ / xi_scale` (default `xi_scale = 0.1`, matching R-INLA's
  `gev.scale.xi`). Per-observation `weights` carry R-INLA's `scale`
  argument as a precision multiplier. Gumbel-limit guard at
  `|ξ| < 1.0e-6` keeps derivatives continuous at `ξ = 0`. Note: the
  GEV density is not globally log-concave; pick `initial_η` well
  inside the support. Implementation:
  [`src/likelihoods/gev.jl`](packages/LatentGaussianModels.jl/src/likelihoods/gev.jl).
- **`POMLikelihood`** (PR-6) — R-INLA `family = "pom"` proportional-
  odds ordinal regression with `K − 1` internal-scale cut points
  carried as likelihood hyperparameters under a single Dirichlet
  prior on the implied class probabilities at `η = 0`. Closed-form
  derivatives for the cumulative-logit link (globally log-concave in
  η). The Jacobian from `θ → α → π` is included in full — R-INLA's
  internal Dirichlet omits this correction (documented in
  `inla.doc("pom")`), so absolute mlik values differ between Julia
  and R-INLA by a fixed θ-independent constant while every posterior
  moment matches. Implementation:
  [`src/likelihoods/pom.jl`](packages/LatentGaussianModels.jl/src/likelihoods/pom.jl).
- **Multinomial via independent-Poisson reformulation** (PR-7,
  ADR-024) — categorical / multinomial regression supported through
  R-INLA's Multinomial-Poisson trick (Baker 1994; Chen 1985), not as
  a new likelihood type. Helpers `multinomial_to_poisson(Y;
  class_names)` and `multinomial_design_matrix(helper, X;
  reference_class)` in
  [`src/multinomial.jl`](packages/LatentGaussianModels.jl/src/multinomial.jl)
  reshape an `n × K` count matrix into long-format
  `(y, row_id, class_id, …)` with row-major layout, and build the
  sparse class-specific covariate block with corner-point β_K = 0
  identifiability. Per-row nuisance intercept attaches as
  `IID(n_rows; τ_init, fix_τ = true)` — `IID` gains both kwargs,
  mirroring R-INLA's `prec = list(initial = ..., fixed = TRUE)`. With
  `fix_τ = true`, `nhyperparameters(c) == 0`. The `fit(::INLA)` fast
  path for `n_hyperparameters(model) == 0` skips θ-mode optimisation
  and grid integration entirely — a single Laplace at
  `θ = Float64[]` is the exact posterior. Vignette at
  [`docs/src/vignettes/multinomial.md`](docs/src/vignettes/multinomial.md).
- **`PCGevtail(λ, interval; xi_scale)`** — penalised-complexity prior
  on the GEV tail (shape) parameter `ξ ∈ [low, high]` (Opitz et al.
  2018), reference `ξ = 0` (Gumbel). Linearised PC distance
  `d(ξ) = ξ`; user-scale density matches `inla.pc.dgevtail`; internal
  scale `θ = ξ / xi_scale` matches `GEVLikelihood`. Defaults
  `λ = 7, interval = (0, 1), xi_scale = 0.1` match R-INLA composed
  with `gev.scale.xi`. Implementation:
  [`src/priors/pc_gevtail.jl`](packages/LatentGaussianModels.jl/src/priors/pc_gevtail.jl).
- **`BetaPrior(a, b)`** — generic Beta(a, b) prior on a bounded-ratio
  user-scale parameter `p ∈ (0, 1)` via the logit transform.
  Distributions.jl-friendly entry point complementing the existing
  `LogitBeta`. Implementation:
  [`src/priors/beta.jl`](packages/LatentGaussianModels.jl/src/priors/beta.jl).
- **Six new R-INLA oracle fixtures**: `synthetic_beta`,
  `synthetic_betabinomial`, `synthetic_studentt`,
  `synthetic_skewnormal`, `synthetic_gev`, `synthetic_pom`,
  `synthetic_multinomial` — one per likelihood family, all under
  [`test/oracle/`](packages/LatentGaussianModels.jl/test/oracle/).

**Phase I-A tail — multivariate IID continued**

- **`IID3D` via LKJ stick-breaking** (PR-1b) — extends
  `IIDND_Sep{N}` to `N = 3` with the canonical-partial-correlation
  parameterisation: 3 log-precisions + 3 atanh CPCs. `Λ = G'G` with
  `G = L⁻¹ · D_τ^{1/2}` where `L` is the LKJ stick-breaking Cholesky
  factor; log NC generalises the N=2 case via `_logcosh` over the
  3 CPCs. Constructor caps at `N ≤ 3` per ADR-022 (separable form).
  `IID3D(n; …)` ergonomic alias ships alongside. Implementation:
  [`src/components/iidnd.jl`](packages/LatentGaussianModels.jl/src/components/iidnd.jl).
- **`PCCor1` PC prior** (PR-1c) — PC prior on a correlation
  `ρ ∈ (-1, 1)` with reference at `ρ = 1`, mirroring R-INLA's
  `pc.cor1` (Sørbye-Rue 2017). The λ calibration root has no closed
  form — solved by log-bisection on the monotone calibration
  function. Routes `log(1 ± ρ)` and `√(1 - ρ)` through
  `softplus(±2θ)` to stay finite at saturation (`|θ| ≳ 19`). Available
  as opt-in `ρprior = PCCor1(...)` for `AR1`. Per ADR-022, `PCCor1` is
  *not* the bivariate-IID correlation prior (which uses `PCCor0`,
  reference at `ρ = 0`); the name reflects R-INLA's reservation of
  `pc.cor1` for the AR(1) lag-1 prior. Implementation:
  [`src/priors/pc_cor1.jl`](packages/LatentGaussianModels.jl/src/priors/pc_cor1.jl).
- **`Generic2`** (PR-1d) — R-INLA `model = "generic2"`: hierarchical
  `(u, v)` joint Gaussian with `v ~ N(0, (τ_v C)⁻¹)` and
  `u | v ~ N(v, (τ_u I)⁻¹)`, block precision
  `Q = [[τ_u I, -τ_u I]; [-τ_u I, τ_u I + τ_v C]]`. Hyperparameter
  ordering matches R-INLA (`θ_1 = log τ_v, θ_2 = log τ_u`).
  Schur-complement-derived
  `½ log|Q|_+ = ½ n log τ_u + ½ (n − rd) log τ_v + ½ log|C|_+`; the
  user-independent `½ log|C|_+` is dropped per the
  F_GENERIC0/F_BYM2 "up to a constant" convention. `scale_model =
  true` applies the Sørbye-Rue geometric-mean scaling to `C`.
  Implementation:
  [`src/components/generic2.jl`](packages/LatentGaussianModels.jl/src/components/generic2.jl).
- **`IID2D` mlik gap closed (Phase F.5 tail)** — the ~8 nat gap that
  v0.1.4 documented as `isfinite`-only was a setup mismatch, not a
  `2diid` normalising-constant adjustment. With R's `y ~ -1 +
  intercept_1 + intercept_2 + …` convention, R-INLA does not
  auto-insert an improper intercept; switching the Julia oracle to
  `Intercept(prec = 1.0e-3, improper = false) × 2` matches. Per-
  intercept offset verified empirically as `½ log(1000) = 3.454`
  nats. Residual gap drops to ~1.5 nats (1.8% relative on
  `mlik_R = -85.59`); assertion tightened to `_rel_iid2d(...) <
  0.05`.

**Phase I-B — measurement error (ADR-023)**

- **`MEB(values; scale, τ_u_prior, τ_u_init)`** (PR-2b) — R-INLA's
  Berkson measurement-error component `model = "meb"`: proper
  `N(values, (τ_u · diag(scale))⁻¹)` prior with R-INLA defaults
  (`GammaPrecision(1, 1.0e-4)` on `log τ_u`, `τ_u_init = log(1000)`).
  β-via-`Copy` on the receiving likelihood per ADR-021/ADR-023 — the
  component owns only the latent block; the β-scaling lives on the
  receiving likelihood. Implementation:
  [`src/components/meb.jl`](packages/LatentGaussianModels.jl/src/components/meb.jl).
- **`MEC(values; scale, τ_*_prior, *_init, fix_*)`** (PR-2c) —
  R-INLA's classical measurement-error component `model = "mec"`:
  proper Gaussian prior `x ~ N(μ̂(θ), Q̂(θ)⁻¹)` with the conjugate-
  Gaussian Berkson tie folded in: `Q̂ = τ_x I + τ_u D`,
  `μ̂ = Q̂⁻¹ (τ_x μ_x 1 + τ_u D · values)`. Three component
  hyperparameters (`log τ_u`, `μ_x`, `log τ_x`) with per-slot
  `fix_*` toggles preserving R-INLA's "all default-fixed" pattern.
  Demonstrates the slope-attenuation phenomenon: naïve OLS shrinks
  by `τ_u / (τ_u + τ_x)`, while `MEC` recovers the truth.
  Implementation:
  [`src/components/mec.jl`](packages/LatentGaussianModels.jl/src/components/mec.jl).
- **`prior_mean(c, θ)` promoted to load-bearing** in the Laplace
  Newton loop. Gradient, log-joint, and log-marginal now use
  `Q (x − μ)` instead of `Q x`; the new `joint_prior_mean(m, θ)`
  helper stacks per-component prior means along
  `m.latent_ranges`. No-op for v0.1 zero-mean components, but
  required for `MEB` (`μ = values`) and `MEC` (`μ` is θ-dependent
  via the conjugate-Gaussian update). Touches
  [`src/inference/laplace.jl`](packages/LatentGaussianModels.jl/src/inference/laplace.jl)
  and
  [`src/model.jl`](packages/LatentGaussianModels.jl/src/model.jl).
- **Measurement-error regression vignette** at
  [`docs/src/vignettes/measurement-error-regression.md`](docs/src/vignettes/measurement-error-regression.md) —
  walks through both flavours, including the OLS-attenuation
  comparison and the structural-model variant of Muff et al. (2015)
  with `τ_x` left free.

**Phase I-C — replicate / group**

- **`Replicate(component, n_replicates)`** (PR-3a) — stacks
  `n_replicates` independent copies of an inner
  `AbstractLatentComponent` sharing one hyperparameter block,
  mirroring R-INLA's `f(idx, model = M, replicate = id)` for
  uniformly-sized panels. Length scales as
  `n_replicates · length(component)`; θ slot count is unchanged
  (the prior lives on a single shared θ block, *not* scaled by
  `n_replicates`). Precision is `blockdiag(Q_inner, …, Q_inner)`;
  `prior_mean` repeats the inner; constraints block-stack onto a
  block-diagonal `A`; `log_NC` is `n_replicates · log_NC(inner)`
  (factorises exactly); `gmrf` rankdef multiplies, so intrinsic
  stacks (Besag, RW) inherit the correct null-space dimension.
  Implementation:
  [`src/components/replicate.jl`](packages/LatentGaussianModels.jl/src/components/replicate.jl).
- **`Group(components | factory, group_id)`** (PR-3b) — generalises
  `Replicate` to non-uniform sizes: a single shared θ block across
  an arbitrary mix of inner components (e.g. unequal-length AR1
  panels for subjects with different visit counts). Vector form is
  the primary API; `Group(factory, group_id)` is the convenience
  overload that counts members per integer label. Inner constructor
  enforces the shared-θ contract by checking that all components
  have the same `nhyperparameters`. Implementation:
  [`src/components/group.jl`](packages/LatentGaussianModels.jl/src/components/group.jl).
- **Replicate(AR1) R-INLA oracle** (PR-3c,
  [`test/oracle/test_synthetic_replicate_ar1.jl`](packages/LatentGaussianModels.jl/test/oracle/test_synthetic_replicate_ar1.jl)) —
  `R = 30` panels of length `n = 20` (600 observations) sharing
  `(τ, ρ)`, fit by both Julia and R-INLA's
  `f(t, model = "ar1", replicate = id)`. Tolerances: τ_x within
  10% relative, ρ within 0.05 absolute, mlik within 2% relative.

### Changed

- **ADR-026: `AbstractMarginalStrategy` type-dispatch refactor**
  (Phase L PR-1) — the symbol-keyed `latent_strategy = :gaussian |
  :simplified_laplace` knob on
  [`INLA(...)`](packages/LatentGaussianModels.jl/src/inference/inla.jl),
  [`posterior_marginal_x(...)`](packages/LatentGaussianModels.jl/src/inference/marginals.jl),
  and `refine_hyperposterior(...)` is promoted to a multiple-dispatch
  type hierarchy: `Gaussian()`, `SimplifiedLaplace()`,
  `FullLaplace(n_grid, span)`. A new strategy adds methods on
  `_integration_mean_shift` (integration-stage hook) and
  `_density_skewness` (per-coordinate marginal hook) instead of
  another `if` arm. A symbol shim `_resolve_marginal_strategy(::Symbol)`
  mirrors `_resolve_scheme` so legacy `:gaussian` /
  `:simplified_laplace` keep working unchanged — no migration impact
  for existing code. The two facets of `SimplifiedLaplace` (mean
  shift in integration, Edgeworth in marginals) dispatch on the
  same type via different methods, preserving their semantic
  independence per ADR-016. Sets the seam Phase L PR-3 plugs
  `FullLaplace` into.
- **`UserComponent` extension guide** at
  [`docs/src/extending.md`](docs/src/extending.md) — documents the
  two extension paths (`UserComponent` vs subtyping
  `AbstractLatentComponent` directly) as first-class alternatives,
  not tiered fallbacks. Subtyping remains the right tool for
  component-specific `prior_mean(c, θ)` overrides (`MEC`-style
  shifted priors), custom `gmrf(c, θ)` factorisations, and
  lazy/structured precision matrices that don't fit
  `SparseMatrixCSC`.

### Fixed

- **`fit(::INLA)` fast path for `n_hyperparameters(model) == 0`**
  (Phase J PR-7) — skips θ-mode optimisation and θ-grid integration
  when the model has no free hyperparameters; a single Laplace at
  `θ = Float64[]` is the exact posterior. Triggered by the
  multinomial-via-Poisson reformulation (`IID(...; fix_τ = true)`
  + class-by-covariate `FixedEffects`), but benefits any future
  zero-hyperparameter model.
- **JET narrowing on `integration_nodes(::Grid, …)`** (Phase K PR-3)
  — hoist `scheme.stdev_corr_pos` / `…_neg` into typed
  `Vector{Float64}` locals before the per-point comprehension so
  the inner closure is inferred without the
  `Union{Nothing, Vector{Float64}}` field type leaking through.
  Functional behaviour unchanged.
- **Phantom plan-doc link cleanup** (Phase L PR-7) — strip references
  to a never-written `phase-i-and-onwards-mighty-emerson.md` from
  the changelog, the Phase I-A IID2D ADR cross-reference list, and
  two vignette pointers. Replaced with "tracked separately as
  post-v0.1 work" or removed where the heading already named the
  phase.

### Known limitations

- **Marginal log-likelihood gap on `weibullsurv`, `lognormalsurv`,
  `gammasurv`, and `coxph` oracles.** Carried forward unchanged from
  v0.1.4. Phase F.5 excavation (2026-05-02) traced this to a
  polynomial-form Laplace approximation in R-INLA's `GMRFLib` that
  differs from Julia's textbook formula at three points: the cubic
  contribution to the centered polynomial (`+⅙ x0³ dddf` vs the
  strict-Taylor `−⅙`), a modified Hessian carrying an `η̂·dddf`
  correction, and `*logdens` evaluated at sample = 0 rather than at
  the posterior mode. Closure requires modifying
  [`src/inference/laplace.jl`](packages/LatentGaussianModels.jl/src/inference/laplace.jl);
  deferred to v0.3 per the Phase Q rolling plan. Fixed-effect and
  hyperparameter posteriors agree tightly with R-INLA on these
  oracles. Oracle tests assert `isfinite(log_marginal)` while the
  gap is being characterised.
- **`FullLaplace` per-coordinate runtime is structurally above the
  replan's aspirational ≤ 5× `SimplifiedLaplace` target.** SL is
  closed-form per grid point (~µs); FL refits a constrained Laplace
  per grid point (~100µs minimum even with rank-1 CHOLMOD updates and
  per-`(θ, i)` caching). The Phase L acceptance gate is correctness
  against the Brunei oracle (PR-5), not the per-coordinate timing
  ratio. Future perf work (rank-1 + caching) can compress the ratio
  by ~2-3× but cannot reach single digits — see
  [`bench/README.md`](packages/LatentGaussianModels.jl/bench/README.md).
- **Brunei FL SD systematic offset.** Modern R-INLA (25.10.19)'s
  post-Laplace VB correction tightens SDs by a few percent on the
  Brunei regime; Julia's pure FL leaves this offset uncorrected, so
  the Brunei oracle widens the SD tolerance accordingly. The mean
  posterior matches at standard tolerance.

### Validated against

R-INLA `25.10.19`, R 4.5.x. Oracle suite expanded from 21 (v0.1.4) to
28 fixtures across
[`test/oracle/`](packages/LatentGaussianModels.jl/test/oracle/).
Fixture generation scripts under
[`scripts/generate-fixtures/`](scripts/generate-fixtures/).

## [v0.1.4] — 2026-05-02

Phase I-A PR-1a. Patch release on `LatentGaussianModels.jl` and the
`INLA.jl` umbrella; `GMRFs.jl`, `INLASPDE.jl`, and `INLASPDERasters.jl`
are unchanged at v0.1.1. First multivariate-IID building block lands:
the bivariate slot for joint-longitudinal-survival random effects,
paired-areal disease mapping, and bivariate meta-analysis. No public
API changed for existing components.

### Added

- **`PCCor0` PC prior on a correlation `ρ ∈ (-1, 1)`** with reference at
  `ρ = 0` (independence). Mirrors R-INLA's `pc.cor0` — used by `2diid`
  and `iid3d`. User-facing parameters `(U, α)` with `P(|ρ| > U) = α`
  give `λ = -log(α) / √(-log(1 - U²))`. Internal scale is
  `θ = atanh(ρ)`; the Jacobian cancels exactly with the
  Kullback-Leibler `|dd/dρ|` factor, leaving a closed-form log-density
  with a Taylor short-circuit at `|ρ|² < 1.0e-7` to avoid the formal
  `0/0`. Implementation:
  [`src/priors/pc_cor0.jl`](packages/LatentGaussianModels.jl/src/priors/pc_cor0.jl);
  regression suite covers symmetry, branch-boundary continuity,
  `λ`-monotonicity, and integral-to-1 over `θ ∈ ℝ`
  ([`test/regression/test_priors.jl`](packages/LatentGaussianModels.jl/test/regression/test_priors.jl)).
- **`IIDND_Sep{N}` family of multivariate IID random effects** with N=2
  shipped (PR-1a). The latent vector is `n·N` slots laid out as N
  consecutive `n`-blocks; joint precision is `Q = Λ ⊗ I_n` with `Λ` the
  inverse of the marginal covariance. For N=2 the constructor
  parameterises via `(τ_1, τ_2, ρ)` on internal scale `(log τ_1, log
  τ_2, atanh ρ)`, with PC priors on the marginal precisions and a
  `PCCor0` on the correlation by default — matches R-INLA's `2diid`
  default exactly. Implementation:
  [`src/components/iidnd.jl`](packages/LatentGaussianModels.jl/src/components/iidnd.jl).
- **`IID2D(n; …)` ergonomic alias** for `IIDND_Sep{2}` with sensible
  default priors (`PCPrecision()` × 2 + `PCCor0()`); accepts a Gaussian
  prior on Fisher-z if the user wants R-INLA's alternate
  `loggamma + atanh-ρ-Gaussian` form. PR-1b territory (`IID3D` +
  Cholesky/LKJ stick-breaking) and PR-1c (`Wishart`/`InvWishart` joint
  prior path) are scoped but not in this release.
- **Argument-validation tests** for `IIDND` reject `n ≤ 0`, `N = 1`,
  `N ≥ 3` (PR-1b territory), and conflicting `hyperprior_corr` /
  `hyperprior_corrs` kwargs
  ([`test/regression/test_iidnd.jl`](packages/LatentGaussianModels.jl/test/regression/test_iidnd.jl)).

### Changed

- **ADR-022 rename `PCCor1` → `PCCor0`** in
  [`plans/decisions.md`](plans/decisions.md). R-INLA's `pc.cor0`
  reserves the reference-at-`ρ = 0` name for the independence-anchored
  prior used by `2diid` / `iid3d`; `pc.cor1` is the
  reference-at-`ρ = 1` companion used by AR(1)'s lag-1 correlation.
  The ADR update was caught and corrected before any code shipped, so
  no migration impact for users — but the wrong name has now been
  burned into PR-1a's public API by the right one.

## [v0.1.3] — 2026-05-02

Phase Q PR-1. Patch release on `LatentGaussianModels.jl` and the
`INLA.jl` umbrella; `GMRFs.jl`, `INLASPDE.jl`, and `INLASPDERasters.jl`
are unchanged at v0.1.1. Closes both Phase Q v0.1 performance regressions
with a single LBFGS-tuning fix; quality numbers unchanged byte-for-byte.

### Changed

- **Default `g_tol = 1.0e-4`** in the outer θ-mode LBFGS for both
  [`INLA`](packages/LatentGaussianModels.jl/src/inference/inla.jl) and
  [`EmpiricalBayes`](packages/LatentGaussianModels.jl/src/inference/empirical_bayes.jl).
  The `Optimization.jl + AutoFiniteDiff` FD-gradient noise floor sits
  near `√eps ≈ 1.0e-8` — exactly Optim.jl's default `g_tol`, so LBFGS
  exhausted the 1000-iteration limit chasing noise. Raising `g_tol` to
  `1.0e-4` recovers the same θ̂ to ≲ 2e-4 in θ-space (well under any
  oracle test tolerance) at a fraction of the work. Users who want the
  prior tolerance can pass `optim_options = (; g_tol = 1.0e-8)` to
  override.
- **Pennsylvania BYM2 wall-clock**: 18.3 s → 0.15 s (119× faster),
  flipping the v0.1.x regression from `0.47×` of R-INLA to **48.8×
  faster** than R-INLA on `bench/oracle_compare.jl`. Phase Q's
  ≤ 1.2× R-INLA acceptance criterion met with 40× headroom.
- **Meuse SPDE wall-clock**: 142 s → 1.32 s (108× faster) under the
  default SuiteSparse backend, flipping the v0.1.x regression from
  27× *slower* than R-INLA to **4.46× faster**. Phase Q's
  ≤ 2× R-INLA-under-Pardiso acceptance criterion met under SuiteSparse
  alone, without needing the `GMRFsPardiso.jl` backend.
- **Quality unchanged.** `bench/oracle_compare_julia.md` Quality table
  is byte-identical to v0.1.2: every `fixed_max_rel`, `hyperpar_max_rel`,
  `mlik_rel` matches to four significant figures across all 11 oracle
  problems.

### Diagnostic

- [`bench/diagnostics/pa_bym2_hessian.jl`](bench/diagnostics/pa_bym2_hessian.jl)
  — eight-stage diagnostic that pinned the regression to the outer
  LBFGS rather than the integration grid or the FD Hessian. Refuted the
  replan-2026-04-28 hypothesis ("wider Hessian at θ̂ → wider grid
  envelope"); FD Hessian was correct, mode-finding was the bottleneck.

## [v0.1.2] — 2026-05-02

Phase F.5 close. Patch release on `LatentGaussianModels.jl` and the
`INLA.jl` umbrella; `GMRFs.jl`, `INLASPDE.jl`, and `INLASPDERasters.jl`
are unchanged at v0.1.1.

### Added

- **`synthetic_baghfalaki` R-INLA oracle** for the joint longitudinal-
  Gaussian + Weibull-PH survival model
  ([`test/oracle/test_synthetic_baghfalaki.jl`](packages/LatentGaussianModels.jl/test/oracle/test_synthetic_baghfalaki.jl)).
  Promotes `test/regression/test_inla_joint_baghfalaki.jl` to oracle
  parity: Julia and R-INLA fit the same dataset via
  `y_resp = list(y_gauss, inla.surv(...))` with `family = c("gaussian",
  "weibullsurv")` and `f(b_surv_idx, copy = "b_long_idx", fixed = FALSE)`.
  Asserts fixed-effects (10%), hyperparameters (20%), β-Copy (15% mean /
  50% sd, wide on the asymmetric posterior), and per-subject `b̂`
  correlation > 0.99. mlik kept as `isfinite` only — joint inherits the
  polynomial-form Laplace gap from the Weibull arm.
- **Joint longitudinal + survival vignette** at
  [`docs/src/vignettes/joint-longitudinal-survival.md`](docs/src/vignettes/joint-longitudinal-survival.md).
  End-to-end Baghfalaki et al. (2024)-style synthetic recovery via
  `StackedMapping`, `Copy`, and the multi-likelihood
  `LatentGaussianModel` (the marquee Phase G PR3 deliverable that was
  missing from the docs site at v0.1.1).

### Changed

- **Phase F.5 calibration excavation.**
  Survival oracle headers
  ([`weibullsurv`](packages/LatentGaussianModels.jl/test/oracle/test_synthetic_weibull_survival.jl),
  [`lognormalsurv`](packages/LatentGaussianModels.jl/test/oracle/test_synthetic_lognormal_survival.jl),
  [`gammasurv`](packages/LatentGaussianModels.jl/test/oracle/test_synthetic_gamma_survival.jl))
  rewritten to record the polynomial-form-Laplace finding that closed
  the 2-week excavation: R-INLA's `GMRFLib` differs from Julia's
  textbook formula at three points (cubic contribution `+⅙ x0³ dddf` vs
  `−⅙`, `η̂·dddf`-corrected Hessian, `*logdens` evaluated at sample = 0
  rather than at the mode). Closure requires modifying
  `src/inference/laplace.jl`; deferred to v0.3 / Phase Q.
- **Documentation: registry policy.** READMEs and `docs/src/index.md`
  drop "not yet on the General registry / submission planned" framing;
  the personal registry at `haavardhvarnes/JuliaRegistry` is the
  documented install path.

## [v0.1.1] — 2026-05-02

Second release line. Multi-likelihood `LatentGaussianModel`, censoring
infrastructure, five new survival likelihoods, the zero-inflated count
pack, the `Copy` component, and the `AbstractObservationMapping` seam.
Julia 1.12+ requirement.

### Added

- **Multi-likelihood `LatentGaussianModel`** (Phase G PR2,
  [`8d66bc9`](https://github.com/HaavardHvarnes/INLA.jl/commit/8d66bc9)).
  A single `LatentGaussianModel` mounts more than one likelihood block
  over a stacked observation vector via `StackedMapping`. Block-diagonal
  observation mappings compose with per-block likelihoods, so joint
  Gauss-Poisson, joint longitudinal-survival, and similar mixed-family
  models share one latent and one hyperparameter posterior. Covered by
  the new `synthetic_joint_gauss_pois` oracle.
- **`AbstractObservationMapping` seam** (Phase G PR1, ADR-017,
  [`fb46f71`](https://github.com/HaavardHvarnes/INLA.jl/commit/fb46f71)).
  The projector from latent `x` to linear predictor `η` is now a typed
  abstraction with three concrete implementations: `IdentityMapping`,
  `LinearProjector` (sparse `A`-matrix; the existing default),
  `StackedMapping` (multi-likelihood block-row stack), plus a
  `KroneckerMapping` stub reserved for Phase M space-time SPDE.
- **`Copy` component** (Phase G PR3, ADR-021,
  [`f20bbfd`](https://github.com/HaavardHvarnes/INLA.jl/commit/f20bbfd) →
  [`2113623`](https://github.com/HaavardHvarnes/INLA.jl/commit/2113623)).
  Joint-effect sharing à la R-INLA's `f(., copy = "name")`. The
  β-scaling lives on the receiving likelihood (not on the projection
  mapping), preserving the separation between observation maps and
  likelihood logic. Backed by an `add_copy_contributions!` hook on
  `AbstractLikelihood`, a closed-form fixed-β oracle, and a joint
  longitudinal + Weibull survival regression test.
- **Censoring infrastructure** (ADR-018,
  [`b97b84f`](https://github.com/HaavardHvarnes/INLA.jl/commit/b97b84f)).
  Per-observation `Censoring` enum (`NONE`, `LEFT`, `RIGHT`,
  `INTERVAL`); survival likelihoods accept a
  `censoring::Vector{Censoring}` field and dispatch internally for
  log-density and η-derivatives.
- **Five new survival likelihoods** (ADR-018):
  - `ExponentialLikelihood`
    ([`b97b84f`](https://github.com/HaavardHvarnes/INLA.jl/commit/b97b84f)).
  - `WeibullLikelihood` — PH parameterisation, shape `α_w` as a
    hyperparameter
    ([`5071990`](https://github.com/HaavardHvarnes/INLA.jl/commit/5071990));
    `PCAlphaW` PC prior (Sørbye–Rue 2017) alongside the
    `loggamma(1, 0.001)` default
    ([`ec18458`](https://github.com/HaavardHvarnes/INLA.jl/commit/ec18458)).
  - `LognormalSurvLikelihood` — AFT parameterisation, precision `τ` on
    `log T`
    ([`1b6f54c`](https://github.com/HaavardHvarnes/INLA.jl/commit/1b6f54c)).
  - `GammaSurvLikelihood` — mean parameterisation, shape `φ`
    ([`5a7c327`](https://github.com/HaavardHvarnes/INLA.jl/commit/5a7c327)).
  - `WeibullCureLikelihood` — Weibull mixture-cure with logistic cure
    fraction
    ([`1501296`](https://github.com/HaavardHvarnes/INLA.jl/commit/1501296)).
- **Cox proportional-hazards via data augmentation** (ADR-018 PR4,
  [`6876b3a`](https://github.com/HaavardHvarnes/INLA.jl/commit/6876b3a)).
  `inla_coxph(time, event)` produces the Holford / Laird-Olivier
  piecewise-exponential-as-Poisson augmentation; `coxph_design` builds
  the matching design matrix.
- **Zero-inflated likelihood pack** (ADR-019,
  [`925d853`](https://github.com/HaavardHvarnes/INLA.jl/commit/925d853)).
  Three R-INLA parameterisations (types 0, 1, 2) × three base count
  families (Poisson, Binomial, NegativeBinomial) = nine new
  likelihoods. ZIP1 oracle vs R-INLA's
  `family = "zeroinflatedpoisson1"`
  ([`b1ab680`](https://github.com/HaavardHvarnes/INLA.jl/commit/b1ab680)).
- **Opt-in simplified-Laplace mean-shift** (ADR-016,
  [`fbe9b50`](https://github.com/HaavardHvarnes/INLA.jl/commit/fbe9b50)).
  `inla(...; latent_strategy = :simplified_laplace)` applies a per-row
  mean-shift correction at the cost of one extra Newton step per
  integration node. The variance correction remains deferred to v0.3
  (Phase Q in the rolling plan). `pennsylvania_bym2` oracle covers the
  new pathway.
- **Survival vignettes.** CoxPH and Weibull survival under the new
  censoring infrastructure, published in
  [`docs/src/vignettes/coxph-weibull-survival.md`](docs/src/vignettes/coxph-weibull-survival.md)
  ([`fd5ac78`](https://github.com/HaavardHvarnes/INLA.jl/commit/fd5ac78)).
- **Seven new R-INLA oracle fixtures**:
  `synthetic_exponential_survival`, `synthetic_weibull_survival`,
  `synthetic_lognormal_survival`, `synthetic_gamma_survival`,
  `synthetic_coxph`, `synthetic_zip1`, plus
  `synthetic_joint_gauss_pois`.

### Changed

- **Drop Julia 1.10 LTS support** (ADR-020,
  [`4b90410`](https://github.com/HaavardHvarnes/INLA.jl/commit/4b90410)).
  All four `src/`-bearing packages now require Julia 1.12+. Back-compat
  shims for `Returns`, `public` markers, and similar 1.11+ features
  have been removed.
- **Project versions bumped to 0.1.1** across `GMRFs.jl`,
  `LatentGaussianModels.jl`, `INLASPDE.jl`, `INLASPDERasters.jl`, and
  the `INLA.jl` umbrella
  ([`700f218`](https://github.com/HaavardHvarnes/INLA.jl/commit/700f218));
  `[compat]` widened for fresh installs.
- **R-INLA fixture regen against `25.10.19`**
  ([`9f98a64`](https://github.com/HaavardHvarnes/INLA.jl/commit/9f98a64)
  + follow-up CI hardening through
  [`44093ab`](https://github.com/HaavardHvarnes/INLA.jl/commit/44093ab)).
  Tolerance comparator replaces the previous byte-level diff to
  accommodate floating-point drift across R-INLA point releases.

### Known limitations

- **Marginal log-likelihood gap on `weibullsurv`, `lognormalsurv`,
  `gammasurv`, and `coxph` oracles.** Phase F.5 excavation
  (2026-05-02) traced
  this to a polynomial-form Laplace approximation in R-INLA's
  `GMRFLib` that differs from Julia's textbook formula at three
  points: the cubic contribution to the centered polynomial
  (`+⅙ x0³ dddf` vs the strict-Taylor `−⅙`), a modified Hessian
  carrying an `η̂·dddf` correction, and `*logdens` evaluated at
  sample = 0 rather than at the posterior mode. Closure requires
  modifying
  [`src/inference/laplace.jl`](packages/LatentGaussianModels.jl/src/inference/laplace.jl);
  deferred to v0.3 per the Phase Q rolling plan. Fixed-effect and
  hyperparameter posteriors agree tightly with R-INLA on these
  oracles. Oracle tests assert `isfinite(log_marginal)` while the
  gap is being characterised.
- **Coxph augmentation `mlik` shifted by `Σ_events log E_{k_last,i}`**
  — the η-independent exposure of the interval the event lands in.
  Cancels in the posterior of `(γ, β)` so it does not affect
  inference. See the algebraic-equivalence regression test
  ([`test/regression/test_coxph_augmentation.jl`](packages/LatentGaussianModels.jl/test/regression/test_coxph_augmentation.jl)).

### Validated against

R-INLA `25.10.19` (CI fixture regen on
[`9f98a64`](https://github.com/HaavardHvarnes/INLA.jl/commit/9f98a64)),
R 4.5.x. Fixture generation scripts under
[`scripts/generate-fixtures/`](scripts/generate-fixtures/).

## [v0.1.0] — 2026-04-28

First tagged release on the user's personal Julia registry. No content
changes versus `v0.1.0-rc1`; release-prep cleanup only.

### Changed

- **Drop `[sources]` blocks from `INLASPDERasters.jl`** to enable
  registration
  ([`06df56a`](https://github.com/HaavardHvarnes/INLA.jl/commit/06df56a)).
- **Version bump** to `v0.1.0` across `GMRFs.jl`,
  `LatentGaussianModels.jl`, `INLASPDE.jl`, and `INLASPDERasters.jl`
  ([`dad9f17`](https://github.com/HaavardHvarnes/INLA.jl/commit/dad9f17)).

## [v0.1.0-rc1] — 2026-04-28

First publicly-usable release line of the Julia INLA ecosystem. The four
`src/`-bearing packages —
[`GMRFs.jl`](packages/GMRFs.jl/),
[`LatentGaussianModels.jl`](packages/LatentGaussianModels.jl/),
[`INLASPDE.jl`](packages/INLASPDE.jl/),
[`INLASPDERasters.jl`](packages/INLASPDERasters.jl/) — cover the
canonical R-INLA datasets within the testing-strategy tolerances.

### Added

- **GMRFs.jl** — sparse Gaussian Markov random field core. Concrete
  types: `IIDGMRF`, `RW1GMRF`, `RW2GMRF`, `AR1GMRF`, `SeasonalGMRF`,
  `BesagGMRF`, `Generic0GMRF`. `GMRFGraph` wraps any sparse adjacency
  for `Graphs.jl` interop. Sampling, log-density, log-determinant,
  marginal variances via selected inversion, sparse-Cholesky
  factor caching (`FactorCache`).
- **LatentGaussianModels.jl** — LGM stack on top of GMRFs.
  Components: `Intercept`, `FixedEffects`, `IID`, `RW1`, `RW2`, `AR1`,
  `Seasonal`, `Besag`, `BYM`, `BYM2`, `Leroux`, `Generic0`,
  `Generic1`. Likelihoods: `Gaussian`, `Poisson`, `Binomial`,
  `NegativeBinomial`, `Gamma` (closed-form gradients/Hessians;
  ForwardDiff fallback for user-defined). Inference strategies:
  `Laplace`, `EmpiricalBayes`, `INLA`. θ-integration schemes: `Grid`,
  `GaussHermite`, `CCD` (`int_strategy = :auto` chooses CCD for
  dim θ > 2, Grid otherwise). Diagnostics: DIC, WAIC, CPO, PIT.
  Hyperpriors: `PCPrecision`, `GammaPrecision`, `LogNormalPrecision`,
  `WeakPrior`, `PCBYM2Phi`, `LogitBeta`. `LogDensityProblems` seam for
  external samplers.
- **INLASPDE.jl** — SPDE–Matérn FEM on triangulated meshes. `SPDE2`
  component for α = 2. `PCMatern` joint PC prior on (range, σ).
  `inla_mesh_2d` constrained-Delaunay mesh generator (DT.jl-native;
  fmesher-equivalent on convex domains). `MeshProjector` A-matrix as
  a `SciMLOperators.AbstractSciMLOperator`.
- **INLASPDERasters.jl** — package scaffolding only; raster ↔ SPDE
  glue (`extract_at_mesh`, `predict_raster`) is planning. Activates
  in v0.2.
- **Oracle test suite.** Eleven R-INLA-derived JLD2 fixtures across
  Scotland and Pennsylvania BYM2, classical BYM, synthetic Gamma /
  Negative Binomial / Generic0 / Generic1 / Seasonal / Leroux /
  disconnected Besag, and Meuse SPDE. R-INLA wall-clock timings
  (`fit$cpu.used`) are stored alongside posteriors so reproductions
  are fully self-contained.
- **`bench/oracle_compare.jl`.** Reproducible parity benchmark over
  all eleven oracle problems; emits a markdown table of relative
  errors and side-by-side wall-clock seconds vs the stored R-INLA
  timing. See [`bench/README.md`](bench/README.md).
- **Documenter site** with three vignettes (Scotland BYM2, Tokyo
  rainfall, Meuse SPDE) and per-package API pages under
  [`docs/src/`](docs/src/).
- **LGMTuring.jl** sub-package providing the NUTS bridge for
  INLA-vs-MCMC triangulation.

### Changed

- **Seasonal log-NC and constraint convention**
  ([`4020589`](https://github.com/HaavardHvarnes/INLA.jl/commit/4020589)).
  `SeasonalGMRF` declares a single sum-to-zero constraint matching
  R-INLA's `model = "seasonal"`. Per-component
  `log_normalizing_constant` uses `rd_eff = period` (not `period − 1`)
  because the constraint hits `range(Q)` rather than `null(Q)`,
  consuming one PD direction. Closes the τ\_seas / mlik gap to R-INLA.
- **BYM log-NC**
  ([`7d1cab7`](https://github.com/HaavardHvarnes/INLA.jl/commit/7d1cab7)).
  `BYM` per-component `log_normalizing_constant` matches R-INLA's
  `extra()` for `F_BYM`: `−¼(2n − K) log(2π)` where `K` is the number
  of connected components. Closes the Scotland BYM mlik gap.
- **Generic0 / Generic1 log-NC**
  ([`3e28604`](https://github.com/HaavardHvarnes/INLA.jl/commit/3e28604)).
  Both match R-INLA's shared `F_GENERIC0` `extra()` branch
  (`inla.c:2986-2987`), with the Gaussian normaliser
  `−½(n − rd) log(2π) + ½(n − rd) θ` applied per component.
- **`Intercept()` is improper by default**
  ([`41c986b`](https://github.com/HaavardHvarnes/INLA.jl/commit/41c986b)),
  matching R-INLA's `prec.intercept = 0`. Closes the constant
  ½ log(prec) shift in BYM2 / BYM joint Gaussian normalising
  constants. Pass `Intercept(prec = …)` for the proper-Normal
  variant.
- **Phase-B feature scope trimmed to MVP**
  ([`ebf8b42`](https://github.com/HaavardHvarnes/INLA.jl/commit/ebf8b42))
  ahead of the rc1 cut.

### Fixed

- **Per-component Sørbye-Rue scaling on disconnected graphs**
  ([`c6547a4`](https://github.com/HaavardHvarnes/INLA.jl/commit/c6547a4)).
  `BesagGMRF` and `BYM2` now scale each connected component
  independently per Freni-Sterrantino, Ventrucci & Rue (2018), and
  emit one sum-to-zero constraint per component. Was the most common
  silent failure mode in disease-mapping models on disconnected
  regions.
- **SPDE2 log-normalizing-constant** for Meuse-class meshes
  ([`bd70f40`](https://github.com/HaavardHvarnes/INLA.jl/commit/bd70f40)).

### Known limitations

Honest list of cases where the rc1 line knowingly diverges from R-INLA:

- **Scotland classical-BYM `τ_b` weakly identified.** Posterior mean
  diverges from R-INLA's by ≈ 60 % at n = 56; the `b`-vs-`u` split
  is not data-identified. Marked `@test_broken` in
  [`test/oracle/test_scotland_bym.jl:101`](packages/LatentGaussianModels.jl/test/oracle/test_scotland_bym.jl).
  Mean of `b + u` and the marginal log-likelihood agree to 1 %.
- **`disconnected_besag` τ posterior mean is heavy-tailed** at n = 12
  (R-INLA's mean ≈ 7587 with sd ≈ 102906; median ≈ 48). The oracle
  test asserts only that the Julia mlik is finite. Use the median or
  smaller fixed-θ grids for tight comparisons here.
- **BYM2 / Leroux φ and ρ are weakly identified by design** at the
  sample sizes in the oracle suite. Reported residual errors of
  ~10–20 % on these are expected, not regressions.
- **Performance regressions vs R-INLA on two cases.**
  `pennsylvania_bym2` runs in ~17.5 s vs R-INLA's 7.4 s (≈ 2.4×
  slower; suspected θ-grid envelope from a wider Hessian at θ̂),
  and `meuse_spde` runs in ~140 s vs R-INLA's 5.3 s (≈ 27× slower;
  R-INLA uses GMRFLib's tuned sparse Cholesky on the mesh-scale
  precision). Every other oracle problem is 10×–1230× faster than
  R-INLA — see the bench harness output for full numbers.

### Validated against

R-INLA `25.x` (see fixture `inla_version` field), R 4.5.x. Fixtures
are regenerated via the scripts under
[`scripts/generate-fixtures/`](scripts/generate-fixtures/).
