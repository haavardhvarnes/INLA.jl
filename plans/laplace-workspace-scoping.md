# Scoping: `LaplaceWorkspace` — and why it's the wrong investment

> 2026-07-05. Scoping exercise for Tier-3 items 13/14 of the 2026-07
> review remediation (`plans/review-2026-07-remediation.md`). Conclusion:
> **do not build the workspace.** The dominant allocation is elsewhere and
> has a much smaller, higher-value fix.

## The question

Items 13 (fixed-sparsity-pattern `joint_precision`) and 14 (route the
hot path through `apply!`/`apply_adjoint!` instead of `as_matrix`; longer
term a reusable `LaplaceWorkspace`) both aim to cut per-iteration
allocation in the inner Newton loop by preallocating buffers and reusing
fixed sparsity patterns across `laplace_mode` calls.

Before designing that, we profiled where `laplace_mode` allocation
actually goes, on a representative constrained disease-mapping model
(625-node connected Besag, Poisson likelihood, sum-to-zero constraint).

## Measured allocation profile (warm, stable across repeats)

`laplace_mode` total: **22.0 MiB**. Breakdown:

| Source | Alloc | Share | What a workspace targets |
|---|---:|---:|---|
| **`joint_precision`** | **18.6 MiB** | **85%** | — (not the buffer/pattern path) |
| Sparse Cholesky (`FactorCache` + `update!`×iters + solves) | ~1.5 MiB | 7% | not addressable by buffers (CHOLMOD internals) |
| Hessian assembly `Q + JᵀDJ` | 0.13 MiB/iter (~0.8 MiB) | 3.4% | ✅ the workspace's main target |
| Vectors (`η`, `∇η`, `∇²η`, `g`, `Δx`) + likelihood derivatives | ~0.05 MiB | <1% | ✅ workspace |
| selected inversion (`marginal_variances`) | 0.5 MiB | — (separate, per θ-point) | — |

**A buffer-reuse / fixed-pattern `LaplaceWorkspace` targets ~4% of
allocation.** It is not worth the plumbing.

## The real source: `scale_factor` recomputed every `precision_matrix` call

`joint_precision → precision_matrix(c::Besag, θ) →
GMRFs.precision_matrix(gmrf(c, θ))` rebuilds a fresh `BesagGMRF` each
call, which recomputes the Sørbye–Rue geometric-mean-variance
`scale_factor`. `per_component_scale_factors`
([`GMRFs.jl/src/gmrf.jl:413`](../packages/GMRFs.jl/src/gmrf.jl)) forms a
**dense** `inv(Qperp)` (`O(n³)`, ~3 MiB for a 625×625) per connected
component, per call.

But `scale_factor` depends only on the **graph** — it is **constant in
θ** (measured: identical value `0.802765` across calls). It is recomputed
on every one of the dozens of `laplace_mode` calls per INLA fit
(θ-optimisation + finite-difference Hessian + integration points).

Measured:

- `scale_factor(g)`: **18.3 MiB/call**, θ-independent.
- `precision_matrix` with a **cached** scale factor (`τ·sf·L`, sparse):
  **53 KiB**.
- **Caching cuts ~18.3 MiB per `precision_matrix` call**, i.e. ~85% of
  `laplace_mode`, multiplied by the call count per fit.

Affects every intrinsic scaled component that routes through
`scale_factor`: **Besag, BYM, BYM2, ICAR, Leroux** — the disease-mapping
workhorses, and the bulk of the R-INLA parity surface. Impact grows with
graph size² (dense inverse), so it is small on the tiny benchmark graphs
(Scotland 56, Pennsylvania 67 areas) but dominant on realistic
hundred-to-thousand-area fields.

## Recommendation

1. **Do not build the `LaplaceWorkspace`** (items 13/14 as framed). The
   4% it targets is not worth the complexity, and it does not touch the
   85%.
2. **Cache the θ-independent `scale_factor`** (and the fixed structure
   matrix) at component construction — compute once in the `Besag` / `BYM`
   / `BYM2` / `Leroux` constructor (it depends only on the graph), store
   it in the struct, and have `precision_matrix` reuse it. `precision_matrix`
   then reduces to a scalar × sparse-Laplacian assembly. Bit-identical
   (the scale factor value is unchanged), oracle-gated by the existing
   constrained regression + R-INLA fixtures.
3. Defer/close items 13/14: the fixed-pattern `joint_precision` (item 13)
   is subsumed once `precision_matrix` is cheap; the lazy-`apply!` path
   (item 14) matters only for large *space-time* mapping (`KroneckerMapping`)
   and is a separate, later concern — it is <1% here.

## Design note for the caching fix (contained)

- The scale factor is a per-graph constant. Cheapest cache: a `Float64`
  field on the component struct, populated in the constructor via
  `scale_factor(graph)` when `scale_model = true` (and `1.0` otherwise).
- `gmrf(c, θ)` / `precision_matrix(c, θ)` pass the cached factor through
  instead of recomputing. This may need a `BesagGMRF` constructor variant
  that accepts a precomputed scale factor (to avoid the recompute inside
  `GMRFs.precision_matrix`), or a `precision_matrix(c, θ)` that assembles
  `τ · sf · L` directly from the cached `sf` and `laplacian_matrix(graph)`.
- Watch: disconnected graphs (per-component factors), `scale_model = false`
  (factor ≡ 1), and the BYM2 mixing where the scaled structure enters the
  `φ`-weighted combination. All are covered by existing regression +
  oracle tests.
- Risk: low. No numerical change (cached value is the same); the only
  behaviour change is speed/allocation. Gate on the full LGM suite.
