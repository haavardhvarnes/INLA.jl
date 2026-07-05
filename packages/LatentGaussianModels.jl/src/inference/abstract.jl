"""
    AbstractInferenceStrategy

Dispatch type for the `fit` entry point. Concrete strategies:

- `Laplace` — fit at fixed `θ`, return the Gaussian approximation
  `x | θ, y ≈ N(x̂, (Q + A' D A)⁻¹)` (`D` from the likelihood Hessian).
- `EmpiricalBayes` — plug-in estimate `θ̂ = argmax π(θ | y)` via the
  outer Laplace log-marginal, then Laplace at `θ̂`.
- `INLA` — the full thing (deferred).
"""
abstract type AbstractInferenceStrategy end

"""
    AbstractInferenceResult

Return type for `fit`.
"""
abstract type AbstractInferenceResult end

"""
    AbstractMarginalStrategy

Dispatch type for the per-coordinate posterior marginal density of
`x_i | y` and, in the INLA integration stage, for the per-θ
approximation of `x_mean` / `x_var`. Mirrors R-INLA's
`control.inla\$strategy`.

Concrete strategies (see ADR-026):

- [`Gaussian`](@ref) — Gaussian centred at the Newton mode, with
  constraint-corrected Laplace marginal variance. R-INLA's
  `strategy = "gaussian"`.
- [`SimplifiedLaplace`](@ref) — Newton mode plus the Rue-Martino mean
  shift in the integration stage; Edgeworth first-order skewness
  correction in the per-coordinate density. Reduces to `Gaussian`
  when the likelihood third derivative `∇³_η log p` is zero. R-INLA's
  `strategy = "simplified.laplace"`.
- [`FullLaplace`](@ref) — per-`x_i` refitted Laplace via constraint
  injection. The reference quality strategy when the latent posterior
  is sharply non-Gaussian. R-INLA's `strategy = "laplace"`. Drives
  both [`posterior_marginal_x`](@ref) and — since the ADR-026 "PR-4"
  follow-up — the integration-stage `x_mean` / `x_var` summaries,
  which are replaced by trapezoid moments of the refitted-Laplace
  mixture after the Gaussian accumulation pass.

The mean-shift facet of `SimplifiedLaplace` (integration-stage) and
its density-skew facet (per-coordinate marginals) are independent —
they happen at different points in the pipeline and are dispatched on
the same type but through different internal hooks. See ADR-016 and
ADR-026.
"""
abstract type AbstractMarginalStrategy end

"""
    Gaussian()

Marginal strategy: per-θ Gaussian centred at the Newton mode, with
constraint-corrected Laplace marginal variance. R-INLA's
`strategy = "gaussian"`. Default for both
[`INLA`](@ref) (integration stage) and
[`posterior_marginal_x`](@ref) (per-coordinate density).

See [`AbstractMarginalStrategy`](@ref).
"""
struct Gaussian <: AbstractMarginalStrategy end

"""
    SimplifiedLaplace()

Marginal strategy: Rue-Martino mean shift in the integration stage,
plus Edgeworth first-order skewness correction in the per-coordinate
density. Reduces to [`Gaussian`](@ref) when the likelihood third
derivative `∇³_η log p` is zero. R-INLA's
`strategy = "simplified.laplace"`.

See [`AbstractMarginalStrategy`](@ref) and ADR-016.
"""
struct SimplifiedLaplace <: AbstractMarginalStrategy end

"""
    FullLaplace(; n_grid = 51, span = 5.0)

Per-`x_i` refitted Laplace marginal strategy via constraint injection.
R-INLA's `strategy = "laplace"`. Definition + helpers in
`inference/full_laplace.jl`; this declaration lives in `abstract.jl`
alongside [`Gaussian`](@ref) and [`SimplifiedLaplace`](@ref) so that
method tables in `inference/marginals.jl` and `inference/inla.jl` can
dispatch on `::FullLaplace` at load time.

`n_grid` and `span` control the per-coordinate refit grid used for the
*integration-stage* summary replacement (`INLA(latent_strategy =
FullLaplace())`): each coordinate's `x_mean[i]` / `x_var[i]` becomes
the trapezoid moment of the refitted-Laplace mixture evaluated on
`n_grid` points spanning `±span` pilot standard deviations. The
per-coordinate density accessor [`posterior_marginal_x`](@ref) keeps
its own `grid_size` / `span` keywords (defaults 75 / 5.0) — the struct
fields do not affect it (ADR-026).

Cost of the integration-stage replacement: up to
`n_x × n_grid × n_design_points` warm-started Newton refits per fit
(tail truncation prunes most grid points). Opt-in, for sharply
non-Gaussian latent posteriors.

See [`AbstractMarginalStrategy`](@ref) and ADR-026.
"""
Base.@kwdef struct FullLaplace <: AbstractMarginalStrategy
    n_grid::Int = 51
    span::Float64 = 5.0

    function FullLaplace(n_grid::Int, span::Float64)
        n_grid ≥ 3 ||
            throw(ArgumentError("FullLaplace: n_grid must be ≥ 3, got $n_grid"))
        span > 0 ||
            throw(ArgumentError("FullLaplace: span must be > 0, got $span"))
        return new(n_grid, span)
    end
end

"""
    _resolve_marginal_strategy(s) -> AbstractMarginalStrategy

Backwards-compatibility shim accepting either an
[`AbstractMarginalStrategy`](@ref) instance (returned as-is) or a symbol
from the legacy whitelist (`:gaussian`, `:simplified_laplace`,
`:full_laplace`). Mirrors `_resolve_scheme(::Symbol, ::Int)` for the
integration-scheme keyword.

Throws `ArgumentError` for unknown symbols.
"""
_resolve_marginal_strategy(s::AbstractMarginalStrategy) = s
function _resolve_marginal_strategy(s::Symbol)
    s === :gaussian && return Gaussian()
    s === :simplified_laplace && return SimplifiedLaplace()
    s === :full_laplace && return FullLaplace()
    throw(ArgumentError("unknown marginal strategy :$s; " *
                        "use Gaussian(), SimplifiedLaplace(), FullLaplace(), " *
                        "or :gaussian / :simplified_laplace / :full_laplace"))
end
