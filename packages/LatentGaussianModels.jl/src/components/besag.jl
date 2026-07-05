"""
    Besag(graph; hyperprior = PCPrecision(), scale_model = true)

Intrinsic CAR (Besag) component on the supplied `graph`
(`AbstractGMRFGraph` or adjacency matrix). One hyperparameter on the
internal scale `θ = log(τ)`.

With `scale_model = true` (default) the Sørbye-Rue (2014) geometric-
mean scaling is applied, matching R-INLA ≥ 17.06. Carries one
sum-to-zero constraint per connected component of `graph`
(Freni-Sterrantino et al. 2018).
"""
# Per-node Sørbye-Rue (2014) scale constants on `graph` — the scaling of the
# connected component containing each node (Freni-Sterrantino et al. 2018),
# or all-ones when unscaled. Depends only on the graph and is θ-independent,
# so intrinsic scaled components (Besag, BYM) compute it once at construction
# and reuse it: `per_component_scale_factors` builds a dense `inv(Qperp)` per
# component, which otherwise dominates `precision_matrix` allocation.
function _sorbye_rue_scale_constants(g::GMRFs.AbstractGMRFGraph, scale_model::Bool)
    scale_model || return ones(Float64, GMRFs.num_nodes(g))
    c_k = GMRFs.per_component_scale_factors(g)
    labels = GMRFs.connected_component_labels(g)
    return Float64[c_k[labels[i]] for i in 1:GMRFs.num_nodes(g)]
end

struct Besag{P <: AbstractHyperPrior, G <: GMRFs.AbstractGMRFGraph} <:
       AbstractLatentComponent
    graph::G
    hyperprior::P
    scale_model::Bool
    scale_constants::Vector{Float64}   # cached per-node Sørbye-Rue constants
end

function Besag(graph::GMRFs.AbstractGMRFGraph;
        hyperprior::AbstractHyperPrior=PCPrecision(),
        scale_model::Bool=true)
    sc = _sorbye_rue_scale_constants(graph, scale_model)
    return Besag(graph, hyperprior, scale_model, sc)
end

Besag(W::AbstractMatrix; kwargs...) = Besag(GMRFs.GMRFGraph(W); kwargs...)

Base.length(c::Besag) = GMRFs.num_nodes(c.graph)
nhyperparameters(::Besag) = 1
initial_hyperparameters(::Besag) = [0.0]

function gmrf(c::Besag, θ)
    return GMRFs.BesagGMRF(c.graph; τ=exp(θ[1]), scale_model=c.scale_model,
        c=c.scale_constants)
end

precision_matrix(c::Besag, θ) = GMRFs.precision_matrix(gmrf(c, θ))
log_hyperprior(c::Besag, θ) = log_prior_density(c.hyperprior, θ[1])
GMRFs.constraints(c::Besag) = GMRFs.sum_to_zero_constraints(c.graph)

# Intrinsic CAR (rank `n - K` where K = # connected components).
# Structural `½ log|R̃|_+` dropped per R-INLA convention.
function log_normalizing_constant(c::Besag, θ)
    n = GMRFs.num_nodes(c.graph)
    K = GMRFs.nconnected_components(c.graph)
    return 0.5 * (n - K) * θ[1]
end
