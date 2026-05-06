"""
    Exceedance(c)

Wrapper requesting per-cell posterior exceedance probabilities
`P(u(s) > c | y)` from the sample-based [`predict_raster`](@ref)
overload.

Probability functionals do not compose linearly under barycentric
averaging, so exceedance rasters cannot be derived from Gaussian-
approximation summaries — the honest path is sample-based: draw
joint posterior samples of the latent SPDE field, project each draw
through the same mesh→raster projector, and reduce with
`mean(η_samples .> c, dims = 2)`.

# Example

```julia
rng = Xoshiro(1234)
ex = predict_raster(rng, model, res, template;
    component = "SPDE2[3]", quantity = Exceedance(log(500.0)),
    n_samples = 2000)
```
"""
struct Exceedance{T <: Real}
    c::T
end
