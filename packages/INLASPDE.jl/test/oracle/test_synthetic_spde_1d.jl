# Synthetic 1D Matérn SPDE oracle — Tier 2 against R-INLA.
#
# Model (matches scripts/generate-fixtures/spde/synthetic_spde_1d.R):
#   y_i = β + u(t_i) + ε_i,    i = 1,…,n,  t_i ∈ [0, 10]
#   u   ~ SPDE-Matérn 1D, α = 2  (ν = 1.5)
#   ε   ~ N(0, 1/τ_obs),  1/τ_obs ~ PC-prec(1, 0.01)
#   β   ~ N(0, 1000)
#   range ~ PC(0.5, 0.5),  σ ~ PC(1, 0.5)
#
# The fixture ships the fmesher mesh (51 vertices uniformly on [0, 10])
# and the fmesher-built spatial projector A; the Julia side rebuilds
# the LGM directly from those, so the oracle isolates SPDE1D + INLA
# integration from any 1D-mesher parity gap.

using JLD2
using SparseArrays
using LatentGaussianModels:
                            LatentGaussianModel, GaussianLikelihood, PCPrecision,
                            Intercept, inla, fixed_effects, hyperparameters

@testset "Synthetic 1D Matérn SPDE — posterior agreement with R-INLA" begin
    fxt = load(joinpath(@__DIR__, "fixtures", "synthetic_spde_1d.jld2"))["fixture"]

    # --- unpack -----------------------------------------------------
    y = Float64.(fxt["input"]["y"])
    points = Float64.(fxt["mesh_1d"]["points"])
    A_field = SparseMatrixCSC{Float64, Int}(fxt["A_field"])

    n_obs = length(y)
    n_v = length(points)

    @test size(A_field) == (n_obs, n_v)

    # --- Julia-side model ------------------------------------------
    # Build a 1D mesh from the same vertices the R-INLA fit used. The
    # FEM matrices then match (up to floating-point) what R-INLA
    # assembled.
    segments = hcat(1:(n_v - 1), 2:n_v)
    spde = SPDE1D(points, segments; α=2,
        pc=PCMatern{1}(
            range_U=0.5, range_α=0.5,
            sigma_U=1.0, sigma_α=0.5
        ))

    # R-INLA's `prec.intercept = 1e-3` is a *proper* N(0, 1000)
    # intercept; opt in via `improper = false` to match.
    intercept = Intercept(prec=1.0e-3, improper=false)

    # Stack projector: x = [α, u(field)]
    A_intercept = ones(n_obs, 1)
    A = hcat(A_intercept, A_field)

    like = GaussianLikelihood(hyperprior=PCPrecision(1.0, 0.01))
    model = LatentGaussianModel(like, (intercept, spde), A)

    res = inla(model, y)

    # --- R-INLA reference ------------------------------------------
    sf_rows = fxt["summary_fixed"]["rownames"]
    sf_mean = Float64.(fxt["summary_fixed"]["mean"])
    sf_sd = Float64.(fxt["summary_fixed"]["sd"])

    sh_rows = fxt["summary_hyperpar"]["rownames"]
    sh_mean = Float64.(fxt["summary_hyperpar"]["mean"])

    # `summary_fixed` in this fixture has a single row "intercept";
    # `summary_fixed.rownames` may come back as a String (single row)
    # or a Vector. Coerce.
    fixed_names = sf_rows isa AbstractVector ? sf_rows : [sf_rows]
    r_intercept = sf_mean[findfirst(==("intercept"), fixed_names)]
    r_sd_intercept = sf_sd[findfirst(==("intercept"), fixed_names)]

    r_prec_noise = sh_mean[findfirst(
        ==("Precision for the Gaussian observations"), sh_rows)]
    r_range = sh_mean[findfirst(==("Range for field"), sh_rows)]
    r_sigma = sh_mean[findfirst(==("Stdev for field"), sh_rows)]

    # --- Julia posterior summaries ---------------------------------
    fe = fixed_effects(model, res)
    j_intercept = fe[1].mean
    j_sd_intercept = fe[1].sd

    θ̂ = res.θ̂
    τ_noise_hat = exp(θ̂[1])                    # noise precision at mode
    ρ_hat, σ_hat = spde_user_scale(spde, θ̂[2:3])

    # --- Tolerances ------------------------------------------------
    # Same band as `meuse_spde` — fixed-effect mean within 0.5×R-SD,
    # SDs within 25%, hyperparameter means within 25–30% (R-INLA
    # CCD vs Julia mode + grid integrate differently). Mlik within 2
    # nats.
    @test abs(j_intercept - r_intercept) ≤ 0.5 * r_sd_intercept
    @test abs(j_sd_intercept - r_sd_intercept) / r_sd_intercept ≤ 0.25

    @test abs(ρ_hat - r_range) / r_range ≤ 0.25
    @test abs(σ_hat - r_sigma) / r_sigma ≤ 0.25
    @test abs(τ_noise_hat - r_prec_noise) / r_prec_noise ≤ 0.30

    r_mlik = Float64(fxt["mlik"][1])
    @test abs(res.log_marginal - r_mlik) ≤ 2.0
end
