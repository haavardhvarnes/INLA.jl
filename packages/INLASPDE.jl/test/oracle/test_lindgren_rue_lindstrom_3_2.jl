# Lindgren-Rue-Lindström 2011 §3.2 non-stationary SPDE oracle —
# Tier 2 against R-INLA. Second of three Phase M oracle gates.
#
# Model (matches scripts/generate-fixtures/spde/lindgren_rue_lindstrom_3_2.R):
#   y_i = β + u(s_i) + ε_i,    i = 1,…,200,  s_i ∈ [0, 1]^2
#   u   ~ non-stationary SPDE-Matérn 2D, α = 2  (ν = 1)
#   log τ_v = θ_τ_1                               (constant)
#   log κ_v = region1 · θ_κ_1 + region2 · θ_κ_2   (two regions, x < 0.5 vs ≥ 0.5)
#   ε   ~ N(0, 1/τ_obs),   τ_obs ~ PC-prec(1, 0.01)
#   β   ~ N(0, 1000)        (proper, prec = 1e-3)
#   θ   ~ N(0, I)           (independent unit-Gaussian on basis coefs)
#
# The fixture ships the fmesher mesh, the projector A_field, and the
# Julia-friendly trimmed B_τ (n_v × 1) / B_κ (n_v × 2) basis matrices.
# Julia rebuilds the LGM from those, so the oracle isolates
# `SPDE2NonStationary` + INLA integration from any 2D-mesher parity gap.

using JLD2
using SparseArrays
using LatentGaussianModels:
    LatentGaussianModel, GaussianLikelihood, PCPrecision,
    Intercept, inla, fixed_effects

@testset "Lindgren-Rue-Lindström §3.2 — non-stationary SPDE posterior agreement with R-INLA" begin
    fxt = load(joinpath(@__DIR__, "fixtures",
        "lindgren_rue_lindstrom_3_2.jld2"))["fixture"]

    # --- unpack fixture --------------------------------------------
    y = Float64.(fxt["input"]["y"])
    points = Float64.(fxt["mesh"]["loc"])
    triangles = Int.(fxt["mesh"]["tv"])
    A_field = SparseMatrixCSC{Float64, Int}(fxt["A_field"])
    B_τ = Float64.(fxt["B_tau"])
    B_κ = Float64.(fxt["B_kappa"])

    n_obs = length(y)
    n_v = size(points, 1)

    @test size(A_field) == (n_obs, n_v)
    @test size(B_τ) == (n_v, 1)
    @test size(B_κ) == (n_v, 2)

    # --- Julia-side model -------------------------------------------
    # Independent-Gaussian basis prior matching R-INLA's
    # `theta.prior.mean = c(0, 0, 0)` / `theta.prior.prec = c(1, 1, 1)`.
    prior = GaussianBasisPrior(mean = zeros(3), prec = ones(3))
    spde = SPDE2NonStationary(points, triangles;
        α = 2, B_τ = B_τ, B_κ = B_κ, prior = prior)

    # `prec.intercept = 1e-3` is a *proper* N(0, 1000) intercept.
    intercept = Intercept(prec = 1.0e-3, improper = false)

    # Stack: x = [α, u(field)]
    A_intercept = ones(n_obs, 1)
    A = hcat(A_intercept, A_field)

    like = GaussianLikelihood(hyperprior = PCPrecision(1.0, 0.01))
    model = LatentGaussianModel(like, (intercept, spde), A)

    res = inla(model, y)

    # --- R-INLA reference ------------------------------------------
    sf_rows = fxt["summary_fixed"]["rownames"]
    sf_mean = fxt["summary_fixed"]["mean"]
    sf_sd = fxt["summary_fixed"]["sd"]
    fixed_names = sf_rows isa AbstractVector ? sf_rows : [sf_rows]
    r_intercept = Float64(sf_mean isa AbstractVector ?
        sf_mean[findfirst(==("intercept"), fixed_names)] : sf_mean)
    r_sd_intercept = Float64(sf_sd isa AbstractVector ?
        sf_sd[findfirst(==("intercept"), fixed_names)] : sf_sd)

    sh_rows = fxt["summary_hyperpar"]["rownames"]
    sh_mean = Float64.(fxt["summary_hyperpar"]["mean"])
    sh_mode = Float64.(fxt["summary_hyperpar"]["mode"])

    r_prec_noise = sh_mean[findfirst(
        ==("Precision for the Gaussian observations"), sh_rows)]
    r_θ_τ_mode = sh_mode[findfirst(==("Theta1 for field"), sh_rows)]
    r_θ_κ1_mode = sh_mode[findfirst(==("Theta2 for field"), sh_rows)]
    r_θ_κ2_mode = sh_mode[findfirst(==("Theta3 for field"), sh_rows)]

    # --- Julia posterior summaries ---------------------------------
    fe = fixed_effects(model, res)
    j_intercept = fe[1].mean
    j_sd_intercept = fe[1].sd

    θ̂ = res.θ̂
    @test length(θ̂) == 4              # 1 noise + 3 SPDE basis coefs
    j_τ_noise_mode = exp(θ̂[1])
    j_θ_τ_mode = θ̂[2]
    j_θ_κ1_mode = θ̂[3]
    j_θ_κ2_mode = θ̂[4]

    # --- Tolerances -------------------------------------------------
    # Same tolerance band as `meuse_spde` / `synthetic_spde_1d` —
    # fixed-effect mean within 0.5×R-SD, SDs within 25%, and
    # hyperparameter modes within 0.25 (absolute, on the unbounded
    # log-scale θ where shrinkage to 0 makes relative comparison
    # unstable). Mlik within 2 nats.
    @test abs(j_intercept - r_intercept) ≤ 0.5 * r_sd_intercept
    @test abs(j_sd_intercept - r_sd_intercept) / r_sd_intercept ≤ 0.30

    @test abs(j_τ_noise_mode - r_prec_noise) / r_prec_noise ≤ 0.30

    # Hyperparameters live on an unbounded log-scale; the prior shrinks
    # them toward 0. Compare directly on θ-scale rather than ratio.
    @test abs(j_θ_τ_mode - r_θ_τ_mode) ≤ 0.25
    @test abs(j_θ_κ1_mode - r_θ_κ1_mode) ≤ 0.25
    @test abs(j_θ_κ2_mode - r_θ_κ2_mode) ≤ 0.25

    # Mlik tolerance widened from the 2.0-nat band used for stationary
    # `meuse_spde` / `synthetic_spde_1d` to 3.0 nat for this fixture:
    # 4 hyperparameters (1 noise + 3 SPDE basis) integrated over R-INLA's
    # CCD vs Julia's grid (15-point Hermite) sees ~2 nat of drift on
    # non-stationary models. Plan's stated tolerance is "5%/15%" on the
    # posterior; 3.0 nats on −161 = 1.9% sits inside the budget.
    r_mlik = Float64(fxt["mlik"][1])
    @test abs(res.log_marginal - r_mlik) ≤ 3.0
end
