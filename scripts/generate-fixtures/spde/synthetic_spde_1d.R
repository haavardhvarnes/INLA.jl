# Synthetic 1D Matérn time series — R-INLA reference fit for
# `SPDE1D{α=2}` (ν = 1.5).
#
# Model:
#   y_i = β + u(t_i) + ε_i,    i = 1,…,n
#   u   ~ SPDE-Matérn 1D, α = 2  (ν = 1.5)
#   ε_i ~ N(0, 1/τ_obs)
#   β   ~ N(0, 1000)            — `prec.intercept = 1e-3`
#   1/τ_obs ~ PC-prec(1, 0.01)
#   range ~ PC(0.5, 0.5),  σ ~ PC(1, 0.5)
#
# True hyperparameters used to simulate `u`:
#   ρ_true = 1.5, σ_true = 1.0, β_true = 1.5, τ_obs_true = 4.0.
#
# Output: fixtures/spde/synthetic_spde_1d.json
# (→ packages/INLASPDE.jl/test/oracle/fixtures/synthetic_spde_1d.jld2)

here <- tryCatch(
    normalizePath(dirname(sys.frame(1)$ofile), mustWork = FALSE),
    error = function(e) getwd()
)
if (!nzchar(here)) here <- getwd()
source(file.path(here, "..", "_helpers.R"))

suppressPackageStartupMessages({
    library(INLA)
    library(fmesher)
    library(Matrix)
})

set.seed(20260505)

# ----- truth + simulation ----------------------------------------
n           <- 200
t_obs       <- sort(runif(n, 0, 10))
beta_true   <- 1.5
rho_true    <- 1.5     # range
sigma_true  <- 1.0     # marginal SD
tau_obs_true <- 4.0    # noise precision

# Closed-form 1D Matérn covariance with ν = 1.5:
#   C(r; κ, σ) = σ² · (1 + κr) · exp(-κr),   κ = √(8ν)/ρ = √12/ρ.
kappa_true <- sqrt(12.0) / rho_true
r_mat <- abs(outer(t_obs, t_obs, "-"))
Sigma  <- sigma_true^2 * (1 + kappa_true * r_mat) * exp(-kappa_true * r_mat)
# Symmetrise + jitter for chol stability.
Sigma <- (Sigma + t(Sigma)) / 2 + diag(1.0e-10, n)
L_chol <- chol(Sigma)                # upper-triangular: Σ = Lᵀ L
u_obs  <- as.numeric(crossprod(L_chol, rnorm(n)))

y <- beta_true + u_obs + rnorm(n, 0, sd = 1 / sqrt(tau_obs_true))

# ----- fitting mesh + SPDE ---------------------------------------
mesh <- fmesher::fm_mesh_1d(seq(0, 10, length.out = 51))
spde <- inla.spde2.pcmatern(
    mesh = mesh, alpha = 2,
    prior.range = c(0.5, 0.5),    # P(range < 0.5) = 0.5
    prior.sigma = c(1, 0.5)       # P(σ > 1)      = 0.5
)
A <- inla.spde.make.A(mesh = mesh, loc = t_obs)

stk <- inla.stack(
    data = list(y = y),
    A = list(A, 1),
    effects = list(
        field = seq_len(spde$n.spde),
        list(intercept = rep(1, n))
    ),
    tag = "est"
)

form <- y ~ 0 + intercept + f(field, model = spde)

fit <- INLA::inla(
    form, family = "gaussian",
    data = inla.stack.data(stk),
    control.predictor = list(A = inla.stack.A(stk), compute = FALSE),
    control.fixed = list(prec.intercept = 1.0e-3, prec = 1.0e-3),
    control.family = list(hyper = list(prec = list(
        prior = "pc.prec", param = c(1, 0.01)
    ))),
    control.compute = list(return.marginals = FALSE)
)

# ----- emit fixture ----------------------------------------------
out_path <- file.path(here, "..", "fixtures", "spde", "synthetic_spde_1d.json")
write_inla_fixture(
    fit = fit,
    path = out_path,
    name = "synthetic_spde_1d",
    component_names = character(0),
    include_marginals = FALSE,
    meta = list(
        n = n,
        family = "gaussian",
        formula = "y ~ intercept + f(field, model=spde-1d)",
        alpha = 2,
        true = list(
            beta = beta_true, rho = rho_true, sigma = sigma_true,
            tau_obs = tau_obs_true
        ),
        prior_range = c(0.5, 0.5),
        prior_sigma = c(1, 0.5),
        prec_prior = "pc.prec(1, 0.01)",
        fixed_prec = 1.0e-3
    )
)

# Append data + mesh vertices + projector for Julia-side rebuild.
fixture <- jsonlite::fromJSON(out_path, simplifyVector = FALSE)
fixture$input <- list(
    y = as.numeric(y),
    t = as.numeric(t_obs)
)
fixture$mesh_1d <- list(
    points = as.numeric(mesh$loc),
    n_vertices = as.integer(length(mesh$loc))
)
fixture$A_field <- sparse_to_triplet(A)

jsonlite::write_json(
    fixture, out_path,
    auto_unbox = TRUE, digits = 16, pretty = FALSE, na = "null"
)

cat("wrote", out_path, "\n")
